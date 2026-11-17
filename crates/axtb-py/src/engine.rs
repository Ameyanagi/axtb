// SPDX-License-Identifier: LGPL-3.0-or-later
//! The private engine of the Python interface (`axtb.compute`, `axtb.Session`): a driver
//! [`Calculator`] behind a mutex, with one committed wavefunction per lineage.
//!
//! The Python layer (`axtb/_qc.py`, `axtb/_session.py`) checks a `qcdata.ProgramInput`, turns its
//! keywords into the driver's flat request (JSON) and its task, and maps the response dict back
//! onto qcdata. This module only runs:
//!
//! - [`compute_once`]: one calculation as `axtb_driver::run` does it (the same threads rule,
//!   [`axtb_driver::task_threads`]), so `axtb.compute` gives the CLI's numbers bit for bit;
//! - [`Engine`]: a session's calculator. Each call names its lineage (the canonical JSON key of
//!   `qcdata.lineage_key`), the seeds to start from when the lineage has no history, and whether
//!   a success becomes the lineage's history. The history is installed with
//!   `Calculator::set_history` (the restart policy decides, as for one calculator per lineage),
//!   and only when the method does not already hold it.
//!
//! Every call releases the GIL while it runs. A driver error raises [`Error`] with the driver's
//! `kind` and `code`, its message and its fields (`z`; `iterations`, `energy_change`,
//! `residual`); a panic raises it with kind and code `panic` and closes the engine.

use axtb_driver::{Calculator, Request, Response, Task, Wavefunction};
use numpy::{AllowTypeChange, PyArray1, PyArrayLike2, PyArrayMethods};
use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::collections::HashMap;
use std::sync::Mutex;

pyo3::create_exception!(
    axtb,
    Error,
    PyException,
    "A failed engine call: `kind` (the driver's: input, unsupported_element, topology, \
     scf_not_converged, numerical, not_applicable; panic, closed), `code` (invalid_input, \
     unsupported_element, topology_failed, scf_not_converged, numerical_failure, not_applicable; \
     panic, session_closed), `message` and `fields` (a dict: `z`; `iterations`, `energy_change`, \
     `residual`)."
);

/// What an engine call failed with, before it becomes a Python exception (built without the
/// GIL).
pub(crate) struct Failure {
    kind: &'static str,
    code: &'static str,
    message: String,
    fields: Vec<(&'static str, FieldValue)>,
}

pub(crate) enum FieldValue {
    Int(usize),
    Float(f64),
}

impl Failure {
    fn of(kind: &'static str, code: &'static str, message: impl Into<String>) -> Self {
        Self { kind, code, message: message.into(), fields: Vec::new() }
    }

    /// A panic's message: the payload of `std::panic::catch_unwind`.
    fn panic(payload: Box<dyn std::any::Any + Send>) -> Self {
        let what = payload
            .downcast_ref::<&str>()
            .map(|s| s.to_string())
            .or_else(|| payload.downcast_ref::<String>().cloned())
            .unwrap_or_default();
        Self::of("panic", "panic", format!("internal error (a bug in axtb; please report it): {what}"))
    }

    pub(crate) fn into_py(self, py: Python<'_>) -> PyErr {
        let err = Error::new_err(self.message.clone());
        let value = err.value(py);
        let fields = PyDict::new(py);
        for (k, v) in &self.fields {
            // setting items of a new dict cannot fail
            let _ = match v {
                FieldValue::Int(i) => fields.set_item(k, i),
                FieldValue::Float(x) => fields.set_item(k, x),
            };
        }
        // setting attributes of a new exception instance cannot fail
        let _ = value.setattr("kind", self.kind);
        let _ = value.setattr("code", self.code);
        let _ = value.setattr("message", self.message);
        let _ = value.setattr("fields", fields);
        err
    }
}

impl From<axtb_driver::Error> for Failure {
    fn from(e: axtb_driver::Error) -> Self {
        let mut f = Failure::of(e.kind(), e.code(), e.to_string());
        match e {
            axtb_driver::Error::UnsupportedElement { z, .. } => f.fields.push(("z", FieldValue::Int(z))),
            axtb_driver::Error::ScfNotConverged { iterations, energy_change, residual } => {
                f.fields.push(("iterations", FieldValue::Int(iterations)));
                f.fields.push(("energy_change", FieldValue::Float(energy_change)));
                f.fields.push(("residual", FieldValue::Float(residual)));
            }
            _ => {}
        }
        f
    }
}

/// `f` with panics caught.
fn guard<T>(f: impl FnOnce() -> Result<T, Failure>) -> Result<T, Failure> {
    std::panic::catch_unwind(std::panic::AssertUnwindSafe(f)).unwrap_or_else(|p| Err(Failure::panic(p)))
}

/// Rows of three from an array-like of shape (`rows`, 3).
fn rows3(a: &PyArrayLike2<'_, f64, AllowTypeChange>, rows: usize, what: &str) -> Result<Vec<[f64; 3]>, Failure> {
    let v = a.as_array();
    if v.shape() != [rows, 3] {
        return Err(Failure::of(
            "input",
            "invalid_input",
            format!("{what} must have shape ({rows}, 3), not {:?}", v.shape()),
        ));
    }
    Ok(v.rows().into_iter().map(|r| [r[0], r[1], r[2]]).collect())
}

/// The driver request of `request` (JSON of the driver's request without `numbers` and
/// `positions`, which come as arrays) for `numbers` at `positions`.
fn parse_request(
    request: &str,
    numbers: Vec<usize>,
    positions: &PyArrayLike2<'_, f64, AllowTypeChange>,
) -> Result<Request, Failure> {
    let nat = numbers.len();
    let positions = rows3(positions, nat, "positions")?;
    let mut req: Request = serde_json::from_str(request)
        .map_err(|e| Failure::of("input", "invalid_input", format!("invalid input: the driver request: {e}")))?;
    req.numbers = numbers;
    req.positions = positions;
    Ok(req)
}

/// The task of `task` (JSON of the driver request's task fields: `gradient`, `hessian`,
/// `hessian_method`, ..., `bond_orders`, `overlap`, `hamiltonian`), with the structured orbitals
/// when `orbitals`.
fn parse_task(task: &str, orbitals: bool) -> Result<Task, Failure> {
    let req: Request = serde_json::from_str(task)
        .map_err(|e| Failure::of("input", "invalid_input", format!("invalid input: the driver task: {e}")))?;
    if !req.numbers.is_empty() || !req.positions.is_empty() {
        return Err(Failure::of("input", "invalid_input", "invalid input: a task has no numbers or positions"));
    }
    let mut t = Task::from_request(&req)?;
    t.orbitals = orbitals;
    Ok(t)
}

/// One calculation as `axtb_driver::run` runs it: the request's calculator on
/// `axtb_driver::task_threads` threads, then the task; `accuracy` sets the SCF accuracy of the
/// calculation after the method is built (`Calculator::set_scf_accuracy`, only when it differs
/// from the request's). Returns the response as `Engine.compute` does, with `threads` (used).
#[pyfunction]
#[pyo3(signature = (request, numbers, positions, task, accuracy=None, orbitals=false))]
pub(crate) fn compute_once<'py>(
    py: Python<'py>,
    request: &str,
    numbers: Vec<usize>,
    positions: PyArrayLike2<'py, f64, AllowTypeChange>,
    task: &str,
    accuracy: Option<f64>,
    orbitals: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let req = parse_request(request, numbers, &positions).map_err(|f| f.into_py(py))?;
    let task = parse_task(task, orbitals).map_err(|f| f.into_py(py))?;
    let run = py.allow_threads(|| {
        guard(|| {
            let threads = axtb_driver::task_threads(&req, &task);
            let mut calc = if req.threads.is_some_and(|t| t != threads) {
                Calculator::new(&Request { threads: Some(threads), ..req.clone() })?
            } else {
                Calculator::new(&req)?
            };
            if let Some(a) = accuracy.filter(|&a| a != req.accuracy.unwrap_or(1.0)) {
                calc.set_scf_accuracy(a)?;
            }
            let r = calc.compute(&req.positions, &task)?;
            Ok((r, threads))
        })
    });
    let (r, threads) = run.map_err(|f| f.into_py(py))?;
    let d = response_dict(py, r)?;
    d.set_item("threads", threads)?;
    Ok(d)
}

/// Which committed state the method holds.
#[derive(Clone, Debug, PartialEq)]
enum Holds {
    /// None: a fresh start (after building, a reset or an error).
    Nothing,
    /// The committed history of this lineage.
    Lineage(String),
    /// A state that is no lineage's history (after a request with `commit=False`).
    Other,
}

struct State {
    calc: Option<Calculator>,
    history: HashMap<String, Wavefunction>,
    holds: Holds,
    /// The SCF accuracy the method runs at now.
    accuracy: f64,
    /// Parsed tasks by their JSON.
    tasks: HashMap<(String, bool), Task>,
}

/// A session's calculator (private to `axtb`): the driver calculator of one request, built
/// once, and the committed wavefunction of each lineage.
#[pyclass(module = "axtb._axtb", frozen)]
pub(crate) struct Engine {
    state: Mutex<State>,
    method: &'static str,
    conventions: &'static str,
    numerics: String,
    natoms: usize,
    threads: usize,
}

impl Engine {
    /// Run `f` on the state without the GIL; a panic closes the engine.
    fn with_state<T: Send>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut State) -> Result<T, Failure> + Send,
    ) -> PyResult<T> {
        py.allow_threads(|| {
            let mut s = self.state.lock().unwrap_or_else(|p| p.into_inner());
            let out = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| f(&mut s)));
            match out {
                Ok(r) => r,
                Err(p) => {
                    s.calc = None;
                    s.history.clear();
                    s.holds = Holds::Nothing;
                    Err(Failure::panic(p))
                }
            }
        })
        .map_err(|f| f.into_py(py))
    }
}

fn closed() -> Failure {
    Failure::of("closed", "session_closed", "the session is closed")
}

#[pymethods]
impl Engine {
    /// Build the calculator of `request` (JSON of the driver's request without `numbers` and
    /// `positions`) for `numbers` at `positions` (bohr, N x 3); its `threads` are the engine's.
    #[new]
    fn new(
        py: Python<'_>,
        request: &str,
        numbers: Vec<usize>,
        positions: PyArrayLike2<'_, f64, AllowTypeChange>,
    ) -> PyResult<Self> {
        let req = parse_request(request, numbers, &positions).map_err(|f| f.into_py(py))?;
        let calc = py.allow_threads(|| guard(|| Ok(Calculator::new(&req)?))).map_err(|f| f.into_py(py))?;
        Ok(Self {
            method: calc.method().name(),
            conventions: calc.conventions(),
            numerics: calc.numerics().to_string(),
            natoms: req.numbers.len(),
            threads: calc.threads(),
            state: Mutex::new(State {
                calc: Some(calc),
                history: HashMap::new(),
                holds: Holds::Nothing,
                accuracy: req.accuracy.unwrap_or(1.0),
                tasks: HashMap::new(),
            }),
        })
    }

    /// The method (canonical driver name).
    #[getter]
    fn method(&self) -> &'static str {
        self.method
    }

    /// The method's conventions (canonical driver name).
    #[getter]
    fn conventions(&self) -> &'static str {
        self.conventions
    }

    /// The numerics identity of the calculator's results.
    #[getter]
    fn numerics(&self) -> &str {
        &self.numerics
    }

    /// The threads of the calculator's pool.
    #[getter]
    fn threads(&self) -> usize {
        self.threads
    }

    /// Run `task` (JSON of the driver request's task fields) at `positions` (bohr, N x 3) in the
    /// lineage `lineage` (its key). The SCF starts from the lineage's committed history, else
    /// from the first of `seeds` (keys) with one, else fresh, as the restart policy says.
    /// `accuracy`: the SCF accuracy of this call (set only when it differs from the last call's);
    /// `orbitals`: the structured orbitals in the response. On success with `commit`, the
    /// method's wavefunction after the call becomes the lineage's history. Returns the response
    /// dict with `history`: the key whose history was installed, or None.
    #[pyo3(signature = (positions, task, lineage, seeds=Vec::new(), commit=true, accuracy=None, orbitals=false))]
    #[allow(clippy::too_many_arguments)]
    fn compute<'py>(
        &self,
        py: Python<'py>,
        positions: PyArrayLike2<'py, f64, AllowTypeChange>,
        task: String,
        lineage: String,
        seeds: Vec<String>,
        commit: bool,
        accuracy: Option<f64>,
        orbitals: bool,
    ) -> PyResult<Bound<'py, PyDict>> {
        let pos = rows3(&positions, self.natoms, "positions").map_err(|f| f.into_py(py))?;
        let (r, from) = self.with_state(py, move |s| {
            let calc = s.calc.as_mut().ok_or_else(closed)?;
            let key = (task, orbitals);
            let task = match s.tasks.get(&key) {
                Some(t) => t.clone(),
                None => {
                    let t = parse_task(&key.0, orbitals)?;
                    s.tasks.insert(key, t.clone());
                    t
                }
            };
            // the history to start from: the lineage's own, else the first seed's
            let from =
                std::iter::once(&lineage).chain(seeds.iter()).find(|k| s.history.contains_key(k.as_str())).cloned();
            let want = from.clone().map_or(Holds::Nothing, Holds::Lineage);
            if want != s.holds {
                calc.set_history(from.as_ref().map(|k| &s.history[k]))?;
                s.holds = want;
            }
            if let Some(a) = accuracy.filter(|&a| a != s.accuracy) {
                calc.set_scf_accuracy(a)?;
                s.accuracy = a;
            }
            match calc.compute(&pos, &task) {
                Ok(r) => {
                    if commit {
                        match calc.wavefunction() {
                            Some(wf) => {
                                s.history.insert(lineage.clone(), wf);
                                s.holds = Holds::Lineage(lineage);
                            }
                            None => s.holds = Holds::Nothing,
                        }
                    } else {
                        s.holds = Holds::Other;
                    }
                    Ok((r, from))
                }
                Err(e) => {
                    // the method holds no wavefunction after an error
                    s.holds = Holds::Nothing;
                    Err(e.into())
                }
            }
        })?;
        let d = response_dict(py, r)?;
        d.set_item("history", from)?;
        d.set_item("threads", self.threads)?;
        Ok(d)
    }

    /// Forget the history of the lineage `key`, or of every lineage (None).
    #[pyo3(signature = (key=None))]
    fn reset(&self, py: Python<'_>, key: Option<String>) -> PyResult<()> {
        self.with_state(py, move |s| {
            match key {
                Some(k) => {
                    s.history.remove(&k);
                }
                None => s.history.clear(),
            }
            Ok(())
        })
    }

    /// The keys of the lineages with history, in no particular order.
    fn lineages(&self, py: Python<'_>) -> PyResult<Vec<String>> {
        self.with_state(py, |s| Ok(s.history.keys().cloned().collect()))
    }

    /// Drop the calculator and every lineage's history; later calls raise kind `closed`.
    fn close(&self, py: Python<'_>) -> PyResult<()> {
        self.with_state(py, |s| {
            s.calc = None;
            s.history.clear();
            s.holds = Holds::Nothing;
            Ok(())
        })
    }

    /// Whether the engine is closed (by `close` or a panic).
    #[getter]
    fn closed(&self, py: Python<'_>) -> PyResult<bool> {
        self.with_state(py, |s| Ok(s.calc.is_none()))
    }
}

fn array1<'py>(py: Python<'py>, v: &[f64]) -> Bound<'py, PyArray1<f64>> {
    PyArray1::from_slice(py, v)
}

/// An (n, K) array of rows.
fn rows<'py, const K: usize>(py: Python<'py>, rows: &[[f64; K]]) -> PyResult<Bound<'py, PyAny>> {
    let flat: Vec<f64> = rows.iter().flatten().copied().collect();
    Ok(PyArray1::from_vec(py, flat).reshape([rows.len(), K])?.into_any())
}

/// A square (n, n) array of row-major `data`.
fn square(py: Python<'_>, data: Vec<f64>) -> PyResult<Bound<'_, PyAny>> {
    let n = (data.len() as f64).sqrt().round() as usize;
    Ok(PyArray1::from_vec(py, data).reshape([n, n])?.into_any())
}

/// Named values as a dict in their order.
fn named<'py>(py: Python<'py>, n: &axtb_driver::Named) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    for (k, v) in &n.0 {
        d.set_item(k, v)?;
    }
    Ok(d)
}

/// A serializable value as Python objects (through JSON: the small structured groups).
fn json<'py, T: serde::Serialize>(py: Python<'py>, v: &T) -> PyResult<Bound<'py, PyAny>> {
    let value = serde_json::to_value(v)
        .map_err(|e| Failure::of("numerical", "numerical_failure", format!("response: {e}")).into_py(py))?;
    Ok(crate::to_py(py, &value)?.into_bound(py))
}

/// The response as the Python layer reads it: the response's fields under their CLI names,
/// `scf`, `gfnff`, `hessian`, `vibrations`, `composite` and `orbitals` as nested dicts, arrays
/// as numpy arrays (`gradient` N x 3, `hessian` 3N x 3N, `atomic_quadrupoles` N x 6 as packed,
/// `normal_modes` modes x N x 3), the binding-only fields included (`subtotals`;
/// `core_electrons`, `shell_populations`, `shell_atoms`, `shell_angular_momenta`,
/// `reference_shell_occupations`, `quadrupole_traceless`; the thermochemistry's
/// `inverted_frequencies`, `dropped_frequencies`, `average_moment`).
fn response_dict(py: Python<'_>, r: Response) -> PyResult<Bound<'_, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("method", r.method)?;
    d.set_item("conventions", r.conventions)?;
    d.set_item("numerics", &r.numerics)?;
    d.set_item("eigensolver", &r.eigensolver)?;
    d.set_item("eigensolver_warning", axtb_driver::eigensolver_warning(&r))?;
    d.set_item("energy", r.energy)?;
    if let Some(g) = &r.gradient {
        d.set_item("gradient", rows(py, g)?)?;
    }
    if let Some(c) = &r.charges {
        d.set_item("charges", array1(py, c))?;
    }
    if let Some(m) = &r.dipole {
        d.set_item("dipole", array1(py, m))?;
    }
    if let Some(t) = &r.terms {
        d.set_item("terms", named(py, t)?)?;
    }
    d.set_item("subtotals", PyList::new(py, &r.subtotals)?)?;
    if let Some(t) = &r.timings {
        d.set_item("timings", named(py, t)?)?;
    }
    d.set_item("seconds", r.seconds)?;
    if let Some(s) = r.scf {
        let x = PyDict::new(py);
        x.set_item("iterations", s.iterations)?;
        x.set_item("guess", s.guess)?;
        x.set_item("mixer_used", s.mixer_used)?;
        x.set_item("scf_attempts", json(py, &s.scf_attempts)?)?;
        x.set_item("homo_lumo_gap", s.homo_lumo_gap)?;
        x.set_item("fermi_levels", (s.fermi_levels[0], s.fermi_levels[1]))?;
        x.set_item("electronic_temperature", s.electronic_temperature)?;
        x.set_item("quadrupole", array1(py, &s.quadrupole))?;
        x.set_item("quadrupole_traceless", array1(py, &s.quadrupole_traceless))?;
        x.set_item("shell_charges", array1(py, &s.shell_charges))?;
        x.set_item("atomic_dipoles", rows(py, &s.atomic_dipoles)?)?;
        x.set_item("atomic_quadrupoles", rows(py, &s.atomic_quadrupoles)?)?;
        x.set_item("orbital_energies_alpha", array1(py, &s.orbital_energies_alpha))?;
        x.set_item("orbital_energies_beta", array1(py, &s.orbital_energies_beta))?;
        x.set_item("occupations_alpha", array1(py, &s.occupations_alpha))?;
        x.set_item("occupations_beta", array1(py, &s.occupations_beta))?;
        if let Some(b) = s.bond_orders {
            x.set_item("bond_orders", square(py, b.data)?)?;
        }
        if let Some(o) = s.overlap {
            x.set_item("overlap", square(py, o.data)?)?;
        }
        if let Some(h) = s.hamiltonian {
            let l = PyList::empty(py);
            for m in h {
                l.append(square(py, m.data)?)?;
            }
            x.set_item("hamiltonian", l)?;
        }
        x.set_item("core_electrons", PyList::new(py, &s.core_electrons)?)?;
        x.set_item("shell_populations", array1(py, &s.shell_populations))?;
        x.set_item("shell_atoms", PyList::new(py, &s.shell_atoms)?)?;
        x.set_item("shell_angular_momenta", PyList::new(py, &s.shell_angular_momenta)?)?;
        x.set_item("reference_shell_occupations", array1(py, &s.reference_shell_occupations))?;
        d.set_item("scf", x)?;
    }
    if let Some(g) = &r.gfnff {
        d.set_item("gfnff", json(py, g)?)?;
    }
    if let Some(h) = r.hessian {
        let x = PyDict::new(py);
        x.set_item("hessian_method", h.hessian_method)?;
        x.set_item("hessian_max_asymmetry", h.hessian_max_asymmetry)?;
        x.set_item("hessian_timings", named(py, &h.hessian_timings)?)?;
        if let Some(dd) = &h.dipole_derivatives {
            x.set_item("dipole_derivatives", rows(py, dd)?)?;
        }
        x.set_item("hessian", square(py, h.hessian.data)?)?;
        d.set_item("hessian", x)?;
    }
    if let Some(v) = r.vibrations {
        let nat = v.masses.len();
        let x = PyDict::new(py);
        x.set_item("vibrational_conventions", v.vibrational_conventions)?;
        x.set_item("frequencies", array1(py, &v.frequencies))?;
        let flat: Vec<f64> = v.normal_modes.iter().flatten().flatten().copied().collect();
        x.set_item("normal_modes", PyArray1::from_vec(py, flat).reshape([v.normal_modes.len(), nat, 3])?)?;
        x.set_item("reduced_masses", array1(py, &v.reduced_masses))?;
        x.set_item("force_constants", array1(py, &v.force_constants))?;
        if let Some(ir) = &v.ir_intensities {
            x.set_item("ir_intensities", array1(py, ir))?;
        }
        x.set_item("n_imaginary", v.n_imaginary)?;
        x.set_item("linear", v.linear)?;
        x.set_item("masses", array1(py, &v.masses))?;
        if let Some(t) = &v.thermo {
            let tx = json(py, t)?.downcast_into::<PyDict>()?;
            tx.set_item("inverted_frequencies", PyList::new(py, &t.inverted_frequencies)?)?;
            tx.set_item("dropped_frequencies", PyList::new(py, &t.dropped_frequencies)?)?;
            tx.set_item("average_moment", t.average_moment)?;
            x.set_item("thermo", tx)?;
        }
        d.set_item("vibrations", x)?;
    }
    if let Some(c) = &r.composite {
        d.set_item("composite", json(py, c)?)?;
    }
    if let Some(o) = r.orbitals {
        d.set_item("orbitals", orbitals_dict(py, o)?)?;
    }
    Ok(d)
}

/// Structured orbitals: `shells` (per atom, per shell: `(l, exponents, coefficients)`, the
/// coefficients of `axtb::Cgto`: x^l-normalized primitives) and `sets` (per spin set: `spin`,
/// `energies`, `occupations` as the driver gives them, `coefficients` nmo x nao in axtb's AO
/// order).
fn orbitals_dict(py: Python<'_>, o: axtb::Orbitals) -> PyResult<Bound<'_, PyDict>> {
    let d = PyDict::new(py);
    let nao = o.nao();
    let shells = PyList::empty(py);
    for atom in &o.shells {
        let a = PyList::empty(py);
        for s in atom {
            a.append((s.l, array1(py, &s.alpha[..s.nprim]), array1(py, &s.coeff[..s.nprim])))?;
        }
        shells.append(a)?;
    }
    d.set_item("shells", shells)?;
    let sets = PyList::empty(py);
    for set in o.sets {
        let x = PyDict::new(py);
        x.set_item("spin", set.spin.label())?;
        x.set_item("energies", array1(py, &set.energies))?;
        x.set_item("occupations", array1(py, &set.occupations))?;
        let nmo = set.energies.len();
        x.set_item("coefficients", PyArray1::from_vec(py, set.coefficients).reshape([nmo, nao])?)?;
        sets.append(x)?;
    }
    d.set_item("sets", sets)?;
    Ok(d)
}

/// The canonical name of `name` of a kind: `method`, `conventions` of `method`, `solvent`
/// (`axtb_driver::canonical`).
#[pyfunction]
#[pyo3(signature = (kind, name, method=None))]
pub(crate) fn canonical(py: Python<'_>, kind: &str, name: &str, method: Option<&str>) -> PyResult<&'static str> {
    axtb_driver::canonical(kind, name, method).map_err(|e| Failure::from(e).into_py(py))
}

/// The atomic numbers `method` has parameters for.
#[pyfunction]
pub(crate) fn supported_elements(py: Python<'_>, method: &str) -> PyResult<Vec<usize>> {
    let m = axtb_driver::Method::parse(method).map_err(|e| Failure::from(e).into_py(py))?;
    Ok(py.allow_threads(|| axtb_driver::supported_elements(m)))
}

/// The numerics identity of `method` (`composite`: the g-xTB solvation composite).
#[pyfunction]
pub(crate) fn numerics(py: Python<'_>, method: &str, composite: bool) -> PyResult<String> {
    let m = axtb_driver::Method::parse(method).map_err(|e| Failure::from(e).into_py(py))?;
    Ok(axtb_driver::numerics(m, composite))
}

/// The electrons a neutral atom brings to the orbitals of `method`, for Z = 1..=118 (index Z - 1;
/// None beyond the method's tables): `axtb_driver::valence_electrons`.
#[pyfunction]
pub(crate) fn valence_electrons(py: Python<'_>, method: &str) -> PyResult<Vec<Option<f64>>> {
    let m = axtb_driver::Method::parse(method).map_err(|e| Failure::from(e).into_py(py))?;
    Ok((1..=118).map(|z| axtb_driver::valence_electrons(m, z)).collect())
}

/// The largest molecule (atoms) whose single point (`hessian` false) or Hessian of `method`
/// runs on one thread (`axtb_driver::single_thread_atoms`).
#[pyfunction]
pub(crate) fn single_thread_atoms(py: Python<'_>, method: &str, hessian: bool) -> PyResult<usize> {
    let m = axtb_driver::Method::parse(method).map_err(|e| Failure::from(e).into_py(py))?;
    Ok(axtb_driver::single_thread_atoms(m, hessian))
}

/// The solvent names of the implicit-solvent model `model` (`alpb` or `gbsa`), as its
/// parameter table names them.
#[pyfunction]
pub(crate) fn solvents(py: Python<'_>, model: &str) -> PyResult<Vec<&'static str>> {
    Ok(axtb_solvation::solvent_names(solvation_model(py, model)?))
}

/// The table name of the solvent `name` for the model `model` (`alpb` or `gbsa`), or None.
#[pyfunction]
pub(crate) fn canonical_solvent(py: Python<'_>, name: &str, model: &str) -> PyResult<Option<&'static str>> {
    Ok(axtb_solvation::canonical_solvent(name, solvation_model(py, model)?))
}

fn solvation_model(py: Python<'_>, model: &str) -> PyResult<axtb_solvation::Model> {
    match axtb_driver::names::key(model).as_str() {
        "alpb" => Ok(axtb_solvation::Model::Alpb),
        "gbsa" => Ok(axtb_solvation::Model::Gbsa),
        _ => Err(Failure::of("input", "invalid_input", format!("unknown solvation model {model:?}")).into_py(py)),
    }
}

/// Why `AXTB_LAPACK` cannot be used, or None (loads the library if no SCF has).
#[pyfunction]
pub(crate) fn lapack_warning(py: Python<'_>) -> Option<String> {
    py.allow_threads(axtb_driver::lapack_warning)
}
