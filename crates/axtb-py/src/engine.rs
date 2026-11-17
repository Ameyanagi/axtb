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
//!   and only when the method does not already hold it. Each call runs on the threads it names,
//!   on the calculator or a fork of it on those threads (`Calculator::fork`, made when first
//!   needed and kept). `compute_many` runs a batch's lineages side by side as their threads fit
//!   into the threads it is given, each lineage's requests in order, and gives what one call at
//!   a time in request order gives, bit for bit.
//!
//! Every call releases the GIL while it runs. A driver error raises [`Error`] with the driver's
//! `kind` and `code`, its message, its fields (`z`; `iterations`, `energy_change`, `residual`)
//! and, for an SCF that did not converge on its `scf_retries`, the ladder's `attempts`; a panic
//! raises it with kind and code `panic` and closes the engine.

use axtb_driver::{Calculator, Request, Response, Task, Wavefunction};
use numpy::{AllowTypeChange, PyArray1, PyArrayLike2, PyArrayMethods};
use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::collections::HashMap;
use std::sync::{Arc, Condvar, Mutex};
use std::time::Instant;

pyo3::create_exception!(
    axtb,
    Error,
    PyException,
    "A failed engine call: `kind` (the driver's: input, unsupported_element, topology, \
     scf_not_converged, numerical, not_applicable; structure, panic, closed), `code` (invalid_input, \
     unsupported_element, topology_failed, scf_not_converged, numerical_failure, not_applicable; \
     invalid_structure, panic, session_closed), `message`, `fields` (a dict: `z`; `iterations`, `energy_change`, \
     `residual`) and `attempts` (the tries of an `scf_retries` ladder that did not converge, as dicts; else None)."
);

/// What an engine call failed with, before it becomes a Python exception (built without the
/// GIL).
pub(crate) struct Failure {
    kind: &'static str,
    code: &'static str,
    message: String,
    fields: Vec<(&'static str, FieldValue)>,
    /// The tries of an `scf_retries` ladder (JSON), if any.
    attempts: Option<serde_json::Value>,
}

pub(crate) enum FieldValue {
    Int(usize),
    Float(f64),
}

impl Failure {
    fn of(kind: &'static str, code: &'static str, message: impl Into<String>) -> Self {
        Self { kind, code, message: message.into(), fields: Vec::new(), attempts: None }
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
        let attempts = self.attempts.as_ref().and_then(|a| crate::to_py(py, a).ok()).unwrap_or_else(|| py.None());
        let _ = value.setattr("attempts", attempts);
        err
    }
}

impl From<axtb_driver::Error> for Failure {
    fn from(e: axtb_driver::Error) -> Self {
        let mut f = Failure::of(e.kind(), e.code(), e.to_string());
        match e {
            axtb_driver::Error::UnsupportedElement { z, .. } => f.fields.push(("z", FieldValue::Int(z))),
            axtb_driver::Error::ScfNotConverged { iterations, energy_change, residual, attempts } => {
                f.fields.push(("iterations", FieldValue::Int(iterations)));
                f.fields.push(("energy_change", FieldValue::Float(energy_change)));
                f.fields.push(("residual", FieldValue::Float(residual)));
                if !attempts.is_empty() {
                    f.attempts = serde_json::to_value(&attempts).ok();
                }
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
            "structure",
            "invalid_structure",
            format!("{what}: {rows} atoms x 3, not shape {:?}", v.shape()),
        ));
    }
    if v.iter().any(|x| !x.is_finite()) {
        return Err(Failure::of("structure", "invalid_structure", format!("{what}: must be finite numbers")));
    }
    Ok(v.rows().into_iter().map(|r| [r[0], r[1], r[2]]).collect())
}

/// The 3 x 3 tensor of a packed one (xx, xy, yy, xz, yz, zz), values copied.
fn tensor(q: &[f64; 6]) -> [[f64; 3]; 3] {
    let [xx, xy, yy, xz, yz, zz] = *q;
    [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]
}

/// An (n, 3, 3) array of packed tensors.
fn tensors<'py>(py: Python<'py>, packed: &[[f64; 6]]) -> PyResult<Bound<'py, PyAny>> {
    let flat: Vec<f64> = packed.iter().flat_map(|q| tensor(q).into_iter().flatten()).collect();
    Ok(PyArray1::from_vec(py, flat).reshape([packed.len(), 3, 3])?.into_any())
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

/// Densities (row-major n x n each) of a list of arrays.
fn densities(d: &[PyArrayLike2<'_, f64, AllowTypeChange>]) -> Vec<Vec<f64>> {
    d.iter().map(|m| m.as_array().iter().copied().collect()).collect()
}

/// GFN2-xTB's start moments: shell charges, atomic dipoles (N x 3), packed quadrupoles (N x 6).
type Moments = (Vec<f64>, Vec<[f64; 3]>, Vec<[f64; 6]>);

fn moments(
    m: &(Vec<f64>, PyArrayLike2<'_, f64, AllowTypeChange>, PyArrayLike2<'_, f64, AllowTypeChange>),
    natoms: usize,
) -> Result<Moments, Failure> {
    let dp = rows3(&m.1, natoms, "atomic dipoles")?;
    let q = m.2.as_array();
    if q.shape() != [natoms, 6] {
        return Err(Failure::of("input", "invalid_input", "atomic quadrupoles: N x 6 (packed)"));
    }
    let qp = q.rows().into_iter().map(|r| std::array::from_fn(|k| r[k])).collect();
    Ok((m.0.clone(), dp, qp))
}

/// Install a given start on a calculator: densities (`Calculator::set_density`) or moments
/// (`Calculator::set_moments`).
fn install_given(calc: &mut Calculator, d: Option<Vec<Vec<f64>>>, m: Option<Moments>) -> Result<(), Failure> {
    if let Some(d) = d {
        calc.set_density(d)?;
    }
    if let Some((q, dp, qp)) = m {
        calc.set_moments(q, dp, qp)?;
    }
    Ok(())
}

/// One calculation as `axtb_driver::run` runs it: the request's calculator on
/// `axtb_driver::task_threads` threads, then the task; `accuracy` sets the SCF accuracy of the
/// calculation after the method is built (`Calculator::set_scf_accuracy`, only when it differs
/// from the request's). Returns the response as `Engine.compute` does, with `threads` (used).
#[pyfunction]
#[pyo3(signature = (request, numbers, positions, task, accuracy=None, orbitals=false, densities=None, moments=None))]
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub(crate) fn compute_once<'py>(
    py: Python<'py>,
    request: &str,
    numbers: Vec<usize>,
    positions: PyArrayLike2<'py, f64, AllowTypeChange>,
    task: &str,
    accuracy: Option<f64>,
    orbitals: bool,
    densities: Option<Vec<PyArrayLike2<'py, f64, AllowTypeChange>>>,
    moments: Option<(Vec<f64>, PyArrayLike2<'py, f64, AllowTypeChange>, PyArrayLike2<'py, f64, AllowTypeChange>)>,
) -> PyResult<Bound<'py, PyDict>> {
    let nat = numbers.len();
    let req = parse_request(request, numbers, &positions).map_err(|f| f.into_py(py))?;
    let task = parse_task(task, orbitals).map_err(|f| f.into_py(py))?;
    let given_d = densities.as_deref().map(self::densities);
    let given_m = moments.as_ref().map(|m| self::moments(m, nat)).transpose().map_err(|f| f.into_py(py))?;
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
            install_given(&mut calc, given_d, given_m)?;
            let r = calc.compute(&req.positions, &task)?;
            Ok((r, threads))
        })
    });
    let (r, threads) = run.map_err(|f| f.into_py(py))?;
    let d = response_dict(py, r)?;
    d.set_item("threads", threads)?;
    Ok(d)
}

/// Which committed state a calculator's method holds.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Holds {
    /// None: a fresh start (after building, a reset or an error).
    Nothing,
    /// The committed history with this id.
    Commit(u64),
    /// A state that is no lineage's history (after a request with `commit=False`, or the given
    /// start).
    Other,
}

/// A lineage's committed history: its id (unique in the engine) and the wavefunction. It does
/// not depend on the threads of the calculator that made it.
type Commit = (u64, Arc<Wavefunction>);

/// A calculator of the engine, what its method holds, and the SCF accuracy it runs at now.
struct Worker {
    calc: Calculator,
    holds: Holds,
    accuracy: f64,
}

impl Worker {
    /// The threads of its pool.
    fn threads(&self) -> usize {
        self.calc.threads()
    }
}

struct State {
    /// The engine's calculator, then its forks on the threads requests asked for (forked when
    /// first needed and kept until the engine closes); empty once closed.
    workers: Vec<Worker>,
    history: HashMap<String, Commit>,
    /// The id of the next commit.
    next: u64,
    /// Parsed tasks by their JSON.
    tasks: HashMap<(String, bool), Task>,
    /// The given start of every request without history (`set_given`), as the method holds it.
    given: Option<Wavefunction>,
}

impl State {
    /// Drop the calculators and every lineage's history.
    fn close(&mut self) {
        self.workers.clear();
        self.history.clear();
    }

    /// The task of `json` (parsed once).
    fn task(&mut self, json: String, orbitals: bool) -> Result<Task, Failure> {
        let key = (json, orbitals);
        if let Some(t) = self.tasks.get(&key) {
            return Ok(t.clone());
        }
        let t = parse_task(&key.0, orbitals)?;
        self.tasks.insert(key, t.clone());
        Ok(t)
    }

    /// One more calculator: a fork of the engine's on `threads` threads, at its SCF accuracy.
    fn fork(&mut self, threads: usize) -> Result<usize, Failure> {
        let first = &self.workers[0];
        let worker = Worker { calc: first.calc.fork(threads)?, holds: Holds::Nothing, accuracy: first.accuracy };
        self.workers.push(worker);
        Ok(self.workers.len() - 1)
    }

    /// The calculator on `threads` threads for a call that starts from `want`: one that holds it,
    /// else the first on those threads, else a new fork.
    fn worker_for(&mut self, threads: usize, want: Holds) -> Result<usize, Failure> {
        let mut on = self.workers.iter().enumerate().filter(|(_, w)| w.threads() == threads);
        let first = on.next();
        match first.into_iter().chain(on).find(|(_, w)| w.holds == want).or(first) {
            Some((i, _)) => Ok(i),
            None => self.fork(threads),
        }
    }
}

/// One request: positions (bohr), task, the keys of its lineage and seeds, whether a success
/// becomes the lineage's history, its SCF accuracy and its threads.
struct Job {
    positions: Vec<[f64; 3]>,
    task: Task,
    lineage: String,
    seeds: Vec<String>,
    commit: bool,
    accuracy: Option<f64>,
    threads: usize,
}

/// A request's response, the key of the lineage whose history it started from, whether it
/// started from the given start, and its wall time on its calculator (s).
struct Done {
    response: Response,
    from: Option<String>,
    given: bool,
    seconds: f64,
}

/// Run `job` on `w`, starting from `from` (a lineage's history), else from `given`, else from
/// the method's own start, as the restart policy says. Returns the response and, for a commit,
/// the wavefunction to commit (its id in `w.holds` is the caller's to set).
fn run(
    w: &mut Worker,
    job: &Job,
    from: Option<&Commit>,
    given: Option<&Wavefunction>,
) -> Result<(Response, Option<Wavefunction>), Failure> {
    let calc = &mut w.calc;
    match (from, given) {
        (None, Some(g)) => {
            // the given start replaces the method's own one for every request without history
            calc.set_wavefunction(g)?;
            w.holds = Holds::Other;
        }
        _ => {
            let want = from.map_or(Holds::Nothing, |c| Holds::Commit(c.0));
            if want != w.holds {
                calc.set_history(from.map(|c| &*c.1))?;
                w.holds = want;
            }
        }
    }
    if let Some(a) = job.accuracy.filter(|&a| a != w.accuracy) {
        calc.set_scf_accuracy(a)?;
        w.accuracy = a;
    }
    match calc.compute(&job.positions, &job.task) {
        Ok(r) => {
            let wf = if job.commit { calc.wavefunction() } else { None };
            w.holds = if job.commit { Holds::Nothing } else { Holds::Other };
            Ok((r, wf))
        }
        Err(e) => {
            // the method holds no wavefunction after an error
            w.holds = Holds::Nothing;
            Err(e.into())
        }
    }
}

/// The lineages of a batch: each request's lineage and seeds (indices into `keys`), and each
/// lineage's requests in order.
struct Lineages {
    keys: Vec<String>,
    lineage: Vec<usize>,
    seeds: Vec<Vec<usize>>,
    requests: Vec<Vec<usize>>,
}

/// Where a request starts: a lineage (its index) and its history, or nothing.
type Start = Option<(usize, Commit)>;

impl Lineages {
    fn new(jobs: &[Job]) -> Self {
        let mut ids: HashMap<String, usize> = HashMap::new();
        let mut keys = Vec::new();
        let mut id = |k: &str| {
            *ids.entry(k.to_string()).or_insert_with(|| {
                keys.push(k.to_string());
                keys.len() - 1
            })
        };
        let lineage: Vec<usize> = jobs.iter().map(|j| id(&j.lineage)).collect();
        let seeds: Vec<Vec<usize>> = jobs.iter().map(|j| j.seeds.iter().map(|k| id(k)).collect()).collect();
        let mut requests = vec![Vec::new(); keys.len()];
        for (i, &l) in lineage.iter().enumerate() {
            requests[l].push(i);
        }
        Self { keys, lineage, seeds, requests }
    }

    /// The requests of lineage `l` before request `i`.
    fn before(&self, l: usize, i: usize) -> &[usize] {
        let r = &self.requests[l];
        &r[..r.partition_point(|&j| j < i)]
    }

    /// Where request `i` starts, as one request at a time in request order would start it, once
    /// that is settled (`None` until then): its lineage's latest commit before it, else the
    /// first of its seeds with one (each seed's history as its requests before `i` leave it:
    /// those must have run), else nothing. A lineage's requests run in order, so a request of a
    /// lineage with history never waits for its seeds.
    fn start(&self, i: usize, p: &Progress, history: &HashMap<String, Commit>) -> Option<Start> {
        let ran = |l: usize| self.before(l, i).last().is_none_or(|&j| p.state[j] == Run::Done);
        // the latest commit of lineage `l` before request `i` (its earlier requests have run)
        let latest = |l: usize| {
            let committed = self.before(l, i).iter().rev().find_map(|&j| p.commits[j].clone());
            committed.or_else(|| history.get(&self.keys[l]).cloned())
        };
        let own = self.lineage[i];
        if !ran(own) {
            return None;
        }
        if let Some(c) = latest(own) {
            return Some(Some((own, c)));
        }
        for &s in &self.seeds[i] {
            if !ran(s) {
                return None;
            }
            if let Some(c) = latest(s) {
                return Some(Some((s, c)));
            }
        }
        Some(None)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Run {
    Waiting,
    Running,
    Done,
}

/// A batch in progress: each request's state, commit and outcome, and the threads in use.
struct Progress {
    state: Vec<Run>,
    commits: Vec<Option<Commit>>,
    out: Vec<Option<Result<Done, Failure>>>,
    /// The threads of the requests running, and the most that may run at once.
    busy: usize,
    total: usize,
    /// The id of the next commit.
    next: u64,
    /// A calculation panicked: the requests not yet started fail (`closed`), and the engine
    /// closes.
    panicked: bool,
}

impl Progress {
    /// Whether a request on `threads` threads fits beside the ones running (alone, it always
    /// does).
    fn fits(&self, threads: usize) -> bool {
        self.busy == 0 || self.busy + threads <= self.total
    }
}

/// The stack of a thread that runs a calculator of `compute_many` (a main thread's).
const WORKER_STACK: usize = 8 << 20;

/// One calculator's part of a batch: until every request on its threads has started, take one
/// that can start and fits into the threads in all (one that starts from the history `w` holds,
/// else the first in request order), run it and publish its outcome and commit.
fn work(
    w: &mut Worker,
    jobs: &[Job],
    lineages: &Lineages,
    history: &HashMap<String, Commit>,
    given: Option<&Wavefunction>,
    sync: &(Mutex<Progress>, Condvar),
) {
    let (lock, ready) = sync;
    let threads = w.threads();
    let mine = |i: &usize| jobs[*i].threads == threads;
    loop {
        let (i, start, panicked) = {
            let mut p = lock.lock().unwrap_or_else(|e| e.into_inner());
            loop {
                let mut waiting = (0..jobs.len()).filter(mine).filter(|&i| p.state[i] == Run::Waiting).peekable();
                if waiting.peek().is_none() {
                    return;
                }
                let mut pick: Option<(usize, Start)> = None;
                if p.fits(threads) {
                    for i in waiting {
                        let Some(start) = lineages.start(i, &p, history) else { continue };
                        let here = start.as_ref().map_or(Holds::Nothing, |(_, c)| Holds::Commit(c.0)) == w.holds;
                        if pick.is_none() || here {
                            pick = Some((i, start));
                        }
                        if here {
                            break;
                        }
                    }
                }
                if let Some((i, start)) = pick {
                    p.state[i] = Run::Running;
                    p.busy += threads;
                    break (i, start, p.panicked);
                }
                p = ready.wait(p).unwrap_or_else(|e| e.into_inner());
            }
        };
        let t0 = Instant::now();
        let from = start.as_ref().map(|(_, c)| c);
        let given = if start.is_none() { given } else { None };
        let (out, panic) = if panicked {
            (Err(closed()), false)
        } else {
            match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| run(w, &jobs[i], from, given))) {
                Ok(r) => (r, false),
                Err(payload) => (Err(Failure::panic(payload)), true),
            }
        };
        let seconds = t0.elapsed().as_secs_f64();
        let mut p = lock.lock().unwrap_or_else(|e| e.into_inner());
        let out = out.map(|(response, wf)| {
            if let Some(wf) = wf {
                let id = p.next;
                p.next += 1;
                p.commits[i] = Some((id, Arc::new(wf)));
                w.holds = Holds::Commit(id);
            }
            let from = start.as_ref().map(|(l, _)| lineages.keys[*l].clone());
            Done { response, from, given: given.is_some(), seconds }
        });
        p.out[i] = Some(out);
        p.state[i] = Run::Done;
        p.busy -= threads;
        p.panicked |= panic;
        drop(p);
        ready.notify_all();
    }
}

/// Run `jobs` as `compute` runs them one at a time in their order, each on its own threads
/// (`Job::threads`), different lineages side by side as long as their threads fit into `total`
/// (taken greedily in request order), one lineage on one calculator at a time; the outcomes in
/// order. The calculators are the engine's and its forks on each request's threads, forked
/// when first needed: per thread count as many as fit into `total` at once and its requests'
/// lineages need. Afterwards each lineage's history is its last commit; a panic closes the
/// engine.
fn batch(s: &mut State, jobs: Vec<Job>, total: usize) -> Vec<Result<Done, Failure>> {
    let lineages = Lineages::new(&jobs);
    let n = jobs.len();
    let mut progress = Progress {
        state: vec![Run::Waiting; n],
        commits: vec![None; n],
        out: (0..n).map(|_| None).collect(),
        busy: 0,
        total,
        next: s.next,
        panicked: false,
    };
    // the lineages of each thread count's requests
    let mut counts: Vec<(usize, Vec<usize>)> = Vec::new();
    for (i, j) in jobs.iter().enumerate() {
        let k = counts.iter().position(|(t, _)| *t == j.threads).unwrap_or_else(|| {
            counts.push((j.threads, Vec::new()));
            counts.len() - 1
        });
        if !counts[k].1.contains(&lineages.lineage[i]) {
            counts[k].1.push(lineages.lineage[i]);
        }
    }
    let mut chosen = Vec::new();
    for (threads, of) in &counts {
        let want = (total / threads).min(of.len()).max(1);
        let mut have: Vec<usize> =
            (0..s.workers.len()).filter(|&i| s.workers[i].threads() == *threads).take(want).collect();
        while have.len() < want {
            match s.fork(*threads) {
                Ok(i) => have.push(i),
                Err(f) => {
                    // no calculator on these threads: their requests fail with the reason
                    if have.is_empty() {
                        for i in (0..n).filter(|&i| jobs[i].threads == *threads) {
                            progress.out[i] = Some(Err(Failure::of(f.kind, f.code, f.message.clone())));
                            progress.state[i] = Run::Done;
                        }
                    }
                    break;
                }
            }
        }
        chosen.extend(have);
    }
    let sync = (Mutex::new(progress), Condvar::new());
    let (history, given) = (&s.history, s.given.as_ref());
    let mut used: Vec<&mut Worker> =
        s.workers.iter_mut().enumerate().filter(|(i, _)| chosen.contains(i)).map(|(_, w)| w).collect();
    if let Some((first, rest)) = used.split_first_mut() {
        std::thread::scope(|scope| {
            for w in rest {
                let (jobs, lineages, sync) = (&jobs, &lineages, &sync);
                // a thread that cannot start leaves its share to the calculators on the same threads
                let _ = std::thread::Builder::new()
                    .stack_size(WORKER_STACK)
                    .spawn_scoped(scope, move || work(w, jobs, lineages, history, given, sync));
            }
            work(first, &jobs, &lineages, history, given, &sync);
        });
    }
    let p = sync.0.into_inner().unwrap_or_else(|e| e.into_inner());
    s.next = p.next;
    if p.panicked {
        s.close();
    } else {
        for (job, c) in jobs.iter().zip(p.commits) {
            if let Some(c) = c {
                s.history.insert(job.lineage.clone(), c);
            }
        }
    }
    p.out.into_iter().map(|o| o.expect("every request ran")).collect()
}

/// A session's calculator (private to `axtb`): the driver calculator of one request, built
/// once, forks of it on the threads requests ask for and to run lineages side by side, and the
/// committed wavefunction of each lineage.
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
                    s.close();
                    Err(Failure::panic(p))
                }
            }
        })
        .map_err(|f| f.into_py(py))
    }

    /// The response dict of a request: the response with `history` (the key whose history it
    /// started from, or None), `given` and `threads`.
    fn response<'py>(
        &self,
        py: Python<'py>,
        r: Response,
        from: Option<String>,
        given: bool,
        threads: usize,
    ) -> PyResult<Bound<'py, PyDict>> {
        let d = response_dict(py, r)?;
        d.set_item("history", from)?;
        d.set_item("given", given)?;
        d.set_item("threads", threads)?;
        Ok(d)
    }
}

fn closed() -> Failure {
    Failure::of("closed", "session_closed", "the session is closed")
}

/// A request of `Engine.compute_many`: positions, task, lineage, seeds, commit, accuracy,
/// orbitals, threads (as `Engine.compute` takes them).
type ManyRequest<'py> =
    (PyArrayLike2<'py, f64, AllowTypeChange>, String, String, Vec<String>, bool, Option<f64>, bool, Option<usize>);

/// A request of `Engine.compute_many` with checked positions, by its index.
type Checked = (usize, Vec<[f64; 3]>, String, String, Vec<String>, bool, Option<f64>, bool, usize);

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
                workers: vec![Worker { calc, holds: Holds::Nothing, accuracy: req.accuracy.unwrap_or(1.0) }],
                history: HashMap::new(),
                next: 0,
                tasks: HashMap::new(),
                given: None,
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

    /// The threads of the engine's own calculator (a call's default).
    #[getter]
    fn threads(&self) -> usize {
        self.threads
    }

    /// The threads of each calculator the engine holds: its own, then the forks made for calls
    /// on other threads and for `compute_many` (empty once closed).
    #[getter]
    fn calculators(&self, py: Python<'_>) -> PyResult<Vec<usize>> {
        self.with_state(py, |s| Ok(s.workers.iter().map(Worker::threads).collect()))
    }

    /// Run `task` (JSON of the driver request's task fields) at `positions` (bohr, N x 3) in the
    /// lineage `lineage` (its key) on `threads` threads (None: the engine's), on a calculator
    /// on those threads (a fork made when first needed). The SCF starts from the lineage's
    /// committed history, else from the first of `seeds` (keys) with one, else fresh, as the
    /// restart policy says. `accuracy`: the SCF accuracy of this call (set only when it differs
    /// from the calculator's last); `orbitals`: the structured orbitals in the response. On
    /// success with `commit`, the method's wavefunction after the call becomes the lineage's
    /// history. Returns the response dict with `history`: the key whose history was installed,
    /// or None, and `threads`.
    #[pyo3(signature = (positions, task, lineage, seeds=Vec::new(), commit=true, accuracy=None, orbitals=false, threads=None))]
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
        threads: Option<usize>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let positions = rows3(&positions, self.natoms, "positions").map_err(|f| f.into_py(py))?;
        let threads = threads.unwrap_or(self.threads);
        let (r, from, given) = self.with_state(py, move |s| {
            if s.workers.is_empty() {
                return Err(closed());
            }
            let task = s.task(task, orbitals)?;
            let job = Job { positions, task, lineage, seeds, commit, accuracy, threads };
            // the history to start from: the lineage's own, else the first seed's
            let from = std::iter::once(&job.lineage)
                .chain(job.seeds.iter())
                .find(|k| s.history.contains_key(k.as_str()))
                .cloned();
            let start = from.as_ref().map(|k| s.history[k].clone());
            let given = start.is_none() && s.given.is_some();
            let i = s.worker_for(threads, start.as_ref().map_or(Holds::Nothing, |c| Holds::Commit(c.0)))?;
            let (r, wf) = run(&mut s.workers[i], &job, start.as_ref(), s.given.as_ref())?;
            if let Some(wf) = wf {
                let id = s.next;
                s.next += 1;
                s.history.insert(job.lineage, (id, Arc::new(wf)));
                s.workers[i].holds = Holds::Commit(id);
            }
            Ok((r, from, given))
        })?;
        self.response(py, r, from, given, threads)
    }

    /// Run `requests` (each `compute`'s arguments in order: positions, task, lineage, seeds,
    /// commit, accuracy, orbitals, threads) as `compute` runs them one at a time in their order,
    /// the same bits, with different lineages side by side as long as their threads fit into
    /// `threads` in all (taken greedily in request order), on the engine's calculator and forks
    /// of it on each request's threads (kept for later calls). Each lineage's requests run in
    /// order; a request of a lineage without history waits for its seeds' earlier requests.
    /// Returns, in order, each request's response dict (as `compute`'s, with `wall`: its time on
    /// its calculator, s) or its `Error`.
    #[pyo3(signature = (requests, threads))]
    fn compute_many<'py>(
        &self,
        py: Python<'py>,
        requests: Vec<ManyRequest<'py>>,
        threads: usize,
    ) -> PyResult<Bound<'py, PyList>> {
        let n = requests.len();
        let mut out: Vec<Option<(Result<Done, Failure>, usize)>> = (0..n).map(|_| None).collect();
        let mut checked: Vec<Checked> = Vec::new();
        for (i, (positions, task, lineage, seeds, commit, accuracy, orbitals, t)) in requests.into_iter().enumerate() {
            let t = t.unwrap_or(self.threads);
            match rows3(&positions, self.natoms, "positions") {
                Ok(x) => checked.push((i, x, task, lineage, seeds, commit, accuracy, orbitals, t)),
                Err(f) => out[i] = Some((Err(f), t)),
            }
        }
        let ran = self.with_state(py, move |s| {
            if s.workers.is_empty() {
                return Ok(checked.into_iter().map(|c| (c.0, Err(closed()), c.8)).collect::<Vec<_>>());
            }
            let (mut index, mut jobs, mut ran) = (Vec::new(), Vec::new(), Vec::new());
            for (i, positions, task, lineage, seeds, commit, accuracy, orbitals, t) in checked {
                match s.task(task, orbitals) {
                    Ok(task) => {
                        index.push((i, t));
                        jobs.push(Job { positions, task, lineage, seeds, commit, accuracy, threads: t });
                    }
                    Err(f) => ran.push((i, Err(f), t)),
                }
            }
            if !jobs.is_empty() {
                ran.extend(index.into_iter().zip(batch(s, jobs, threads.max(1))).map(|((i, t), r)| (i, r, t)));
            }
            Ok(ran)
        })?;
        for (i, r, t) in ran {
            out[i] = Some((r, t));
        }
        let list = PyList::empty(py);
        for r in out {
            match r.expect("every request has an outcome") {
                (Ok(d), t) => {
                    let x = self.response(py, d.response, d.from, d.given, t)?;
                    x.set_item("wall", d.seconds)?;
                    list.append(x)?;
                }
                (Err(f), _) => list.append(f.into_py(py).into_value(py))?,
            }
        }
        Ok(list)
    }

    /// Start every request without history from a given start: densities (g-xTB, one n x n per
    /// spin channel, a closed shell's the total density) or GFN2-xTB's moments (shell charges,
    /// atomic dipoles N x 3, packed quadrupoles N x 6), derived at the engine's current positions
    /// (the opening geometry) and checked by the method.
    #[pyo3(signature = (densities=None, moments=None))]
    #[allow(clippy::type_complexity)]
    fn set_given<'py>(
        &self,
        py: Python<'py>,
        densities: Option<Vec<PyArrayLike2<'py, f64, AllowTypeChange>>>,
        moments: Option<(Vec<f64>, PyArrayLike2<'py, f64, AllowTypeChange>, PyArrayLike2<'py, f64, AllowTypeChange>)>,
    ) -> PyResult<()> {
        let d = densities.as_deref().map(self::densities);
        let m = moments.as_ref().map(|m| self::moments(m, self.natoms)).transpose().map_err(|f| f.into_py(py))?;
        self.with_state(py, move |s| {
            let w = s.workers.first_mut().ok_or_else(closed)?;
            install_given(&mut w.calc, d, m)?;
            s.given = w.calc.wavefunction();
            for w in &mut s.workers {
                w.calc.set_history(None)?;
                w.holds = Holds::Nothing;
            }
            Ok(())
        })
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

    /// Drop the calculators and every lineage's history; later calls raise kind `closed`.
    fn close(&self, py: Python<'_>) -> PyResult<()> {
        self.with_state(py, |s| {
            s.close();
            Ok(())
        })
    }

    /// Whether the engine is closed (by `close` or a panic).
    #[getter]
    fn closed(&self, py: Python<'_>) -> PyResult<bool> {
        self.with_state(py, |s| Ok(s.workers.is_empty()))
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
/// `normal_modes` modes x N x 3; the traceless molecular quadrupole and the atomic quadrupoles
/// also unpacked: `quadrupole_tensor` 3 x 3, `atomic_quadrupole_tensors` N x 3 x 3), the
/// binding-only fields included (`subtotals`;
/// `core_electrons`, `shell_populations`, `shell_atoms`, `shell_angular_momenta`,
/// `reference_shell_occupations`, `quadrupole_traceless`; the thermochemistry's
/// `inverted_frequencies`, `dropped_frequencies`, `average_moment`), and the warnings the
/// Python layer adds (`eigensolver_warning`, `retry_warnings`).
fn response_dict(py: Python<'_>, r: Response) -> PyResult<Bound<'_, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("method", r.method)?;
    d.set_item("conventions", r.conventions)?;
    d.set_item("numerics", &r.numerics)?;
    d.set_item("eigensolver", &r.eigensolver)?;
    d.set_item("eigensolver_warning", axtb_driver::eigensolver_warning(&r))?;
    d.set_item("retry_warnings", axtb_driver::retry_warnings(&r))?;
    match &r.scf_retry {
        Some(x) => d.set_item("scf_retry", json(py, x)?)?,
        None => d.set_item("scf_retry", py.None())?,
    }
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
        x.set_item("quadrupole_tensor", rows(py, &tensor(&s.quadrupole_traceless))?)?;
        x.set_item("shell_charges", array1(py, &s.shell_charges))?;
        x.set_item("atomic_dipoles", rows(py, &s.atomic_dipoles)?)?;
        x.set_item("atomic_quadrupoles", rows(py, &s.atomic_quadrupoles)?)?;
        x.set_item("atomic_quadrupole_tensors", tensors(py, &s.atomic_quadrupoles)?)?;
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

/// The shells of `method`'s basis for `numbers` from its tables: a list of (atom, l, reference
/// occupation, exponents (empty for GFN2-xTB)); `axtb_driver::shell_layout`.
#[pyfunction]
#[allow(clippy::type_complexity)]
pub(crate) fn shell_layout(
    py: Python<'_>,
    method: &str,
    numbers: Vec<usize>,
) -> PyResult<Option<Vec<(usize, usize, f64, Vec<f64>)>>> {
    let m = axtb_driver::Method::parse(method).map_err(|e| Failure::from(e).into_py(py))?;
    Ok(axtb_driver::shell_layout(m, &numbers)
        .map(|v| v.into_iter().map(|s| (s.atom, s.l, s.reference_occupation, s.exponents)).collect()))
}

/// The largest molecule (atoms) whose single point (`hessian` false) or Hessian of `method`
/// runs on one thread (`axtb_driver::single_thread_atoms`).
#[pyfunction]
pub(crate) fn single_thread_atoms(py: Python<'_>, method: &str, hessian: bool) -> PyResult<usize> {
    let m = axtb_driver::Method::parse(method).map_err(|e| Failure::from(e).into_py(py))?;
    Ok(axtb_driver::single_thread_atoms(m, hessian))
}

/// The threads a plan states for one energy or gradient (`hessian` false) or Hessian of `method`
/// for `atoms` atoms (`axtb_driver::planned_threads`).
#[pyfunction]
pub(crate) fn planned_threads(py: Python<'_>, method: &str, hessian: bool, atoms: usize) -> PyResult<usize> {
    let m = axtb_driver::Method::parse(method).map_err(|e| Failure::from(e).into_py(py))?;
    Ok(axtb_driver::planned_threads(m, hessian, atoms))
}

/// The plan's thread bands (`axtb_driver::THREAD_BANDS`): (method, hessian, [(from atoms,
/// threads)]), for the generated docs.
#[pyfunction]
#[allow(clippy::type_complexity)]
pub(crate) fn thread_bands() -> Vec<(&'static str, bool, Vec<(usize, usize)>)> {
    axtb_driver::THREAD_BANDS.iter().map(|(m, h, b)| (m.name(), *h, b.to_vec())).collect()
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
