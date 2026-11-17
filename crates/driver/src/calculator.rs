// SPDX-License-Identifier: LGPL-3.0-or-later
//! A method bound to one molecule: the request's method built behind axtb-base's [`Method`]
//! trait (GFN2-xTB, g-xTB and its solvation composite, GFN-FF), and the response of any of them.

use crate::error::{Error, Result};
use crate::names::{self, Method};
use crate::request::Request;
use crate::response::{
    self, Attempt, Composite, Details, Gfnff, HessianPart, Matrix, Named, Part, Parts, Response, Retried, Scf,
    Topology, Vibrations,
};
use crate::task::Task;
use axtb::method::{
    Fingerprint, HessianOptions, Method as MethodTrait, Restart, Retry, RetryLadder, RetryOn, Rung, ScfData,
    ScfSettings, SinglePoint, TermKind, Try, Want, WavefunctionData,
};
use axtb::ToMolden;
use std::time::Instant;

pub use axtb::method::Wavefunction;

/// A method bound to one molecule (atomic numbers, charge, unpaired electrons) with its
/// settings, built from a [`Request`]: [`Calculator::compute`] runs a [`Task`] at any
/// geometry and starts each SCF from the previous geometry's wavefunction (as the method's
/// conventions or the request's `guess` say).
///
/// With `scf_retries`, a calculation that fails (an SCF that does not converge or breaks down)
/// runs again at each rung in order, each from the method's own start, until one succeeds: the
/// whole calculation, so that a Hessian (its reference and every displaced or perturbed SCF) and
/// the single point of its response come from one rung. The g-xTB solvation composite climbs
/// them per part instead ([`axtb_gxtb::CompositeOptions::scf_retries`]). Afterwards the method has
/// its own settings again, and the wavefunction of the rung that converged is the next SCF's
/// history (a session still tries its own settings first).
pub struct Calculator {
    method: Box<dyn MethodTrait>,
    kind: Method,
    /// The numerics identity of its responses ([`crate::numerics`]).
    numerics: String,
    numbers: Vec<usize>,
    charge: f64,
    uhf: Option<usize>,
    threads: usize,
    /// The `scf_retries` the calculator climbs itself (empty: none, or the composite's own).
    retries: Vec<Retry>,
    /// How many `scf_retries` the request gave (the composite's included).
    given_retries: usize,
}

/// A calculation's Hessian and single point, and with `scf_retries` the rung that gave them
/// and every try.
type Climbed = (Option<axtb::Hessian>, SinglePoint, Option<(Retried, Vec<Attempt>)>);

fn input(msg: impl Into<String>) -> Error {
    Error::Input(msg.into())
}

/// Refuse the settings in `bad` for `method`.
fn refuse(method: Method, bad: &[&str]) -> Result<()> {
    if bad.is_empty() {
        Ok(())
    } else {
        Err(input(format!("{method} does not support {}", bad.join(", "))))
    }
}

/// The implicit-solvent setting of a request, if any.
fn solvation(req: &Request) -> Result<Option<axtb::SolvationSpec>> {
    let state = req.solvation_state.as_deref().map(names::solvation_state).transpose()?;
    let (model, solvent) = match (&req.alpb, &req.gbsa) {
        (Some(_), Some(_)) => return Err(input("give either alpb or gbsa, not both")),
        (Some(s), None) => (axtb::SolvationModel::Alpb, s),
        (None, Some(s)) => (axtb::SolvationModel::Gbsa, s),
        (None, None) => {
            if state.is_some() {
                return Err(input("solvation_state is a setting of the implicit solvent (alpb or gbsa)"));
            }
            return Ok(None);
        }
    };
    Ok(Some(axtb::SolvationSpec { model, solvent: solvent.clone(), state: state.unwrap_or_default() }))
}

/// The hint for a g-xTB method given an implicit solvent alone: g-xTB has no solvent
/// parameters; the solvation composite is the way.
pub fn gxtb_solvent_hint(method: Method) -> String {
    format!(
        "{method} has no implicit-solvent parameters (the g-xTB program stops with --alpb/--gbsa). For g-xTB in \
         solvent use the solvation composite E = E(g-xTB) + [E(GFN2-xTB + ALPB/GBSA) - E(GFN2-xTB)]: add \
         solvation_correction \"gfn2\" (CLI: --solvation-correction gfn2) to the solvent (alpb or gbsa)"
    )
}

impl Calculator {
    /// Set up the method of `req` for its molecule at `req.positions` (bohr): its settings are
    /// checked (unknown names, settings of other methods and impossible values are
    /// [`Error::Input`]); the fields of what to compute are not read (see
    /// [`Task::from_request`]).
    pub fn new(req: &Request) -> Result<Self> {
        let kind = Method::parse(req.method.as_deref().unwrap_or("gfn2"))?;
        let nat = req.numbers.len();
        if req.positions.len() != nat {
            return Err(input(format!("positions must have one row per atom ({nat}), not {}", req.positions.len())));
        }
        if req.uhf.is_some_and(|u| u < 0) {
            return Err(input(format!("uhf (unpaired electrons) must not be negative, not {}", req.uhf.unwrap())));
        }
        if req.threads == Some(0) {
            return Err(input("threads must be at least 1"));
        }
        if req.max_iter == Some(0) {
            return Err(input("max_iter must be at least 1"));
        }
        let retries: Vec<Retry> = req.scf_retries.iter().flatten().map(|&r| r.into()).collect();
        if req.scf_retries.is_some() && kind != Method::Gfnff {
            Retry::check(&retries, req.etemp)?;
        }
        if kind != Method::Gfnff && (req.fragment_charges.is_some() || req.topology_positions.is_some()) {
            return Err(input("fragment_charges and topology_positions are GFN-FF settings"));
        }
        if !kind.is_gxtb() && (req.solvation_correction.is_some() || req.correction_conventions.is_some()) {
            return Err(input(
                "solvation_correction and correction_conventions are settings of g-xTB in solvent; GFN2-xTB has its \
                 own implicit solvent (alpb or gbsa alone)",
            ));
        }
        let threads = req.threads.unwrap_or(1);
        let uhf = req.uhf.map(|u| u as usize);
        let guess = req.guess.as_deref().map(names::guess).transpose()?;
        // `restart` sets the policy; `guess` then names a start alone
        let restart = req.restart.as_deref().map(names::restart).transpose()?;
        if restart.is_some() && matches!(guess, Some(names::Guess::Previous | names::Guess::Fresh)) {
            return Err(input(
                "with restart, guess names where an SCF starts without a previous state (\"eeq\" or \"sad\" for \
                 GFN2-xTB, \"eeq-bc\" for g-xTB), not \"previous\" or \"fresh\"",
            ));
        }
        let method: Box<dyn MethodTrait> = match kind {
            Method::Gfn2 => {
                let conventions =
                    req.conventions.as_deref().map(names::gfn2_conventions).transpose()?.unwrap_or_default();
                let mut opts = axtb::Options::default();
                if let Some(a) = req.accuracy {
                    opts.scf.accuracy = a;
                }
                opts.scf.etemp = req.etemp;
                opts.scf.max_iter = req.max_iter;
                opts.mixing.mixer = req.mixer.as_deref().map(names::mixer).transpose()?;
                match guess {
                    Some(names::Guess::Sad) => opts.initial_guess = axtb::InitialGuess::Sad,
                    Some(names::Guess::EeqBc) => {
                        return Err(input(
                            "GFN2-xTB starts from EEQ charges or neutral atoms: guess \"eeq\" or \"sad\"",
                        ))
                    }
                    _ => {}
                }
                // a named start at every geometry is the method's own start at every geometry,
                // unless `restart` says otherwise
                opts.scf.restart = restart.or(guess.map(|g| match g {
                    names::Guess::Previous => Restart::Previous,
                    _ => Restart::Fresh,
                }));
                opts.threads = threads;
                opts.trace = req.scf_trace;
                opts.dispersion = req.dispersion.unwrap_or(true);
                opts.conventions = conventions;
                opts.solvation = solvation(req)?;
                let calc = axtb::Calculator::new(
                    axtb::Hamiltonian::Gfn2,
                    &req.numbers,
                    &req.positions,
                    req.charge,
                    uhf,
                    opts,
                )?;
                Box::new(calc)
            }
            Method::Gfnff => {
                let mut bad = Vec::new();
                if req.alpb.is_some() || req.gbsa.is_some() || req.solvation_state.is_some() {
                    bad.push("implicit solvation");
                }
                if req.dispersion == Some(false) {
                    bad.push("switching off dispersion");
                }
                if req.accuracy.is_some()
                    || req.etemp.is_some()
                    || req.max_iter.is_some()
                    || req.scf_retries.is_some()
                    || req.mixer.is_some()
                    || req.guess.is_some()
                    || req.restart.is_some()
                    || req.scf_trace
                {
                    bad.push("SCF settings (accuracy, etemp, max_iter, scf_retries, mixer, guess, restart, scf_trace)");
                }
                refuse(kind, &bad)?;
                let conventions = names::gfnff_conventions(req.conventions.as_deref().unwrap_or("axtb"))?;
                let setup = match &req.topology_positions {
                    Some(t) if t.len() != nat => {
                        return Err(input(format!("topology_positions must have one row per atom ({nat})")))
                    }
                    Some(t) => t.as_slice(),
                    None => &req.positions,
                };
                let opts = axtb_gfnff::Options {
                    threads,
                    fragment_charges: req.fragment_charges.clone(),
                    conventions,
                    ..Default::default()
                };
                Box::new(axtb_gfnff::Calculator::new(&req.numbers, setup, req.charge, opts)?)
            }
            Method::GxtbV1_1 | Method::GxtbV2_0_1 => gxtb_method(req, kind, threads, uhf, guess, restart)?,
        };
        let composite = kind.is_gxtb() && req.solvation_correction.is_some();
        let numerics = crate::numerics(kind, composite);
        let given_retries = retries.len();
        let retries = if composite { Vec::new() } else { retries };
        Ok(Self {
            method,
            kind,
            numerics,
            numbers: req.numbers.clone(),
            charge: req.charge,
            uhf,
            threads,
            retries,
            given_retries,
        })
    }

    /// The method.
    pub fn method(&self) -> Method {
        self.kind
    }

    /// The numerics identity of its results ([`crate::numerics`]).
    pub fn numerics(&self) -> &str {
        &self.numerics
    }

    /// The canonical name of the method's conventions.
    pub fn conventions(&self) -> &'static str {
        self.method.info().conventions
    }

    /// Atomic numbers of the molecule.
    pub fn numbers(&self) -> &[usize] {
        &self.numbers
    }

    /// Threads of the calculator.
    pub fn threads(&self) -> usize {
        self.threads
    }

    /// A copy of the calculator on a pool of its own of `threads` threads, to compute beside
    /// this one (`axtb.Session.compute_many`'s lineages side by side): the same settings (the
    /// SCF accuracy as it is now) and the method's setup shared where it is immutable
    /// (`Method::fork`), with its own SCF state and work buffers and no wavefunction. On as many
    /// threads as this calculator and from the same start, it computes what this one computes,
    /// bit for bit.
    pub fn fork(&self, threads: usize) -> Result<Self> {
        let pool = axtb::Threads::new(threads)?;
        Ok(Self {
            method: self.method.fork(&pool)?,
            kind: self.kind,
            numerics: self.numerics.clone(),
            numbers: self.numbers.clone(),
            charge: self.charge,
            uhf: self.uhf,
            threads,
            retries: self.retries.clone(),
            given_retries: self.given_retries,
        })
    }

    /// Run `task` at `positions` (bohr). With a Hessian, the response is that of the single
    /// point at `positions` (energy, gradient, charges, orbitals; g-xTB's SCF from its own
    /// start, as without the Hessian) plus the Hessian, its dipole derivatives and its
    /// vibrational analysis.
    pub fn compute(&mut self, positions: &[[f64; 3]], task: &Task) -> Result<Response> {
        let t0 = Instant::now();
        let nat = self.numbers.len();
        if positions.len() != nat {
            return Err(input(format!("positions must have one row per atom ({nat}), not {}", positions.len())));
        }
        let want = Want {
            gradient: task.gradient || task.hessian.is_some(),
            orbitals: task.molden.is_some() || task.orbitals,
            bond_orders: task.bond_orders,
            overlap: task.overlap,
            hamiltonian: task.hamiltonian,
        };
        let hessian =
            task.hessian.as_ref().map(|h| HessianOptions { method: h.method, step: h.step, accuracy: h.accuracy });
        task.check(self.method.info())?;
        self.method.check(&want, hessian.as_ref())?;
        let (hessian, sp, retried) = if !self.retries.is_empty() {
            self.climb(positions, &want, hessian.as_ref())?
        } else if self.given_retries > 0 {
            self.composite_tries(positions, &want, hessian.as_ref())?
        } else {
            let hessian = hessian.map(|o| self.method.hessian(positions, &o)).transpose()?;
            (hessian, self.method.singlepoint(positions, &want)?, None)
        };
        let info = self.method.info();
        let mut out = Response {
            program: response::PROGRAM,
            version: response::VERSION,
            method: Some(info.name),
            conventions: info.conventions,
            numerics: self.numerics.clone(),
            ..Default::default()
        };
        let orbitals = fill(&mut out, sp, hessian, self.given_retries > 0);
        if let Some((retried, attempts)) = retried {
            out.scf_retry = Some(retried);
            if let Some(s) = out.scf.as_mut() {
                s.scf_attempts = Some(attempts);
            }
        }
        check_finite(&out)?;
        if let Some(m) = &task.molden {
            let orbitals =
                orbitals.as_ref().ok_or_else(|| Error::Numerical("the calculation returned no orbitals".into()))?;
            let title = format!("axtb {} {}", response::VERSION, self.kind);
            let text = orbitals.to_molden(m.functions, &title)?;
            if let Some(path) = &m.file {
                std::fs::write(path, &text).map_err(|e| input(format!("cannot write the Molden file {path}: {e}")))?;
            }
            out.molden_file = m.file.clone();
            out.molden = m.text.then_some(text);
        }
        if task.orbitals {
            out.orbitals = orbitals;
        }
        if let (Some(h), Some(hp)) = (&task.hessian, &out.hessian) {
            if let Some(v) = &h.vibrations {
                let opts = v.options(&info.analysis);
                let a = axtb::vibrations::analyze(
                    &self.numbers,
                    positions,
                    &hp.hessian.data,
                    hp.dipole_derivatives.as_deref(),
                    &opts,
                )?;
                out.vibrations = Some(Vibrations::new(a, out.energy));
            }
        }
        out.seconds = Some(t0.elapsed().as_secs_f64());
        Ok(out)
    }

    /// The electronic temperature (K) and iteration limit of `s` for this method, its defaults
    /// where `s` has none.
    fn effective(&self, s: &ScfSettings) -> (f64, usize) {
        let (etemp, max_iter) = match self.kind {
            Method::Gfn2 => (axtb::ETEMP, axtb::MAX_ITER),
            Method::GxtbV1_1 => (axtb_gxtb::composite::GXTB_ETEMP, axtb_gxtb::options::MAX_ITER_V1_1),
            Method::GxtbV2_0_1 => (axtb_gxtb::composite::GXTB_ETEMP, axtb_gxtb::scf::MAX_ITER),
            Method::Gfnff => (0.0, 0),
        };
        (s.etemp.unwrap_or(etemp), s.max_iter.unwrap_or(max_iter))
    }

    /// The calculation of `compute` on the `scf_retries` ladder: the Hessian (if any) and the
    /// single point at the method's own settings, then at each rung from the method's own start
    /// until both succeed; returns them with the rung and every try.
    fn climb(&mut self, positions: &[[f64; 3]], want: &Want, hessian: Option<&HessianOptions>) -> Result<Climbed> {
        let base = self.method.scf_settings().cloned().expect("an SCF method has SCF settings");
        let first = self.effective(&base);
        let mut rungs = vec![Rung::default()];
        rungs.extend(Retry::rungs(first, &self.retries));
        let ladder = RetryLadder { rungs, on: RetryOn::Failed };
        let mut attempts = Vec::new();
        let m: &mut dyn MethodTrait = self.method.as_mut();
        let r = ladder.climb_methods(&mut [m], 0, |_, m, i, rung| {
            let (etemp, max_iter) = (rung.etemp.unwrap_or(first.0), rung.max_iter.unwrap_or(first.1));
            let tried = Attempt { retry: Some(i), etemp: Some(etemp), max_iter: Some(max_iter), ..Attempt::default() };
            let r = hessian
                .map(|o| m.hessian(positions, o))
                .transpose()
                .and_then(|h| Ok((h, m.singlepoint(positions, want)?)));
            match &r {
                Ok((_, sp)) => {
                    let scf = sp.scf().expect("an SCF method's single point has its SCF");
                    let one = Attempt { start: Some(scf.start), ..tried };
                    if scf.attempts.is_empty() {
                        attempts.push(Attempt {
                            mixer: scf.mixer_used,
                            iterations: Some(scf.iterations),
                            converged: true,
                            ..one
                        });
                    } else {
                        attempts.extend(scf.attempts.iter().map(|a| Attempt {
                            retry: one.retry,
                            etemp: one.etemp,
                            max_iter: one.max_iter,
                            start: one.start,
                            ..Attempt::mixer(a)
                        }));
                    }
                }
                Err(e) => attempts.push(failed_try(tried, e)),
            }
            r
        });
        match r {
            Ok((mut v, i)) => {
                let (h, sp) = v.pop().expect("one method's result");
                let label = ladder.rungs[i].label.clone();
                Ok((h, sp, Some((Retried { rung: i, of: self.given_retries, label }, attempts))))
            }
            Err(e) => Err(with_attempts(e.into(), attempts)),
        }
    }

    /// The calculation of `compute` with the composite's own `scf_retries` ladders: every try of
    /// its parts (its Hessian's, then its single point's), and the highest rung of a part.
    fn composite_tries(
        &mut self,
        positions: &[[f64; 3]],
        want: &Want,
        hessian: Option<&HessianOptions>,
    ) -> Result<Climbed> {
        let mut attempts = Vec::new();
        let h = hessian.map(|o| self.method.hessian(positions, o)).transpose();
        if hessian.is_some() {
            attempts.extend(self.method.tries().iter().map(|t| Attempt { hessian: Some(true), ..try_attempt(t) }));
        }
        let h = h.map_err(|e| with_attempts(e.into(), attempts.clone()))?;
        let sp = self.method.singlepoint(positions, want);
        attempts.extend(self.method.tries().iter().map(try_attempt));
        let sp = sp.map_err(|e| with_attempts(e.into(), attempts.clone()))?;
        let rung = attempts.iter().filter(|a| a.converged).filter_map(|a| a.retry).max().unwrap_or(0);
        Ok((h, sp, Some((Retried { rung, of: self.given_retries, label: None }, attempts))))
    }

    /// The converged wavefunction of the last SCF (`None` before the first, after an error and
    /// for GFN-FF), for another calculator of the same method and molecule
    /// ([`Calculator::set_wavefunction`]).
    pub fn wavefunction(&self) -> Option<Wavefunction> {
        self.method.wavefunction()
    }

    /// Start the next SCF from `wf` (from [`Calculator::wavefunction`] of a calculator for the
    /// same method and molecule), whatever the `guess` setting says.
    pub fn set_wavefunction(&mut self, wf: &Wavefunction) -> Result<()> {
        Ok(self.method.set_wavefunction(wf)?)
    }

    /// Start the next GFN2-xTB SCF from given moments: shell charges, atomic dipoles and
    /// quadrupoles (e.g. those of a previous response, or another program's).
    pub fn set_moments(
        &mut self,
        shell_charges: Vec<f64>,
        dipoles: Vec<[f64; 3]>,
        quadrupoles: Vec<[f64; 6]>,
    ) -> Result<()> {
        let data =
            ScfData { shell_charges, atomic_dipoles: dipoles, atomic_quadrupoles: quadrupoles, ..ScfData::default() };
        let fingerprint = Fingerprint::new(self.method.info().name, &self.numbers, self.charge, self.uhf);
        Ok(self.method.set_wavefunction(&Wavefunction::new(fingerprint, WavefunctionData::Scf(data)))?)
    }

    /// Start the next SCF from given density matrices (g-xTB and the composite's g-xTB part:
    /// row-major n x n in axtb's AO order, one per spin channel, a closed shell's the total
    /// density), whatever the restart policy says: e.g. C n C^T of a given wavefunction's
    /// orbitals. Checked against the molecule (channels, size, electrons to 1e-6); GFN2-xTB starts
    /// from moments ([`Calculator::set_moments`]), GFN-FF from nothing.
    pub fn set_density(&mut self, densities: Vec<Vec<f64>>) -> Result<()> {
        Ok(self.method.set_density(densities)?)
    }

    /// Install `wf` as the method's own previous wavefunction (`None`: none): unlike
    /// [`Calculator::set_wavefunction`], the restart policy decides whether the next SCF starts
    /// from it. A session that keeps one wavefunction per lineage installs it before each call of
    /// that lineage and computes what one calculator per lineage computes.
    pub fn set_history(&mut self, wf: Option<&Wavefunction>) -> Result<()> {
        Ok(self.method.set_history(wf)?)
    }

    /// The SCF accuracy of later calls: the SCF thresholds (GFN2-xTB keeps the integral screening
    /// it was built with); every part's for the composite. GFN-FF has none.
    pub fn set_scf_accuracy(&mut self, accuracy: f64) -> Result<()> {
        Ok(self.method.set_scf_accuracy(accuracy)?)
    }

    /// Forget the previous wavefunction(s): the next SCF starts from the method's own guess.
    /// `gxtb_only`: of g-xTB only (the composite's g-xTB part; for every other method nothing).
    pub fn reset_wavefunction(&mut self, gxtb_only: bool) {
        if !gxtb_only {
            return self.method.reset_wavefunction();
        }
        if !self.kind.is_gxtb() {
            return;
        }
        // a composite: every part's but the base's back (WavefunctionData::Parts, base first)
        let kept = self.method.wavefunction();
        self.method.reset_wavefunction();
        if let Some(wf) = kept {
            if let WavefunctionData::Parts(parts) = wf.data() {
                let mut parts = parts.clone();
                parts[0] = None;
                let rest = Wavefunction::new(wf.fingerprint().clone(), WavefunctionData::Parts(parts));
                self.method.set_wavefunction(&rest).expect("a composite takes its own parts back");
            }
        }
    }
}

/// The g-xTB method of a request: plain, or the solvation composite with
/// `solvation_correction`.
fn gxtb_method(
    req: &Request,
    kind: Method,
    threads: usize,
    uhf: Option<usize>,
    guess: Option<names::Guess>,
    restart: Option<Restart>,
) -> Result<Box<dyn MethodTrait>> {
    let version = kind.gxtb_version().expect("a g-xTB method");
    let conventions = req.conventions.as_deref().map(names::gxtb_conventions).transpose()?.unwrap_or_default();
    if req.scf_trace {
        return Err(input(format!("{kind} does not support scf_trace (GFN2-xTB only)")));
    }
    if let Some(c) = &req.solvation_correction {
        if req.dispersion == Some(false) {
            return Err(input("the g-xTB solvation composite does not support switching off dispersion"));
        }
        let correction = names::correction(c)?;
        let solvation = solvation(req)?
            .ok_or_else(|| input("solvation_correction needs a solvent: alpb or gbsa (e.g. alpb \"water\")"))?;
        let mut o = axtb_gxtb::CompositeOptions {
            version,
            correction,
            solvation,
            scf: ScfSettings {
                accuracy: req.accuracy.unwrap_or(1.0),
                etemp: req.etemp,
                max_iter: req.max_iter,
                restart: match guess {
                    None => restart,
                    Some(names::Guess::Previous) => Some(Restart::Previous),
                    Some(names::Guess::Fresh) => Some(Restart::Fresh),
                    Some(_) => {
                        return Err(input(
                            "the g-xTB solvation composite starts each part from its method's own start: guess \
                         \"previous\" or \"fresh\"",
                        ))
                    }
                },
            },
            threads,
            conventions,
            ..Default::default()
        };
        o.correction_conventions = req.correction_conventions.as_deref().map(names::gfn2_conventions).transpose()?;
        o.mixer = req.mixer.as_deref().map(names::mixer).transpose()?;
        o.scf_retries = req.scf_retries.as_ref().map(|r| r.iter().map(|&x| x.into()).collect());
        return Ok(Box::new(axtb_gxtb::CompositeCalculator::new(&req.numbers, &req.positions, req.charge, uhf, o)?));
    }
    if req.alpb.is_some() || req.gbsa.is_some() || req.solvation_state.is_some() || req.correction_conventions.is_some()
    {
        return Err(input(gxtb_solvent_hint(kind)));
    }
    if req.dispersion == Some(false) {
        return Err(input(format!("{kind} does not support switching off dispersion")));
    }
    if req.mixer.is_some() {
        return Err(input(format!("{kind} does not support choosing the SCF mixer")));
    }
    let mut o = axtb_gxtb::Options { threads, conventions, ..Default::default() };
    if let Some(a) = req.accuracy {
        o.scf.accuracy = a;
    }
    o.scf.etemp = req.etemp;
    o.scf.max_iter = req.max_iter;
    o.scf.restart = match guess {
        None => restart,
        Some(names::Guess::EeqBc) if restart.is_some() => restart,
        Some(names::Guess::Previous) => Some(Restart::Previous),
        Some(names::Guess::Fresh | names::Guess::EeqBc) => Some(Restart::Fresh),
        Some(names::Guess::Eeq | names::Guess::Sad) => {
            return Err(input(format!("{kind} starts from the program's EEQ-BC charges: guess \"eeq-bc\"")))
        }
    };
    Ok(Box::new(axtb_gxtb::Calculator::new(version, &req.numbers, &req.positions, req.charge, uhf, o)?))
}

/// A Hessian in the response.
fn hessian_part(h: axtb::Hessian) -> HessianPart {
    HessianPart {
        hessian: Matrix::square(h.hessian),
        hessian_method: names::hessian_method_label(h.method),
        hessian_max_asymmetry: h.max_asymmetry,
        hessian_timings: Named::new(h.timings),
        dipole_derivatives: h.dipole_derivatives,
    }
}

/// The attempt of a failed try `tried`: its error's iterations and code.
fn failed_try(tried: Attempt, e: &axtb::Error) -> Attempt {
    let iterations = match e {
        axtb::Error::ScfNotConverged { iterations, .. } => Some(*iterations),
        _ => None,
    };
    Attempt { iterations, converged: false, error: Some(Error::from(e.clone()).code()), ..tried }
}

/// A composite's try as an attempt.
fn try_attempt(t: &Try) -> Attempt {
    let tried = Attempt {
        part: t.part,
        retry: Some(t.rung),
        etemp: Some(t.etemp),
        max_iter: Some(t.max_iter),
        start: t.start,
        mixer: t.mixer,
        ..Attempt::default()
    };
    match &t.error {
        None => Attempt { iterations: t.iterations, converged: true, ..tried },
        Some(e) => failed_try(tried, e),
    }
}

/// `e` with the tries of a ladder (an SCF that did not converge carries them).
fn with_attempts(e: Error, attempts: Vec<Attempt>) -> Error {
    match e {
        Error::ScfNotConverged { iterations, energy_change, residual, .. } => {
            Error::ScfNotConverged { iterations, energy_change, residual, attempts }
        }
        other => other,
    }
}

/// Fill a single point and its Hessian into `out`: the energy, gradient, charges, dipole,
/// terms and timings of every method, the SCF group of an SCF method, GFN-FF's and a
/// composite's groups (`rungs`: with each part's rung on its `scf_retries`); returns the
/// orbitals (when kept).
fn fill(
    out: &mut Response,
    mut sp: SinglePoint,
    hessian: Option<axtb::Hessian>,
    rungs: bool,
) -> Option<axtb::Orbitals> {
    if let Some(report) = sp.extras.composite.take() {
        out.composite = Some(composite_response(*report, hessian.as_ref().map(|h| h.parts.as_slice()), rungs));
    }
    if let Some(g) = sp.extras.gfnff.take() {
        let t = g.topology;
        out.gfnff = Some(Gfnff {
            n_hbonds: g.n_hbonds[0] + g.n_hbonds[1],
            n_xbonds: g.n_xbonds,
            topology: Topology {
                bonds: t.bonds,
                pi_bond_orders: t.pi_bond_orders,
                bond_types: t.bond_types,
                hybridization: t.hybridization,
                pi_atoms: t.pi_atoms,
                fragments: t.fragments,
                fragment_charges: t.fragment_charges,
                topology_charges: t.topology_charges,
                n_angles: t.n_angles,
                n_torsions: t.n_torsions,
            },
        });
    }
    out.hessian = hessian.map(hessian_part);
    let nat = sp.charges.len();
    out.energy = Some(sp.energy);
    out.gradient = sp.gradient;
    out.dipole = Some(sp.dipole);
    out.terms = Some(Named::new(sp.terms.iter().map(|t| (t.name, t.value)).collect()));
    out.subtotals = sp.terms.iter().filter(|t| t.kind == TermKind::Subtotal).map(|t| t.name).collect();
    out.timings = Some(Named::new(sp.timings));
    out.charges = Some(sp.charges);
    let e = sp.electronic?;
    let [a, b] = axtb::properties::alpha_beta(&e.levels);
    let s = e.scf;
    out.eigensolver = Some(s.eigensolver);
    out.scf = Some(Scf {
        iterations: s.iterations,
        guess: s.start,
        mixer_used: s.mixer_used,
        scf_attempts: (!s.attempts.is_empty()).then(|| s.attempts.iter().map(Attempt::mixer).collect()),
        homo_lumo_gap: e.homo_lumo_gap,
        fermi_levels: e.fermi_levels,
        electronic_temperature: e.electronic_temperature,
        quadrupole: e.quadrupole,
        shell_charges: e.shell_charges,
        atomic_dipoles: e.atomic_dipoles,
        atomic_quadrupoles: e.atomic_quadrupoles,
        orbital_energies_alpha: a.energies,
        orbital_energies_beta: b.energies,
        occupations_alpha: a.occupations,
        occupations_beta: b.occupations,
        bond_orders: e.bond_orders.map(|bo| Matrix { rows: nat, cols: nat, data: bo }),
        overlap: e.overlap.map(Matrix::square),
        hamiltonian: e.hamiltonian.map(|h| h.into_iter().map(Matrix::square).collect()),
        core_electrons: e.core_electrons,
        shell_populations: e.shell_populations,
        shell_atoms: e.shell_atoms,
        shell_angular_momenta: e.shell_angular_momenta,
        reference_shell_occupations: e.reference_shell_occupations,
        quadrupole_traceless: e.quadrupole_traceless,
    });
    e.orbitals
}

/// Refuse a non-finite energy, gradient or Hessian (JSON cannot hold them; they would be
/// `null`).
fn check_finite(out: &Response) -> Result<()> {
    let bad = |what: &str| Err(Error::Numerical(format!("the calculation gave a non-finite {what}")));
    if out.energy.is_some_and(|e| !e.is_finite()) {
        return bad("energy");
    }
    if out.gradient.as_ref().is_some_and(|g| g.iter().flatten().any(|x| !x.is_finite())) {
        return bad("gradient");
    }
    if out.hessian.as_ref().is_some_and(|h| h.hessian.data.iter().any(|x| !x.is_finite())) {
        return bad("Hessian");
    }
    Ok(())
}

/// A composite's report as the response's `composite` group (`hessian`: the parts of its
/// Hessian; `rungs`: with each part's rung on its `scf_retries`).
fn composite_response(
    r: axtb::method::CompositeReport,
    hessian: Option<&[axtb::method::PartReport]>,
    rungs: bool,
) -> Composite {
    let part = |p: &axtb::method::PartReport| {
        let scf = p.scf.as_ref();
        let attempts = scf.map_or(&[][..], |s| &s.attempts[..]);
        let part = Part {
            method: p.label.clone(),
            conventions: p.conventions,
            energy: p.energy,
            iterations: scf.map_or(0, |s| s.iterations),
            mixer_used: scf.and_then(|s| s.mixer_used),
            scf_attempts: attempts.iter().map(Attempt::mixer).collect(),
            etemp: p.etemp,
            max_iter: p.max_iter,
            retry: p.retry.clone(),
            rung: rungs.then_some(p.rung),
            seconds: p.seconds,
            hessian_method: p.hessian_method.map(names::hessian_method_label),
        };
        (p.key, part)
    };
    Composite {
        composite: Details(r.description),
        components: Named(r.components),
        parts: Parts(r.parts.iter().map(part).collect()),
        hessian_parts: hessian.map(|h| Parts(h.iter().map(part).collect())),
    }
}
