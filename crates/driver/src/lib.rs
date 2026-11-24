// SPDX-License-Identifier: LGPL-3.0-or-later
//! axtb-driver: the one request/response layer of the axtb front ends.
//!
//! The CLI fills a [`Request`] from its arguments or a JSON request, Python from keyword
//! arguments; [`run`] checks it, dispatches it to GFN2-xTB (`axtb`), g-xTB and its solvation
//! composite (`axtb-gxtb`) or GFN-FF (`axtb-gfnff`), and returns a
//! [`Response`] with the same fields for every front end (the CLI writes it as JSON, Python
//! returns it as a dict). A stateful [`Calculator`] (under Python's `axtb.Session`) builds the
//! request's method behind axtb-base's `Method` trait and runs a [`Task`] at one geometry
//! after another, each SCF from the previous geometry's [`Wavefunction`] as the method's
//! restart policy says; its forks ([`Calculator::fork`]) share the method's setup and compute
//! beside it (a session's lineages side by side).
//!
//! Rules every front end gets from here:
//! - one alias rule for names ([`names`]): case and anything but letters and digits are
//!   ignored; responses report the canonical names;
//! - unknown fields, settings of another method and settings without their calculation (a
//!   Hessian's step without a Hessian) are [`Error::Input`], never ignored;
//! - one error type with kinds (input, unsupported element, topology, SCF not converged,
//!   numerical), [`Error::kind`];
//! - at least one thread, no clamps; a one-off task ([`run`]) of a molecule small enough that
//!   threads only cost runs on one ([`task_threads`]: the same results); axtb's internal cap
//!   on the threads one calculation uses by method, calculation and size ([`thread_cap`],
//!   applied by Python; the CLI runs on the threads a request gives, default 1); warnings about
//!   more
//!   than 8 threads for an SCF method's single point ([`threads_warning`]), about unpaired
//!   electrons for GFN-FF ([`uhf_warning`]) and about a run-time LAPACK that was not used
//!   ([`eigensolver_warning`]) are the front ends' to print;
//! - XYZ files in Ångström are converted with the factor of the method's conventions
//!   ([`read_xyz`]);
//! - every response names the numerics identity of what it ran ([`numerics`]).
//!
//! # Numerics epochs
//!
//! [`numerics`] gives the identity of the code that produced a result, per method (afh-stack
//! decision 29): `"<method>:<epoch>"`, e.g. `"gfn2:1"`, `"gxtb-v1.1:1"`, `"gxtb-v2.0.1:1"`,
//! `"gfnff:1"`, and for the g-xTB solvation composite its parts joined, `"gxtb-v2.0.1:1+gfn2:1"`.
//! A method's epoch changes with every commit whose results for an unchanged request of that
//! method may change (classed R, round-off, or M, a model or default change, in the changelog),
//! and never with a commit that keeps every result bit for bit (B). One epoch covers every
//! conventions of a method: a change in any of them bumps it. GFN2-xTB in solvent is GFN2-xTB's.
//! The vibrational analysis belongs to every method's numbers (a change there bumps every
//! method), and the analysis of a given Hessian has its own epoch ([`ANALYSIS_NUMERICS`]). Thread
//! counts and a run-time LAPACK move results only by round-off and are not part of it. The
//! epochs start at 1 with alpha 1 of 0.2; their history is in the changelog.

mod calculator;
mod error;
pub mod names;
mod request;
mod response;
mod task;

pub use calculator::{gxtb_solvent_hint, Calculator, Wavefunction};
pub use error::{Error, Result};
pub use names::{
    canonical, canonical_method, shell_layout, supported_elements, valence_electrons, Method, ShellLayout,
};
pub use request::{Request, ScfRetry, Temperatures};
pub use response::{
    Attempt, Composite, Contribution, Details, Gfnff, HessianPart, Matrix, Named, Part, Parts, Response, Retried, Scf,
    Thermo, ThermoPoint, Topology, Vibrations, PROGRAM, VERSION,
};
pub use task::{HessianTask, MoldenTask, Task, ThermoSettings, VibSettings};

use axtb::units::{reference::MCTC_AATOAU, AATOAU};
use std::time::Instant;

/// The numerics epoch of each method ([`numerics`], "Numerics epochs" above): bumped by every
/// commit classed R or M that changes that method's results, never by a B commit.
pub const NUMERICS_EPOCHS: [(Method, u32); 4] =
    [(Method::Gfn2, 1), (Method::GxtbV1_1, 1), (Method::GxtbV2_0_1, 1), (Method::Gfnff, 1)];

/// The numerics identity of the analysis of a given Hessian (`analyze_hessian`).
pub const ANALYSIS_NUMERICS: &str = "analysis:1";

/// The numerics identity of `method` (afh-stack decision 29; "Numerics epochs" above):
/// `"<canonical name>:<epoch>"`, the same for every conventions of the method; with `composite`,
/// the g-xTB solvation composite (a g-xTB method with `solvation_correction`), its parts joined by
/// `+`, the g-xTB part first: `"gxtb-v2.0.1:1+gfn2:1"` (`composite` means nothing for the other
/// methods: GFN2-xTB in solvent is GFN2-xTB's).
pub fn numerics(method: Method, composite: bool) -> String {
    let one = |m: Method| {
        let epoch = NUMERICS_EPOCHS.iter().find(|(x, _)| *x == m).map(|(_, e)| *e).expect("every method has an epoch");
        format!("{}:{epoch}", m.name())
    };
    if composite && method.is_gxtb() {
        format!("{}+{}", one(method), one(Method::Gfn2))
    } else {
        one(method)
    }
}

/// Run a request: the calculation it describes, or with `analyze_hessian` the vibrational
/// analysis of the given Hessian ([`analyze`]). `seconds` of the response includes the
/// setup. The calculation runs on [`task_threads`] threads.
pub fn run(req: &Request) -> Result<Response> {
    // a method this build leaves out is refused before any other validation
    if let Some(method) = req.method.as_deref().and_then(|name| Method::parse(name).ok()) {
        method.check_available()?;
    }
    let t0 = Instant::now();
    let mut out = if req.analyze_hessian.is_some() {
        analyze(req)?
    } else {
        let task = Task::from_request(req)?;
        let threads = task_threads(req, &task);
        let mut calc = if req.threads.is_some_and(|t| t != threads) {
            Calculator::new(&Request { threads: Some(threads), ..req.clone() })?
        } else {
            Calculator::new(req)?
        };
        calc.compute(&req.positions, &task)?
    };
    out.seconds = Some(t0.elapsed().as_secs_f64());
    Ok(out)
}

/// The threads one task of `req` runs on in [`run`]: the requested ones, or one for a molecule
/// of at most [`single_thread_atoms`] atoms, where starting the pool and synchronizing its
/// parallel loops cost more than they save. Systems that small (far below 200 basis functions)
/// give the same bits on any thread count. A [`Calculator`] of its own (under Python's
/// `axtb.Session`) keeps the requested threads: its tasks vary.
pub fn task_threads(req: &Request, task: &Task) -> usize {
    let threads = req.threads.unwrap_or(1);
    if threads <= 1 {
        // (zero threads is the calculator's input error)
        return threads;
    }
    match Method::parse(req.method.as_deref().unwrap_or("gfn2")) {
        Ok(method) if req.numbers.len() <= single_thread_atoms(method, task.hessian.is_some()) => 1,
        _ => threads,
    }
}

/// The largest molecule (atoms) whose single point (`hessian` false) or Hessian of `method`
/// runs on one thread in [`run`]: up to this size one thread was faster than 8 for every
/// molecule measured (laptop, i9-11900H; energy + gradient: water, methane, ammonium, ethanol,
/// glycine, benzene, caffeine, ferrocene, taxol and alkanes up to 152 atoms for GFN-FF; the
/// method's default Hessian). g-xTB's Hessians gain from threads at every size.
pub fn single_thread_atoms(method: Method, hessian: bool) -> usize {
    match (method, hessian) {
        (Method::Gfn2, false) => 10,
        (Method::Gfn2, true) => 6,
        (Method::GxtbV2_0_1, false) => 9,
        (Method::GxtbV1_1, false) => 6,
        (Method::GxtbV1_1 | Method::GxtbV2_0_1, true) => 0,
        (Method::Gfnff, false) => 100,
        (Method::Gfnff, true) => 10,
    }
}

/// Thread caps ([`THREAD_CAPS`]): (from this many atoms, at most this many threads), from 0
/// atoms on.
pub type Caps = &'static [(usize, usize)];

/// axtb's internal cap on the threads of one energy or gradient, per method, by size: at most
/// this many of the threads a calculation is given, where more do not make it faster (the
/// plateau: the fewest threads within about 10 % of the fastest wall time measured), and never
/// more than 8, since beyond 8 the SCF eigensolver only gets slower ([`thread_cap`]). Measured
/// with the release CLI 0.2.0a1 on an AMD Threadripper PRO 9985WX node, 1 to 16 threads (mb33,
/// beni-acetone, caffeine, gh106, taxol, (H2O)100: 16 to 300 atoms; the CLI's whole process): 1
/// thread up to about 40 atoms (more were slower, or gained within the noise of a few
/// milliseconds), 4 at taxol (113 atoms: g-xTB 2.0.1 218 ms on 1 thread, 126 ms on 4, 134 ms on
/// 8; a 119-atom gradient 245 ms on 1, 140 ms on 8, 187 ms on 64), 8 from (H2O)100 (300 atoms:
/// g-xTB 2.0.1 912 ms on 1, 408 ms on 8). GFN-FF on an Intel i9-11900H laptop: no gain to 300
/// atoms, 4 from crambin (640 atoms: 79 ms on 1, 58 ms on 4, 57 ms on 8, 78 ms on 16).
pub const THREAD_CAPS: [(Method, Caps); 4] = [
    (Method::Gfn2, &[(0, 1), (60, 4), (300, 8)]),
    (Method::GxtbV2_0_1, &[(0, 1), (60, 4), (300, 8)]),
    (Method::GxtbV1_1, &[(0, 1), (60, 4), (300, 8)]),
    (Method::Gfnff, &[(0, 1), (600, 4)]),
];

/// The most threads one energy or gradient (`hessian` false) or Hessian of `method` for a
/// molecule of `atoms` atoms uses of the threads it is given (`None`: all of them): 1 at or
/// below [`single_thread_atoms`]; an energy or gradient its [`THREAD_CAPS`] band; a Hessian, whose
/// displacements and perturbations run side by side, no cap. Python's `axtb.compute` and
/// `axtb.Session` apply it; the CLI runs on the threads a request gives ([`task_threads`]).
pub fn thread_cap(method: Method, hessian: bool, atoms: usize) -> Option<usize> {
    if atoms <= single_thread_atoms(method, hessian) {
        return Some(1);
    }
    if hessian {
        return None;
    }
    let (_, caps) = THREAD_CAPS.iter().find(|(m, _)| *m == method).expect("every method has thread caps");
    caps.iter().rev().find(|(from, _)| atoms >= *from).map(|&(_, n)| n)
}

/// The vibrational analysis of a given Hessian (`analyze_hessian`, 3N x 3N, Eh/bohr²) at
/// `positions`, with optional `dipole_derivatives` (IR intensities) and `energy` (total
/// enthalpies and free energies), as `xtb thermo` does. The calculation's fields are refused;
/// the conventions are `analysis_conventions` (default `axtb`).
pub fn analyze(req: &Request) -> Result<Response> {
    let t0 = Instant::now();
    let Some(h) = &req.analyze_hessian else {
        return Err(Error::Input("an analysis needs the Hessian (analyze_hessian)".into()));
    };
    let calculation = [
        ("method", req.method.is_some()),
        ("conventions (use analysis_conventions)", req.conventions.is_some()),
        ("charge", req.charge != 0.0),
        ("uhf", req.uhf.is_some()),
        ("accuracy", req.accuracy.is_some()),
        ("etemp", req.etemp.is_some()),
        ("max_iter", req.max_iter.is_some()),
        ("scf_retries", req.scf_retries.is_some()),
        ("threads", req.threads.is_some()),
        ("dispersion", req.dispersion.is_some()),
        ("mixer", req.mixer.is_some()),
        ("scf_trace", req.scf_trace),
        ("alpb", req.alpb.is_some()),
        ("gbsa", req.gbsa.is_some()),
        ("solvation_state", req.solvation_state.is_some()),
        ("solvation_correction", req.solvation_correction.is_some()),
        ("correction_conventions", req.correction_conventions.is_some()),
        ("guess", req.guess.is_some()),
        ("restart", req.restart.is_some()),
        ("fragment_charges", req.fragment_charges.is_some()),
        ("topology_positions", req.topology_positions.is_some()),
        ("gradient", req.gradient),
        ("hessian", req.hessian),
        ("hessian_method", req.hessian_method.is_some()),
        ("hessian_step", req.hessian_step.is_some()),
        ("hessian_accuracy", req.hessian_accuracy.is_some()),
        ("molden", req.molden || req.molden_file.is_some() || req.molden_cartesian),
        ("bond_orders", req.bond_orders),
        ("overlap", req.overlap),
        ("hamiltonian", req.hamiltonian),
        ("vibrations", req.vibrations.is_some()),
    ];
    let given: Vec<&str> = calculation.iter().filter(|(_, g)| *g).map(|(k, _)| *k).collect();
    if !given.is_empty() {
        return Err(Error::Input(format!(
            "an analysis of a given Hessian (analyze_hessian) takes no calculation settings: {}",
            given.join(", ")
        )));
    }
    let nat = req.numbers.len();
    let n3 = 3 * nat;
    if req.positions.len() != nat {
        return Err(Error::Input(format!("positions must have one row per atom ({nat}), not {}", req.positions.len())));
    }
    if h.len() != n3 || h.iter().any(|r| r.len() != n3) {
        return Err(Error::Input(format!("analyze_hessian must be {n3} x {n3}")));
    }
    if req.dipole_derivatives.as_ref().is_some_and(|d| d.len() != n3) {
        return Err(Error::Input(format!("dipole_derivatives must be {n3} x 3")));
    }
    let flat: Vec<f64> = h.iter().flatten().copied().collect();
    let settings = VibSettings::from_request(req)?;
    let defaults = axtb::method::AnalysisDefaults { conventions: Default::default(), xtb_frequency_scale: None };
    let opts = settings.options(&defaults);
    let v = axtb::vibrations::analyze(&req.numbers, &req.positions, &flat, req.dipole_derivatives.as_deref(), &opts)?;
    Ok(Response {
        program: PROGRAM,
        version: VERSION,
        method: None,
        conventions: v.conventions.label(),
        numerics: ANALYSIS_NUMERICS.into(),
        energy: req.energy,
        vibrations: Some(Vibrations::new(v, req.energy)),
        seconds: Some(t0.elapsed().as_secs_f64()),
        ..Default::default()
    })
}

/// Atomic numbers and positions (bohr) of an XYZ file (Ångström) for `method` in
/// `conventions` (default `axtb`), converted as the reproduced program reads XYZ files:
/// mctc-lib's Ångström-to-bohr factor for GFN2-xTB in the `xtb`, `xtb-6.7.1` and `tblite`
/// conventions, GFN-FF in `xtb` and `xtb-6.7.1` and g-xTB 2.0.1 in `gxtb`; the g-xTB 1.1.0
/// program's for g-xTB v1.1; CODATA 2022 otherwise.
pub fn read_xyz(text: &str, method: Method, conventions: Option<&str>) -> Result<(Vec<usize>, Vec<[f64; 3]>)> {
    let aatoau = match method {
        Method::Gfn2 => match conventions.map(names::gfn2_conventions).transpose()?.unwrap_or_default() {
            axtb::Conventions::Axtb => AATOAU,
            _ => MCTC_AATOAU,
        },
        Method::Gfnff => {
            if names::gfnff_conventions(conventions.unwrap_or("axtb"))?.is_xtb() {
                MCTC_AATOAU
            } else {
                AATOAU
            }
        }
        #[cfg(feature = "gxtb")]
        Method::GxtbV1_1 => {
            let mol = axtb_gxtb::read_xyz(text, 0.0, None)?;
            return Ok((mol.numbers, mol.positions));
        }
        #[cfg(feature = "gxtb")]
        Method::GxtbV2_0_1 => {
            let c = conventions.map(names::gxtb_conventions).transpose()?.unwrap_or_default();
            let mol = axtb_gxtb::model::read_xyz(text, 0.0, None, c)?;
            return Ok((mol.numbers, mol.positions));
        }
        #[cfg(not(feature = "gxtb"))]
        Method::GxtbV1_1 | Method::GxtbV2_0_1 => return Err(Error::gxtb_unavailable()),
    };
    Ok(axtb::io::read_xyz(text, aatoau)?)
}

/// The warning for one single point of an SCF method on `threads` threads, more than help
/// (`axtb::threads`); Hessians, which run their displacements in parallel, and GFN-FF (no
/// SCF) get none.
pub fn threads_warning(method: Method, threads: usize, hessian: bool) -> Option<String> {
    if method.is_scf() && !hessian {
        axtb::threads::scf_threads_warning(threads)
    } else {
        None
    }
}

/// The warnings of a response whose result came from a rung of `scf_retries`: "converged on
/// retry 2 of 2: etemp 3000 K, max_iter 1000"; the composite: one per part whose single point
/// converged on a rung ("the gxtb part converged on ..."), and one per part whose Hessian
/// converged on another rung than its single point ("the gxtb part's Hessian converged on ...").
pub fn retry_warnings(response: &Response) -> Vec<String> {
    let Some(r) = &response.scf_retry else { return Vec::new() };
    let Some(c) = &response.composite else {
        return r.label.iter().map(|l| format!("converged on {l}")).collect();
    };
    fn fired(p: &Part) -> Option<&str> {
        p.rung.is_some_and(|i| i > 0).then_some(p.retry.as_deref()).flatten()
    }
    let mut out: Vec<String> =
        c.parts.0.iter().filter_map(|(k, p)| fired(p).map(|l| format!("the {k} part converged on {l}"))).collect();
    for (k, p) in c.hessian_parts.iter().flat_map(|h| &h.0) {
        if let Some(l) = fired(p).filter(|_| c.parts.get(k).is_none_or(|s| s.rung != p.rung)) {
            out.push(format!("the {k} part's Hessian converged on {l}"));
        }
    }
    out
}

/// The warning for unpaired electrons given to GFN-FF, which has no electronic structure, so
/// that `uhf` has no effect (a positive value is accepted, so that a radical's multiplicity
/// from a structure file passes through, but is not silently ignored).
pub fn uhf_warning(method: Method, uhf: Option<i64>) -> Option<String> {
    (method == Method::Gfnff && uhf.is_some_and(|u| u > 0))
        .then(|| "GFN-FF has no electronic structure: uhf (unpaired electrons) has no effect".to_string())
}

/// The warning for a response whose SCF ran faer although a run-time LAPACK was asked for
/// (`AXTB_LAPACK` refused, or the library failed in a call; `SinglePoint::eigensolver`).
pub fn eigensolver_warning(response: &Response) -> Option<String> {
    let e = response.eigensolver.as_deref()?;
    e.strip_prefix("faer (")
        .and_then(|s| s.strip_suffix(')'))
        .map(|why| format!("the SCF eigensolver ran faer instead of the run-time LAPACK: {why}"))
}

/// The warning when `AXTB_LAPACK` names a library that cannot be used, for a failed request
/// (a successful one: [`eigensolver_warning`]). Loads the library if no SCF has.
pub fn lapack_warning() -> Option<String> {
    axtb::linalg::runtime_lapack().err().map(|why| format!("AXTB_LAPACK not used, faer instead: {why}"))
}

/// Why this build cannot run on this CPU, or `None`. The x86-64 builds assume x86-64-v3
/// (`.cargo/config.toml`: AVX, AVX2, FMA, BMI1, BMI2, F16C, LZCNT, MOVBE, XSAVE), which CPUs
/// before Intel Haswell (2013) and AMD Excavator (2015) and many low-end Intel Atom, Celeron and
/// Pentium (N-series) CPUs lack; there the first such instruction would stop the process with
/// SIGILL. The CLI calls this first thing; the Python package checks `/proc/cpuinfo` in
/// `axtb/__init__.py` before it loads the extension module, with the same message. It asks the
/// CPU (`cpuid`, and `xgetbv` for the operating system's AVX support) instead of
/// `is_x86_feature_detected!`, which a build for these features answers at compile time, and it
/// is out of line and free of vector instructions (checked in the disassembly).
#[inline(never)]
pub fn cpu_mismatch() -> Option<&'static std::ffi::CStr> {
    #[cfg(all(target_arch = "x86_64", any(target_feature = "avx2", target_feature = "fma")))]
    if !x86_64_v3() {
        return Some(CPU_MISMATCH);
    }
    None
}

/// Whether the CPU and the operating system support x86-64-v3.
#[cfg(all(target_arch = "x86_64", any(target_feature = "avx2", target_feature = "fma")))]
#[inline(never)]
fn x86_64_v3() -> bool {
    use std::arch::x86_64::{__cpuid, __cpuid_count, _xgetbv};
    let bit = |reg: u32, b: u32| reg & (1 << b) != 0;
    // `__cpuid` is unsafe in Rust 1.89 (the minimum supported version) and safe in newer ones
    #[allow(unused_unsafe)]
    // SAFETY: cpuid exists on every x86-64 CPU; leaves above the maximum are not asked
    let (max, max_ext) = unsafe { (__cpuid(0).eax, __cpuid(0x8000_0000).eax) };
    if max < 7 || max_ext < 0x8000_0001 {
        return false;
    }
    #[allow(unused_unsafe)]
    // SAFETY: leaves 1, 7 and 0x80000001 exist (checked above)
    let (l1, l7, ext) = unsafe { (__cpuid(1), __cpuid_count(7, 0), __cpuid(0x8000_0001)) };
    // leaf 1 ECX: FMA 12, MOVBE 22, XSAVE 26, OSXSAVE 27, AVX 28, F16C 29
    if ![12, 22, 26, 27, 28, 29].iter().all(|&b| bit(l1.ecx, b)) {
        return false;
    }
    // SAFETY: OSXSAVE is set, so the OS enabled xgetbv
    let xcr0 = unsafe { _xgetbv(0) };
    // the OS saves the SSE and AVX registers; leaf 7 EBX: BMI1 3, AVX2 5, BMI2 8; LZCNT
    xcr0 & 0b110 == 0b110 && [3, 5, 8].iter().all(|&b| bit(l7.ebx, b)) && bit(ext.ecx, 5)
}

/// The message of [`cpu_mismatch`].
#[cfg(all(target_arch = "x86_64", any(target_feature = "avx2", target_feature = "fma")))]
const CPU_MISMATCH: &std::ffi::CStr = c"this axtb was built for x86-64-v3 CPUs (AVX2, FMA, BMI2: Intel since Haswell (2013), AMD \
since Excavator (2015); not many low-end Intel Atom, Celeron and Pentium N-series CPUs), and this CPU (or its operating \
system) does not support them. Build axtb from source for any x86-64 CPU with RUSTFLAGS=\"-C target-cpu=x86-64\" (CLI: \
cargo install axtb-cli; Python: pip install --no-binary axtb axtb)";
