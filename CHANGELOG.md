# Changelog

All notable changes to axtb. Versions follow [Semantic Versioning](https://semver.org); the
Rust crates, the CLI and the Python package share one version.

## 0.2.0 (2026-11-19)

The final release of the 0.2.0 alphas, for the final qc* releases: the alphas' entries below
make up its changes since 0.1.2. Since 0.2.0-alpha.3:

- The Python package requires `qcdata>=0.24.0,<0.25` (the final 0.24.0) and its `test` extra
  `qccompute[testing]>=0.23.0,<0.24`; the README installs without `--pre`.
- No result changes: the numerics epochs stay `gfn2:1`, `gxtb-v1.1:1`, `gxtb-v2.0.1:1` and
  `gfnff:1`.
- qccompute 0.23.0 lowers a `threads` keyword above its runner's count with a `RuntimeWarning`;
  axtb still reports `resource_capped` and runs on the lowered count (tested).

## 0.2.0-alpha.3 (2026-11-17)

### Threads: one count per calculation, axtb's internal cap within it

- Python runs a calculation on the count it is given: the `threads` keyword (at most the grant,
  with `resource_capped`), else the runner's grant, else (a direct call) this process's CPUs.
  Within it axtb uses at most its internal cap, where more threads stop making a calculation
  faster: SCF energies and gradients at most 8 threads at any size and fewer for small molecules
  (GFN2-xTB, g-xTB 2.0.1 and v1.1: 1 below 60 atoms, 4 from 60, 8 from 300; GFN-FF: 1 below 600,
  4 from 600); Hessians all they are given; at or below the one-thread size 1, as before. The
  `threads` keyword is a maximum too. Before, a direct call without the keyword ran on the plan's
  bands and a grant was used whole.
- The plan states no thread count: `axtb.program.plan` and `Session.plan` no longer set
  `resources.threads` (qccompute chooses the count; a runner that still reads the field takes 1).
  Their record gives `thread_cap`, the most threads the calculation uses (None: no cap).
- Sessions share their count like a list: a single `compute` runs alone (concurrent callers wait
  their turn) on the whole count, and `compute_many` runs min(lineages, count) lineages at once,
  each request on `count // k` at most its cap, instead of per-calctype bands packed into the
  grant. A session's count is its `threads` keyword, else the grant, else this process's CPUs. A
  frame equals a one-at-a-time `compute` on the same thread count, bit for bit; a large molecule's
  bits can depend on k (round-off), and the `threads` keyword fixes them.
- `Session.plan()`'s `resources`, `record.threads_used`, `threads_by_calctype` and
  `lineages_at_once` are gone; `resolved.settings.threads` is the session's count.
- Rust: axtb-driver's `THREAD_BANDS` and `planned_threads` are `THREAD_CAPS` (energies and
  gradients) and `thread_cap(method, hessian, atoms) -> Option<usize>`. The CLI is unchanged: it
  runs on a request's threads (default 1) with the one-thread rule.
- Numerics epochs unchanged: no result of a given thread count changes (golden set 188 of 188).

## 0.2.0-alpha.2 (2026-11-17)

### Sessions run lineages side by side

- `axtb.Session.compute_many` runs different lineages side by side on forks of the session's
  calculator, as long as their threads fit into the session's threads in all (the grant; a direct
  call: this process's CPUs), taken in request order; `plan()["record"]` gives
  `lineages_at_once` for the input's calctype. Each lineage's requests run in request order; a
  request whose lineage has no history waits for its seeds' earlier requests in the batch. The
  frames are what one-at-a-time `compute` calls in request order give, bit for bit, failures in
  place. The forks share the method's immutable setup, run on one thread count each, have their
  own SCF state and work buffers, and are kept until the session closes (one calculation's working memory each: 1.9 MiB for a caffeine g-xTB 2.0.1
  gradient, 24 MiB for taxol's). Measured on the laptop, 12 caffeine g-xTB 2.0.1 gradient
  lineages per call in a session of 8 threads: 31 ms per call, against 160 ms one request at a
  time on 1 thread (5.1x) and 93 ms one at a time on 8 threads, as 0.2.0a1 ran them.
- A session's requests, single `compute` calls included, run on their own threads (the
  `threads` keyword, else the plan's for the request's calctype and the molecule's size, at most
  the grant), no longer on the whole grant: a gradient on the gradient's threads, a Hessian on
  the Hessian's, each on a fork on those threads. `plan()` gives the input calctype's threads,
  the most a calculation takes as `resources.threads`, and `record.threads_by_calctype`. A single
  `compute` stays serialized by the session's lock.
- Rust: `Method::fork` (axtb-base; GFN2-xTB, both g-xTB versions, the composite and GFN-FF
  implement it), `Composite::forked`, `Part::fork`, axtb-driver `Calculator::fork`. Breaking:
  `method::Guard` is an `Arc<dyn Fn>` (was a `Box<dyn FnMut>`), so that forks share it.
- Numerics epochs unchanged: no result changes (golden set 188 of 188).

### Plan threads from measurements

- The plan states a thread count for every calculation by method, calctype and size, never the
  whole grant (`resources.threads` is no longer None for a Hessian): the largest count whose
  parallel efficiency stayed at least 0.5 in thread-scaling runs (1-16 threads on an AMD
  Threadripper PRO 9985WX node, 3-300 atoms; larger energies and gradients on a laptop).
  Energies and gradients: 1 below 100 atoms, then 2, and 4 from 600 (g-xTB v1.1: from 300;
  GFN-FF: 1 below 600, then 2), where the plan said 8 from 300 (GFN-FF from 1000). Hessians:
  g-xTB 2.0.1 2/4/8/16 from 0/10/20/60 atoms, v1.1 4/16 from 0/10, GFN2-xTB 1/4/8 from 0/60/150,
  GFN-FF 1/2 from 0/100. The `threads` keyword still wins; a direct `axtb.compute` runs on the
  plan's threads, under a runner on the grant. The table is `axtb_driver::THREAD_BANDS`
  (`planned_threads`), and docs/python.md's Threads section is generated from it.
- The CLI keeps its default of 1 thread and its `task_threads` (a request's threads, 1 for a
  molecule so small that threads only cost): the golden set's cases depend on both.

- `scf.retries` (Python), `scf_retries` (the JSON request), `--scf-retry etemp=K,max_iter=N` (the
  CLI, repeated): an opt-in ladder of rungs `{etemp, max_iter}` (either key optional; a missing one
  keeps the calculation's value) for GFN2-xTB and g-xTB, both versions, every conventions. After
  a calculation whose SCF fails (does not converge, or breaks down), each rung runs it again as a
  fresh SCF from the method's own start, in order, until one converges. Default none: nothing
  changes without it. GFN-FF refuses it (`no_effect`).
  - A ladder that does not fire changes no number; a result from a rung is the calculation at the
    rung's settings, bit for bit (the composite's single point after its Hessian, which starts its
    GFN2 parts fresh at the rung: within the SCF thresholds).
  - The floor: with `scf.etemp` the ladder starts there, and a rung below it is refused
    (`no_effect`, `scf.retries[i].etemp`).
  - Hessians run whole at one rung: a try that fails anywhere (the reference SCF, a displaced
    SCF, the response) reruns the whole Hessian and its single point at the next rung, so no
    rung mixes temperatures inside one Hessian and no displacement climbs alone. (A displacement
    that fails at the reference's settings moves the Hessian to the next rung rather than failing
    it: one Hessian at one temperature instead of an error to retry by hand.)
  - Sessions: a request starts from history as usual, then fresh at its own settings
    (`restart_fallback`), then the rungs; a frame from a rung commits, and the next frame tries the
    calculation's own settings first. A sticky floor is `scf.etemp` on later inputs.
  - The composite: a given ladder replaces the built-in retries (`GXTB_RETRIES`, `PAIR_RETRIES`)
    for every part, each part group climbing it after its first try; the single point after a
    composite Hessian starts each part at its Hessian's rung. Without it the composite is
    unchanged, bit for bit.
  - Recording: `resolved.settings.scf.retries` echoes the ladder (it is in the key);
    `record.contingent["scf.retries"]` has `fired` and the `rung` (the composite: and the `parts`
    that climbed); `record.scf.attempts` lists every try (`retry`, `etemp`, `max_iter`,
    `converged`, `iterations`, `start` or `error`; the composite's `part`), in failed outputs and
    frames too; the warning `scf_retry` names the rung ("converged on retry 2 of 2: etemp 3000 K,
    max_iter 1000"), per part for the composite. The response's `scf_retry` (`rung`, `of`, `label`)
    and `scf_attempts`; the CLI prints the warning on stderr.
  - When every rung fails: kind `convergence`, code `scf_not_converged`, the last rung's
    `iterations`, `energy_change` and `residual`, and every try (`attempts` of the CLI's JSON
    error, `record.scf.attempts` in Python). A last rung whose SCF broke down is
    `numerical_failure`.
- Rust: axtb-base `method::Retry` (`check`, `rungs`), `method::Try`, `Rung::apply`,
  `RetryLadder::climb_methods` (the composite's rung machinery, shared), `Method::tries`,
  `Composite::with_pinned_rungs` and `PartReport::rung`; axtb-gxtb
  `CompositeOptions::scf_retries`; axtb-driver `Request::scf_retries` (`ScfRetry`),
  `Response::scf_retry` (`Retried`), `retry_warnings`, `Part::rung`, and `Attempt` with the
  fields of a try (`mixer` and `iterations` now optional); `Error::ScfNotConverged` carries the
  tries (`attempts`).
- Numerics epochs unchanged: every existing result is bit for bit the same (the golden set,
  177 of 177); 11 golden cases for the ladder (stretched water in each method: a gradient, a
  Hessian, a trajectory; the composite's Hessian; every rung failing).

### Docs and tests

- docs/python.md states how far a session's warm start moves a frame from a fresh `compute`, per
  method and conventions (none where `scf.restart` is off, the default for g-xTB in its programs'
  conventions), that sessions work in forked children, and how results differ across CPU types
  (round-off; measured Intel against AMD).
- docs/python.md: `results.hessian` is the full Cartesian second derivative (3N x 3N,
  symmetrized), never projected; translations and rotations are projected out only in
  `vibrational_analysis` (checked at non-stationary geometries for every Hessian method).
- `crates/axtb-py/tests/test_fork.py`: a forked child computes, whether or not its parent did.

## 0.2.0-alpha.1 (2026-11-17)

The first alpha of 0.2. The Python package is new: axtb is a program of the qc* stack (qcdata
0.24 and qccompute 0.23, both in alpha), with `axtb.compute`, `axtb.Session` and the program
qccompute finds by its entry point. The 0.1 Python API and the qccodec codec are removed, with
no shims; the table below maps every 0.1 call. The CLI, its JSON and the Rust crates keep their
0.1 form apart from the additions below; every number is 0.1.1's except where an entry says
otherwise (the numerics epochs at the end start here). Python 0.2.0a1 needs Python 3.11 or newer
(one abi3 wheel) and `qcdata>=0.24.0a5,<0.25`; the crates are published as `0.2.0-alpha.1` and
require each other at exactly that version.

### Python: axtb.compute, axtb.Session, axtb.program

- `axtb.compute(inp, ctx=None)`: one `energy`, `gradient` or `hessian` calculation of a
  `qcdata.ProgramInput` (a Hessian with its vibrational analysis and thermochemistry), the
  CLI's numbers bit for bit (`scripts/golden_python.py` checks every case of the golden set).
  Failures raise `qcdata.ProgramError` with the failed output attached.
- `axtb.Session(inp, *, calctypes=None, ctx=None)`: one molecule at many geometries on one
  engine (qccompute's session protocol): one committed wavefunction per lineage, `seed`,
  `commit=False`, `compute_many`, `reset`; interleaved lineages compute what separate sessions
  compute, bit for bit; thread-safe, the GIL released. About 30 µs of Python per step, 100 µs
  per `compute` (water).
- `axtb.program`, registered under the `qccompute.programs` entry points (no qccompute import):
  capabilities, `settings_schema`, `record_schema`, `plan`, `compute`, `open_session`,
  `numerics`; `qccompute.testing.check_program` passes for every method, and `python -m
  qccompute.serve axtb` serves it.
- Keywords are grouped settings with qcdata's shared names, checked strictly before anything
  runs, in qcdata's order (`unknown_key` with a did-you-mean from the 0.1 names, `invalid_value`,
  `no_effect`, `not_applicable`, `missing_setting`, the structure, `unsupported`).
  `resolved.settings` echoes every setting with its canonical name and default; docs/python.md
  lists them (generated from the table).
- Results fill qcdata's typed fields, never `extras`: energy components that add up to the
  energy, the Buckingham quadrupole, Mulliken shell populations, 3 x 3 atomic multipoles, the
  orbital levels on every SCF result and on request the orbitals with their basis set, the AO
  overlap and Hamiltonian, bond orders, `VibrationalAnalysis` and one `Thermochemistry` per
  temperature. axtb's detail goes in `record` (the SCF's start, iterations and mixer, the
  composite's parts, GFN-FF's bonds and topology, the accuracy estimate), phase timings in
  `execution.timings`.
- `accuracy` targets (energy, gradient, Hessian) choose the SCF thresholds through an error model
  fitted to measurements (`tests/accuracy_calibration.py`, `bench/results/accuracy.json`; the
  constants in docs/python.md); `record.accuracy` holds the estimate.
- `ProgramInput.given` starts the SCF from a given result: g-xTB from the density of given
  orbitals, GFN2-xTB from given shell populations and atomic multipoles; misfits are refused.
- Threads: the `threads` keyword, at most a runner's grant (`resource_capped`), else the grant,
  else the plan's (one thread for molecules too small for more, at most 8 for SCF energies and
  gradients). The 0.1 `RuntimeWarning`s are gone: warnings are `ProgramWarning`s in the output.
- The multiplicity's parity must be that of the electrons the method's orbitals hold
  (valence-only methods leave the cores out); GFN-FF runs singlets with an even electron count
  only (`not_modeled`; 0.1 ignored `uhf` with a warning).

Migrating from 0.1 (`s` a `qcdata.Structure`, positions in bohr; docs/python.md has every key):

| 0.1 | 0.2 |
|---|---|
| `axtb.singlepoint(numbers, positions, method=M, gradient=True, **settings)` | `axtb.compute(ProgramInput(program="axtb", calctype="gradient", structure=s, model={"method": M}, keywords={...}))` |
| `r["energy"]`, `r["gradient"]`, `r["charges"]`, `r["dipole"]` | `out.results.energy`, `.gradient`, `.atomic_charges`, `.dipole_moment` |
| `r["terms"]` (with the subtotal `electronic`) | `out.results.energy_components` (a partition), `out.record["energy_subtotals"]` |
| `r["iterations"]`, `r["guess"]`, `r["mixer_used"]`, `r["scf_attempts"]` | `out.results.scf_iterations`, `out.record["scf"]` (`start`, `iterations`, `mixer_used`, `attempts`) |
| `r["shell_charges"]`, `r["atomic_dipoles"]`, `r["atomic_quadrupoles"]` (packed) | `out.results.shell_populations` (with `shell_atoms`, `shell_angular_momenta`; charge = reference occupation − population), `.atomic_dipoles`, `.atomic_quadrupoles` (3 x 3) |
| `r["quadrupole"]` (the conventions' form, packed) | `out.results.quadrupole_moment`: the traceless 3 x 3 tensor in every conventions (the CLI keeps the conventions' form) |
| `r["orbital_energies_alpha"]`, `r["occupations_alpha"]`, ..., `r["fermi_levels"]`, `r["electronic_temperature"]` | `out.results.wavefunction.orbitals`: `energies_a/_b`, `occupations_a/_b`, `fermi_level_a/_b`, `electronic_temperature`, `core_electrons` |
| `r["homo_lumo_gap"]` | from the orbital levels (the CLI keeps `homo_lumo_gap`) |
| `bond_orders=True` → `r["bond_orders"]` | keyword `output.bond_orders` → `out.results.bond_orders` |
| `molden=True` → `r["molden"]` | the CLI's `--molden`; or `output.orbitals` (orbitals with their basis set) and qcinf |
| `axtb.Calculator(M, numbers, positions, **settings)`, `.singlepoint(x, gradient=True)` | `axtb.Session(inp)`, `.compute(x, "gradient")` → a `qcdata.Frame` |
| `.hessian(x, step=, accuracy=, method=)` | calctype `hessian`; keywords `hessian.step`, `hessian.accuracy_scale`, `hessian.method` |
| `h["frequencies"]`, `["normal_modes"]`, `["ir_intensities"]`, `["reduced_masses"]`, `["force_constants"]`, `["masses"]` | `out.results.vibrational_analysis`: `frequencies`, `normal_modes_cartesian`, `ir_intensities`, `reduced_masses`, `force_constants`, `masses` |
| `h["thermo"]` (`zpve`, `enthalpy`, `g_rrho`, `free_energy`, `points`) | `out.results.thermochemistry`, one entry per temperature (`zero_point_energy`; `total.enthalpy` = E + H, `total.gibbs` = E + G) |
| `.guess_state()`, `.set_guess_state(state)` | lineages and `seed` in a session; `ProgramInput.given` across inputs |
| `.set_guess(shell_charges, atomic_dipoles, atomic_quadrupoles)` | `ProgramInput.given` with `shell_populations` and the atomic multipoles |
| `.reset_guess()` | `Session.reset()` |
| `guess="previous"`, `"fresh"` | `scf.restart`: `previous`, `off` |
| `guess="eeq"`, `"sad"`, `"eeq-bc"` (that start at every geometry) | `scf.guess` (the start without history; the restart policy stays) |
| `axtb.vibrational_analysis(numbers, positions, hessian, ...)` | the CLI's `analyze_hessian` request (qcinf after alpha 1) |
| `axtb.canonical_method(name)` | `axtb.program.plan(inp)["resolved"]`, `out.record["method"]` |
| `axtb.gxtb_solvent_hint(method)` | the message of the `missing_setting` error |
| `axtb.InputError` (a `ValueError`) | `ProgramError`, kind `input` (codes above; `unsupported_element` with `details.z`, `topology_failed`) |
| `axtb.ScfNotConvergedError` (`iterations`, `energy_change`, `residual`) | `ProgramError`, kind `convergence`, code `scf_not_converged`, the same in `error.details` |
| `RuntimeError` (numerical, panic) | kind `numerical`, `internal` |
| `RuntimeWarning` (threads, LAPACK) | `out.warnings` (`eigensolver_fallback`); no threads advice |
| `charge=`, `uhf=` | `structure.charge`, `structure.multiplicity` (`uhf` + 1) |
| `accuracy`, `etemp`, `max_iter`, `mixer` (`auto`) | `accuracy_scale`, `scf.etemp`, `scf.max_iter`, `scf.mixer` (`broyden-then-diis`); `accuracy` is now the error targets |
| `alpb=S`, `gbsa=S`, `solvation_state`, `solvation_correction`, `correction_conventions` | `solvent.name` (with `solvent.model: gbsa`), `solvent.state`, `solvent.correction`, `solvent.correction_conventions` |
| `fragment_charges`, `topology_positions` | `topology.fragment_charges`, `topology.geometry` |
| `temperature`, `pressure`, `frequency_scale`, `imaginary_cutoff`, `rotor_cutoff`, `symmetry_number`, `symmetry_threshold`, `masses`, `analysis_conventions`, `vibrations=False`, `thermo=False` | `freq.temperature`, `freq.pressure`, `freq.scale`, `freq.imaginary_cutoff`, `freq.rotor_cutoff`, `freq.symmetry_number`, `freq.symmetry_tol`, `freq.masses`, `freq.conventions`, `freq.analysis: off`, `freq.thermo: off` |
| `dispersion=False`, `scf_trace` | the CLI's (validation) |
| conventions `xtb`, `tblite`, `gxtb` | canonical `xtb-6.8.0`, `tblite-0.8.0`, `gxtb-1.1.0` / `gxtb-2.0.1` (the short names stay aliases) |
| `pip install "axtb[qc]"`: the qccodec codec, `propagate_wfn` | `qccompute.compute`, `qccompute.session` through the entry point; warm starts through sessions |

For qcopt: once `find_program("axtb")` reports sessions, every axtb task goes to its
`SessionEvaluator`; `AxtbEvaluator` and its `axtb.Calculator` routing go, the allow-lists take
the dotted 0.2 keys (`scf.*`, `solvent.*`, `hessian.*`, `freq.*`, `output.*`; `accuracy_scale`
stays a surface key), one session with calctypes {gradient, hessian} serves a shared spec (a
gradient-only spec with `hessian.*` is `not_applicable`), a Hessian spec that needs no analysis
sets `freq.analysis: off`, answers read `frame.results` and `scf_iterations`, and IRC branches can
seed from the TS lineage. For mepd: its xtb-style arguments map to `scf.etemp`, `accuracy_scale`
and `scf.max_iter`, its retry ladder to `{"scf": {"etemp": ..., "max_iter": 1000}}`,
`conventions="gxtb"` stays an alias, and its Molden files come from the CLI (`--molden`) or from
`output.orbitals` once qcinf writes them.

### Rust, the CLI and the driver

- For the Python engine (no CLI change): `axtb_driver::valence_electrons(method, z)` (the
  electrons a method's orbitals hold per element: valence-only methods leave the cores out) and
  `shell_layout(method, numbers)` (each shell's atom, angular momentum, reference occupation and
  exponents); `Calculator::set_density` with `Method::set_density` (g-xTB v1.1 and 2.0.1, and the
  composite's g-xTB part, start the next SCF from given density matrices; GFN2-xTB and GFN-FF
  refuse).
- Driver calls for sessions and the Python binding: `Calculator::set_history` (with the new
  `Method::set_history`) installs a wavefunction as the method's own previous one, so the restart
  policy still decides whether the next SCF starts from it (`set_wavefunction` marks it given); a
  session that keeps one wavefunction per lineage computes what one calculator per lineage
  computes, bit for bit. `Request::restart` (`previous` or `fresh`) sets the restart policy beside
  `guess`, which then names the start alone (GFN2-xTB `sad` with `previous`); the old `guess`
  values keep working. `Calculator::set_scf_accuracy` (with `Method::set_scf_accuracy`, every part's
  for the composite) sets the SCF thresholds of later calls; GFN2-xTB keeps the integral screening
  it was built with. `Error::code()` gives a stable code for every error variant
  (`invalid_input`, `unsupported_element`, `topology_failed`, `scf_not_converged`,
  `numerical_failure`, `not_applicable`; the CLI's JSON keeps `kind`); `axtb_driver::canonical`
  gives the canonical name of a method, a method's conventions or a solvent, and
  `supported_elements` the elements a method has parameters for (GFN2-xTB H-Rn, GFN-FF and g-xTB
  2.0.1 H-Lr, v1.1 H-Ce, Yb-Ac and U); `axtb_solvation::canonical_solvent`.
- For the Python binding (in Rust and the driver's response, not in the CLI's JSON): the
  orbitals with their basis, structured, without Molden text (`Task::orbitals`,
  `Response::orbitals`); per atom the core electrons outside the valence-only orbitals, Z less
  the reference shell occupations (`Electronic::core_electrons`; Zn: 28 in GFN2-xTB and g-xTB);
  Mulliken shell populations Σ(P S) of the total density with each shell's atom, angular
  momentum and reference occupation (`Electronic::shell_populations`, `shell_atoms`,
  `shell_angular_momenta`, `reference_shell_occupations`); the traceless (Buckingham) molecular
  quadrupole in every conventions beside the conventions' form (`quadrupole_traceless`); the
  thermochemistry's average moment of inertia B_av (kg m²) and the imaginary frequencies it
  inverts and leaves out (`Thermo::average_moment`, `inverted_frequencies`,
  `dropped_frequencies`). `axtb_base::method::core_electrons`, `units::ELECTRON_MASS_KG`. No
  number and no CLI output changes.
- GFN2-xTB's energy terms partition the energy: a new piece `band`, the band energy Tr(H0 P),
  which xtb reports only inside the electronic energy, so that the pieces (every term but the
  subtotal `electronic`) add up to the energy, as g-xTB's and GFN-FF's do. `dispersion` and
  `solvation` stay whole pieces (their charge-independent and self-consistent parts; the latter
  are also inside `electronic`, a subtotal). The CLI's `terms` gain the key `band`; no number
  changes. The driver's response names the subtotals for the Python binding
  (`Response::subtotals`, not in the CLI's JSON).
- Numerics epochs (afh-stack decision 29): every response names the identity of the code that
  produced it, `numerics`, and `axtb_driver::numerics(method, composite)` gives it before a run.
  The rule: each method has an epoch, `"<method>:<epoch>"` (`"gfn2:1"`, `"gxtb-v1.1:1"`,
  `"gxtb-v2.0.1:1"`, `"gfnff:1"`), one epoch per method across all its conventions; the g-xTB
  solvation composite joins its parts' (`"gxtb-v2.0.1:1+gfn2:1"`), and GFN2-xTB in solvent is
  GFN2-xTB's. A method's epoch changes with every commit that may change its results for an
  unchanged request (classed R, round-off, or M, a model or default change, below) and never
  with one that keeps every result bit for bit (B). The vibrational analysis is part of every
  method's numbers; the analysis of a given Hessian has its own identity, `"analysis:1"`. Thread
  counts and a run-time LAPACK are not part of it (round-off). The history is in "Numerics
  epochs" at the end of this file.
- AO matrices on request: `overlap` (the SCF's S, with its screening; g-xTB: of the SCF basis)
  and `hamiltonian` (the matrix whose generalized eigenvectors with S are the returned orbitals:
  GFN2-xTB's last Fock matrix, one for both spins; g-xTB's unmixed last build per spin channel),
  which brings the overlap. Request fields `overlap` and `hamiltonian` (CLI `--overlap`,
  `--hamiltonian`), response fields of the same names, in Rust `Want::overlap`,
  `Want::hamiltonian`, `Electronic::overlap`, `Electronic::hamiltonian`; axtb's AO order (atom by
  atom, shell by shell, m = -l..l, p as y, z, x). Asking for them changes no other number.
  g-xTB 2.0.1 in the program's conventions refuses the Hamiltonian (its SCF diagonalizes a
  DIIS-mixed matrix last) and GFN-FF both, before anything runs, with the new error kind
  `not_applicable` (`axtb_base::Error::NotApplicable`; Python: `ProgramError` kind `unsupported`,
  code `not_applicable`).
  GFN-FF's refusals of orbitals (a Molden file) and bond orders take the same kind (they were
  `input`).
- Every SCF result reports the Fermi levels of the alpha and beta electrons (`fermi_levels`, Eh;
  `null` where the occupations are the Aufbau filling, the same level twice for a closed shell)
  and the electronic temperature the occupations used (`electronic_temperature`, K; 0 for the
  Aufbau filling, e.g. GFN2-xTB's program conventions up to 0.1 K). In Rust:
  `Electronic::fermi_levels` and `Electronic::electronic_temperature`. New keys in the CLI's JSON
  (Python: `Orbitals.fermi_level_a/_b`, `electronic_temperature`); nothing else changes.
- g-xTB 2.0.1 in the axtb conventions ends its SCF with one unmixed iteration: once the
  program's convergence test is met after iteration k, the Fock matrix of P_k's state is
  diagonalized as built (not DIIS-mixed), and the resulting density P_{k+1} gives the charges,
  moments, energy and gradient, as the g-xTB v1.1 programs end their SCFs (`Gxtb201Rule`). The
  returned orbitals are then the fresh eigenvectors of one Fock build, and the results move toward
  those of a tightly converged SCF: on water, caffeine, taxol, ferrocene, [Cu(NH3)4]²⁺, [ZnCl4]²⁻
  and (H2O)100, energies by up to 9.7e-8 Eh (taxol), gradients by up to 3.4e-5 Eh/bohr, charges by
  up to 2.4e-4, each closer to the `--acc 1e-4` result. One iteration more (9 → 10); an energy and
  gradient takes 2-4 % more time on one thread and about 6 % more on eight (taxol, (H2O)200). It
  applies wherever 2.0.1 runs in the axtb conventions (single points, Hessians and their
  reference SCFs, restarts, the solvation composite's g-xTB part). The unmixed iteration is
  axtb's, not the program's: it runs even when the program's test is met in the last iteration
  `max_iter` allows (the reported iterations count it), so that the axtb conventions converge
  where the program does. The program conventions (`gxtb`) keep the program's ending, bit for
  bit.
- Asking for an output changes no other number. In the axtb conventions GFN2-xTB's SCF
  computed every eigenvector in every iteration when the orbitals (a Molden file), a Hessian's
  response or open-shell bond orders were asked for, and the occupied ones alone otherwise; g-xTB
  computed every orbital again after its SCF. So energies, gradients and charges with and without
  a Molden file differed by round-off (GFN2-xTB), and g-xTB's gradients did. Now every SCF takes
  the same path whatever is asked for, and the orbitals without occupation are computed once
  after it, from its last Fock matrix (`GenEigSolver::solve_rest`, which replaces
  `solve_all_with`); the occupied orbitals returned are those of the SCF, bit for bit. GFN2-xTB's
  open-shell bond orders use the orbitals of the SCF's density. The program conventions, which use
  the full eigensolver throughout, are unchanged. Results that kept every orbital move by
  round-off: GFN2-xTB Hessians in the axtb conventions (energies ≤ 3e-14 Eh, gradients ≤ 2e-15
  Eh/bohr, Hessians ≤ 3e-13 Eh/bohr², frequencies ≤ 1.5e-8 cm⁻¹) and g-xTB 2.0.1 Hessians
  (≤ 5.7e-13 Eh/bohr², ≤ 4.7e-8 cm⁻¹).
- License: every crate is LGPL-3.0-or-later, `axtb-gxtb` (g-xTB v1.1 and 2.0.1 with their
  parameters) included, and so are `axtb-driver`, the `axtb` program and the Python package and
  its wheels; in 0.1.x they were GPL-3.0-only as built by default. g-xTB itself stays GPL-3.0:
  its authors, Thomas Froitzheim and Stefan Grimme, gave explicit permission to use g-xTB,
  including its parameters (1.1 and 2.0.1), in axtb and to distribute the implementation under
  LGPL-3.0-or-later. [NOTICE](NOTICE) quotes the permission; [LICENSE](LICENSE) lists the
  crates. Every published crate carries `COPYING`, `COPYING.LESSER` and `NOTICE`, and every Rust
  and Python source file starts with `SPDX-License-Identifier: LGPL-3.0-or-later`
  (`scripts/check.sh fast` checks it).
- The Cargo feature `gxtb` of `axtb-driver`, `axtb-cli` and `axtb-py` is gone: it existed only
  for an LGPL build without g-xTB, and g-xTB is now always built in. `--no-default-features` no
  longer leaves g-xTB out; a dependency that names the feature (`features = ["gxtb"]`) must
  drop it.

## 0.1.2 (2026-11-17)

- g-xTB 2.0.1 analytic Hessians: the response terms in a symmetric form
  (`AnalyticOptions::symmetric`, on by default) whose error is quadratic in the response
  solver's tolerance, and a default tolerance of 1e-6 (was 1e-10). Over 15 cases (closed and
  open shells, 0-5000 K, both conventions) the Hessians equal the 0.1.x ones to ≤ 2.2e-11
  Eh/bohr² (vibrational frequencies ≤ 2.5e-7 cm⁻¹), and they take 1.33-1.37x less time for taxol
  (1 and 8 threads) and (H2O)100 (8 threads) on the laptop (`tests/gxtb2_hessian_symmetric.py`).
  Memory: + 2 · 3N · n_occ · n_virt · 4 bytes (taxol +53 MB, (H2O)100 +0.67 GB). The reported
  `hessian_max_asymmetry` of these Hessians is now ~1e-7, removed by the symmetrization.
- Fixes found by the engineering review (no validated energy, gradient or Hessian changes):
  - g-xTB v1.1 results report the electronic entropy term `entropy` (−TS) among their terms and
    in `EnergyTerms`, as 2.0.1 does; it was missing and read 0 at a finite temperature. The
    totals were right.
  - `guess="sad"` starts GFN2-xTB SCFs from neutral atoms; it used to start them from EEQ
    charges without saying so. A named start of another method (`"eeq-bc"` for GFN2-xTB,
    `"eeq"` or `"sad"` for g-xTB, any named start for the solvation composite) is refused.
  - Unpaired electrons (`uhf` > 0) for GFN-FF, which has no electronic structure, give a
    warning instead of being ignored silently.
  - The documented default of `max_iter` is 1000 for g-xTB v1.1 (250 for GFN2-xTB and
    g-xTB 2.0.1); the docs said 250 for all.
  - The g-xTB 2.0.1 response and analytic Hessians reuse the calculator's setup of the
    geometry instead of building it again (results unchanged).
  - The D4 three-body term no longer gives NaN for a pair whose C6 is exactly zero (no element
    of the shipped tables has one).
- The Python package's `qc` extra requires `qcdata<0.24` and `qccodec<0.22`: qcdata 0.24 removes the
  `SinglePointData` fields the qccodec codec fills for Hessians (`freqs_wavenumber`,
  `normal_modes_cartesian`, `gibbs_free_energy`). axtb 0.2 fills `vibrational_analysis` and
  `thermochemistry` instead.
- Every validation gate of `scripts/check.sh programs` now fails when a result exceeds its
  documented limit (the numbers of docs/validation.md and the specs; `tests/gate.py`), or a
  calculation fails, except the known differences the docs explain. 17 of the 24 gates used to
  print their numbers and always pass, and `vibrational-*` and `solvation-composite` checked
  limits 10 to 7000 times looser than their documented results. `scripts/check.sh list` prints
  each gate's limits.

## 0.1.1 (2026-11-10)

- The minimum supported Rust is 1.89 for every crate (`rust-version`; 0.1.0 declared 1.99, which
  made Cargo refuse the crates in projects that support older compilers). CI checks the build on
  1.89. Nothing computed changes: the results and timings in this repository are those of 0.1.0.

## 0.1.0 (2026-11-07)

The first release: the Rust crates `axtb`, `axtb-d4`, `axtb-solvation`, `axtb-gfnff`,
`axtb-gxtb`, `axtb-driver` and `axtb-cli` (the `axtb` program, JSON in and out) on crates.io,
the Python package `axtb` on PyPI, and its codec for qccodec/qccompute. The README describes
the methods, conventions and settings; [docs/validation.md](docs/validation.md) has the
agreement with each program and the test sets, and [bench/results](bench/results) the result
files they come from.

### Methods, and the program each convention reproduces

Every method has axtb's default conventions (`axtb`) and the conventions of the program it
reproduces, which give that program's numbers to the last digits, quirks included.

- **GFN2-xTB** (`axtb`): energies, analytic gradients, Hessians (analytic, from the
  coupled-perturbed SCC response, by default; response and numerical), self-consistent D4
  dispersion (`axtb-d4`, bit-identical to dftd4 4.3.1), multipole electrostatics.
  Conventions `xtb`: the xtb 6.8.0 program (4b8c987); `xtb-6.7.1`: xtb 6.7.1 as the g-xTB 2.0.1
  build of xtb has it (26dd68d; the 6.7.1 release edcfbbe differs in some GFN2-xTB gradients,
  fixed upstream before that build); `tblite`: tblite 0.8.0. In `xtb` and `xtb-6.7.1` the SCF
  is the program's own SCC loop and Broyden mixer.
- **ALPB/GBSA implicit solvation** for GFN2-xTB (`axtb-solvation`): 25 ALPB and 14 GBSA
  solvents, the three reference states; xtb 6.8.0 and 6.7.1, and tblite 0.8.0 in `tblite`
  conventions.
- **g-xTB v1.1** (`axtb-gxtb`, method `gxtb-v1.1`), the public preprint version: energies,
  analytic gradients (which the program does not have), semi-numerical Hessians, finite
  electronic temperature. Conventions `gxtb`: the g-xTB 1.1.0 program (grimme-lab/g-xtb
  release v1.1.0).
- **g-xTB 2.0.1** (method `gxtb-v2.0.1`): energies, analytic gradients, analytic (default),
  response and numerical Hessians, finite electronic temperature, all 103 elements H–Lr.
  Conventions `gxtb`: the g-xTB 2.0.1 program (`xtb --gxtb` of the grimme-lab/g-xtb 2.0.1
  release). Its source is not public; [docs/gxtb-v2.0.1-spec.md](docs/gxtb-v2.0.1-spec.md)
  says how the method was reproduced from the released program.
- **g-xTB in implicit solvent**: the composite E_gxtb + [E_GFN2+ALPB/GBSA − E_GFN2] for
  energies, gradients and Hessians in one process (`CompositeCalculator`), reproducing the
  three xtb program runs it replaces.
- **GFN-FF** (`axtb-gfnff`, method `gfnff`): topology for every element H–Lr, energies,
  analytic gradients, analytic (default) and numerical Hessians with dipole derivatives.
  Conventions `xtb`: xtb 6.8.0 (`--gfnff`); `xtb-6.7.1`: xtb 6.7.1.
- **Vibrational analysis** of any Hessian: frequencies, normal modes, IR intensities and the
  quasi-RRHO thermochemistry. Conventions `xtb`: `xtb --hess` and `xtb thermo` of 6.8.0;
  `xtb-6.7.1`: those of 6.7.1.
- Also: Molden files of the converged orbitals; restarts from the previous geometry's
  wavefunction; Mayer–Wiberg bond orders, orbital energies, HOMO–LUMO gap, molecular
  quadrupole.

### Default conventions versus the program modes

The default `axtb` conventions use CODATA 2022 constants and correct what the programs get
wrong; the program modes keep it.

- GFN2-xTB: double-precision tables, one coordination number for energy and gradient, smooth
  long-range cutoffs, the solvation fixes of grimme-lab/xtb#1504 and #1505, and DIIS when
  Broyden mixing does not converge. For ordinary molecules the modes differ by about 1e-8 Eh.
- g-xTB v1.1: the exact derivative of the energy as the Fock matrix, so the SCF is stationary
  (the program's Fock matrix leaves out two potential terms).
- g-xTB 2.0.1: an integer electron count and smooth cutoffs; the same energy expression and
  SCF as the program.
- Vibrational analysis: standard atomic masses, rotations and translations projected in
  mass-weighted coordinates, always 3N−6 (3N−5) modes, the point group's exact symmetry
  number.

### Licenses

- `axtb`, `axtb-d4`, `axtb-solvation` and `axtb-gfnff`: LGPL-3.0-or-later.
- `axtb-gxtb`: GPL-3.0-only, because it embeds the g-xTB parameters (grimme-lab/g-xtb,
  GPL-3.0).
- `axtb-driver`, `axtb-cli` and the Python package build g-xTB in by default and are
  GPL-3.0-only as built. With `--no-default-features` they leave g-xTB out, and the build is
  LGPL-3.0-or-later: `cargo build -p axtb-cli --no-default-features`,
  `maturin build --no-default-features -m crates/axtb-py/Cargo.toml`.
- [LICENSE](LICENSE) explains the split; [NOTICE](NOTICE) credits the upstream work.

### Requirements

- Rust 1.89 or newer to build.
- x86-64 Linux is the developed and tested platform. The builds from this repository, the
  release binary and the wheels target x86-64-v3 (AVX2, FMA; Intel Haswell and AMD Excavator
  or newer). The CLI stops with an error and `import axtb` raises `ImportError` on a CPU
  without these instructions; `RUSTFLAGS="-C target-cpu=x86-64"` builds for any x86-64 CPU,
  with bit-identical results.
- Wheels: manylinux_2_28 (glibc 2.28 or newer), Python 3.10 or newer (abi3). The release
  binary needs glibc 2.28 or newer.

### Known limits

- No periodic systems; molecules only.
- No GFN1-xTB.
- No implicit solvation for GFN-FF.
- g-xTB has no implicit solvent of its own: only through the composite with GFN2-xTB's
  ALPB/GBSA.
- No GPU: the eigensolver runs on the CPU (faer, or a LAPACK library loaded at run time).
- No PDB residue charges.

### Program issues found and reported upstream

- [grimme-lab/xtb#1504](https://github.com/grimme-lab/xtb/issues/1504): a stale variable in
  the Born-radius integral makes ALPB/GBSA energies depend on the atom order (fixed in xtb
  6.8.0); the same in tblite,
  [tblite/tblite#430](https://github.com/tblite/tblite/issues/430).
- [grimme-lab/xtb#1505](https://github.com/grimme-lab/xtb/issues/1505): the ALPB/GBSA
  `reference` state used the wrong solvent density (fixed in xtb 6.8.0).
- [grimme-lab/xtb#1506](https://github.com/grimme-lab/xtb/issues/1506): single-precision
  literals in the D4 coordination number, and a small-gap Broyden damping that never ran
  (fixed and removed in xtb 6.8.0).
- [grimme-lab/xtb#1515](https://github.com/grimme-lab/xtb/issues/1515): GFN-FF allocates a
  hydrogen-bond array that grows as N³ in every energy evaluation (still in xtb 6.8.0).
- [grimme-lab/xtb#1516](https://github.com/grimme-lab/xtb/issues/1516): GFN-FF's carbonyl
  hydrogen-bond energy is not smooth where C=O···H is linear, and xtb 6.7.1's `--hess`
  crashes there (still in xtb 6.8.0).
- [sarah-quinones/gemm#53](https://github.com/sarah-quinones/gemm/issues/53) and
  [sarah-quinones/gemm-x64-v2#8](https://github.com/sarah-quinones/gemm-x64-v2/issues/8): the
  faer dependencies gemm-common and private-gemm-x86 open every CPU's cache files at
  start-up. This repository builds patched copies (`third_party/`); crates.io builds of the
  axtb crates do not get the patch.

docs/validation.md lists the other program behaviours found on the way, which are documented
but not reported.

## Numerics epochs

Each method's epoch (`axtb_driver::numerics`; the rule in the alpha 1 entry above) and the
commits that changed it. All start at 1 with alpha 1 of 0.2, whose round-off and model changes
before the epochs existed are listed in the entries above.

| method | epoch | since | what changed |
|---|---|---|---|
| `gfn2` | 1 | 0.2.0 alpha 1 | start |
| `gxtb-v1.1` | 1 | 0.2.0 alpha 1 | start |
| `gxtb-v2.0.1` | 1 | 0.2.0 alpha 1 | start (includes the unmixed final iteration of the axtb conventions) |
| `gfnff` | 1 | 0.2.0 alpha 1 | start |
| `analysis` | 1 | 0.2.0 alpha 1 | start (the analysis of a given Hessian) |
