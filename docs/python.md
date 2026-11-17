# axtb's Python interface

The reference of `axtb.compute`, `axtb.Session` and the program registered with qccompute
(axtb 0.2). axtb is a program of the qc* stack: a `qcdata.ProgramInput` goes in, a
`qcdata.ProgramOutput` comes out, and qccompute runs it like any other program without
importing anything of axtb's beyond the entry point. Units are atomic (bohr, hartree), as
qcdata's are.

```python
import axtb
from qcdata import ProgramInput, Structure

water = Structure(symbols=["O", "H", "H"], geometry=[[0, 0, 0.22], [0, 1.43, -0.89], [0, -1.43, -0.89]])
inp = ProgramInput(program="axtb", calctype="gradient", structure=water,
                   model={"method": "gxtb-v2.0.1"}, keywords={"threads": 4})
out = axtb.compute(inp)
out.results.energy, out.results.gradient, out.record["scf"]["iterations"]

with axtb.Session(inp) as s:                       # one molecule, many geometries
    for x in trajectory:                           # n_atoms x 3, bohr
        frame = s.compute(x)                       # warm-started from the previous geometry
        frame.results.energy
```

The public names are `compute`, `Session`, `program` and `__version__`; everything else
(`axtb._axtb`, the modules with a leading underscore) is private and changes without notice.

## compute(inp, ctx=None)

One calculation of a `ProgramInput`: calctype `energy`, `gradient` or `hessian` (the Hessian
with its vibrational analysis and thermochemistry unless `freq.analysis` is `off`). It returns
the successful `ProgramOutput`; a refused input or a failed calculation raises
`qcdata.ProgramError` whose `output` is the failed `ProgramOutput` (its `error`, `resolved`,
`warnings`, `record` and `execution` as far as the calculation got). `ctx` is a runner's
`qccompute.program.RunContext`; axtb reads only `ctx.threads`.

The calculation is the CLI's for the same settings, bit for bit: the planning turns the input
into the driver's flat request (from the settings given, never from defaults, which the
driver applies) and the driver runs it as `axtb --json` does. `tests/flat.py` maps a 0.1 flat
request onto a `ProgramInput`, and `scripts/golden_python.py` checks every case of the CLI's
golden set through `axtb.compute` bit for bit.

The input's fields:

| field | axtb |
|---|---|
| `program` | `"axtb"` |
| `calctype` | `energy`, `gradient`, `hessian`; anything else is refused (`unsupported_calctype`) |
| `structure` | symbols, geometry (bohr), `charge` (an integer), `multiplicity` (default by parity; its parity must be that of the electrons the method's orbitals hold: GFN2-xTB and g-xTB are valence-only, so cores do not count; GFN-FF refuses a multiplicity above 1, `not_modeled`) |
| `model.method` | `gfn2` (GFN2-xTB), `gxtb-v2.0.1`, `gxtb-v1.1`, `gfnff` (aliases: case and anything but letters and digits are ignored, `GFN2-xTB`, `g-xTB 2.0.1`); `model.basis` must be empty |
| `keywords` | the settings, [below](#keywords), among them `accuracy`: targets for the energy, gradient and Hessian errors ([Accuracy](#accuracy)) |
| `given` | a converged result to start the SCF from ([Given starts](#given-starts)) |
| `files`, `cmdline_args`, `subprograms`, `structures` | refused (`not_applicable`) |

## Session(inp, *, calctypes=None, ctx=None)

One molecule at many geometries on one engine, as qccompute's session protocol describes. The
input pins the atoms, charge, multiplicity, method and keywords for the session's lifetime; its
geometry and calctype are the defaults of a request. `calctypes` declares what the session may
run (each with its lower derivatives; None: all three); `hessian.*` and `freq.*` keywords need
`hessian` among them. The engine builds the method once (GFN-FF its topology too) and keeps one
committed wavefunction per lineage. A session is a context manager; `close()` drops the engine
and every lineage's history.

`Session.compute(positions=None, calctype=None, *, accuracy=None, lineage=None, outputs=(),
seed=(), commit=True)` evaluates one geometry and returns a `qcdata.Frame` (`structure`,
`results`, `warnings`, `record`, `execution`):

- `positions`: n_atoms x 3 in bohr (None: the input's geometry); `calctype`: one the session
  allows (None: the input's).
- `lineage`: any JSON value naming the trajectory the request continues (None: the default
  lineage). A request starts from its lineage's committed wavefunction, else from the first of
  its `seed` lineages that has one, else fresh, as the method's restart policy says
  (`scf.restart`). Its record says `warm_start` and `from_lineage`.
- `commit`: whether a success becomes the lineage's history (`commit=False` for probes, line
  searches and finite differences that must not move the trajectory's start).
- `accuracy`: this request's targets (a bare number targets the calctype's highest derivative);
  `outputs`: optional outputs of this request (`orbitals`, `overlap`, `hamiltonian`,
  `bond_orders`; `ir` with a Hessian), beside the keywords' `output.*`.
- A refused or failed request raises `ProgramError` with the failed frame as `exc.frame`; the
  session stays usable and the lineage's history is untouched (a panic closes the session).

`compute_many(requests)` evaluates `qcdata.Request`s, different lineages side by side:

- **Order.** Each lineage's requests run in request order, each reading its lineage's history as
  committed by the lineage's requests before it. A request whose lineage has no history waits for
  its `seed` lineages' earlier requests in the batch and reads their history as those leave it,
  so no result depends on the scheduling.
- **Side by side.** Each request runs on its own threads, those of its calctype and size
  ([Threads](#threads)), and different lineages run at once on forks of the session's calculator
  as long as their threads fit into the session's threads in all, taken in request order. A
  session of 8 threads runs 8 caffeine g-xTB 2.0.1 gradients at once (1 thread each), and a
  caffeine Hessian (8 threads) alone. `plan()["record"]["lineages_at_once"]` says how many
  requests of the input's calctype run at once. The method's setup is built once and shared
  where it is immutable (parameters, basis, GFN-FF's topology); each fork runs on one thread
  count and has its own SCF state and work buffers. A lineage's history does not depend on the
  threads, so its requests may run on different forks. The forks are made when a call first
  needs them and kept until the session closes.
- **Memory.** One calculator per lineage that runs at once, each forked one about the working
  memory of one calculation: per fork 1.9 MiB for a caffeine g-xTB 2.0.1 gradient, 24 MiB for
  taxol's (the session's own calculator takes 1.8 and 26 MiB). Measured on the laptop with 12
  lineages per call in a session of 8 threads, 8 at once on 1 thread each: caffeine 31 ms per
  call against 160 ms one request at a time (5.1x) and 93 ms one at a time on 8 threads each;
  taxol 0.73 s against 3.6 s (4.9x).
- **Results.** The frames are what one-at-a-time `compute` calls in request order give, bit for
  bit (`crates/axtb-py/tests/test_session.py`, with interleaved and seeded lineages and mixed
  calctypes): every request runs on its calctype's threads, side by side or not. They are not
  bit for bit what the same calls on all of the session's threads would give: a calculation's
  round-off depends on its thread count.
- **Failures** stay in place and raise `qcdata.ProgramErrors` at the end, with every frame in
  request order.

A single `compute` runs on its calctype's threads too, serialized with the session's other
calls by its lock; only `compute_many` runs lineages side by side. `reset(lineage=...)` forgets
every lineage's history (no argument) or one lineage's; `lineages()` lists those with history;
`plan()` returns the settings pinned at opening (the threads of the input's calctype), the
resources (the most threads one calculation takes) and the record's static part, with
`threads_by_calctype`.

Interleaved lineages give what separate sequential sessions give, bit for bit, and a session's
frames equal the CLI's trajectory cases of the golden set bit for bit. A session's methods are
thread-safe (a lock serializes them; the GIL is released while the engine runs). A session step
costs about 30 µs of Python on top of the calculation, `compute` about 100 µs (water).

**How far a warm start moves a result.** A frame that starts from its lineage's history converges
to the same SCF solution as a fresh `compute` at that geometry, within the SCF's thresholds, so
the two differ by up to the accuracy model's bound at the session's factor
([Accuracy](#accuracy)). With `scf.restart: off` no frame warm-starts, and every frame equals
`compute` bit for bit. That is the default for g-xTB in its programs' conventions (`gxtb-2.0.1`,
`gxtb-1.1.0`), as the g-xTB programs start every SCF fresh. Measured: a 6-geometry path of
caffeine, taxol and [Cu(NH3)4]²⁺ (steps of 0.01 bohr per atom), the largest difference of any
frame from a fresh `compute`:

| method, conventions | `accuracy_scale` 1: energy (Eh), gradient (Eh/bohr) | 0.01: energy, gradient |
|---|---|---|
| g-xTB 2.0.1, `axtb` | 3.4e-8, 4.4e-5 | 7.5e-11, 1.6e-6 |
| g-xTB 2.0.1, `gxtb-2.0.1` | 0, 0 (no warm start) | 0, 0 |
| g-xTB v1.1, `axtb` | 1.3e-8, 1.0e-5 | 2.9e-8, 5.2e-6 ([Cu(NH3)4]²⁺; closed shells ≤ 5.3e-7) |
| g-xTB v1.1, `gxtb-1.1.0` | 0, 0 (no warm start) | 0, 0 |
| GFN2-xTB, `axtb`, `xtb-6.8.0`, `tblite-0.8.0` | 8.5e-9, 1.2e-5 | 1.0e-12, 9.2e-8 |

A session can be forked: a child process (`os.fork`) opens sessions and computes whether or not
the parent computed before the fork. axtb starts no thread that outlives a call; threads a
process has at a fork are other libraries' (numpy's OpenBLAS pool, for one)
(`crates/axtb-py/tests/test_fork.py`).

## Given starts

`ProgramInput.given` takes a converged result (`qcdata.SinglePointData`, e.g. the `results` of
an earlier output of the same molecule) as the SCF's start, in place of the method's own guess:

- g-xTB (both versions, and the solvation composite's g-xTB part) starts from the density
  P = C n Cᵀ of the given orbitals (`wavefunction.orbitals` with their coefficients and basis
  set: ask for `output.orbitals`), per spin channel;
- GFN2-xTB starts from the given shell populations and atomic multipoles
  (`shell_populations` with `shell_atoms` and `shell_angular_momenta`, `atomic_dipoles`,
  `atomic_quadrupoles`), which every SCF result of axtb has (another method's on the same
  shells is a guess too).

The used fields must fit: the method's shells on these atoms, its basis set (centers, angular
momenta, exponents), the spin channels of the molecule; the other fields of `given` are left
alone. A misfit, or a given result with nothing usable, is refused before the calculation
(`not_applicable`). `resolved.settings.scf.guess` is `given` and the record says
`scf.start = "given"`; in a session every request without history starts from it. A converged
frame given back reproduces its result within the accuracy estimate in at most 3 iterations
(measured: 2-3, from 9-12).

## SCF retries

`scf.retries` (default none) is a ladder of rungs tried in order after a calculation whose SCF
fails: it does not converge within `scf.max_iter`, or breaks down (a failed mixer, a lost Fermi
level). Each rung `{"etemp": K, "max_iter": n}` runs the calculation again as a fresh SCF from the
method's own start with those settings; a missing key keeps the calculation's value. It applies
to GFN2-xTB and g-xTB (both versions, every conventions) and the solvation composite; GFN-FF has
no SCF (`no_effect`). mepd's ladders, for example:

```python
keywords = {"scf": {"retries": [{"etemp": 1000, "max_iter": 1000}, {"etemp": 3000, "max_iter": 1000}]}}
# GFN2-xTB in solvent, the first try at the defaults
keywords = {"solvent": {"name": "water"},
            "scf": {"retries": [{"etemp": 1000, "max_iter": 1000}, {"max_iter": 1000},
                                {"etemp": 5000, "max_iter": 1000}]}}
```

- **Inside each try** the method does what it does without a ladder: GFN2-xTB's Broyden-then-DIIS
  fallback in the axtb conventions, and in a session an SCF from history that does not converge
  runs once more fresh (`record.scf.restart_fallback`). The rungs come after that.
- **A ladder that does not fire changes no number**, and a result from a rung is the calculation
  at that rung's settings, bit for bit (both are a fresh SCF with the same settings). The
  composite's single point after its Hessian is the exception: it starts its GFN2 parts fresh at
  the rung, not from their Hessians' references, so it agrees within the SCF thresholds (its
  Hessian is bit for bit).
- **The floor.** With `scf.etemp` the ladder starts at that temperature, and a rung below it is
  refused (`no_effect`, path `scf.retries[i].etemp`), so a path never drops below its floor.
- **Hessians** run whole at one rung. When a try fails anywhere (the reference SCF, a displaced
  SCF of a numerical Hessian, the response of an analytic one), the whole Hessian and the single
  point of its output run again at the next rung. Every SCF of a result is at the rung that
  converged, and no displacement climbs alone. A Hessian whose reference converges but one of
  whose displaced or perturbed SCFs fails has no result at the reference's settings; rerunning it
  whole at the next rung gives one Hessian at one temperature instead of an error to retry by
  hand. The composite's own retries already climb its parts' Hessians whole.
- **Sessions.** A request starts from its lineage's history as usual; when that fails, a fresh SCF
  at the calculation's own settings, then the rungs. A frame from a rung commits normally: its
  state is the next frame's history, and the next frame tries the calculation's own settings
  first. To keep a lineage at a rung's temperature, set `scf.etemp` on its later inputs (it is
  in the key, so those results are stored and reused like any other).
- **The composite.** A given ladder replaces the built-in retries (the g-xTB part's at 1000 and
  3000 K, the GFN2 pair's three rungs) for every part: the g-xTB part and the GFN2 pair each climb
  the rungs in order after their first try (a missing key: the part's own value). The single point
  after a composite Hessian starts each part at the rung its Hessian converged on.
  `record.composite.definition` shows `retries: false` and the ladder (`scf_retries`), and each
  part its `rung` and `retry` label. Without `scf.retries` the composite is unchanged.
- **The record.** `resolved.settings.scf.retries` echoes the ladder (it is part of the input's
  key). `record.contingent["scf.retries"]` says whether a rung gave the result and which (`rung`:
  0 the calculation's own settings, i the i-th rung; the composite: the highest of its parts' and
  the `parts` that climbed). `record.scf.attempts` lists every try: `retry` (its rung), `etemp`,
  `max_iter`, `converged`, `iterations`; a converged try its `start` (GFN2-xTB: one entry per mixer,
  with `mixer`), a failed one its `error` code (`scf_not_converged`, `numerical_failure`); the
  composite's tries carry their `part`, and `hessian` for its Hessian's. Every rung starts from the
  method's own start. A result from a rung has the warning `scf_retry`, e.g. "converged on retry 2
  of 2: etemp 3000 K, max_iter 1000" (the composite one per part: "the gxtb part converged on
  ...").
- **When every rung fails**, the error is the last rung's: kind `convergence`, code
  `scf_not_converged`, `details` with the last rung's `iterations`, `energy_change` and
  `residual`, and the failed output's (or frame's) `record.scf.attempts` lists every try. A last
  rung whose SCF broke down raises `numerical_failure` instead.

The CLI takes the ladder as `--scf-retry etemp=1000,max_iter=1000`, repeated, and the JSON
request as `scf_retries`.

## Keywords

`ProgramInput.keywords` hold axtb's settings, grouped as qcdata's conventions say
(`group.key`). Shared keys carry qcdata's names; the others are axtb's own. Names follow the
alias rule (case and anything but letters and digits are ignored) and `resolved.settings`
echoes the canonical spellings. Methods: G2 GFN2-xTB, X1 g-xTB v1.1, X2 g-xTB 2.0.1, XC the
g-xTB solvation composite, FF GFN-FF. Calctypes: E energy, G gradient, H hessian.

The checks are strict and run in qcdata's order before anything is computed; the first error
raises (kind `input`):

1. `unknown_key`: not in the table, with a did-you-mean from axtb 0.1's names and close
   spellings;
2. `invalid_value`: a wrong type, range or shape, or an unknown name;
3. `no_effect`: a setting the method, the conventions, another setting or the calculation
   makes irrelevant (an SCF setting for GFN-FF, a solvent's reference state without a solvent,
   a Hessian step for an analytic Hessian);
4. `not_applicable`: a setting of a calctype that will not run;
5. `missing_setting`: a value that needs another key (g-xTB in a solvent needs
   `solvent.correction`);
6. the structure: `invalid_structure`, `unsupported_element`, `not_modeled`;
7. kind `unsupported`: a shared setting axtb cannot honour (`unsupported_setting`), an output
   it cannot give (`not_applicable`).

`axtb.program.settings_schema(calctype)` gives the keywords of a calctype as a JSON schema.

<!-- keywords: generated by scripts/gen_python_docs.py -->

| key | tier | unit | default | methods | calctypes | description |
|---|---|---|---|---|---|---|
| `threads` | shared |  | the runner's grant (a plan without one: deferred), else the plan's | G2, X1, X2, XC, FF | EGH | CPU threads of the calculation, at most the grant (more: the grant and the warning resource_capped). |
| `accuracy` | shared | Eh, Eh/bohr, Eh/bohr² | none | G2, X1, X2, XC, FF | EGH | Accuracy targets (energy, gradient, hessian), turned into the SCF accuracy factor by axtb's error model and recorded under record.accuracy; GFN-FF meets every target (no thresholds). |
| `accuracy_scale` | axtb |  | 1 | G2, X1, X2, XC | EGH | Multiplier of every SCF threshold, and of GFN2-xTB's integral screening (xtb's --acc); a floor under the accuracy targets. |
| `conventions` | axtb |  | axtb | G2, X1, X2, XC, FF | EGH | Whose numbers to reproduce: axtb (CODATA 2022, fixes), or a reference program: GFN2-xTB xtb-6.8.0, xtb-6.7.1, tblite-0.8.0; g-xTB gxtb-1.1.0, gxtb-2.0.1; GFN-FF xtb-6.8.0, xtb-6.7.1 (aliases xtb, tblite, gxtb). |
| `scf.max_iter` | shared |  | GFN2-xTB, g-xTB 2.0.1 250; v1.1 1000 | G2, X1, X2, XC | EGH | Iteration limit of one SCF. |
| `scf.etemp` | axtb | K | GFN2-xTB 300; g-xTB 0; the composite each part's (g-xTB 0, GFN2 300; not echoed) | G2, X1, X2, XC | EGH | Electronic temperature of the Fermi smearing; given for the composite, it sets every part (and turns off the g-xTB part's retries at 1000 and 3000 K). |
| `scf.retries` | axtb |  | none | G2, X1, X2, XC | EGH | SCF retry ladder: after a calculation whose SCF fails, a fresh SCF at each rung in order until one converges; a rung is {"etemp": K, "max_iter": n}, either key optional (missing: the calculation's). No rung below a given scf.etemp. A Hessian runs whole at one rung. A result from a rung has the warning scf_retry (contingent). The composite: replaces its built-in retries, for every part. |
| `scf.mixer` | axtb |  | axtb conventions broyden-then-diis; the programs' broyden | G2, XC | EGH | Mixer of GFN2-xTB's SCF (and of the composite's GFN2 parts): broyden-then-diis falls back to DIIS when Broyden fails (contingent). |
| `scf.guess` | shared |  | GFN2-xTB eeq (tblite-0.8.0 sad); g-xTB eeq-bc | G2, X1, X2 | EGH | Where an SCF without history starts: GFN2-xTB eeq (charges) or sad (neutral atoms), g-xTB eeq-bc; given: from the input's given wavefunction. |
| `scf.restart` | axtb |  | previous; g-xTB in its program's conventions off | G2, X1, X2, XC | EGH | Sessions: whether an SCF starts from its lineage's history (previous) or fresh (off). |
| `scf.e_tol` | shared | | | refused: kind `unsupported` | | axtb's SCF thresholds follow accuracy_scale and accuracy, not explicit tolerances |
| `scf.g_tol` | shared | | | refused: kind `unsupported` | | axtb's SCF thresholds follow accuracy_scale and accuracy, not explicit tolerances |
| `solvent.name` | shared |  | none | G2, XC | EGH | Implicit solvent (ALPB or GBSA parameters of GFN2-xTB). g-xTB has none: with solvent.correction gfn2 it runs the composite E(g-xTB) + [E(GFN2-xTB + solvent) - E(GFN2-xTB)]. |
| `solvent.model` | shared |  | alpb | G2, XC | EGH | The continuum model of solvent.name. |
| `solvent.eps` | shared | | | refused: kind `unsupported` | | axtb's ALPB and GBSA models are parametrized per solvent name (solvent.name), not by a permittivity |
| `solvent.state` | axtb |  | gsolv | G2, XC | EGH | Reference state of the solvation free energy. |
| `solvent.correction` | axtb |  | none | XC | EGH | g-xTB in solvent: the composite with GFN2-xTB's solvation correction. |
| `solvent.correction_conventions` | axtb |  | axtb; xtb-6.7.1 with g-xTB's program conventions | XC | EGH | Conventions of the composite's GFN2-xTB parts. |
| `topology.geometry` | axtb | bohr | the input's geometry (a session: the opening one; not echoed) | FF | EGH | The positions GFN-FF sets its topology up from. |
| `topology.fragment_charges` | axtb | e | the charge on atom 1's fragment (in record.gfnff; not echoed) | FF | EGH | Charges of the molecular fragments, in the order of their first atoms, summing to the charge. |
| `hessian.method` | axtb |  | analytic; g-xTB v1.1 numerical; the composite each part's | G2, X1, X2, XC, FF | H | analytic second derivatives, response (central differences of the response gradient) or numerical (central differences of analytic gradients). g-xTB v1.1: numerical, response (axtb conventions); GFN-FF: analytic, numerical. |
| `hessian.step` | axtb | bohr | numerical 0.005; response GFN2-xTB 5e-4, g-xTB 2e-4 | G2, X1, X2, XC, FF | H | Finite-difference step of the numerical and response Hessians. |
| `hessian.accuracy_scale` | axtb |  | numerical 0.01; analytic and response GFN2-xTB 1e-3, g-xTB 1e-5 | G2, X1, X2, XC | H | SCF accuracy factor of the numerical Hessian's displacements, or of the analytic and response Hessians' reference; never looser than the calculation's. |
| `freq.analysis` | axtb |  | harmonic | G2, X1, X2, XC, FF | H | The harmonic analysis of the Hessian (vibrational_analysis), and its thermochemistry. |
| `freq.thermo` | axtb |  | qrrho | G2, X1, X2, XC, FF | H | Thermochemistry: Grimme's quasi-RRHO (xtb's), or none. |
| `freq.conventions` | axtb |  | the calculation's (axtb; xtb-6.8.0 or xtb-6.7.1 in those; g-xTB's programs xtb-6.7.1) | G2, X1, X2, XC, FF | H | Constants, masses and projection of the analysis. |
| `freq.masses` | axtb | u | the analysis conventions' table (not echoed) | G2, X1, X2, XC, FF | H | Atomic masses of the analysis. |
| `freq.temperature` | shared | K | 298.15 | G2, X1, X2, XC, FF | H | Temperature, or several: one thermochemistry entry each, in order. |
| `freq.pressure` | shared | Pa | 101325 | G2, X1, X2, XC, FF | H | Pressure of the ideal gas. |
| `freq.scale` | shared |  | 1 (GFN-FF in the xtb analysis conventions 1.03) | G2, X1, X2, XC, FF | H | Scale factor of the frequencies in the thermochemistry. |
| `freq.imaginary_cutoff` | shared | cm⁻¹ | 20 | G2, X1, X2, XC, FF | H | Imaginary modes below this magnitude count as real in the thermochemistry (xtb's ithr). |
| `freq.rotor_cutoff` | axtb | cm⁻¹ | 50 | G2, X1, X2, XC, FF | H | The quasi-RRHO rotor cutoff (xtb's sthr). |
| `freq.symmetry_number` | axtb |  | from the point group (not echoed) | G2, X1, X2, XC, FF | H | Rotational symmetry number. |
| `freq.symmetry_tol` | axtb | bohr | 0.1 | G2, X1, X2, XC, FF | H | Tolerance of the point-group detection (xtb's desy). |
| `output.orbitals` | shared |  | false | G2, X1, X2, XC | EGH | The orbitals' coefficients with their basis set (wavefunction.orbitals). |
| `output.overlap` | shared |  | false | G2, X1, X2, XC | EGH | The AO overlap matrix (wavefunction.overlap), with the basis set. |
| `output.hamiltonian` | shared |  | false | G2, X1, X2, XC | EGH | The Hamiltonian whose eigenvectors are the orbitals (wavefunction.hamiltonian_a/_b), with the overlap and the orbitals; not g-xTB 2.0.1 in its program's conventions. |
| `output.core_hamiltonian` | shared | | | refused: kind `unsupported` | | tight binding's H0 is model data, not a core Hamiltonian |
| `output.bond_orders` | shared |  | false | G2, X1, X2, XC | EGH | Mayer bond orders (bond_orders). |
| `output.molden` | shared | | | refused: kind `unsupported` | | request output.orbitals and write the Molden file with qcinf |
| `output.density` | axtb | | | refused: kind `unsupported` | | the SCF density is Orbitals.density() of output.orbitals |
| `memory_gib` | shared | | | refused: kind `unsupported` | | axtb has no memory setting |
| `backend` | shared | | | refused: kind `unsupported` | | axtb runs on the CPU only |
| `device_memory_gib` | shared | | | refused: kind `unsupported` | | axtb runs on the CPU only |

<!-- end keywords -->

## Outputs

`results` is a `qcdata.SinglePointData`:

| field | axtb |
|---|---|
| `energy`, `gradient`, `hessian` | Eh, Eh/bohr (n_atoms x 3), Eh/bohr² (3N x 3N), as the calctype asks; the Hessian is the full Cartesian second derivative, symmetrized and never projected: translations and rotations are projected out only in `vibrational_analysis` |
| `energy_components` | the energy's pieces, which add up to the energy (GFN2-xTB: `repulsion`, `dispersion`, `isotropic_es`, `third_order`, `anisotropic`, `entropy`, `solvation`, `band`; g-xTB's and GFN-FF's own terms; the composite's `solvation_correction`); the subtotal `electronic`, which overlaps them, is in `record.energy_subtotals` |
| `dipole_moment`, `quadrupole_moment` | e·bohr; e·bohr², the traceless (Buckingham) tensor about the origin of the input's coordinates, in every conventions |
| `atomic_charges`, `atomic_charges_kind` | Mulliken charges (`mulliken`); GFN-FF's EEQ charges (`eeq`) |
| `shell_populations`, `shell_atoms`, `shell_angular_momenta` | Mulliken shell populations Σ(PS) of the total density, with each shell's atom and angular momentum (SCF methods) |
| `atomic_dipoles`, `atomic_quadrupoles` | the SCF's atomic multipoles (GFN2-xTB's moments, g-xTB's), e·bohr and e·bohr² (traceless) |
| `bond_orders`, `bond_orders_kind` | with `output.bond_orders`: Mayer-Wiberg bond orders (n_atoms x n_atoms, xtb's `wbo`), `mayer` |
| `scf_iterations`, `calcinfo_*` | the SCF iterations; atoms, basis functions, orbitals, alpha and beta electrons |
| `wavefunction` | the orbital levels always (energies, occupations per spin in 0..1, the electronic temperature and Fermi levels, `core_electrons` per atom outside the valence-only orbitals); with `output.orbitals` the coefficients and the basis set; with `output.overlap` the AO overlap S; with `output.hamiltonian` the matrix whose generalized eigenvectors with S are the orbitals (with the overlap and orbitals). GFN2-xTB is restricted (one set of orbitals, open shells included), g-xTB restricted for a closed shell and unrestricted for an open one; the composite's are its g-xTB part's |
| `vibrational_analysis` | hessian: masses (u), frequencies (cm⁻¹, imaginary negative), Cartesian normal modes, reduced masses, force constants, IR intensities and dipole derivatives (SCF methods and GFN-FF), the rotor type, rotational constants, point group and symmetry number |
| `thermochemistry` | hessian: one entry per `freq.temperature`, in order: Grimme's quasi-RRHO entropy (`qrrho_grimme`, rotor cutoff `freq.rotor_cutoff`, α = 4, B_av) and RRHO enthalpy; each term's enthalpy, entropy and heat capacity, the total's Gibbs energy (total = E + H; the vibrational term includes the zero-point energy) |
| `provenance` | `axtb`, its version (PEP 440) and, from a git build, the commit |

The basis set of `wavefunction` is axtb's: one spherical shell per tight-binding shell, atom by
atom, m = -l..l for every l (`AOConvention(spherical_order="increasing_m")`; p as y, z, x).
Smeared occupations of a spin without electrons have no Fermi level, which qcdata's `Orbitals`
need: such a wavefunction is left out with the warning `orbitals_omitted`.

`resolved.settings` echoes every setting with its canonical spelling and the defaults that
applied, `resolved.reasons` why (`threads`: the keyword, the grant, or the plan).
`execution.timings` holds the phases of the single point, the Hessian's (`hessian:...`) and
`total` (s); `execution.threads` the threads the calculation used.

### record

`record` is axtb's detail, diagnostic and not stable across versions
(`axtb.program.record_schema(calctype)`):

| key | content |
|---|---|
| `method`, `conventions`, `numerics` | the canonical method and conventions, the numerics identity ([Numerics](#numerics)) |
| `eigensolver` | `faer`, or the run-time LAPACK library (`AXTB_LAPACK`) |
| `scf` | `start` (`eeq`, `sad`, `eeq-bc`, `previous`, `given`), `iterations`, `mixer_used` and `attempts` (GFN2-xTB's mixers; with `scf.retries` every try of the ladder, [SCF retries](#scf-retries), in a failed output too), `restart_fallback` (a session's SCF from history failed and ran again fresh) |
| `energy_subtotals` | terms that overlap the energy components (the electronic energy) |
| `contingent` | contingent settings and whether each fired (`scf.mixer`, `scf.retries` with its `rung`, `composite.retries`, `scf.restart`, a session's `hessian.*` and `freq.*`) |
| `accuracy` | with a target: the target, the thresholds used and the model's estimate at them |
| `hessian` | the Hessian's method and largest asymmetry before symmetrization |
| `composite` | the solvation composite: its definition, the energies of its parts and their SCFs (`wavefunction.part` names the part the orbitals are) |
| `gfnff` | GFN-FF: hydrogen and halogen bonds, the fragment charges, and (compute only) the topology |
| `warm_start`, `from_lineage` | session frames: whether the SCF started from history, and whose |
| `error` | a failure's fields: `z` of an unsupported element; `iterations`, `energy_change`, `residual` of an SCF that did not converge |

### Warnings

`warnings` are `qcdata.ProgramWarning`s with a code: `resource_capped` (the `threads` keyword
above the grant), `accuracy_limited` (a target below what axtb's thresholds reach), `scf_retry`
(a result from a rung of `scf.retries`, naming the rung's settings),
`orbitals_omitted` (no Fermi level for qcdata's `Orbitals`), `eigensolver_fallback`
(`AXTB_LAPACK` was set, but faer solved). Settings problems are errors, not warnings, and axtb
writes no text in process (`logs` is None).

## Errors

Every refusal and failure raises `qcdata.ProgramError`; `exc.error` is a `qcdata.ErrorInfo`
with axtb's `kind` and `code`, and `exc.output` (compute) or `exc.frame` (a session) the failed
output. The input's checks come first ([Keywords](#keywords)); the calculation's errors map from
the driver's:

| kind | code | when |
|---|---|---|
| `input` | `unknown_key`, `invalid_value`, `no_effect`, `not_applicable`, `missing_setting` | the keywords, `given` and the other fields |
| `input` | `invalid_structure` | non-finite positions, a wrong shape, a multiplicity of the wrong parity |
| `input` | `unsupported_element` | an element without parameters (`details.z`) |
| `input` | `not_modeled` | GFN-FF with a multiplicity above 1 |
| `input` | `topology_failed` | a structure GFN-FF cannot set up |
| `input` | `session_closed` | a request to a closed session |
| `unsupported` | `unsupported_calctype`, `unsupported_setting`, `not_applicable` | a calctype, shared setting or output axtb does not have |
| `convergence` | `scf_not_converged` | the SCF did not converge (`details`: `iterations`, `energy_change`, `residual`; with `scf.retries` every rung failed, the details are the last rung's and `record.scf.attempts` lists every try) |
| `numerical` | `numerical_failure` | a numerical failure (a singular matrix, a failed eigensolver) |
| `internal` | `panic` | a bug: please report it (a session closes) |

## Threads

A calculation runs on exactly its threads and never spawns more. Its threads are:

- **`compute` under a runner:** the `threads` keyword, at most the grant (above it: the grant
  and the warning `resource_capped`), or else the grant.
- **A direct `compute`:** the keyword, or else the threads the plan states.
- **A session's requests** (`compute` and `compute_many` alike): the keyword, at most the grant
  (`resource_capped` above it), or else the plan's threads for the request's calctype and the
  molecule's size, at most the grant. `compute_many` runs lineages side by side as long as
  their threads fit into the session's threads in all: the grant, and without one (a direct
  call) this process's CPUs (`os.sched_getaffinity`). So `axtb.Session(inp).compute_many(...)`
  uses the whole machine when it has many lineages to run; a program that runs several sessions
  at once passes each a `ctx` with its share of the threads (or `threads` to qccompute's
  `Runner.session`), or they oversubscribe the CPUs.

A molecule at or below the size where threads only cost runs on one thread whatever it asked
for, with the same bits (GFN2-xTB single points up to 10 atoms, Hessians up to 6; g-xTB 2.0.1
single points up to 9, v1.1 up to 6; GFN-FF single points up to 100 and Hessians up to 10).

The plan (`axtb.program.plan(inp)` gives it as `resources.threads`; a session's `plan()` per
calctype) states a thread count for every calculation, never the whole grant: the
largest count whose parallel efficiency T1/(n Tn) stayed at least 0.5 when measured, by method,
calctype and size (atoms). The bands, from `axtb_driver::THREAD_BANDS`:

<!-- threads: generated by scripts/gen_python_docs.py -->

| method | calculation | threads by size (atoms) |
|---|---|---|
| GFN2-xTB | energy, gradient | 1 (below 100), 2 (100-599), 4 (from 600) |
| g-xTB 2.0.1 | energy, gradient | 1 (below 100), 2 (100-599), 4 (from 600) |
| g-xTB v1.1 | energy, gradient | 1 (below 100), 2 (100-299), 4 (from 300) |
| GFN-FF | energy, gradient | 1 (below 600), 2 (from 600) |
| GFN2-xTB | Hessian | 1 (below 60), 4 (60-149), 8 (from 150) |
| g-xTB 2.0.1 | Hessian | 2 (below 10), 4 (10-19), 8 (20-59), 16 (from 60) |
| g-xTB v1.1 | Hessian | 4 (below 10), 16 (from 10) |
| GFN-FF | Hessian | 1 (below 100), 2 (from 100) |

<!-- end threads -->

Measured with the release CLI 0.2.0a1 on an AMD Threadripper PRO 9985WX node, 1 to 16 threads,
on water, an SN2 complex, mb33, beni-acetone, caffeine, gh106, taxol and (H2O)100 (3 to 300
atoms), the whole process timed; energies and gradients above 300 atoms on an Intel i9-11900H
laptop ((H2O)200; GFN-FF: crambin, (H2O)333). Energies and gradients gain little from threads
below about 100 atoms: a caffeine g-xTB 2.0.1 gradient takes 13.2 ms on 1 thread and 11.9 ms on
8; taxol's 218 ms on 1 and 156 ms on 2. Hessians parallelize over their displacements and
perturbations: a taxol g-xTB 2.0.1 Hessian takes 36 s on 1 thread and 3.6 s on 16. A calculation
that should use more threads asks with the `threads` keyword.

## Accuracy

The `accuracy` keyword (a `qcdata.Accuracy`, or its fields) sets error targets for the energy (Eh), the
gradient's largest component (Eh/bohr) and the Hessian (Eh/bohr²); a session's request may set
its own. axtb's SCF thresholds all scale with one factor (`accuracy_scale`), and an error model
fitted to measurements turns the targets into that factor: per method, conventions kind (axtb
or the reference program's) and shell (closed or open), the upper envelope err ≤ k acc^p over
molecules and fresh and warm starts. A target's factor is the largest measured factor (10, 3, 1,
0.3, 0.1, 0.03, 0.01, 0.001) whose bound, times a safety factor of 2, meets every target in every
model that applies: the composite its own fit and its two methods', an open shell its fit and the
closed shells'. Between measured factors the errors plateau, so an estimate rests on a
measurement at its factor or above. An explicit `accuracy_scale` is a floor under the factor.
g-xTB v1.1's open shells, whose errors do not fall monotonically, use a table of the worst error
at each measured factor. A Hessian target sets `hessian.accuracy_scale` the same way (factors
0.1 to 1e-6), never looser than the method's default or an explicit value. GFN-FF has no
thresholds: it meets every target (estimate 0).

`record.accuracy` holds the target, the thresholds used and the model's estimate at them: a
bound fitted to measurements, not a measurement of the run. A target below what the floor
reaches runs at the floor and warns `accuracy_limited`. The constants come from
`tests/accuracy_calibration.py` (data in `bench/results/accuracy.json`):

<!-- accuracy: generated by scripts/gen_python_docs.py -->

Envelopes err <= k acc^p of the energy (Eh) and the gradient's largest component (Eh/bohr), per method, conventions kind (`axtb`, or `program`: the reference program's) and shell:

| fit | energy k | energy p | gradient k | gradient p |
|---|---|---|---|---|
| composite axtb closed | 2.6e-07 | 1.25 | 4.32e-05 | 0.75 |
| composite axtb open | 6.51e-08 | 1.25 | 1.97e-05 | 0.75 |
| gfn2 axtb closed | 3.04e-08 | 1.75 | 1.36e-05 | 1 |
| gfn2 axtb open | 1.27e-09 | 2 | 6.87e-06 | 0.75 |
| gfn2 program closed | 3.04e-08 | 1.75 | 9.45e-06 | 1 |
| gfn2 program open | 1.28e-09 | 2 | 4.76e-06 | 1 |
| gxtb-v1.1 axtb closed | 1.59e-07 | 0.75 | 1.12e-05 | 0.5 |
| gxtb-v1.1 program closed | 8.27e-07 | 0.75 | 0.000109 | 0.5 |
| gxtb-v2.0.1 axtb closed | 1.21e-07 | 1.25 | 4.38e-05 | 0.75 |
| gxtb-v2.0.1 axtb open | 6.34e-08 | 1.25 | 2.15e-05 | 0.75 |
| gxtb-v2.0.1 program closed | 1.35e-07 | 1.25 | 7.8e-05 | 0.75 |
| gxtb-v2.0.1 program open | 2.19e-07 | 1.25 | 3.47e-05 | 1 |

Tables (no envelope fits): the largest error measured at each factor.

| acc | gxtb-v1.1 axtb open energy | gxtb-v1.1 axtb open gradient | gxtb-v1.1 program open energy | gxtb-v1.1 program open gradient |
|---|---|---|---|---|
| 10 | 8.28e-07 | 0.00043 | 2.9e-07 | 7.34e-05 |
| 3 | 4.09e-08 | 2.01e-05 | 5.64e-08 | 1.69e-05 |
| 1 | 3.04e-08 | 1.21e-05 | 2.72e-08 | 4.09e-06 |
| 0.3 | 2.43e-08 | 9.83e-06 | 2.17e-08 | 3.66e-06 |
| 0.1 | 8.46e-09 | 1.81e-06 | 1.12e-08 | 1.83e-06 |
| 0.03 | 8.46e-09 | 1.81e-06 | 1.12e-08 | 1.83e-06 |
| 0.01 | 8.43e-09 | 1.72e-06 | 7.19e-09 | 1.59e-06 |
| 0.001 | 1.78e-10 | 3.04e-07 | 5.64e-09 | 7.99e-07 |

Hessians: err <= k h^p (Eh/bohr², the largest element) at `hessian.accuracy_scale` h:

| method and Hessian method | k | p |
|---|---|---|
| gfn2 analytic | 0.001 | 1 |
| gfn2 response | 0.001 | 1 |
| gxtb-v1.1 numerical | 0.00301 | 0.5 |
| gxtb-v2.0.1 analytic | 2.24e-05 | 0.5 |
| gxtb-v2.0.1 response | 2.24e-05 | 0.5 |

Factors measured: 10, 3, 1, 0.3, 0.1, 0.03, 0.01, 0.001 (SCF), 0.1, 0.01, 0.001, 0.0001, 1e-05, 1e-06 (Hessian). Safety factor 2, factor range 0.001 to 10. Data: `bench/results/accuracy.json`.

<!-- end accuracy -->

## Numerics

Every output names the identity of the code that produced its numbers, `record.numerics`, and
`axtb.program.numerics(spec)` gives it before a run, for a store's cache key: the method's
epoch (`gfn2:1`, `gxtb-v1.1:1`, `gxtb-v2.0.1:1`, `gfnff:1`), one per method across all its
conventions; the solvation composite joins its parts' (`gxtb-v2.0.1:1+gfn2:1`). An epoch
changes with every axtb change that may change that method's results for an unchanged input,
and never with one that keeps every result bit for bit; thread counts and a run-time LAPACK are
not part of it. The CHANGELOG lists the epochs.

**Across machines.** A build repeats its results bit for bit on one machine and thread count. On
another CPU type its numbers differ by round-off, because the matrix kernels pick their blocking and
instructions at run time. The same binary run on an Intel i9-11900H and an AMD Threadripper PRO
9985WX, over the 177 cases of the CLI's golden set:
- **Bit-identical:** 92 cases: most GFN2-xTB and GFN-FF cases, almost no g-xTB case.
- **The rest:** energies ≤ 4.8e-13 Eh, gradients ≤ 1.4e-14 Eh/bohr, charges ≤ 8.8e-14,
  Hessians ≤ 2.5e-12 Eh/bohr² and frequencies ≤ 6.2e-8 cm⁻¹ apart, with the same SCF
  iterations.
- **Exceptions:** an SCF converging slowly or at the edge of its test can take another number
  of iterations. There the two machines agree only within the SCF's thresholds (the accuracy
  model's bound). The golden set's stretched-water retry case took 36 and 39 iterations.
  Ni(CO)4 at 5000 K with g-xTB v1.1 in its program's conventions took 14 and 13 on another
  machine, 2.2e-10 Eh apart.

So a stored result made on another CPU type agrees with a rerun to round-off as a rule, and
within the accuracy model's bound always.

## qccompute

axtb registers `axtb.program` under the `qccompute.programs` entry points, so qccompute finds
it by name, without axtb importing qccompute:

```python
from qccompute import compute, session
out = compute(inp)                                    # the runner grants the threads
with session(inp) as s:                               # qccompute's session protocol
    frame = s.compute(positions)
```

**A session's threads under qccompute.** qccompute grants a session what its input asks
for (the plan's threads for the input's calctype: 1 for a small molecule's gradient), unless
`session(inp, threads=N)` says more. So `compute_many` fans lineages out only when the session
is opened with `threads`. A NEB step of 12 caffeine images (g-xTB 2.0.1 gradients, warm) took
130 ms with `session(inp)` and 20 ms with `session(inp, threads=8)` on an 8-core laptop, with
the same frames bit for bit.

`axtb.program` declares energies, gradients and Hessians of every method, sessions, in-process
and thread-safe; `python -m qccompute.serve axtb` serves it to workers. `capabilities()`,
`settings_schema(calctype)`, `record_schema(calctype)`, `plan(inp)`, `compute(inp, ctx)`,
`open_session(inp, ctx, calctypes)` and `numerics(spec)` are the protocol's methods.
