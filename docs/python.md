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

`compute_many(requests)` evaluates `qcdata.Request`s in order on the one engine, each reading its
lineage's history as committed by the requests before it; failures stay in place and raise
`qcdata.ProgramErrors` at the end, with every frame. `reset(lineage=...)` forgets every
lineage's history (no argument) or one lineage's; `lineages()` lists those with history;
`plan()` returns the settings pinned at opening, the threads and the record's static part.

Interleaved lineages give what separate sequential sessions give, bit for bit, and a session's
frames equal the CLI's trajectory cases of the golden set bit for bit. A session's methods are
thread-safe (a lock serializes them; the GIL is released while the engine runs). A session step
costs about 30 µs of Python on top of the calculation, `compute` about 100 µs (water).

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
| `energy`, `gradient`, `hessian` | Eh, Eh/bohr (n_atoms x 3), Eh/bohr² (3N x 3N), as the calctype asks |
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
| `scf` | `start` (`eeq`, `sad`, `eeq-bc`, `previous`, `given`), `iterations`, `mixer_used` and `attempts` (GFN2-xTB's ladder), `restart_fallback` (a session's SCF from history failed and ran again fresh) |
| `energy_subtotals` | terms that overlap the energy components (the electronic energy) |
| `contingent` | contingent settings and whether each fired (`scf.mixer`, `composite.retries`) |
| `accuracy` | with a target: the target, the thresholds used and the model's estimate at them |
| `hessian` | the Hessian's method and largest asymmetry before symmetrization |
| `composite` | the solvation composite: its definition, the energies of its parts and their SCFs (`wavefunction.part` names the part the orbitals are) |
| `gfnff` | GFN-FF: hydrogen and halogen bonds, the fragment charges, and (compute only) the topology |
| `warm_start`, `from_lineage` | session frames: whether the SCF started from history, and whose |
| `error` | a failure's fields: `z` of an unsupported element; `iterations`, `energy_change`, `residual` of an SCF that did not converge |

### Warnings

`warnings` are `qcdata.ProgramWarning`s with a code: `resource_capped` (the `threads` keyword
above the grant), `accuracy_limited` (a target below what axtb's thresholds reach),
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
| `convergence` | `scf_not_converged` | the SCF did not converge (`details`: `iterations`, `energy_change`, `residual`) |
| `numerical` | `numerical_failure` | a numerical failure (a singular matrix, a failed eigensolver) |
| `internal` | `panic` | a bug: please report it (a session closes) |

## Threads

A calculation runs on exactly its threads and never spawns more. Under a runner they are the
`threads` keyword, at most the grant (above it: the grant and the warning `resource_capped`),
or else the grant. A direct call without the keyword runs on the threads the plan states:

- 1 for a molecule at or below the size where threads only cost (the results are the same
  bits on any thread count there);
- energies and gradients: 1 below 300 atoms, 8 from 300 (GFN-FF: 1 below 1000, 8 from 1000);
- Hessians above the one-thread size: every CPU of the process (the plan asks for the whole
  grant).

`axtb.program.plan(inp)` gives the same number as `resources.threads` (None: the whole grant).
A session chooses its threads at opening, for its largest calctype.

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

## qccompute

axtb registers `axtb.program` under the `qccompute.programs` entry points, so qccompute finds
it by name, without axtb importing qccompute:

```python
from qccompute import compute, session
out = compute(inp)                                    # the runner grants the threads
with session(inp) as s:                               # qccompute's session protocol
    frame = s.compute(positions)
```

`axtb.program` declares energies, gradients and Hessians of every method, sessions, in-process
and thread-safe; `python -m qccompute.serve axtb` serves it to workers. `capabilities()`,
`settings_schema(calctype)`, `record_schema(calctype)`, `plan(inp)`, `compute(inp, ctx)`,
`open_session(inp, ctx, calctypes)` and `numerics(spec)` are the protocol's methods.
