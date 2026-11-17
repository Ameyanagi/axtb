# Changelog

All notable changes to axtb. Versions follow [Semantic Versioning](https://semver.org); the
Rust crates, the CLI and the Python package share one version.

## Unreleased

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
  `not_applicable` (`axtb_base::Error::NotApplicable`; Python: `axtb.InputError` with that kind).
- Every SCF result reports the Fermi levels of the alpha and beta electrons (`fermi_levels`, Eh;
  `null` where the occupations are the Aufbau filling, the same level twice for a closed shell)
  and the electronic temperature the occupations used (`electronic_temperature`, K; 0 for the
  Aufbau filling, e.g. GFN2-xTB's program conventions up to 0.1 K). In Rust:
  `Electronic::fermi_levels` and `Electronic::electronic_temperature`. New keys in the CLI's JSON
  and Python's dict; nothing else changes.
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
- The Python package's `qc` extra requires `qcdata<0.24`: qcdata 0.24 removes the
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
