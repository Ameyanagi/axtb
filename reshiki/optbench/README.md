<!-- SPDX-License-Identifier: LGPL-3.0-or-later -->

# ReShiki optimizer benchmark

This standalone Rust package runs ReShiki's geometry optimizer against axtb's
GFN2-xTB and GFN-FF and COSMolKit's MMFF94. The Python tools run official
`xtb --opt` and compare evaluation counts, timings, convergence, and final
minima. The harness is LGPL-3.0-or-later as part of the axtb fork. It calls
ReShiki's MIT-licensed optimizer as a library. The optimizer crate also offers
Apache-2.0 as a license option.

## Build and prepare

Run these commands from `reshiki/optbench`. The separate Cargo workspace
carries the axtb Linux GEMM cache patches and builds the driver without g-xTB.

```sh
export CARGO_TARGET_DIR="$(cd ../../.. && pwd -P)/optbench-target"
cargo build --release
bench="$CARGO_TARGET_DIR/release/reshiki-optbench"
"$bench" prepare --manifest cases.toml --out starts
```

Use the Cargo path override to work against a local ReShiki checkout. The same
`--config` option works with `cargo build` and `cargo clippy`.

```sh
cargo run --release \
  --config 'patch."https://github.com/Ameyanagi/ReShiki".reshiki-optimize.path="/abs/ReShiki/crates/optimize"' \
  -- prepare --manifest cases.toml --out starts
```

`cases.toml` contains 77 cases in six sets. Baker starts use the manifest's
SMILES and newly embedded coordinates in place of the original Z-matrices.
[Baker, J. Comput. Chem. 14, 1085–1100 (1993)](https://doi.org/10.1002/jcc.540140910)

For SMILES, `prepare` makes explicit hydrogens, sanitizes, and embeds one
ETKDGv3 conformer with seed 42. It writes the raw `etkdg` structure and a
`mmff` structure relaxed by COSMolKit's MMFF94 minimizer. A failed embedding
skips the case and prints the reason. Unavailable MMFF94 typing leaves a copy
of the ETKDG structure as the `mmff` start for xTB runs. The XYZ metadata marks
the fallback and the MMFF runner records an error for that case. The XYZ
element order and an indexed connectivity fingerprint are checked when
rebuilding MMFF from the original SMILES. The COSMolKit `cap-transforms`
feature allows the runner to attach those coordinates to the reconstructed
molecule directly.

Use a fresh `--out` directory after changing the manifest. A repeated
preparation overwrites successful starts in place and leaves older files
on disk for cases absent from the new output.

XYZ cases produce `given` and `rough` starts. The latter adds deterministic
Gaussian noise with σ = 0.10 Å per coordinate, seeded by the case name. The
idealized geometries for the four added metal complexes are documented in the
XYZ comment lines in `structures/`.

## Run benchmarks

```sh
for potential in gfn2 gfnff mmff; do
    "$bench" run --starts starts --potential "$potential" --configs anc,lbfgs \
      --sets baker --cases water,ethanol,caffeine --jobs 3 --out "results/$potential.jsonl"
done
```

Omit the filters for the full sweep or add `--starts-kind etkdg,mmff,given,rough`
to choose specific starts. `--cases` accepts names and `set/name` pairs.
Caffeine occurs in two sets. Use one output directory per repeated experiment
with the same identities, such as a restart or accuracy sweep, so the JSONL
identities and saved geometries stay distinct. `run` appends and flushes one
row per job, including optimizer and potential failures. Geometry paths are
relative to the JSONL file's directory. The worker pool runs `--jobs` jobs
with single-threaded potential calls. Run with `--jobs 1` for timings without
competition among jobs.

A configuration is a method followed by optional `+` modifiers.

| Config element | Meaning |
|---|---|
| `anc`, `lbfgs` | ANC-RFO or Cartesian L-BFGS |
| `internal` | Translation–rotation internal coordinates (TRIC), available with `--features internal` and an optimizer that exports `Method::Internal`. |
| `hess=model` | Model Hessian, the optimizer default |
| `hess=self`, `hess=gfnff` | Analytic Hessian of the target potential or a separate GFN-FF calculator |
| `hess=auto<N>` | Use the target's analytic Hessian for up to N atoms where available, for example `hess=auto30`. |
| `update=bfgs\|damped\|bofill` | Hessian update |
| `accel=none\|gdiis\|gediis` | Extrapolation scheme |
| `adaptive` | Enable the optimizer's adaptive accuracy schedule. |
| `tight`, `notightfinal` | Use tight optimization thresholds or turn off the final tight-accuracy verification. |
| `preopt=gfnff` | Run the same config on GFN-FF first, then use its last accepted geometry for the target stage. |

The default build enables `knobs` for the published optimizer branch. The
`internal` feature enables `knobs` and requires the upstream TRIC variant.
Use `--no-default-features` for the original ANC-RFO and L-BFGS API.
MMFF94 operates on typed SMILES cases.
Its analytic Hessian is unavailable, so use `hess=model`, `hess=auto<N>` for
the model fallback, or `hess=gfnff`.

`--no-restart` resets the GFN2 wavefunction before every potential call.
`--acc-coarse 30 --acc-normal 1 --acc-tight 0.2` are the default GFN2 SCF
accuracy multipliers. GFN-FF and MMFF use fixed numerical accuracy.

### Recorded metrics

- Energy and maximum Cartesian gradient use Eh and Eh/bohr. `evaluations`
  counts all optimizer energy-plus-gradient calls, including rejected trials
  and failures. `trace` records the initial and accepted geometries and any
  accepted re-evaluation at a new accuracy. Each sample is `[evaluations,
  seconds, energy, max_gradient, accuracy, gradient_norm]`, using cumulative
  counts and times. Both gradient measures are in Eh/bohr.
  `scf_iterations_total` counts SCF iterations visible in the driver response
  or SCF error.
- `potential_seconds` covers calculator setup plus gradient calls. The
  `setup_seconds` and `gradient_seconds` fields give those components.
  `hessian_seconds` covers the driver's analytic Hessian task with vibrational
  analysis disabled, including any single point the driver performs as part of
  that task. Optimizer overhead is the remaining wall time. File I/O and job
  queue time are outside the run timer.
- Hessian cost divides Hessian time by the run's median measured gradient-call
  time. The Hessian count follows the optimizer API, which includes failed
  calls and availability probes. `preopt` and `target` hold stage counts and
  times. Top-level counts and times sum the stages. The top-level trace uses
  target energies with the preoptimization time and evaluation offsets. For
  GFN2 targets, GFN-FF gradient time divided by the target stage's median
  gradient time gives `preopt_gfn2_gradient_equivalents`. The total also includes
  the target gradient evaluation count. The CSV includes those equivalents.

## Official xtb and tight references

```sh
uv run --no-project --with numpy xtb_opt.py --starts starts --potential gfn2 \
  --xtb /path/to/xtb --jobs 3 --sets baker --cases water,ethanol,caffeine \
  --out results/xtb-gfn2.jsonl
```

`$AXTB_XTB680` supplies the default executable path. Choose `--potential gfnff`
for official GFN-FF and `--level tight` for a tighter comparison. The script
isolates each job in a fresh temporary directory and saves its output log and
validated atom-ordered final geometry.

For official xtb 6.8.0 ANCopt, `evaluations` counts the printed `CYCLE` blocks.
The JSONL also preserves the separately reported optimizer iteration count.
GFN2 SCC totals include the initial and final singlepoints when all tables
are available. The potential and overhead times come from ANCopt's printed
timers. Potential time covers the optimizer singlepoints, overhead includes
the model Hessian, and process wall time covers the whole executable invocation.
Each complete cycle contributes its printed energy and gradient norm to the
trace, with seconds set to zero, maximum gradient set to `null`, and accuracy
set to `normal`. The final gradient norm is saved in `gradient_norm`.
[xtb 6.8.0 main program](https://github.com/grimme-lab/xtb/blob/v6.8.0/src/prog/main.F90)
[ANCopt loop](https://github.com/grimme-lab/xtb/blob/v6.8.0/src/optimizer.f90)
[Geometry optimization driver](https://github.com/grimme-lab/xtb/blob/v6.8.0/src/geoopt_driver.f90)

```sh
"$bench" reference --starts starts --potential gfn2 --results results/gfn2.jsonl \
  --xtb results/xtb-gfn2.jsonl --out results/gfn2-refs.jsonl
"$bench" reference --starts starts --potential gfnff --results results/gfnff.jsonl \
  --out results/gfnff-refs.jsonl
"$bench" reference --starts starts --potential mmff --results results/mmff.jsonl \
  --out results/mmff-refs.jsonl
```

Omit `--xtb` if no official run is available. For each case and start present
in the supplied results, the reference command scores the candidate
geometries on the requested harness potential and refines the lowest at
`Level::Tight` with up to 3000 steps. GFN2 uses tight SCF throughout. GFN-FF
uses the topology determined at the original start for the benchmark and its
reference. The lowest scored or accepted energy and its geometry are retained
even if refinement stalls, with the status of the tight optimization in the
reference row. This is a local reference seeded by the available benchmark
results. Official xtb's reported energies are retained in its own result rows
and can differ from axtb's conventions.

## Analyze

```sh
uv run --no-project --with numpy --with matplotlib analyze.py \
  results/gfn2.jsonl results/gfnff.jsonl results/mmff.jsonl \
  --refs results/gfn2-refs.jsonl results/gfnff-refs.jsonl results/mmff-refs.jsonl \
  --xtb results/xtb-gfn2.jsonl --baseline anc --out report
```

Omit `--xtb` when benchmarking without the official executable. Use
`--configs anc,lbfgs` to select plotted configs, and
`--molecules water,ethanol,caffeine` to choose up to four convergence panels.
The outputs are `summary.md`, `per_molecule.csv`, and three 150-dpi PNGs for
evaluation counts, energy above the reference versus calls and time, and the
evaluation-count CDF. Paired ratios use runs with the same set, case, start,
and potential. A run reaches the same minimum if its Kabsch-aligned heavy-atom
RMSD is below 0.1 Å or its absolute energy difference is below 10⁻⁶ Eh.
`xtb_criterion_evaluations` gives the cumulative evaluation count at the
first qualifying pair of consecutive samples at the same accuracy. The
absolute energy change must be below 5 × 10⁻⁶ Eh and the second sample's
gradient norm must be below 10⁻³ Eh/bohr. The analyzer also accepts older
five-field traces, whose missing gradient norms leave this metric unavailable.
The summary reports both native and criterion evaluation ratios using the
official run's recorded evaluation count as the denominator.

```sh
uv run --no-project --with numpy --with matplotlib python -m unittest discover -v
uv run --no-project --with ruff ruff check .
uv run --no-project --with ruff ruff format --check .
cargo fmt --check
cargo clippy --all-targets -- -D warnings
```
