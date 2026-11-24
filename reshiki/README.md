<!-- SPDX-License-Identifier: LGPL-3.0-or-later -->
# ReShiki validation

Run from the axtb checkout root with Python 3.10 or newer. The validation and platform tools use NumPy. The full
suite covers the 25 [ReShiki fixtures](structures/README.md), the 13 top-level upstream XYZ files, permanganate,
and singlet para-benzyne. The case IDs are the XYZ stems. [cases.py](cases.py) defines the method expectations,
Hessian subset, and tolerance sources.

## Reference validation on Linux x64

```sh
export AXTB_BIN=/absolute/path/to/axtb                 # default is target/release/axtb in this checkout
export AXTB_XTB680=/absolute/path/to/xtb-6.8.0/bin/xtb
export AXTB_GXTB2=/absolute/path/to/g-xtb-2.0.1/bin/xtb
export RESHIKI_TBLITE_PYTHON=/absolute/path/to/tblite-venv/bin/python
uv run --no-project --with numpy python reshiki/validate.py --jobs 4 --out reshiki/results
uv run --no-project --with numpy python reshiki/validate.py --only water --methods gfn2 --out reshiki/water-results
```

`AXTB_XTB680` and `AXTB_GXTB2` accept paths or command names. Preflight requires xtb `6.8.0 (4b8c987)`, the
g-xTB 2.0.1 build identifying as xtb `6.7.1 (26dd68d)`, and tblite `0.8.0`. It checks all three references
for filtered runs too. `RESHIKI_TBLITE_PYTHON` defaults to the validator's interpreter. When that interpreter
can install tblite, `uv run --no-project --with numpy --with 'tblite==0.8.0' python reshiki/validate.py`
also works without setting `RESHIKI_TBLITE_PYTHON`. Set `AXTB_SCRATCH` to put reference scratch directories
on a specific filesystem. The default output directory is `reshiki/results`.

`validate.py` writes `results.json` with per-check deviations, limits, sources, SCF iteration counts,
eigensolver identities, and diagnostics. `summary.md` contains the program and host provenance, per-check
maxima, every non-pass, and totals. Setup failures, out-of-tolerance comparisons, and unexpected calculation
or parsing failures exit 1. The declared UF₆ refusals, documented `gh106` SCF failures, and in-tolerance
long-molecule differences exit 0 when all other checks pass.

UF₆ is declared `both_error` for GFN2. GFN-FF supports actinide parameters.
See [What xtb 6.8.0 changed in GFN-FF](../docs/validation.md#what-xtb-680-changed-in-gfn-ff).

### Commands and settings

Every axtb calculation runs `"$AXTB_BIN" --json` with atomic numbers and bohr coordinates on stdin. Calculation
requests set `threads: 1` and carry the case charge and UHF. `gxtb` in a script argument maps to
`method: "gxtb-v2.0.1"`. Both LAPACK variables are removed from the axtb environment, and successful SCF
responses must report `eigensolver: "faer"`. Successful responses must identify the requested method and
conventions, including the Hessian method for Hessian calculations. `--version` records each CLI's identity.

Each reference command runs in separate scratch with `OMP_NUM_THREADS=1` and `MKL_NUM_THREADS=1`. xtb and
g-xTB run with `XTBPATH` unset and `XTBHOME` set to their own scratch directory. `coord` is a Turbomole
file in bohr. `Q` and `U` below stand for the case charge and unpaired electrons. GFN-FF uses
only `Q`. tblite runs as `"$RESHIKI_TBLITE_PYTHON" -c <child code>` with the input on stdin, using
`Calculator("GFN2-xTB", ...)` at accuracy `1e-6` and electronic temperature 300 K.

| Check | Reference command | axtb request |
|---|---|---|
| GFN2 `xtb_values` | `$AXTB_XTB680 coord --gfn 2 --grad --json --acc 1e-4 --chrg Q --uhf U` | `gradient: true`, `conventions: "xtb"`, `accuracy: 1e-6` |
| GFN2 `tblite_values` | tblite child as above | `gradient: true`, `conventions: "tblite"`, `accuracy: 1e-6` |
| GFN2 `default_vs_xtb` | axtb JSON with `conventions: "xtb"`, `accuracy: 1e-6` | `gradient: true`, default conventions at `1e-6`, same bohr coordinates |
| GFN2 `molden` | `$AXTB_XTB680 coord --gfn 2 --molden --norestart --acc 1e-4 --chrg Q --uhf U` | `molden: true`, `conventions: "xtb"`, `accuracy: 1e-4` |
| GFN2 Hessians | `$AXTB_XTB680 coord --gfn 2 --hess --acc 0.001 --chrg Q --uhf U` | Analytic Hessian with `conventions` and `analysis_conventions` both `"xtb"`, `hessian_accuracy: 1e-6` |
| g-xTB `gxtb_values` | `$AXTB_GXTB2 coord --gxtb --grad --json --norestart --acc 1e-4 --chrg Q --uhf U` | `gradient: true`, `conventions: "gxtb"`, `accuracy: 1e-4`, `max_iter: 250` |
| g-xTB `molden` | `$AXTB_GXTB2 coord --gxtb --molden --norestart --acc 1.0 --chrg Q --uhf U` | `molden: true`, `conventions: "gxtb"`, `accuracy: 1.0` |
| g-xTB Hessians | `$AXTB_GXTB2 coord --gxtb --hess --acc 1e-4 --chrg Q --uhf U` | Analytic Hessian with `conventions: "gxtb"` and `analysis_conventions: "xtb-6.7.1"` |
| GFN-FF `xtb_values` | `$AXTB_XTB680 coord --gfnff --grad --json --chrg Q` | `gradient: true`, `conventions: "xtb"` |
| GFN-FF `molden` | axtb refusal check | `molden: true`, expects `not_applicable` |
| GFN-FF Hessians | `$AXTB_XTB680 coord --gfnff --hess --chrg Q` | Analytic Hessian with `conventions` and `analysis_conventions` both `"xtb"` |

Value checks read the program's total energy, `gradient`, and charge output (`charges`, or full-precision EEQ
charges from `gfnff_lists.json`). GFN2 reads the dipole from `xtbout.json`. g-xTB charge deviations use the
reference charge file or printed population table and are required and ungated. SCF iteration counts are also
ungated. Molden checks compare sorted `Ene=` values within each `Spin=` channel. GFN2 and closed-shell g-xTB
may duplicate a restricted channel. Open-shell g-xTB requires both spin channels in both files. The reference
writer may omit trailing orbitals above 2 Eh.

Each Hessian run produces two checks on the same geometry. `hessian_analytic` compares the vibrational modes
of axtb's analytic Hessian against the reference printout. `hessian_same` sends the reference's `hessian`
matrix to a separate `axtb --json` request with `analyze_hessian`, the version-specific `analysis_conventions`,
and its Hessian-frame geometry. The g-xTB build's small pre-Hessian rotation is applied to that geometry.
The checks expect 3N−5 vibrational modes for linear molecules and 3N−6 for nonlinear molecules. A deficit in
`g98.out` must match the axtb modes at |ν| ≤ 0.1 cm⁻¹, which are dropped before comparison. GFN-FF's
`vibspectrum` must list all 3N modes. Its leading rigid-body modes are removed before comparison.

The GFN2 value limits come from [the reference matrix](../tests/reference_matrix/run.py). Beyond 25 bohr, the
default-versus-xtb check uses the documented real-space cutoff differences of `9.4e-6 Eh` and
`3.1e-6 Eh/bohr` for energy and gradient. The cutoff check compares at the sources' two significant figures
and retains the raw deviations. Charge and dipole limits stay at the matrix values. The GFN2
analytic Hessian limit comes from [the matrix summary](../tests/reference_matrix/summary.md). Vibrational
same-Hessian limits come from [the vibrational test](../tests/vibrational_check.py) and the GFN-FF table.
The g-xTB single-point gates use tight accuracy. The reproduction test measured the cited energy bound at
default accuracy and the cited gradient bound at tight accuracy. [g-xTB reproduction test](../tests/gxtb2_reproduction.py)
The other limits come from the Molden, vibrational analysis, cutoff, and GFN-FF sections of
[upstream validation](../docs/validation.md). The GFN2 Molden limit was measured with xtb 6.7.1 upstream.
ReShiki applies that limit to xtb 6.8.0 at the requested accuracy. The GFN2 analytic Hessian check measures
matched modes with reference |ν| > 50 cm⁻¹, as in the upstream matrix. In `cases.py`, each entry gives its
source and units.

## Prepare Hessian geometries

```sh
uv run --no-project --with numpy python reshiki/optimize_structures.py
```

The one-off preparation needs the two pinned xtb binaries. It runs the following for each Hessian case and
stores converged coordinates as `reshiki/structures/opt/<id>.<method>.xyz` with charge, UHF, and program provenance.

```sh
"$AXTB_XTB680" coord --gfn 2 --opt tight --chrg Q --uhf U
"$AXTB_XTB680" coord --gfnff --opt tight --chrg Q
"$AXTB_GXTB2" coord --gxtb --opt tight --chrg Q --uhf U
```

Use `--overwrite` to replace prepared geometries. Both Hessian workflows choose the original geometry if the
method-specific file is absent. The analytic frequency bounds were documented on upstream cases, so the new
set can expose differences, especially on unoptimized geometries. Single-point values and Molden checks always
use the original XYZ.

## Compare platforms

Build the fork's axtb CLI on each host and run the dump. Copy the Linux x64 baseline and use the same checkout
and optimized XYZ set on each host. A quick smoke run can add `--only water`.

```sh
uv run --no-project --with numpy python reshiki/platform_dump.py --out linux-x64.json
uv run --no-project --with numpy python reshiki/platform_dump.py --out other.json
uv run --no-project --with numpy python reshiki/platform_compare.py linux-x64.json other.json --md comparison.md
```

The dump calls only `"$AXTB_BIN" --json` with default conventions and one thread, plus `"$AXTB_BIN" --version`
and `rustc --version` for provenance. It keeps gradients, energy, charges, iterations, analytic vibrational
frequencies on the Hessian subset, and complete axtb diagnostics. The comparator gates energy at `1e-8 Eh`,
gradient at `1e-7 Eh/bohr`, and sorted vibrational frequencies at `0.1 cm⁻¹`. These are the lead's chosen
platform gates. Both dumps must contain the selected runs and matching input fingerprints for each pair.
Hessian pairs must have equal frequency counts. At least one numeric pair is required, and the only accepted
error pair is `unsupported_element` for UF₆ with GFN2. Charges and iterations are reported without limits.
The comparator exits 1 on a violation. Regenerate older dumps that lack input fingerprints.

## Regenerate the XYZ inputs and check the scripts

```sh
uv run --no-project --with rdkit --with numpy python reshiki/make_structures.py
uv run --no-project python -m py_compile reshiki/*.py
uvx ruff check reshiki
uvx ruff format --check reshiki
```
