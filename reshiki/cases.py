# SPDX-License-Identifier: LGPL-3.0-or-later
"""Structures, expected outcomes, and provenance for the reference gates."""

import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
STRUCTURES = HERE / "structures"
METHODS = ("gfn2", "gfnff", "gxtb")
HESSIAN_CASES = {
    "aspirin",
    "paracetamol",
    "acetate",
    "methyl-radical",
    "co2",
    "hcn",
    "nickel-tetracarbonyl",
    "water-dimer",
    "caffeine",
}
ELEMENT_SYMBOLS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
    "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
    "Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db "
    "Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
)
ATOMIC_NUMBERS = {symbol: z for z, symbol in enumerate(ELEMENT_SYMBOLS.split(), start=1)}
# Å per bohr, from crates/base/src/units.rs (mctc-lib derived and qcref respectively).
BOHR_2018 = 0.5291772109044924
BOHR_2022 = 0.529177210544
CUTOFF_BOHR = 25.0

MATRIX = "tests/reference_matrix/run.py"
DOCS = "docs/validation.md"
VIB = "tests/vibrational_check.py LIMITS"
GXTB = "tests/gxtb2_reproduction.py LIMITS"

EXPECTED_ERRORS = {("uf6", "gfn2"): "unsupported_element"}
EXPECTED_BOTH_FAIL = {
    ("gh106", "gfn2", "tblite_values"): {
        "source": 'tests/validate_tblite.py KNOWN_DIFFERENCES["gh106"]: tblite 0.8.0 does not converge',
        "diagnostic": r"^(?:TBLiteRuntimeError: )?SCF not converged(?: in \d+ cycles)?$",
    },
    ("gh106", "gxtb", "gxtb_values"): {
        "source": "tests/structures/README.md: gh106 is a pathological SCF case",
        "diagnostic": r"^\s*(?:-\d+-\s+)?tblite_calculator_singlepoint: SCF not converged(?: in \d+ cycles)?\s*$",
    },
}

VALUE_QUANTITIES = ("energy", "gradient", "charges")
GFN2_VALUE_QUANTITIES = (*VALUE_QUANTITIES, "dipole")
REQUIRED_QUANTITIES = {
    ("gfn2", "xtb_values"): GFN2_VALUE_QUANTITIES,
    ("gfn2", "tblite_values"): GFN2_VALUE_QUANTITIES,
    ("gfn2", "default_vs_xtb"): GFN2_VALUE_QUANTITIES,
    ("gfn2", "default_vs_xtb_long"): ("energy", "gradient"),
    ("gfnff", "xtb_values"): VALUE_QUANTITIES,
    ("gxtb", "gxtb_values"): VALUE_QUANTITIES,
    ("gfn2", "molden"): ("orbital_energies",),
    ("gfnff", "molden"): (),
    ("gxtb", "molden"): ("orbital_energies",),
    **{(method, check): ("frequencies",) for method in METHODS for check in ("hessian_analytic", "hessian_same")},
    ("gfnff", "hessian_numerical"): ("frequencies",),
}

# Each value is (limit, source). None means the quantity is reported without a gate.
TOLERANCES = {
    ("gfn2", "xtb_values", "energy"): (1e-7, f'{MATRIX} TOL["xtb"]["E"] = 1e-7 Eh'),
    ("gfn2", "xtb_values", "gradient"): (1e-7, f'{MATRIX} TOL["xtb"]["G"] = 1e-7 Eh/bohr'),
    ("gfn2", "xtb_values", "charges"): (1e-6, f'{MATRIX} TOL["xtb"]["q"] = 1e-6 e'),
    ("gfn2", "xtb_values", "dipole"): (1e-5, f'{MATRIX} TOL["xtb"]["dip"] = 1e-5 e·bohr'),
    ("gfn2", "tblite_values", "energy"): (1e-9, f'{MATRIX} TOL["tblite"]["E"] = 1e-9 Eh'),
    ("gfn2", "tblite_values", "gradient"): (1e-8, f'{MATRIX} TOL["tblite"]["G"] = 1e-8 Eh/bohr'),
    ("gfn2", "tblite_values", "charges"): (1e-7, f'{MATRIX} TOL["tblite"]["q"] = 1e-7 e'),
    ("gfn2", "tblite_values", "dipole"): (1e-6, f'{MATRIX} TOL["tblite"]["dip"] = 1e-6 e·bohr'),
    ("gfn2", "default_vs_xtb", "energy"): (1e-7, f'{MATRIX} TOL["axtb_vs_xtbmode"]["E"] = 1e-7 Eh'),
    ("gfn2", "default_vs_xtb", "gradient"): (1e-6, f'{MATRIX} TOL["axtb_vs_xtbmode"]["G"] = 1e-6 Eh/bohr'),
    ("gfn2", "default_vs_xtb", "charges"): (1e-5, f'{MATRIX} TOL["axtb_vs_xtbmode"]["q"] = 1e-5 e'),
    ("gfn2", "default_vs_xtb", "dipole"): (1e-4, f'{MATRIX} TOL["axtb_vs_xtbmode"]["dip"] = 1e-4 e·bohr'),
    ("gfn2", "default_vs_xtb_long", "energy"): (
        9.4e-6,
        "tests/reference_matrix/summary.md real-space cutoffs, 9.4e-6 Eh for taxol",
    ),
    ("gfn2", "default_vs_xtb_long", "gradient"): (3.1e-6, f"{DOCS} What xtb 6.8.0 changed, 3.1e-6 Eh/bohr"),
    ("gfn2", "molden", "orbital_energies"): (
        7.0e-10,
        f"{DOCS} Molden GFN2-xTB at 1e-4 (both), 7.0e-10 Eh measured with xtb 6.7.1, applied here to 6.8.0",
    ),
    ("gfn2", "hessian_analytic", "frequencies"): (
        6.6,
        (
            "tests/reference_matrix/summary.md hessian_analytic/analytic_xtb max |Δν| 6.6 cm⁻¹ over modes with "
            f"|ν_xtb| > 50 cm⁻¹ ({MATRIX} check_analytic_xtb)"
        ),
    ),
    ("gfn2", "hessian_same", "frequencies"): (5.0e-5, f'{VIB} ["xtb"]["GFN2 same Hessian: frequencies"] 5.0e-5 cm⁻¹'),
    ("gxtb", "gxtb_values", "energy"): (
        2.1e-12,
        f'{GXTB} ["max_abs_dE"] 2.1e-12 Eh at default accuracy, applied here at --acc 1e-4',
    ),
    ("gxtb", "gxtb_values", "gradient"): (6e-11, f'{GXTB} ["max_dgrad_tight"] 6e-11 Eh/bohr at --acc 1e-4'),
    ("gxtb", "gxtb_values", "charges"): (None, "No documented g-xTB 2.0.1 charge gate"),
    ("gxtb", "molden", "orbital_energies"): (3.3e-12, f"{DOCS} Molden g-xTB 2.0.1, 3.3e-12 Eh"),
    ("gxtb", "hessian_analytic", "frequencies"): (
        0.35,
        f"{DOCS} Vibrational analysis, g-xTB 2.0.1 analytic vs the program at --acc 1e-4, 0.35 cm⁻¹",
    ),
    ("gxtb", "hessian_same", "frequencies"): (
        5.0e-5,
        f'{VIB} ["xtb-6.7.1"]["g-xTB program Hessian: frequencies"] 5.0e-5 cm⁻¹',
    ),
    ("gfnff", "xtb_values", "energy"): (7.0e-13, f"{DOCS} GFN-FF vs xtb 6.8.0 energies, 7.0e-13 Eh"),
    ("gfnff", "xtb_values", "gradient"): (4.3e-12, f"{DOCS} GFN-FF vs xtb 6.8.0 gradients, 4.3e-12 Eh/bohr"),
    ("gfnff", "xtb_values", "charges"): (2.3e-12, f"{DOCS} GFN-FF vs xtb 6.8.0 EEQ charges, 2.3e-12 e"),
    ("gfnff", "hessian_analytic", "frequencies"): (0.27, f"{DOCS} GFN-FF analytic Hessian, 0.27 cm⁻¹"),
    ("gfnff", "hessian_same", "frequencies"): (
        0.005,
        'tests/gfnff_reproduction.py LIMITS["xtb"]["hessian"]["program_hessian frequencies_maxdiff_cm1"] 0.005 cm⁻¹',
    ),
    ("gfnff", "hessian_numerical", "frequencies"): (
        0.005,
        f"{DOCS} GFN-FF table: axtb numerical Hessian vs xtb 6.8.0 --hess, frequencies 0.005 cm⁻¹",
    ),
}


def vibrational_mode_count(positions_bohr):
    nat = len(positions_bohr)
    if nat <= 1:
        return 0
    # xtb 6.8.0's is_linear tests the smallest unit-mass moment about the centroid.
    centered = positions_bohr - positions_bohr.mean(axis=0)
    singular_values = np.linalg.svd(centered, compute_uv=False)
    linear = np.sum(singular_values[1:] ** 2) < 1e-4
    return 3 * nat - (5 if linear else 6)


def read_structure(path):
    lines = path.read_text().splitlines()
    rows = [line.split() for line in lines[2:] if line.strip()]
    if len(rows) != int(lines[0]) or not rows:
        raise ValueError(f"XYZ atom count does not match the geometry in {path}")
    symbols = [row[0].capitalize() for row in rows]
    return symbols, np.array([row[1:4] for row in rows], dtype=float)


def read_state(path):
    with path.open() as stream:
        next(stream)
        comment = next(stream)
    charge = re.search(r"\bcharge\s*[:=]\s*([+-]?\d+)", comment, re.IGNORECASE)
    uhf = re.search(r"\buhf\s*=\s*(\d+)", comment, re.IGNORECASE)
    symbols, _ = read_structure(path)
    q = int(charge[1]) if charge else 0
    u = int(uhf[1]) if uhf else (sum(ATOMIC_NUMBERS[s] for s in symbols) - q) % 2
    return q, u


def load_cases(only=None):
    tests = ROOT / "tests"
    paths = [
        *STRUCTURES.glob("*.xyz"),
        *(tests / "structures").glob("*.xyz"),
        tests / "scf_hard/tm-MnO4-.xyz",
        tests / "scf_hard/gap-p-benzyne.xyz",
    ]
    cases = []
    for path in sorted(paths, key=lambda path: path.stem.casefold()):
        if only and only.casefold() not in path.stem.casefold():
            continue
        charge, uhf = read_state(path)
        cases.append(
            {
                "id": path.stem,
                "path": path,
                "charge": charge,
                "uhf": uhf,
                "methods": {m: "both_error" if (path.stem, m) in EXPECTED_ERRORS else "match" for m in METHODS},
                "hessian": path.stem in HESSIAN_CASES,
                "tolerances": TOLERANCES,
            }
        )
    return cases


def optimized_path(case, method):
    return STRUCTURES / "opt" / f"{case['id']}.{method}.xyz"


def hessian_path(case, method):
    path = optimized_path(case, method)
    return path if path.is_file() else case["path"]


def tolerances(case, method, check, long=False):
    specs = case["tolerances"]
    missing = [q for q in REQUIRED_QUANTITIES.get((method, check), ()) if (method, check, q) not in specs]
    if missing:
        raise ValueError(f"Missing tolerance entries for {method}/{check}: {', '.join(missing)}")
    result = {
        q: {"limit": limit, "source": source, **({"status": "not_gated"} if limit is None else {})}
        for (m, c, q), (limit, source) in specs.items()
        if m == method and c == check
    }
    if long and check == "default_vs_xtb":
        for quantity, spec in tolerances(case, method, "default_vs_xtb_long").items():
            result[quantity] = {**spec, "comparison": "two significant figures as published"}
    return result
