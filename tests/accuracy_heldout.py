#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Accuracy targets on molecules outside the calibration set: every measured error at or below
its target.

    tests/accuracy_heldout.py [--methods gfn2,gxtb-v2.0.1,gxtb-v1.1] [--extra]
                              [--calibration FILE] [--out FILE]

The accuracy model (axtb/_accuracy.py) is fitted to tests/accuracy_calibration.py's 20 molecules.
This check runs it on 10 others, from crates/gxtb/tests/reference: acetate, ammonium, benzene,
triplet O2, ethanol, glycine, hydroxide, iodomethane, mesyl chloride and methane. For each
molecule, method (GFN2-xTB, g-xTB 2.0.1 and v1.1, axtb conventions) and target (gradient 1e-4,
1e-5 and 1e-6 Eh/bohr, energy a tenth of it), `axtb.compute` runs a gradient with the `accuracy`
keyword and is compared with the same molecule at `accuracy_scale` 1e-4: |dE| against the energy
target and the largest |dg| component against the gradient target (90 runs). It passes when every
error is at or below its target. It also prints each error over the model's estimate (its bound
at the factor chosen, record.accuracy.estimate) and the worst of those ratios per method.

`--extra` adds Ni(CO)4, CuCl4 2- (doublet) and SiH4/PH3, reported apart so that the 90 runs stay
comparable. `--calibration FILE` replaces the model's constants for the run: a JSON file holding
a CALIBRATION dict (as axtb._accuracy.CALIBRATION), e.g. dumped from another version, or a
calibration's data file (bench/results/accuracy.json), whose fit is turned into constants as
`tests/accuracy_calibration.py --fit FILE --write` would. Run it after every recalibration (a
manual step: it needs the axtb package, not the reference programs). Exit status 1 when an error
exceeds its target.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "crates/gxtb/tests/reference"
MOLECULES = (
    "acetate",
    "ammonium",
    "benzene",
    "dioxygen_triplet",
    "ethanol",
    "glycine",
    "hydroxide",
    "iodomethane",
    "mesyl_chloride",
    "methane",
)
EXTRA = ("nickel_tetracarbonyl", "tetrachlorocuprate", "silane_phosphine")
METHODS = ("gfn2", "gxtb-v2.0.1", "gxtb-v1.1")
GRADIENT_TARGETS = (1e-4, 1e-5, 1e-6)
REFERENCE_SCALE = 1e-4


def structure(name: str):
    """A reference molecule with its charge and multiplicity."""
    from qcdata import Structure

    d = REFERENCE / name
    inp = json.loads((d / "input.json").read_text())
    return Structure.open(d / "mol.xyz").model_copy(
        update={"charge": int(inp["charge"]), "multiplicity": int(inp["uhf"]) + 1}
    )


def gradient(method: str, mol, keywords: dict):
    import axtb
    from qcdata import ProgramInput

    i = ProgramInput(
        program="axtb",
        calctype="gradient",
        structure=mol,
        model={"method": method},
        keywords={"threads": 1, **keywords},
    )
    return axtb.compute(i)


def constants(path: Path) -> dict:
    """The CALIBRATION dict of a file: one itself, or the constants of a calibration's fit."""
    data = json.loads(path.read_text())
    if "rows" not in data:
        return data
    sys.path.insert(0, str(ROOT / "tests"))
    import accuracy_calibration as cal

    return cal.constants(data["fit"], data)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--methods", default=",".join(METHODS))
    ap.add_argument("--extra", action="store_true", help="also Ni(CO)4, CuCl4 2- and SiH4/PH3, reported apart")
    ap.add_argument("--calibration", default="", help="the model's constants for this run (see above)")
    ap.add_argument("--out", default="", help="write every run as JSON")
    a = ap.parse_args()
    from axtb import _accuracy

    if a.calibration:
        _accuracy.CALIBRATION = constants(Path(a.calibration))
    sets = {"held-out": MOLECULES, **({"extra": EXTRA} if a.extra else {})}
    rows, failed = [], 0
    for method in a.methods.split(","):
        for group, names in sets.items():
            for name in names:
                mol = structure(name)
                ref = gradient(method, mol, {"accuracy_scale": REFERENCE_SCALE}).results
                for g in GRADIENT_TARGETS:
                    target = {"energy": g / 10, "gradient": g}
                    out = gradient(method, mol, {"accuracy": target})
                    r = out.results
                    de = abs(r.energy - ref.energy)
                    dg = float(np.abs(np.asarray(r.gradient) - np.asarray(ref.gradient)).max())
                    acc = out.record["accuracy"]
                    est = acc["estimate"]
                    ok = de <= target["energy"] and dg <= target["gradient"]
                    failed += not ok
                    row = {
                        "group": group,
                        "method": method,
                        "molecule": name,
                        "target": target,
                        "accuracy_scale": acc["thresholds"]["accuracy_scale"],
                        "iterations": r.scf_iterations,
                        "dE": de,
                        "dG": dg,
                        "estimate": est,
                        "dE/estimate": de / est["energy"] if est.get("energy") else None,
                        "dG/estimate": dg / est["gradient"] if est.get("gradient") else None,
                        "ok": ok,
                    }
                    rows.append(row)
                    print(
                        f"{group:8} {method:12} {name:22} g {g:.0e}: factor {row['accuracy_scale']:<6g} "
                        f"it {r.scf_iterations:3}  dE {de:.1e} dG {dg:.1e}  dE/est {row['dE/estimate'] or 0:.2f} "
                        f"dG/est {row['dG/estimate'] or 0:.2f}" + ("" if ok else "  OVER TARGET"),
                        flush=True,
                    )
    print("summary (runs, over target, worst error/target, worst error/estimate):")
    for method in a.methods.split(","):
        for group in sets:
            rs = [r for r in rows if r["method"] == method and r["group"] == group]
            worst_t = max(max(r["dE"] / r["target"]["energy"], r["dG"] / r["target"]["gradient"]) for r in rs)
            worst_e = max(max(r["dE/estimate"] or 0.0, r["dG/estimate"] or 0.0) for r in rs)
            over = sum(not r["ok"] for r in rs)
            print(f"  {method:12} {group:8} {len(rs):3} runs, {over} over target, {worst_t:.3f}, {worst_e:.3f}")
    if a.out:
        Path(a.out).write_text(json.dumps(rows, indent=1))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
