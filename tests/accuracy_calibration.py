#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""The calibration of axtb's accuracy model: how the errors of energies, gradients and Hessians
fall with the SCF accuracy factor, per method and conventions (design 7.4).

    tests/accuracy_calibration.py [--jobs N] [--out bench/results/accuracy.json] [--quick]
    tests/accuracy_calibration.py --fit bench/results/accuracy.json     print the constants

For each molecule, method and conventions, the SCF of a gradient runs at accuracy factors 10 ...
0.001 and is compared with one at the method's floor (an SCF converged as tightly as it goes):
max |dE| and the largest |dg| component. Each runs fresh, and warm: converged at the structure,
then started from that state at a 0.02-bohr displacement (sessions mostly warm-start). The
factor is set as a request's accuracy target sets it (the method built at accuracy_scale 1, then
`Calculator::set_scf_accuracy`), so GFN2-xTB's integral screening stays the default's and the
errors are the SCF's. Hessians: GFN2-xTB and g-xTB 2.0.1 analytic and response, v1.1 numerical,
on five molecules, at hessian.accuracy_scale 0.1 ... 1e-6, against a tight one (max |dH_ij|).

The fit (`fit`): per method, conventions kind (axtb or the program's) and quantity, the power law
err <= k acc^p as the upper envelope over molecules and starts (p from the slope, k the largest
err / acc^p), and the smallest factor still worth setting (the floor, where errors stop falling).
g-xTB v1.1's open shells do not fall monotonically: they get a table of the worst error at each
factor instead. The constants go into crates/axtb-py/python/axtb/_accuracy.py.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
S = ROOT / "tests/structures"
BOHR = 0.529177210544
ACCS = (10.0, 3.0, 1.0, 0.3, 0.1, 0.03, 0.01, 0.001)
HACCS = (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
# the factor of each method's reference SCF (its floor: tighter changes nothing but round-off)
REF = {"gfn2": 1e-6, "gxtb-v2.0.1": 1e-4, "gxtb-v1.1": 1e-5}
HREF = {"gfn2": 1e-8, "gxtb-v2.0.1": 1e-7, "gxtb-v1.1": 1e-6}
SYMBOLS = (
    "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
    "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
    "Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn"
).split()

# name: (file, charge, unpaired electrons)
MOLECULES = {
    "water": ("water.xyz", 0, 0),
    "caffeine": ("caffeine.xyz", 0, 0),
    "taxol": ("taxol.xyz", 0, 0),
    "ferrocene": ("ferrocene.xyz", 0, 0),
    "cu-nh3-4": ("cu-nh3-4.xyz", 2, 1),
    "o2-triplet": ("../../crates/gxtb/tests/reference_v2/reference/dioxygen_triplet/mol.xyz", 0, 2),
    "zncl4": ("zncl4.xyz", -2, 0),
    "water100": ("large/water100.xyz", 0, 0),
    "c50h102": ("large/c50h102.xyz", 0, 0),
    "cf3br-benzene": ("cf3br-benzene.xyz", 0, 0),
    "s229-cl": ("s229-cl.xyz", 0, 0),
    "h2o-3a": ("../scf_hard/scan-H2O-3.0.xyz", 0, 0),
    "f2-3a": ("../scf_hard/scan-F2-3.0.xyz", 0, 0),
    "f2-3.5a": ("../scf_hard/scan-F2-3.5.xyz", 0, 0),
    "mbheavy112": ("../scf_hard/random-mbheavy-112.xyz", 0, None),
    "anthracene2": ("anthracene2.xyz", 0, 0),
    "beni-acetone": ("beni-acetone.xyz", 0, 0),
    "gh106": ("gh106.xyz", 0, 0),
    "mb33": ("mb33.xyz", 0, None),
    "tyr-xab": ("tyr-xab.xyz", 0, 0),
}
QUICK = ("water", "caffeine", "cu-nh3-4", "o2-triplet", "beni-acetone")
HESSIAN_MOLECULES = ("water", "caffeine", "beni-acetone", "cu-nh3-4", "o2-triplet")
# (label, method, conventions (driver name), extra request fields)
SETUPS = (
    ("gfn2 axtb", "gfn2", "axtb", {}),
    ("gfn2 xtb-6.8.0", "gfn2", "xtb", {}),
    ("gxtb-v2.0.1 axtb", "gxtb-v2.0.1", "axtb", {}),
    ("gxtb-v2.0.1 gxtb-2.0.1", "gxtb-v2.0.1", "gxtb", {}),
    ("gxtb-v1.1 axtb", "gxtb-v1.1", "axtb", {}),
    ("gxtb-v1.1 gxtb-1.1.0", "gxtb-v1.1", "gxtb", {}),
    ("composite axtb", "gxtb-v2.0.1", "axtb", {"alpb": "water", "solvation_correction": "gfn2"}),
)
HESSIAN_SETUPS = (
    ("gfn2 analytic", "gfn2", "analytic"),
    ("gfn2 response", "gfn2", "response"),
    ("gxtb-v2.0.1 analytic", "gxtb-v2.0.1", "analytic"),
    ("gxtb-v2.0.1 response", "gxtb-v2.0.1", "response"),
    ("gxtb-v1.1 numerical", "gxtb-v1.1", "numerical"),
)


def molecule(name):
    f, charge, uhf = MOLECULES[name]
    lines = (S / f).read_text().split("\n")
    n = int(lines[0])
    atoms = [line.split() for line in lines[2 : 2 + n]]
    numbers = [SYMBOLS.index(a[0].capitalize()) for a in atoms]
    pos = np.array([[float(x) / BOHR for x in a[1:4]] for a in atoms])
    return numbers, pos, charge, uhf


def displaced(pos):
    """The structure moved by up to 0.02 bohr per coordinate (a deterministic pattern)."""
    n = len(pos)
    d = np.array([[((7 * i + 3 * c) % 5) - 2 for c in range(3)] for i in range(n)], dtype=float)
    return pos + 0.01 * d


def request(method, conventions, extra, charge, uhf):
    req = {
        "numbers": [],
        "positions": [],
        "method": method,
        "conventions": conventions,
        "charge": float(charge),
        "threads": 1,
        **extra,
    }
    if uhf is not None:
        req["uhf"] = uhf
    return json.dumps(req)


def errors(r, ref):
    return abs(r["energy"] - ref["energy"]), float(np.abs(r["gradient"] - ref["gradient"]).max())


def run_scf(job):
    """One molecule and setup: fresh and warm errors at every factor, the iterations."""
    from axtb import _axtb

    name, (label, method, conventions, extra) = job
    numbers, pos, charge, uhf = molecule(name)
    if method == "gxtb-v1.1" and any(z in range(59, 70) or z in (90, 91) or z > 92 for z in numbers):
        return []
    req = request(method, conventions, extra, charge, uhf)
    task = json.dumps({"numbers": [], "positions": [], "gradient": True})
    out = []
    t0 = time.perf_counter()
    try:
        ref = _axtb.compute_once(req, numbers, pos, task, accuracy=REF[method])
        x1 = displaced(pos)
        ref1 = _axtb.compute_once(req, numbers, x1, task, accuracy=REF[method])
    except _axtb.Error as e:
        return [{"molecule": name, "setup": label, "error": f"reference: {e.message}"}]
    open_shell = (sum(numbers) - charge) % 2 == 1 or (uhf or 0) > 0
    for acc in ACCS:
        row = {"molecule": name, "setup": label, "acc": acc, "natoms": len(numbers), "open_shell": open_shell}
        try:
            r = _axtb.compute_once(req, numbers, pos, task, accuracy=acc)
            row["dE"], row["dG"] = errors(r, ref)
            row["iterations"] = r["scf"]["iterations"] if "scf" in r else None
        except _axtb.Error as e:
            row["fresh_error"] = e.message
        try:
            eng = _axtb.Engine(req, numbers, pos)
            eng.compute(pos, task, "null", accuracy=acc)
            w = eng.compute(x1, task, "null", accuracy=acc)
            row["warm_dE"], row["warm_dG"] = errors(w, ref1)
            row["warm_iterations"] = w["scf"]["iterations"] if "scf" in w else None
            row["warm_start"] = w["scf"]["guess"] if "scf" in w else None
        except _axtb.Error as e:
            row["warm_error"] = e.message
        out.append(row)
    for row in out:
        row["seconds"] = (time.perf_counter() - t0) / len(out)
    return out


def run_hessian(job):
    from axtb import _axtb

    name, (label, method, hmethod) = job
    numbers, pos, charge, uhf = molecule(name)
    req = request(method, "axtb", {}, charge, uhf)

    def task(h):
        return json.dumps(
            {
                "numbers": [],
                "positions": [],
                "hessian": True,
                "hessian_method": hmethod,
                "hessian_accuracy": h,
                "vibrations": False,
            }
        )

    try:
        ref = _axtb.compute_once(req, numbers, pos, task(HREF[method]))["hessian"]["hessian"]
    except _axtb.Error as e:
        return [{"molecule": name, "setup": label, "error": f"reference: {e.message}"}]
    out = []
    for h in HACCS:
        row = {"molecule": name, "setup": label, "hacc": h}
        try:
            r = _axtb.compute_once(req, numbers, pos, task(h))["hessian"]["hessian"]
            row["dH"] = float(np.abs(r - ref).max())
        except _axtb.Error as e:
            row["error"] = e.message
        out.append(row)
    return out


# --------------------------------------------------------------------------- the fit

FLOOR_ABS = {"dE": 1e-13, "dG": 1e-10, "dH": 1e-9}


def kind_of(label):
    """(method, conventions kind) of a setup label."""
    method, conv = label.split(" ", 1)
    if method == "composite":
        return "composite", "axtb"
    return method, "axtb" if conv == "axtb" else "program"


def envelope(points, p):
    """The largest err / acc^p over (acc, err) points."""
    return max(e / a**p for a, e in points)


def fit_power(points):
    """p from the slope of log err on log acc (the errors above their noise floor), k the upper
    envelope at that p; p rounded to 0.25 and kept in [0.5, 2.5]."""
    pts = [(a, e) for a, e in points if e > 0]
    if len({a for a, _ in pts}) < 3:
        return None
    la = np.log10([a for a, _ in pts])
    le = np.log10([e for _, e in pts])
    p = float(np.polyfit(la, le, 1)[0])
    p = min(2.5, max(0.5, round(p * 4) / 4))
    return p, envelope(pts, p)


def other_states(rows):
    """The (molecule, setup, start) series that do not converge to the reference's state, which no
    threshold model bounds: a gradient error above 1e-6 Eh/bohr at the smallest factor, or one that
    jumps tenfold above 1e-6 at a smaller factor than a smaller error (a stretched bond's SCF
    reaching another state, the random heavy molecule's slow SCF). Reported, not fitted."""
    series = {}
    for r in rows:
        if "acc" in r:
            for start, f in (("fresh", "dG"), ("warm", "warm_dG")):
                if f in r:
                    series.setdefault((r["molecule"], r["setup"], start), {})[r["acc"]] = r[f]
    out = set()
    for key, errs in series.items():
        accs = sorted(errs, reverse=True)
        if errs[accs[-1]] > 1e-6:
            out.add(key)
            continue
        for i, a in enumerate(accs):
            if any(errs[b] > 1e-6 and errs[b] > 10 * errs[a] for b in accs[i + 1 :]):
                out.add(key)
                break
    return out


def fit(rows):
    """The model constants of a calibration's rows."""
    out = {"scf": {}, "hessian": {}, "tables": {}}
    excluded = other_states(rows)
    out["excluded"] = sorted(" ".join(k) for k in excluded)
    groups = {}
    for r in rows:
        if "acc" not in r:
            continue
        key = kind_of(r["setup"])
        shell = "open" if r["open_shell"] else "closed"
        for q, fields in (("energy", ("dE", "warm_dE")), ("gradient", ("dG", "warm_dG"))):
            for f in fields:
                start = "warm" if f.startswith("warm") else "fresh"
                if f in r and (r["molecule"], r["setup"], start) not in excluded:
                    groups.setdefault((key, shell, q), []).append((r["acc"], r[f], r["molecule"]))
    for ((method, conv), shell, q), pts in sorted(groups.items()):
        noise = FLOOR_ABS["dE" if q == "energy" else "dG"]
        usable = [(a, max(e, noise)) for a, e, _ in pts]
        name = f"{method} {conv} {shell}"
        if method == "gxtb-v1.1" and shell == "open":
            table = {}
            for a, e in usable:
                table[a] = max(table.get(a, 0.0), e)
            # the worst error at each factor, made monotone (an error at a factor bounds every
            # larger factor's)
            accs = sorted(table)
            worst = 0.0
            mono = {}
            for a in accs:
                worst = max(worst, table[a])
                mono[a] = worst
            out["tables"].setdefault(name, {})[q] = sorted(mono.items(), reverse=True)
            continue
        above = [(a, e) for a, e in usable if e > 10 * noise]
        res = fit_power(above) or (1.0 if q == "gradient" else 2.0, envelope(usable, 1.0 if q == "gradient" else 2.0))
        p, k = res
        # the floor: the smallest factor whose worst error is still above ten times the noise
        worst = {}
        for a, e in usable:
            worst[a] = max(worst.get(a, 0.0), e)
        floor = min((a for a in worst if worst[a] > 10 * noise), default=min(worst))
        out["scf"].setdefault(name, {})[q] = {"k": k, "p": p, "floor": floor, "points": len(usable)}
    hgroups = {}
    for r in rows:
        if "hacc" in r and "dH" in r:
            hgroups.setdefault(r["setup"], []).append((r["hacc"], max(r["dH"], FLOOR_ABS["dH"])))
    for label, pts in sorted(hgroups.items()):
        above = [(a, e) for a, e in pts if e > 10 * FLOOR_ABS["dH"]]
        p, k = fit_power(above) or (1.0, envelope(pts, 1.0))
        worst = {}
        for a, e in pts:
            worst[a] = max(worst.get(a, 0.0), e)
        floor = min((a for a in worst if worst[a] > 10 * FLOOR_ABS["dH"]), default=min(worst))
        out["hessian"][label] = {"k": k, "p": p, "floor": floor}
    return out


def up(x, digits=3):
    """`x` rounded up to `digits` significant digits (a rounded envelope stays an envelope)."""
    if x <= 0:
        return x
    e = math.floor(math.log10(x)) - digits + 1
    return math.ceil(x / 10**e - 1e-9) * 10**e


def write_constants(f, source, data):
    """Write the fitted constants into the CALIBRATION block of axtb/_accuracy.py."""
    cal = {
        "source": "bench/results/accuracy.json",
        "axtb": data.get("axtb"),
        "commit": data.get("commit"),
        "acc_floor": min(ACCS),
        "accs": list(ACCS),
        "hessian_accs": list(HACCS),
        "scf": {
            name: {q: {"k": float(f"{up(v['k']):.3g}"), "p": v["p"]} for q, v in qs.items()}
            for name, qs in f["scf"].items()
        },
        "tables": {
            name: {q: [[a, float(f"{up(e):.3g}")] for a, e in t] for q, t in qs.items()}
            for name, qs in f["tables"].items()
        },
        "hessian": {
            name: {"k": float(f"{up(v['k']):.3g}"), "p": v["p"], "floor": v["floor"]}
            for name, v in f["hessian"].items()
        },
    }
    path = ROOT / "crates/axtb-py/python/axtb/_accuracy.py"
    text = path.read_text()
    begin = "# BEGIN CALIBRATION (written by tests/accuracy_calibration.py --write; do not edit by hand)\n"
    end = "# END CALIBRATION"
    head, rest = text.split(begin, 1)
    _, tail = rest.split(end, 1)
    body = "CALIBRATION: dict[str, Any] = " + json.dumps(cal, indent=4).replace("null", "None") + "\n"
    path.write_text(head + begin + body + end + tail)
    print(f"wrote the constants of {source} into {path.relative_to(ROOT)}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--out", default=str(ROOT / "bench/results/accuracy.json"))
    ap.add_argument("--quick", action="store_true", help="five molecules")
    ap.add_argument("--fit", default="", help="fit the rows of a calibration file and print the constants")
    ap.add_argument("--write", action="store_true", help="with --fit: write the constants into axtb/_accuracy.py")
    a = ap.parse_args()
    if a.fit:
        data = json.loads(Path(a.fit).read_text())
        f = fit(data["rows"])
        if a.write:
            write_constants(f, a.fit, data)
        print(json.dumps(f, indent=1))
        return
    names = QUICK if a.quick else tuple(MOLECULES)
    jobs = [(n, s) for n in names for s in SETUPS]
    hjobs = [(n, s) for n in HESSIAN_MOLECULES for s in HESSIAN_SETUPS]
    # the large molecules first, so the pool ends with the short jobs
    jobs.sort(key=lambda j: -len(molecule(j[0])[0]))
    t0 = time.perf_counter()
    rows = []
    with ProcessPoolExecutor(a.jobs) as pool:
        futures = [pool.submit(run_scf, j) for j in jobs] + [pool.submit(run_hessian, j) for j in hjobs]
        for i, f in enumerate(futures):
            got = f.result()
            rows.extend(got)
            print(f"{i + 1}/{len(futures)} {time.perf_counter() - t0:.0f} s", file=sys.stderr, flush=True)
    import axtb

    doc = {
        "axtb": axtb.__version__,
        "commit": getattr(axtb._qc, "COMMIT", None),
        "accs": ACCS,
        "hessian_accs": HACCS,
        "reference": REF,
        "hessian_reference": HREF,
        "seconds": time.perf_counter() - t0,
        "jobs": a.jobs,
        "rows": rows,
        "fit": fit(rows),
    }
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(doc, indent=1, default=float))
    print(f"wrote {a.out}: {len(rows)} rows in {doc['seconds']:.0f} s")
    print(json.dumps(doc["fit"], indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
