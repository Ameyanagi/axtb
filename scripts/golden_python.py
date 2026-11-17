#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""The Python interface against the golden outputs of the CLI, bit for bit.

    scripts/golden_python.py compare DIR [--only a,b] [--quick] [--interleave]

DIR is a recording of scripts/golden.py from the same build and machine (golden files are
machine-specific). Every corpus case (scripts/golden.py's) is mapped to a `qcdata.ProgramInput`
(tests/flat.py: the flat request's settings as keywords, ``threads`` explicit) and run through
``axtb.compute``; the trajectory cases run through ``axtb.Session`` on the default lineage (the
request's ``guess`` as ``scf.restart``). Cases without a qcdata form are listed as skipped with
the reason (the analysis of a given Hessian, Molden text, switching off dispersion, a negative
``uhf``). Cases with ``AXTB_LAPACK`` run in a child process with the variable set (it is read
once per process).

Compared exactly (numpy's array_equal on float64; JSON floats read back exactly):

- energy, gradient, Hessian; atomic charges, dipole; the energy components and record's
  subtotals against ``terms``; SCF iterations; the molecular quadrupole in the conventions whose
  form is traceless; bond orders;
- frequencies, normal modes, reduced masses, force constants, IR intensities; every
  thermochemistry entry's enthalpy, entropy, Cp and Gibbs energy (E + G(RRHO)), the zero-point
  energy and each contribution;
- orbital energies and occupations, atomic dipoles, atomic quadrupoles (unpacked) where the
  output has the wavefunction;
- errors: the qcdata kind of the CLI's kind, and the message of an error of the driver.

``--quick`` runs the water and caffeine cases (CI). ``--interleave`` also runs every trajectory
pair as two lineages alternating in one session, which must equal the sequential frames. Exit
status 1 on any difference.
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "tests"))
import flat  # noqa: E402
from golden import corpus  # noqa: E402

# the qcdata kind of each CLI error kind
KINDS = {
    "input": "input",
    "unsupported_element": "input",
    "topology": "input",
    "scf_not_converged": "convergence",
    "numerical": "numerical",
    "not_applicable": "unsupported",
}
# conventions whose molecular quadrupole is tblite's form, not the traceless tensor
TBLITE_QUADRUPOLE = {("gfn2", "tblite"), ("gxtb-v2.0.1", "gxtb")}


class Diff(Exception):
    pass


def same(name, a, b):
    """Raise Diff unless a and b are equal bit for bit (arrays of float64)."""
    if a is None and b is None:
        return
    if a is None or b is None:
        raise Diff(f"{name}: one is missing ({a is None}, {b is None})")
    x, y = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if x.shape != y.shape:
        x = x.reshape(-1) if x.size == y.size else x
        y = y.reshape(-1) if x.size == y.size else y
    if x.shape != y.shape or not np.array_equal(x, y):
        d = np.abs(x - y).max() if x.shape == y.shape else "shape"
        raise Diff(f"{name}: differs (max |d| {d})")


def unpack(q):
    xx, xy, yy, xz, yz, zz = q
    return [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]


def compare_results(o, res, rec, method, conventions):
    """Compare a successful CLI response `o` with qcdata results `res` and record `rec`."""
    same("energy", res.energy, o["energy"])
    if "gradient" in o:
        same("gradient", res.gradient, o["gradient"])
    same("charges", res.atomic_charges, o["charges"])
    same("dipole", res.dipole_moment, o["dipole"])
    terms = dict(res.energy_components)
    terms.update(rec.get("energy_subtotals", {}))
    if set(terms) != set(o["terms"]):
        raise Diff(f"terms: {sorted(terms)} vs {sorted(o['terms'])}")
    for k, x in o["terms"].items():
        same(f"terms.{k}", terms[k], x)
    if "iterations" in o:
        if res.scf_iterations != o["iterations"]:
            raise Diff(f"iterations: {res.scf_iterations} vs {o['iterations']}")
        if (method, conventions) not in TBLITE_QUADRUPOLE:
            same("quadrupole", res.quadrupole_moment, unpack(o["quadrupole"]))
    if "bond_orders" in o:
        same("bond_orders", res.bond_orders, o["bond_orders"])
    if "hessian" in o:
        same("hessian", res.hessian, o["hessian"])
    if "frequencies" in o:
        va = res.vibrational_analysis
        same("frequencies", va.frequencies, o["frequencies"])
        same("normal_modes", va.normal_modes_cartesian, o["normal_modes"])
        same("reduced_masses", va.reduced_masses, o["reduced_masses"])
        same("force_constants", va.force_constants, o["force_constants"])
        same("ir_intensities", va.ir_intensities, o.get("ir_intensities"))
        same("masses", va.masses, o["masses"])
        t = o.get("thermo")
        if t is not None:
            if len(res.thermochemistry) != len(t["points"]):
                raise Diff("thermochemistry: entry count")
            for i, (e, pt) in enumerate(zip(res.thermochemistry, t["points"], strict=True)):
                same(f"thermo[{i}].temperature", e.temperature, pt["temperature"])
                same(f"thermo[{i}].enthalpy", e.total.enthalpy, pt["total_enthalpy"])
                same(f"thermo[{i}].entropy", e.total.entropy, pt["entropy"])
                same(f"thermo[{i}].cp", e.total.cp, pt["heat_capacity"])
                same(f"thermo[{i}].gibbs", e.total.gibbs, pt["free_energy"])
                same(f"thermo[{i}].zpve", e.zero_point_energy, t["zpve"])
                for part, key in (("translational", "translation"), ("rotational", "rotation")):
                    c = pt[key]
                    term = getattr(e, part)
                    same(
                        f"thermo[{i}].{part}",
                        [term.enthalpy, term.entropy, term.cp],
                        [c["enthalpy"], c["entropy"], c["heat_capacity"]],
                    )
                vib = pt["vibration"]
                same(
                    f"thermo[{i}].vibrational",
                    [e.vibrational.enthalpy, e.vibrational.entropy, e.vibrational.cp],
                    [t["zpve"] + vib["enthalpy"], vib["entropy"], vib["heat_capacity"]],
                )
    wf = res.wavefunction
    if wf is not None and "orbital_energies_alpha" in o:
        orb = wf.orbitals
        same("occupations_alpha", orb.occupations_a, o["occupations_alpha"])
        same("occupations_beta", orb.occupations_b, o["occupations_beta"])
        same("orbital_energies_alpha", orb.energies_a, o["orbital_energies_alpha"])
        if orb.energies_b is not None:
            same("orbital_energies_beta", orb.energies_b, o["orbital_energies_beta"])
        else:
            same("orbital_energies_beta", orb.energies_a, o["orbital_energies_beta"])
    if "atomic_dipoles" in o and res.atomic_dipoles is not None:
        same("atomic_dipoles", res.atomic_dipoles, o["atomic_dipoles"])
        same("atomic_quadrupoles", res.atomic_quadrupoles, [unpack(q) for q in o["atomic_quadrupoles"]])


def compare_error(gold_error, exc):
    """Compare the CLI's error with a ProgramError of the Python interface."""
    want = KINDS[gold_error["kind"]]
    got = exc.error.kind.value
    if got != want:
        raise Diff(f"error kind: {got} ({exc.error.code}) vs {want} ({gold_error['kind']}): {exc.error.message}")
    driver = exc.output is not None and exc.output.resolved is not None
    if driver and exc.error.message != gold_error["message"]:
        raise Diff(f"error message: {exc.error.message!r} vs {gold_error['message']!r}")
    return "driver" if driver else f"planning: {exc.error.code}"


def run_compute(req, gold):
    import axtb
    from qcdata import ProgramError

    inp = flat.program_input(req)
    try:
        method = axtb._keywords.canonical_method(inp.model.method)
    except ProgramError:
        method = inp.model.method
    conventions = req.get("conventions") or "axtb"
    try:
        out = axtb.compute(inp)
    except ProgramError as e:
        if gold["exit_code"] == 0:
            raise Diff(f"failed: {e.error.kind.value} {e.error.code}: {e.error.message}") from None
        return compare_error(gold["output"]["error"], e)
    if gold["exit_code"] != 0:
        raise Diff(f"succeeded where the CLI failed: {gold['output'].get('error')}")
    compare_results(gold["output"], out.results, out.record, method, flat_conventions(method, conventions))
    return ""


def flat_conventions(method, conventions):
    """The driver's name of the conventions (as TBLITE_QUADRUPOLE names them)."""
    key = "".join(c for c in conventions.lower() if c.isalnum())
    return {"xtb680": "xtb", "xtb671": "xtb671", "tblite080": "tblite", "gxtb110": "gxtb", "gxtb201": "gxtb"}.get(
        key, key
    )


def run_trajectory(req, gold, interleave):
    """A trajectory case through a session on the default lineage (and interleaved with its
    reverse in two lineages)."""
    import axtb
    from qcdata import ProgramError

    if not hasattr(axtb, "Session"):
        raise flat.NoQcdataForm("a trajectory: needs axtb.Session")
    traj = req["trajectory"]
    base = {k: v for k, v in req.items() if k != "trajectory"}
    inp = flat.program_input(base)
    method = axtb._keywords.canonical_method(inp.model.method)
    conventions = flat_conventions(method, req.get("conventions") or "axtb")
    outs = gold["output"] if isinstance(gold["output"], list) else [gold["output"]]
    with axtb.Session(inp) as s:
        for i, x in enumerate(traj):
            o = outs[i] if i < len(outs) else None
            try:
                f = s.compute(np.asarray(x), inp.calctype)
            except ProgramError as e:
                if o is None or "error" not in o:
                    raise Diff(f"step {i} failed: {e.error.message}") from None
                compare_error(o["error"], e)
                break
            if o is None or "error" in o:
                raise Diff(f"step {i} succeeded where the CLI failed")
            try:
                compare_results(o, f.results, f.record, method, conventions)
                guess = f.record["scf"]["start"] if "scf" in f.record else None
                if o.get("guess") is not None and guess != o["guess"]:
                    raise Diff(f"start {guess} vs {o['guess']}")
            except Diff as d:
                raise Diff(f"step {i}: {d}") from None
    if interleave and gold["exit_code"] == 0:
        rev = [list(reversed(traj)), traj]
        seq = []
        for t in rev:
            with axtb.Session(inp) as s:
                seq.append([s.compute(np.asarray(x), inp.calctype) for x in t])
        with axtb.Session(inp) as s:
            inter = [[], []]
            for a, b in zip(*rev, strict=True):
                inter[0].append(s.compute(np.asarray(a), inp.calctype, lineage="a"))
                inter[1].append(s.compute(np.asarray(b), inp.calctype, lineage=["b", 1]))
        for k in range(2):
            for i, (f, g) in enumerate(zip(seq[k], inter[k], strict=True)):
                try:
                    same("energy", f.results.energy, g.results.energy)
                    same("gradient", f.results.gradient, g.results.gradient)
                    same("charges", f.results.atomic_charges, g.results.atomic_charges)
                    if f.results.scf_iterations != g.results.scf_iterations:
                        raise Diff("iterations")
                except Diff as d:
                    raise Diff(f"interleaved lineage {k} step {i}: {d}") from None
        return "interleaved too"
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["compare"])
    ap.add_argument("dir")
    ap.add_argument("--only", default="", help="comma-separated case-name substrings")
    ap.add_argument("--exact", default="", help="one case by name (the child process of an AXTB_LAPACK case)")
    ap.add_argument("--quick", action="store_true", help="the water and caffeine cases")
    ap.add_argument("--interleave", action="store_true")
    a = ap.parse_args()
    d = Path(a.dir)
    cases = corpus()
    if a.exact:
        cases = [c for c in cases if c[0] == a.exact]
    if a.only:
        keys = a.only.split(",")
        cases = [c for c in cases if any(k in c[0] for k in keys)]
    if a.quick:
        cases = [c for c in cases if "water" in c[0] or "caffeine" in c[0]]
    counts = {"same": 0, "skip": 0, "diff": 0}
    for name, req, env in cases:
        gold = json.loads((d / f"{name}.json").read_text())
        if env and not a.exact:
            p = subprocess.run(
                [sys.executable, __file__, "compare", str(d), "--exact", name]
                + (["--interleave"] if a.interleave else []),
                env={**os.environ, **env},
                capture_output=True,
                text=True,
            )
            line = (p.stdout.strip().splitlines() or [f"DIFF {name}: child failed: {p.stderr[-400:]}"])[0]
            print(line, flush=True)
            counts["same" if line.startswith("same") else "skip" if line.startswith("skip") else "diff"] += 1
            continue
        try:
            if "trajectory" in req:
                how = run_trajectory(req, gold, a.interleave)
            else:
                how = run_compute(req, gold)
        except flat.NoQcdataForm as e:
            print(f"skip {name}: {e}", flush=True)
            counts["skip"] += 1
            continue
        except Diff as e:
            print(f"DIFF {name}: {e}", flush=True)
            counts["diff"] += 1
            continue
        print(f"same {name}" + (f" ({how})" if how else ""), flush=True)
        counts["same"] += 1
    if not a.exact:
        print(
            f"{counts['same']} of {len(cases)} identical, {counts['skip']} skipped (no qcdata form), "
            f"{counts['diff']} differences"
        )
    sys.exit(1 if counts["diff"] else 0)


if __name__ == "__main__":
    main()
