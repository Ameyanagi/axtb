# SPDX-License-Identifier: LGPL-3.0-or-later
"""A flat 0.1 request (the CLI's JSON) as a `qcdata.ProgramInput` for ``axtb.compute`` and
``axtb.Session``: the tests', benchmarks' and scripts/golden_python.py's mapping. Private to
the repository, not an API.

`program_input(req)` maps every key that has a qcdata form (design 12.3):

- numbers and positions -> the structure (bohr); ``uhf`` -> the multiplicity (default by
  parity, as the driver does); ``charge``;
- ``gradient`` / ``hessian`` -> the calctype;
- each 0.1 setting -> its keyword (``etemp`` -> ``scf.etemp``, ``alpb`` -> ``solvent.name``,
  ``hessian_method`` -> ``hessian.method``, ``accuracy`` -> ``accuracy_scale``, ...);
- ``guess``: ``previous`` -> ``scf.restart: previous``, ``fresh`` -> ``off``, a start ->
  ``scf.guess``;
- ``threads`` stays explicit (default 1, the CLI's), so the bits are comparable.

A request with a key that has no qcdata form (``molden``, ``molden_file``, ``dispersion``,
``scf_trace``, ``analyze_hessian``, a negative ``uhf``) raises `NoQcdataForm`.
"""

from __future__ import annotations

from typing import Any

SYMBOLS = (
    "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y "
    "Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os "
    "Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr"
).split()

# 0.1 key -> keyword path
KEYWORD = {
    "accuracy": "accuracy_scale",
    "conventions": "conventions",
    "etemp": "scf.etemp",
    "max_iter": "scf.max_iter",
    "mixer": "scf.mixer",
    "restart": "scf.restart",
    "solvation_state": "solvent.state",
    "solvation_correction": "solvent.correction",
    "correction_conventions": "solvent.correction_conventions",
    "fragment_charges": "topology.fragment_charges",
    "topology_positions": "topology.geometry",
    "hessian_method": "hessian.method",
    "hessian_step": "hessian.step",
    "hessian_accuracy": "hessian.accuracy_scale",
    "analysis_conventions": "freq.conventions",
    "masses": "freq.masses",
    "temperature": "freq.temperature",
    "pressure": "freq.pressure",
    "rotor_cutoff": "freq.rotor_cutoff",
    "imaginary_cutoff": "freq.imaginary_cutoff",
    "frequency_scale": "freq.scale",
    "symmetry_number": "freq.symmetry_number",
    "symmetry_threshold": "freq.symmetry_tol",
}
NO_FORM = (
    "molden",
    "molden_file",
    "molden_cartesian",
    "dispersion",
    "scf_trace",
    "analyze_hessian",
    "dipole_derivatives",
    "energy",
)


class NoQcdataForm(ValueError):
    """The request has a key without a qcdata form."""


def _put(kw: dict[str, Any], path: str, value: Any) -> None:
    group, _, leaf = path.rpartition(".")
    (kw.setdefault(group, {}) if group else kw)[leaf] = value


def keywords(req: dict[str, Any]) -> dict[str, Any]:
    """The keywords of a flat request's settings (numbers, positions, charge, uhf, method and the
    calctype flags aside)."""
    # the CLI's flags default to false; any other key given counts (dispersion: false above all)
    flags = ("molden", "molden_cartesian", "scf_trace")
    bad = [k for k in NO_FORM if req.get(k) not in ((None, False) if k in flags else (None,))]
    if bad:
        raise NoQcdataForm(f"no qcdata form: {', '.join(bad)}")
    kw: dict[str, Any] = {"threads": req.get("threads", 1)}
    for k, path in KEYWORD.items():
        if req.get(k) is not None:
            _put(kw, path, req[k])
    for model in ("alpb", "gbsa"):
        if req.get(model) is not None:
            _put(kw, "solvent.name", req[model])
            if model == "gbsa":
                _put(kw, "solvent.model", "gbsa")
    guess = req.get("guess")
    if guess is not None:
        if guess in ("previous", "restart"):
            _put(kw, "scf.restart", "previous")
        elif guess == "fresh":
            _put(kw, "scf.restart", "off")
        else:
            _put(kw, "scf.guess", guess)
    if req.get("vibrations") is False:
        _put(kw, "freq.analysis", "off")
    if req.get("thermo") is False:
        _put(kw, "freq.thermo", "off")
    for out in ("bond_orders", "overlap", "hamiltonian"):
        if req.get(out):
            _put(kw, f"output.{out}", True)
    return kw


def calctype(req: dict[str, Any]) -> str:
    return "hessian" if req.get("hessian") else ("gradient" if req.get("gradient") else "energy")


def structure(req: dict[str, Any]):
    """The `qcdata.Structure` of a flat request (positions in bohr)."""
    from qcdata import Structure

    numbers = req["numbers"]
    charge = req.get("charge", 0)
    uhf = req.get("uhf")
    if uhf is not None and uhf < 0:
        raise NoQcdataForm("no qcdata form: a negative uhf (multiplicity < 1)")
    if charge != int(charge):
        raise NoQcdataForm("no qcdata form: a fractional charge")
    nel = sum(numbers) - int(charge)
    multiplicity = (uhf if uhf is not None else nel % 2) + 1
    return Structure(
        symbols=[SYMBOLS[z] for z in numbers], geometry=req["positions"], charge=int(charge), multiplicity=multiplicity
    )


def program_input(req: dict[str, Any]):
    """The `qcdata.ProgramInput` of a flat request."""
    from qcdata import ProgramInput

    return ProgramInput(
        program="axtb",
        calctype=calctype(req),
        structure=structure(req),
        model={"method": req.get("method") or "gfn2"},
        keywords=keywords(req),
    )
