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

import json
import os
import subprocess
from pathlib import Path
from typing import Any

import numpy as np

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
    "scf_retries": "scf.retries",
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


def _electrons(numbers, method: str, charge: float) -> int:
    """The electrons the method's orbitals hold, as the driver counts them for the default uhf
    (valence-only methods leave the cores out); all of them where the table has no entry."""
    from axtb import _axtb

    try:
        v = _axtb.valence_electrons(method)
        counts = [v[z - 1] if 0 < z <= len(v) else None for z in numbers]
    except Exception:  # noqa: BLE001 - an unknown method: the input's refusal follows
        counts = [None]
    if None in counts:
        return sum(numbers) - int(charge)
    return int(round(sum(counts))) - int(charge)


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
    multiplicity = (uhf if uhf is not None else _electrons(numbers, req.get("method") or "gfn2", charge) % 2) + 1
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


# ----------------------------------------------------------------------------- the old calls

AXTB_BIN = os.environ.get("AXTB_BIN", str(Path(__file__).resolve().parents[1] / "target/release/axtb"))


def _arrays(out: dict[str, Any]) -> dict[str, Any]:
    """The CLI's numeric lists as numpy arrays (as axtb 0.1's Python results had them)."""
    for k, v in out.items():
        if isinstance(v, list) and v and k != "fermi_levels":
            try:
                a = np.asarray(v, dtype=float)
            except (TypeError, ValueError):
                continue
            out[k] = a
        elif isinstance(v, dict) and k in ("thermo",):
            _arrays(v)
    return out


def cli(req: dict[str, Any], env: dict[str, str] | None = None) -> dict[str, Any]:
    """The CLI's JSON response to a flat request (the interface of what has no qcdata form: Molden
    text, the analysis of a given Hessian, switching off dispersion), numeric lists as arrays.
    Raises `Failed` (`InputFailed` for input) with the CLI's error."""
    body = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in req.items()}
    p = subprocess.run(
        [AXTB_BIN, "--json"], input=json.dumps(body), capture_output=True, text=True, env={**os.environ, **(env or {})}
    )
    try:
        out = json.loads(p.stdout)
    except json.JSONDecodeError:
        raise Failed(p.stderr.strip() or f"the CLI exited with {p.returncode}", "numerical") from None
    if p.returncode != 0 or "error" in out:
        err = out.get("error", {})
        kind = err.get("kind", "numerical")
        cls = InputFailed if kind in INPUT_KINDS else Failed
        raise cls(err.get("message", p.stderr), kind, {k: v for k, v in err.items() if k not in ("kind", "message")})
    return _arrays(out)


def cli_singlepoint(numbers, positions, **settings) -> dict[str, Any]:
    """axtb 0.1's `singlepoint` through the CLI (what Python no longer gives: the conventions' own
    quadrupole, Molden text)."""
    return cli(_request(numbers, positions, settings))


def analysis(numbers, positions, hessian, **settings) -> dict[str, Any]:
    """The CLI's vibrational analysis of a given Hessian (`analyze_hessian`; qcinf's after alpha 1)."""
    req = {
        "numbers": [int(z) for z in numbers],
        "positions": np.asarray(positions, dtype=float).tolist(),
        "analyze_hessian": np.asarray(hessian, dtype=float).tolist(),
    }
    for k, v in settings.items():
        if v is not None:
            req["analysis_conventions" if k == "conventions" else k] = v
    return cli(req)


def _pack(t: np.ndarray) -> np.ndarray:
    """N x 3 x 3 (or 3 x 3) tensors packed as xx, xy, yy, xz, yz, zz."""
    t = np.asarray(t, dtype=float)
    return np.stack([t[..., 0, 0], t[..., 0, 1], t[..., 1, 1], t[..., 0, 2], t[..., 1, 2], t[..., 2, 2]], axis=-1)


def as_dict(
    results,
    record: dict[str, Any],
    timings: dict[str, float] | None = None,
    numbers=None,
    settings: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """A result of axtb.compute or axtb.Session as axtb 0.1's dict (the CLI's keys), for the scripts
    and tests written against it."""
    from axtb import _axtb

    r = results
    out: dict[str, Any] = {
        "method": record.get("method"),
        "conventions": record.get("conventions"),
        "energy": r.energy,
        "charges": np.asarray(r.atomic_charges),
        "dipole": np.asarray(r.dipole_moment),
        "terms": {**r.energy_components, **record.get("energy_subtotals", {})},
        "timings": {k: v for k, v in (timings or {}).items() if k != "total" and not k.startswith("hessian:")},
    }
    if r.gradient is not None:
        out["gradient"] = np.asarray(r.gradient)
    scf = record.get("scf")
    if scf is not None:
        out["iterations"] = scf["iterations"]
        out["guess"] = scf["start"]
        out["mixer_used"] = scf.get("mixer_used")
        out["scf_attempts"] = scf.get("attempts")
        out["eigensolver"] = record.get("eigensolver")
    if r.shell_populations is not None and record.get("method") and numbers is not None:
        layout = _axtb.shell_layout(record["method"], list(numbers))
        if layout is not None and len(layout) == len(r.shell_populations):
            out["shell_charges"] = np.array([x[2] for x in layout]) - np.asarray(r.shell_populations)
    if r.atomic_dipoles is not None:
        out["atomic_dipoles"] = np.asarray(r.atomic_dipoles)
        out["atomic_quadrupoles"] = _pack(r.atomic_quadrupoles)
    if r.quadrupole_moment is not None:
        out["quadrupole"] = _pack(r.quadrupole_moment)
    wf = r.wavefunction
    if wf is not None:
        o = wf.orbitals
        out["orbital_energies_alpha"] = np.asarray(o.energies_a)
        out["orbital_energies_beta"] = np.asarray(o.energies_b if o.energies_b is not None else o.energies_a)
        out["occupations_alpha"] = np.asarray(o.occupations_a)
        out["occupations_beta"] = np.asarray(o.occupations_b)
        out["electronic_temperature"] = o.electronic_temperature
        out["fermi_levels"] = [o.fermi_level_a, o.fermi_level_b]
        out["homo_lumo_gap"] = _gap(o, record.get("conventions") not in (None, "axtb"))
    if r.bond_orders is not None:
        out["bond_orders"] = np.asarray(r.bond_orders)
    if r.hessian is not None:
        out["hessian"] = np.asarray(r.hessian)
        h = record.get("hessian", {})
        out["hessian_method"] = h.get("method")
        out["hessian_max_asymmetry"] = h.get("max_asymmetry")
    va = r.vibrational_analysis
    if va is not None:
        out["frequencies"] = np.asarray(va.frequencies)
        out["normal_modes"] = np.asarray(va.normal_modes_cartesian)
        out["reduced_masses"] = np.asarray(va.reduced_masses)
        out["force_constants"] = np.asarray(va.force_constants)
        out["masses"] = np.asarray(va.masses)
        out["vibrational_conventions"] = ((settings or {}).get("freq") or {}).get("conventions")
        out["linear"] = va.rotor_type is not None and va.rotor_type.value == "linear"
        out["n_imaginary"] = va.n_imaginary(20.0)
        if va.ir_intensities is not None:
            out["ir_intensities"] = np.asarray(va.ir_intensities)
            out["dipole_derivatives"] = np.asarray(va.dipole_derivatives)
    if r.thermochemistry:
        freq = (settings or {}).get("freq") or {}
        last = _thermo_point(r.thermochemistry[-1])
        out["thermo"] = {
            **last,
            "zpve": r.thermochemistry[-1].zero_point_energy,
            "pressure": r.thermochemistry[-1].pressure,
            "rotor_cutoff": r.thermochemistry[-1].qrrho_frequency,
            "frequency_scale": r.thermochemistry[-1].frequency_scale_factor,
            "imaginary_cutoff": freq.get("imaginary_cutoff"),
            "point_group": va.point_group if va is not None else None,
            "symmetry_number": r.thermochemistry[-1].symmetry_number,
            "linear": out.get("linear"),
            "rotational_constants": [x or 0.0 for x in (va.rotational_constants or [])] if va is not None else None,
            "points": [_thermo_point(t) for t in r.thermochemistry],
        }
    comp = record.get("composite")
    if comp is not None:
        out["composite"] = comp["definition"]
        out["components"] = comp["components"]
        out["parts"] = comp["parts"]
        if "hessian_parts" in comp:
            out["hessian_parts"] = comp["hessian_parts"]
    g = record.get("gfnff")
    if g is not None:
        out["n_hbonds"], out["n_xbonds"] = g["n_hbonds"], g["n_xbonds"]
        if "topology" in g:  # numeric lists as arrays, as 0.1's
            out["topology"] = _arrays(dict(g["topology"]))
    return out


def _thermo_point(t) -> dict[str, Any]:
    """A qcdata `Thermochemistry` entry as one of axtb 0.1's thermo points (energies relative to
    the electronic energy; the vibrational term without the zero-point energy)."""
    e, zpve = t.electronic_energy, t.zero_point_energy
    enthalpy = t.total.enthalpy - e

    def term(x, less=0.0):
        return {"enthalpy": x.enthalpy - less, "entropy": x.entropy, "heat_capacity": x.cp}

    return {
        "temperature": t.temperature,
        "enthalpy": enthalpy,
        "thermal_enthalpy": enthalpy - zpve,
        "entropy": t.total.entropy,
        "ts": t.temperature * t.total.entropy,
        "g_rrho": t.total.gibbs - e,
        "heat_capacity": t.total.cp,
        "translation": term(t.translational),
        "rotation": term(t.rotational),
        "vibration": term(t.vibrational, zpve),
        "free_energy": t.total.gibbs,
        "total_enthalpy": t.total.enthalpy,
    }


def _request(numbers, positions, settings: dict[str, Any]) -> dict[str, Any]:
    req = {"numbers": [int(z) for z in numbers], "positions": np.asarray(positions, dtype=float).tolist()}
    for k, v in settings.items():
        if v is not None:
            req[k] = v
    if "charge" in req:
        req["charge"] = float(req["charge"])
    return req


INPUT_KINDS = ("input", "unsupported_element", "topology", "not_applicable")


class Failed(RuntimeError):
    """A refused or failed calculation of the old calls (a `RuntimeError`, as axtb 0.1's
    numerical and SCF errors were): ``kind`` is axtb 0.1's kind (``input``,
    ``unsupported_element``, ``topology``, ``scf_not_converged``, ``numerical``,
    ``not_applicable``), the error's details are attributes, and ``error`` is the
    `qcdata.ErrorInfo` (None from the CLI)."""

    CODES = {
        "unsupported_element": "unsupported_element",
        "topology_failed": "topology",
        "scf_not_converged": "scf_not_converged",
        "numerical_failure": "numerical",
        "not_applicable": "not_applicable",
    }

    def __init__(self, message: str, kind: str, details: dict[str, Any] | None = None, error: Any = None):
        super().__init__(message)
        self.kind = kind
        self.error = error
        for k, v in (details or {}).items():
            setattr(self, k, v)


class InputFailed(Failed, ValueError):
    """A refused input (a `ValueError` too, as axtb 0.1's ``InputError``)."""


def _guarded(f, *args, **kwargs):
    from qcdata import ProgramError

    try:
        return f(*args, **kwargs)
    except ProgramError as e:
        kind = Failed.CODES.get(e.error.code or "", "input")
        cls = InputFailed if kind in INPUT_KINDS else Failed
        raise cls(e.error.message, kind, getattr(e.error, "details", None), e.error) from e


def _gap(o, alpha_rule: bool) -> float | None:
    """The HOMO-LUMO gap of the orbital levels, as the driver's (crates/axtb/src/properties.rs):
    one set of orbitals, or per spin (the alpha channel's under a program's conventions, else the
    lowest LUMO less the highest HOMO)."""

    def homo(nel: float) -> int | None:
        n = np.floor(nel)
        k = n + 1 if nel - n > 0.5 else n
        return int(k) - 1 if k >= 1 else None

    def channel(e, nel):
        h = homo(nel)
        return (e[h], e[h + 1]) if h is not None and h + 1 < len(e) else None

    na, nb = float(np.sum(o.occupations_a)), float(np.sum(o.occupations_b))
    ea = np.asarray(o.energies_a)
    if o.energies_b is None:
        c = channel(ea, max(na, nb))
        return None if c is None else c[1] - c[0]
    eb = np.asarray(o.energies_b)
    ga = channel(ea, na)
    if alpha_rule:
        return None if ga is None else ga[1] - ga[0]
    gb = channel(eb, nb)
    if ga is not None and gb is not None:
        return min(ga[1], gb[1]) - max(ga[0], gb[0])
    if ga is None and gb is None:
        return None
    h, lumo = ga if ga is not None else gb
    other = ea if na < 1 else eb
    return min(lumo, other[0]) - h


def run(req: dict[str, Any]) -> dict[str, Any]:
    """A flat request through axtb.compute, as axtb 0.1's dict (Molden text and the other keys
    without a qcdata form: through the CLI)."""
    import axtb

    try:
        inp = _guarded(program_input, req)
    except NoQcdataForm:
        return cli(req)
    out = _guarded(axtb.compute, inp)
    return as_dict(out.results, out.record, out.execution.timings, inp.structure.atomic_numbers, out.resolved.settings)


def singlepoint(numbers, positions, **settings) -> dict[str, Any]:
    """axtb 0.1's `singlepoint` (its keyword arguments) through axtb.compute."""
    molden = settings.pop("molden", None)
    req = _request(numbers, positions, settings)
    if molden:
        req["molden"] = True
        if molden == "cartesian":
            req["molden_cartesian"] = True
    return run(req)


class Calculator:
    """axtb 0.1's `Calculator` on `axtb.Session`s: one per Hessian setting, each on its default
    lineage (warm starts between calls as before; a Hessian's session apart from the single
    points'). `set_guess` starts the next calls from given moments (`given`)."""

    def __init__(self, method, numbers, positions, **settings):
        self._settings = {"method": method, **settings}
        guess = settings.get("guess")
        if guess not in (None, "previous", "restart", "fresh") and settings.get("restart") is None:
            # 0.1: a named start meant that start at every geometry (0.2 keeps the restart policy)
            self._settings["restart"] = "off"
        self._numbers = [int(z) for z in numbers]
        self._positions = np.asarray(positions, dtype=float)
        self._sessions: dict[Any, Any] = {}
        self._given = None
        _guarded(self._session, None)  # refusals at construction, as 0.1's

    def _session(self, hessian: dict[str, Any] | None):
        import axtb

        key = None if hessian is None else json.dumps(hessian, sort_keys=True)
        if key not in self._sessions:
            req = _request(self._numbers, self._positions, {**self._settings, **(hessian or {})})
            if hessian is not None:
                req["hessian"] = True
            inp = program_input(req)
            if self._given is not None:
                inp = inp.model_copy(update={"given": self._given})
            self._sessions[key] = (axtb.Session(inp), list(inp.structure.atomic_numbers))
        return self._sessions[key]

    def _call(self, positions, calctype, hessian=None, outputs=()):
        session, numbers = _guarded(self._session, hessian)
        f = _guarded(session.compute, np.asarray(positions, dtype=float), calctype, outputs=outputs)
        return as_dict(f.results, f.record, f.execution.timings, numbers, session.plan()["resolved"]["settings"])

    def singlepoint(self, positions, gradient=False, bond_orders=False, molden=None):
        if molden:
            raise NoQcdataForm("Molden text comes from the CLI (flat.cli)")
        return self._call(
            positions, "gradient" if gradient else "energy", outputs=("bond_orders",) if bond_orders else ()
        )

    def hessian(self, positions, step=None, accuracy=None, method=None, **analysis):
        h = {
            k: v
            for k, v in (("hessian_step", step), ("hessian_accuracy", accuracy), ("hessian_method", method))
            if v is not None
        }
        if analysis.pop("vibrations", True) is False:
            h["vibrations"] = False
        if analysis.pop("thermo", True) is False:
            h["thermo"] = False
        h.update({k: v for k, v in analysis.items() if v is not None})
        return self._call(positions, "hessian", hessian=h)

    def set_guess(self, shell_charges, atomic_dipoles, atomic_quadrupoles):
        """Start the next calls from given GFN2-xTB moments (packed quadrupoles, as the CLI's)."""
        from qcdata import Provenance, SinglePointData

        from axtb import _axtb

        if _axtb.canonical("method", self._settings["method"]) != "gfn2":
            raise InputFailed("set_guess: given moments are GFN2-only (a session takes any given start)", "input")
        layout = _axtb.shell_layout(self._settings["method"], self._numbers)
        q = np.asarray(atomic_quadrupoles, dtype=float)
        tensors = np.array([[[x[0], x[1], x[3]], [x[1], x[2], x[4]], [x[3], x[4], x[5]]] for x in q])
        self._given = SinglePointData(
            provenance=Provenance(program="axtb"),
            shell_populations=[ref - s for (_, _, ref, _), s in zip(layout, shell_charges, strict=True)],
            shell_atoms=[a for a, *_ in layout],
            shell_angular_momenta=[ell for _, ell, *_ in layout],
            atomic_dipoles=np.asarray(atomic_dipoles, dtype=float),
            atomic_quadrupoles=tensors,
        )
        for s, _ in self._sessions.values():
            s.close()
        self._sessions.clear()

    def reset_guess(self):
        for s, _ in self._sessions.values():
            s.reset()
