# SPDX-License-Identifier: LGPL-3.0-or-later
"""``axtb.compute``: one calculation of a `qcdata.ProgramInput`, and the mapping of the
driver's response onto qcdata (`SinglePointData`, ``record``, ``execution``).

The input is planned (`axtb._keywords.plan_input`), its threads chosen (`threads`), and the
driver runs it as the CLI does (`axtb._axtb.compute_once`: the same calculator, threads rule and
numbers, bit for bit). The response becomes qcdata's typed fields (`results`), axtb's own detail
goes in ``record`` (`record`), phase timings in ``execution.timings``. A refused input or a
failed calculation raises `qcdata.ProgramError` with the whole failed output attached.
"""

from __future__ import annotations

import json
import os
import time
from collections.abc import Mapping
from typing import Any

import numpy as np
from qcdata import (
    AOConvention,
    BasisSet,
    ExecutionInfo,
    OrbitalKind,
    Orbitals,
    OrbitalSpin,
    ProgramError,
    ProgramInput,
    ProgramOutput,
    ProgramWarning,
    Provenance,
    Resolved,
    Shell,
    SinglePointData,
    Thermochemistry,
    ThermoTerm,
    VibrationalAnalysis,
    Wavefunction,
)

from . import _axtb
from ._errors import from_driver
from ._models import build
from ._keywords import Planned, plan_input, plan_threads

__all__ = ["COMMIT", "VERSION", "compute", "provenance", "record", "results", "static_results", "task", "threads"]


def _version() -> str:
    try:
        from importlib.metadata import version

        return version("axtb")
    except Exception:  # noqa: BLE001 - not installed (a build tree): the crate's version
        return str(_axtb.__version__)


def _commit() -> str | None:
    try:
        from ._build import COMMIT  # written by build.rs

        return COMMIT
    except ImportError:
        return None


#: axtb's version as pip reports it (PEP 440: "0.2.0a1"), and the git commit of the build.
VERSION = _version()
COMMIT = _commit()
#: The output type of a calculation. Outputs, results and their parts are built from trusted
#: values of the right shapes without validation (`model_construct`: validation costs more than
#: a small molecule's whole calculation); the tests validate them again.
_OUTPUT = ProgramOutput[ProgramInput, SinglePointData]


_build = build


def provenance() -> Provenance:
    """axtb's `Provenance`: program, version and (from a git checkout) commit."""
    return _build(Provenance, program="axtb", program_version=VERSION, program_commit=COMMIT)


def _cpus() -> int:
    try:
        return len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return os.cpu_count() or 1


def threads(p: Planned, largest: str, ctx: Any) -> tuple[int, list[ProgramWarning], str]:
    """The threads of a calculation (or a session, by its largest calctype), the warnings and
    the reason: the ``threads`` keyword, at most the grant (``resource_capped`` above it); without
    the keyword the grant; without a grant (a direct call) the plan's threads (None: this
    machine's CPUs)."""
    grant = getattr(ctx, "threads", None) if ctx is not None else None
    if p.threads is not None:
        if grant is not None and p.threads > grant:
            w = ProgramWarning(
                code="resource_capped",
                message=f"threads {p.threads} capped to the {grant} granted",
                path="threads",
            )
            return grant, [w], f"the keyword ({p.threads}) capped to the grant"
        return p.threads, [], "the keyword"
    if grant is not None:
        return grant, [], "the runner's grant"
    planned = plan_threads(p, largest)
    if planned is None:
        return _cpus(), [], "a direct call: this machine's CPUs (the plan asks for all)"
    return planned, [], f"a direct call: the plan's {planned}"


def implied(outputs: frozenset[str] | set[str]) -> frozenset[str]:
    """The outputs with what they imply (decision 25a): ``hamiltonian`` brings ``overlap`` and
    ``orbitals``; ``overlap`` alone brings the basis set and no coefficients."""
    out = set(outputs)
    if "hamiltonian" in out:
        out |= {"overlap", "orbitals"}
    return frozenset(out)


def needs_orbitals(outputs: frozenset[str]) -> bool:
    """Whether the engine must return the structured orbitals (for the coefficients or the basis
    set of an AO matrix)."""
    return bool(outputs & {"orbitals", "overlap", "hamiltonian"})


def task(p: Planned, calctype: str, outputs: frozenset[str]) -> str:
    """The driver task (JSON) of a calctype with outputs."""
    t: dict[str, Any] = {"numbers": [], "positions": []}
    if calctype == "gradient":
        t["gradient"] = True
    elif calctype == "hessian":
        t["hessian"] = True
        t.update(p.task)
    if "bond_orders" in outputs:
        t["bond_orders"] = True
    if "overlap" in outputs:
        t["overlap"] = True
    if "hamiltonian" in outputs:
        t["hamiltonian"] = True
    return json.dumps(t, separators=(",", ":"))


def _thermo(r: Mapping[str, Any], p: Planned) -> tuple[dict[str, Any], list[Thermochemistry]]:
    """The VibrationalAnalysis fields of the thermochemistry, and its entries per temperature."""
    v = r["vibrations"]
    t = v.get("thermo")
    if t is None:
        return {}, []
    energy = r["energy"]
    rot = [c if c != 0.0 else None for c in t["rotational_constants"]]
    fields = {
        "rotational_constants": rot,
        "point_group": t["point_group"],
        "symmetry_number": t["symmetry_number"],
        "symmetry_number_source": "input" if "freq.symmetry_number" in p.given else "detected",
    }
    zpve = t["zpve"]
    entries = []
    for pt in t["points"]:
        tr, ro, vi = pt["translation"], pt["rotation"], pt["vibration"]
        entries.append(
            Thermochemistry(
                temperature=pt["temperature"],
                pressure=t["pressure"],
                electronic_energy=energy,
                zero_point_energy=zpve,
                translational=ThermoTerm(enthalpy=tr["enthalpy"], entropy=tr["entropy"], cp=tr["heat_capacity"]),
                rotational=ThermoTerm(enthalpy=ro["enthalpy"], entropy=ro["entropy"], cp=ro["heat_capacity"]),
                # qcdata's vibrational term includes the zero-point energy, axtb's does not
                vibrational=ThermoTerm(enthalpy=zpve + vi["enthalpy"], entropy=vi["entropy"], cp=vi["heat_capacity"]),
                total=ThermoTerm(
                    enthalpy=pt["total_enthalpy"],
                    entropy=pt["entropy"],
                    cp=pt["heat_capacity"],
                    gibbs=pt["free_energy"],
                ),
                symmetry_number=t["symmetry_number"],
                electronic_degeneracy=1.0,
                entropy_model="qrrho_grimme",
                enthalpy_model="rrho",
                qrrho_frequency=t["rotor_cutoff"],
                qrrho_alpha=4.0,
                qrrho_bav=t["average_moment"],
                frequency_scale_factor=t["frequency_scale"],
                excluded_imaginary=list(t["dropped_frequencies"]),
                inverted_imaginary=list(t["inverted_frequencies"]),
            )
        )
    return fields, entries


def static_results(r: Mapping[str, Any], p: Planned) -> dict[str, Any]:
    """The fields of `results` that one molecule and method keep at every geometry (a session
    computes them once): provenance, the atom, basis, orbital and electron counts, the charges'
    kind."""
    d: dict[str, Any] = {
        "provenance": provenance(),
        "calcinfo_natoms": len(p.numbers),
        "atomic_charges_kind": "mulliken" if p.scf else "eeq",
    }
    scf = r.get("scf")
    if scf is not None:
        d["calcinfo_nmo"] = len(scf["orbital_energies_alpha"])
        d["calcinfo_nbasis"] = sum(2 * ell + 1 for ell in scf["shell_angular_momenta"])
        for spin, name in (("alpha", "calcinfo_nalpha"), ("beta", "calcinfo_nbeta")):
            n = float(np.sum(scf[f"occupations_{spin}"]))
            if abs(n - round(n)) <= 1e-6:
                d[name] = int(round(n))
    return d


#: axtb's AO order: atom by atom, shell by shell, m = -l..l for every l (p as y, z, x); every shell
#: spherical, each function normalized as qcdata's default says (the x^l component's norm)
AO_CONVENTION = AOConvention(spherical_order="increasing_m")


def basis(o: Mapping[str, Any]) -> BasisSet:
    """The basis set of the engine's structured orbitals: one spherical shell per axtb shell, atom
    by atom, shell by shell. axtb's contraction coefficients multiply unnormalized primitives and
    normalize the x^l component, which is qcdata's `Shell` convention."""
    shells = [
        build(Shell, center=atom, angular_momentum=int(ell), pure=True, exponents=a, coefficients=c)
        for atom, per_atom in enumerate(o["shells"])
        for ell, a, c in per_atom
    ]
    return build(BasisSet, shells=shells, ao_convention=AO_CONVENTION)


def wavefunction(r: Mapping[str, Any], p: Planned, outputs: frozenset[str]) -> Wavefunction | None:
    """The `Wavefunction` of an SCF response: the orbital levels always (energies, occupations per
    spin in 0..1, the electronic temperature and Fermi levels, the core electrons each atom keeps
    outside the valence-only orbitals); the coefficients (``orbitals``), the overlap and the
    Hamiltonian on request, with the basis set they share. GFN2-xTB is restricted (one set of
    orbitals occupied per spin, open shells included), g-xTB restricted for a closed shell and
    unrestricted for an open one; the composite's are its g-xTB part's.

    None where qcdata cannot hold the levels: smeared occupations of a spin without electrons have
    no Fermi level (a warning says so)."""
    scf = r["scf"]
    restricted = p.method == "gfn2" or p.multiplicity == 1
    temperature = float(scf["electronic_temperature"])
    levels = scf["fermi_levels"]
    if temperature > 0 and None in levels:
        return None
    orb: dict[str, Any] = {
        "spin": OrbitalSpin.restricted if restricted else OrbitalSpin.unrestricted,
        "kind": OrbitalKind.canonical,
        "energies_a": scf["orbital_energies_alpha"],
        "occupations_a": scf["occupations_alpha"],
        "occupations_b": scf["occupations_beta"],
        "electronic_temperature": temperature,
        "fermi_level_a": levels[0] if temperature > 0 else None,
        "fermi_level_b": levels[1] if temperature > 0 else None,
        "core_electrons": scf["core_electrons"],
    }
    if not restricted:
        orb["energies_b"] = scf["orbital_energies_beta"]
    wf: dict[str, Any] = {}
    if needs_orbitals(outputs):
        o = r["orbitals"]
        wf["basis"] = basis(o)
        if "orbitals" in outputs:
            sets = o["sets"]
            orb["coefficients_a"] = sets[0]["coefficients"].T
            if not restricted:
                orb["coefficients_b"] = sets[1]["coefficients"].T
    if "overlap" in outputs:
        wf["overlap"] = scf["overlap"]
    if "hamiltonian" in outputs:
        h = scf["hamiltonian"]
        wf["hamiltonian_a"] = h[0]
        if not restricted:
            wf["hamiltonian_b"] = h[1]
    return build(Wavefunction, orbitals=build(Orbitals, **orb), **wf)


def results(
    r: Mapping[str, Any],
    p: Planned,
    calctype: str,
    outputs: frozenset[str],
    static: Mapping[str, Any] | None = None,
) -> SinglePointData:
    """The `SinglePointData` of a response (`static`: its `static_results`, if known)."""
    nat = len(p.numbers)
    d: dict[str, Any] = dict(static if static is not None else static_results(r, p))
    d["provenance"] = provenance()  # one per result (models are mutable)
    d["energy"] = r["energy"]
    d["dipole_moment"] = r["dipole"].tolist()
    d["atomic_charges"] = r["charges"]
    subtotals = set(r["subtotals"])
    d["energy_components"] = {k: x for k, x in r["terms"].items() if k not in subtotals}
    if calctype != "energy":
        d["gradient"] = r["gradient"]
    scf = r.get("scf")
    if scf is not None:
        d["scf_iterations"] = scf["iterations"]
        d["quadrupole_moment"] = scf["quadrupole_tensor"]
        if "bond_orders" in scf and "bond_orders" in outputs:
            d["bond_orders"] = scf["bond_orders"]
            d["bond_orders_kind"] = "mayer"
        d["shell_populations"] = scf["shell_populations"]
        d["shell_atoms"] = scf["shell_atoms"]
        d["shell_angular_momenta"] = scf["shell_angular_momenta"]
        d["atomic_dipoles"] = scf["atomic_dipoles"]
        d["atomic_quadrupoles"] = scf["atomic_quadrupole_tensors"]
        wf = wavefunction(r, p, outputs)
        if wf is not None:
            d["wavefunction"] = wf
    h = r.get("hessian")
    if h is not None:
        d["hessian"] = h["hessian"]
    v = r.get("vibrations")
    if v is not None:
        fields, entries = _thermo(r, p)
        va: dict[str, Any] = {
            "masses": v["masses"],
            "projection": "translations_rotations",
            "frequencies": v["frequencies"],
            "normal_modes_cartesian": v["normal_modes"],
            "reduced_masses": v["reduced_masses"],
            "force_constants": v["force_constants"],
            "rotor_type": "atom" if nat == 1 else ("linear" if v["linear"] else "nonlinear"),
            **fields,
        }
        if "ir_intensities" in v:
            va["ir_intensities"] = v["ir_intensities"]
            va["dipole_derivatives"] = h["dipole_derivatives"]
        d["vibrational_analysis"] = VibrationalAnalysis(**va)
        d["thermochemistry"] = entries
    # trusted values of the right shapes: built without validation (tests revalidate them)
    return _build(SinglePointData, **d)


def record(r: Mapping[str, Any] | None, p: Planned, *, session: bool = False) -> dict[str, Any]:
    """axtb's ``record`` of a response (`None`: a failure before the driver answered)."""
    rec: dict[str, Any] = {
        "method": p.method,
        "conventions": p.conventions,
        "numerics": r["numerics"] if r else _axtb.numerics(p.method, p.composite),
    }
    if r is None:
        return rec
    if r.get("eigensolver"):
        rec["eigensolver"] = r["eigensolver"]
    scf = r.get("scf")
    if scf is not None:
        s: dict[str, Any] = {"start": scf["guess"], "iterations": scf["iterations"]}
        if scf["mixer_used"] is not None:
            s["mixer_used"] = scf["mixer_used"]
        if scf["scf_attempts"]:
            s["attempts"] = scf["scf_attempts"]
        rec["scf"] = s
    subtotals = {k: r["terms"][k] for k in r["subtotals"]}
    if subtotals:
        rec["energy_subtotals"] = subtotals
    contingent: dict[str, Any] = {}
    if "scf.mixer" in p.contingent and scf is not None:
        contingent["scf.mixer"] = {"fired": scf["mixer_used"] == "diis-fallback"}
    comp = r.get("composite")
    if comp is not None:
        rec["composite"] = {
            "definition": comp["composite"],
            "components": comp["components"],
            "parts": comp["parts"],
        }
        if comp.get("hessian_parts") is not None:
            rec["composite"]["hessian_parts"] = comp["hessian_parts"]
        rec["wavefunction"] = {"part": "gxtb"}
        retried = [k for k, part in comp["parts"].items() if part.get("retry")]
        contingent["composite.retries"] = {"fired": bool(retried), "parts": retried}
    if contingent:
        rec["contingent"] = contingent
    h = r.get("hessian")
    if h is not None:
        rec["hessian"] = {"method": h["hessian_method"], "max_asymmetry": h["hessian_max_asymmetry"]}
    g = r.get("gfnff")
    if g is not None:
        topo = g["topology"]
        rec["gfnff"] = {
            "n_hbonds": g["n_hbonds"],
            "n_xbonds": g["n_xbonds"],
            "fragment_charges": topo["fragment_charges"],
        }
        if not session:
            rec["gfnff"]["topology"] = topo
    return rec


def timings(r: Mapping[str, Any], total: float) -> dict[str, float]:
    """``execution.timings``: the single point's phases, the Hessian's (``hessian:``), ``total``."""
    t = {k: float(x) for k, x in (r.get("timings") or {}).items()}
    h = r.get("hessian")
    if h is not None:
        t.update({f"hessian:{k}": float(x) for k, x in h["hessian_timings"].items()})
    t["total"] = total
    return t


def warnings_of(r: Mapping[str, Any], res: SinglePointData | None = None) -> list[ProgramWarning]:
    """The warnings of a response: a run-time LAPACK that was not used; orbital levels qcdata
    cannot hold (smeared occupations of a spin without electrons have no Fermi level)."""
    out = []
    if r.get("eigensolver_warning"):
        out.append(build(ProgramWarning, code="eigensolver_fallback", message=r["eigensolver_warning"]))
    if res is not None and r.get("scf") is not None and res.wavefunction is None:
        out.append(
            build(
                ProgramWarning,
                code="orbitals_omitted",
                message="smeared occupations of a spin without electrons have no Fermi level, which qcdata's "
                "Orbitals need: the wavefunction is left out",
            )
        )
    return out


def failed_output(
    inp: ProgramInput,
    e: ProgramError,
    *,
    p: Planned | None = None,
    resolved: Resolved | None = None,
    warnings: list[ProgramWarning] | None = None,
    rec: dict[str, Any] | None = None,
    execution: ExecutionInfo | None = None,
) -> ProgramOutput[Any, Any]:
    """The failed output of an input: empty results with axtb's provenance, the error, and what
    planning learnt (``resolved``, ``record``)."""
    fields: dict[str, Any] = {"record": rec or {}, "warnings": warnings or []}
    if resolved is not None:
        fields["resolved"] = resolved
    if execution is not None:
        fields["execution"] = execution
    return ProgramOutput.failure(inp, e.error, results=SinglePointData(provenance=provenance()), **fields)


def compute(inp: ProgramInput, ctx: Any = None) -> ProgramOutput[ProgramInput, SinglePointData]:
    """One calculation: ``energy``, ``gradient`` or ``hessian`` (with the vibrational analysis
    and thermochemistry unless ``freq.analysis`` is ``off``).

    Args:
        inp: The input: structure, ``model.method``, keywords (see ``docs/python.md``).
        ctx: A runner's `qccompute.program.RunContext`; only ``ctx.threads`` is read.

    Returns:
        The successful output.

    Raises:
        ProgramError: The input is refused or the calculation failed; ``exc.output`` is the
            failed output.
        TypeError: ``inp`` is not a `qcdata.ProgramInput`.
    """
    start = time.perf_counter()
    if not isinstance(inp, ProgramInput):
        raise TypeError(f"axtb.compute takes a qcdata.ProgramInput, not {type(inp).__name__}")
    try:
        p = plan_input(inp)
    except ProgramError as e:
        e.output = failed_output(inp, e)
        raise
    calctype = inp.calctype.value
    n, warnings, why = threads(p, calctype, ctx)
    resolved = _build(Resolved, settings={**p.settings, "threads": n}, reasons={**p.reasons, "threads": why})
    request = json.dumps({**p.request, "threads": n}, separators=(",", ":"))
    try:
        outputs = implied(p.outputs)
        r = _axtb.compute_once(
            request, p.numbers, p.positions, task(p, calctype, outputs), orbitals=needs_orbitals(outputs)
        )
    except _axtb.Error as x:
        error, fields = from_driver(x)
        rec = record(None, p)
        rec["error"] = fields
        out_error = ProgramError(error)
        execution = ExecutionInfo(threads=n, timings={"total": time.perf_counter() - start})
        out_error.output = failed_output(
            inp, out_error, p=p, resolved=resolved, warnings=warnings, rec=rec, execution=execution
        )
        raise out_error from None
    res = results(r, p, calctype, outputs)
    rec = record(r, p)
    out = _build(
        _OUTPUT,
        input_data=inp,
        success=True,
        results=res,
        resolved=resolved,
        warnings=warnings + warnings_of(r, res),
        record=rec,
        execution=_build(
            ExecutionInfo, devices=[], threads=r["threads"], timings=timings(r, time.perf_counter() - start)
        ),
    )
    return out
