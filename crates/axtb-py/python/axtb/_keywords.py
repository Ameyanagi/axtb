# SPDX-License-Identifier: LGPL-3.0-or-later
"""axtb's keywords: the table, the strict checks and the plan of an input.

The table (`KEYS`) is data: one entry per setting, with its tier, unit, concern, the methods
and calculation types it applies to, its default and its description. It generates
`settings_schema` and the keyword reference in ``docs/python.md`` (`markdown_table`; a test
fails when the file differs). `plan_input` checks a `qcdata.ProgramInput` against it in the
order of qcdata's conventions (D41; the first error raises):

1. ``unknown_key``: not in the table, with a did-you-mean from the 0.1 names and from close
   spellings; a group that is not an object is ``invalid_value``;
2. ``invalid_value``: a wrong type, range or shape, or an unknown name;
3. ``no_effect``: made irrelevant by the method, the conventions, another setting or the call;
4. ``not_applicable``: a setting of a calctype that will not run; ``files``,
   ``cmdline_args``, ``subprograms``, ``structures``; a ``given`` with nothing usable;
5. ``missing_setting``: a value that needs another key;
6. the structure: ``invalid_structure`` (non-finite positions), ``unsupported_element``,
   ``invalid_structure`` (a multiplicity whose parity is not that of the electrons the method's
   orbitals hold: valence-only methods leave cores out, `axtb._axtb.valence_electrons`),
   ``not_modeled`` (GFN-FF with multiplicity > 1);
7. kind ``unsupported``: a shared setting axtb cannot honour (``unsupported_setting``), an
   output it cannot give (``not_applicable``).

An unsupported calctype is refused before all of these (``unsupported_calctype``): the other
checks need a calctype axtb computes.

Names follow the alias rule (case and anything but ASCII letters and digits are ignored) and
``resolved`` echoes the canonical spellings. The plan builds the driver's flat request from
the settings given (never from defaults, which the driver applies itself: an explicit
electronic temperature, for example, turns off the composite's g-xTB retry ladder), so a
calculation runs exactly what the CLI runs for the same request.
"""

from __future__ import annotations

import dataclasses
import difflib
import json
import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from qcdata import Accuracy, ProgramError, ProgramInput

from . import _axtb
from ._errors import refuse

__all__ = [
    "CALCTYPES",
    "KEYS",
    "METHODS",
    "Planned",
    "alias",
    "canonical_method",
    "markdown_table",
    "plan_input",
    "settings_schema",
]

METHODS = ("gfn2", "gfnff", "gxtb-v1.1", "gxtb-v2.0.1")
CALCTYPES = ("energy", "gradient", "hessian")
_SCF = frozenset({"gfn2", "gxtb-v1.1", "gxtb-v2.0.1", "composite"})
_ALL = frozenset({"gfn2", "gxtb-v1.1", "gxtb-v2.0.1", "gfnff", "composite"})
_GXTB = frozenset({"gxtb-v1.1", "gxtb-v2.0.1"})
_EGH = frozenset(CALCTYPES)
_H = frozenset({"hessian"})


def alias(name: str) -> str:
    """The key of a name under the alias rule: its ASCII letters and digits, lower case."""
    return "".join(c for c in name.lower() if c.isascii() and c.isalnum())


# --------------------------------------------------------------------------- names


#: Conventions per method: the canonical name and the driver's name.
CONVENTIONS: dict[str, dict[str, str]] = {
    "gfn2": {"axtb": "axtb", "xtb-6.8.0": "xtb", "xtb-6.7.1": "xtb671", "tblite-0.8.0": "tblite"},
    "gfnff": {"axtb": "axtb", "xtb-6.8.0": "xtb", "xtb-6.7.1": "xtb671"},
    "gxtb-v1.1": {"axtb": "axtb", "gxtb-1.1.0": "gxtb"},
    "gxtb-v2.0.1": {"axtb": "axtb", "gxtb-2.0.1": "gxtb"},
}
#: The fixed aliases of the versioned conventions names.
_CONVENTION_ALIASES = {"xtb": "xtb-6.8.0", "tblite": "tblite-0.8.0", "gxtb": None}
#: The vibrational analysis' conventions: canonical name and driver name.
ANALYSIS_CONVENTIONS = {"axtb": "axtb", "xtb-6.8.0": "xtb", "xtb-6.7.1": "xtb671"}
#: The analysis conventions of a calculation's conventions (MethodInfo::analysis).
_ANALYSIS_OF = {
    "xtb-6.8.0": "xtb-6.8.0",
    "xtb-6.7.1": "xtb-6.7.1",
    "gxtb-1.1.0": "xtb-6.7.1",
    "gxtb-2.0.1": "xtb-6.7.1",
}


def _enum(value: Any, allowed: Iterable[str], path: str, aliases: Mapping[str, str] | None = None) -> str:
    """The canonical spelling of `value` among `allowed` (alias rule), else ``invalid_value``."""
    allowed = list(allowed)
    if not isinstance(value, str):
        raise refuse("input", "invalid_value", f"{path}: must be one of {_quoted(allowed)}, not {value!r}")
    k = alias(value)
    for a in allowed:
        if alias(a) == k:
            return a
    for a, target in (aliases or {}).items():
        if alias(a) == k:
            return target
    close = difflib.get_close_matches(value, allowed, n=1, cutoff=0.6)
    hint = f"; did you mean {close[0]!r}?" if close else ""
    raise refuse("input", "invalid_value", f"{path}: unknown value {value!r}; use {_quoted(allowed)}{hint}")


def _quoted(names: Iterable[str]) -> str:
    names = [f'"{n}"' for n in names]
    return names[0] if len(names) == 1 else ", ".join(names[:-1]) + " or " + names[-1]


def canonical_method(name: Any) -> str:
    """The canonical method of a ``model.method`` (``GFN2-xTB`` -> ``gfn2``), else
    ``invalid_value`` (an unversioned g-xTB names the versions)."""
    if not isinstance(name, str) or not name:
        raise refuse("input", "invalid_value", f"model.method: a method name, not {name!r}")
    if alias(name) in ("gxtb", "gxtbv2", "gxtb2", "gxtbv1", "gxtb1"):
        raise refuse(
            "input",
            "invalid_value",
            f'model.method: {name!r} names no version of g-xTB; use "gxtb-v2.0.1" (g-xTB 2.0.1) or "gxtb-v1.1" '
            "(the preprint's version 1.1)",
        )
    try:
        return str(_axtb.canonical("method", name))
    except _axtb.Error:
        raise refuse(
            "input",
            "invalid_value",
            f'model.method: unknown method {name!r}; use "gfn2", "gfnff", "gxtb-v1.1" or "gxtb-v2.0.1"',
        ) from None


def canonical_conventions(value: Any, method: str, path: str = "conventions") -> str:
    """The canonical conventions name of `method`."""
    table = CONVENTIONS[method]
    aliases = {
        "xtb": "xtb-6.8.0",
        "tblite": "tblite-0.8.0",
        "gxtb": "gxtb-1.1.0" if method == "gxtb-v1.1" else "gxtb-2.0.1",
    }
    aliases = {k: v for k, v in aliases.items() if v in table}
    return _enum(value, table, path, aliases)


# --------------------------------------------------------------------------- the table


@dataclass(frozen=True)
class Key:
    """One setting of the table.

    Attributes:
        path: ``group.key`` or a top-level key.
        shared: Whether qcdata's conventions name it (else axtb's own).
        tier: ``common`` or ``expert`` (``x-tier``).
        concern: ``model``, ``approximation``, ``convergence``, ``algorithm``, ``resource`` or
            ``output`` (``x-concern``).
        unit: Its unit (``x-unit``), if any.
        schema: The JSON Schema of its value.
        methods: The methods it applies to (``composite``: g-xTB with ``solvent.correction``).
        calctypes: The calctypes it applies to.
        default: Its default, as the reference shows it.
        description: What it does.
        driver: The driver request's field(s) it sets.
        unsupported: For a shared setting axtb cannot honour: why (kind ``unsupported``).
    """

    path: str
    shared: bool
    tier: str
    concern: str
    unit: str | None
    schema: dict[str, Any]
    methods: frozenset[str]
    calctypes: frozenset[str]
    default: str
    description: str
    driver: str = ""
    unsupported: str | None = None


_NUM = {"type": "number"}
_POS = {"type": "number", "exclusiveMinimum": 0}
_NONNEG = {"type": "number", "minimum": 0}
_INT1 = {"type": "integer", "minimum": 1}
_BOOL = {"type": "boolean"}


def _enum_schema(*values: str) -> dict[str, Any]:
    return {"type": "string", "enum": list(values)}


_ACCURACY_SCHEMA = {
    "type": "object",
    "additionalProperties": False,
    "minProperties": 1,
    "properties": {q: {"type": "number", "exclusiveMinimum": 0} for q in ("energy", "gradient", "hessian")},
}
_GEOMETRY_SCHEMA = {"type": "array", "items": {"type": "array", "items": _NUM, "minItems": 3, "maxItems": 3}}
_RETRIES_SCHEMA = {
    "type": "array",
    "minItems": 1,
    "items": {
        "type": "object",
        "additionalProperties": False,
        "minProperties": 1,
        "properties": {"etemp": {**_NONNEG, "x-unit": "K"}, "max_iter": _INT1},
    },
}
_UNSUPPORTED_SCHEMA: dict[str, Any] = {"not": {}}

KEYS: tuple[Key, ...] = (
    Key(
        "threads",
        True,
        "common",
        "resource",
        None,
        _INT1,
        _ALL,
        _EGH,
        "the runner's grant (a plan without one: deferred), else the plan's",
        "CPU threads of the calculation, at most the grant (more: the grant and the warning resource_capped).",
        "threads",
    ),
    Key(
        "accuracy",
        True,
        "common",
        "convergence",
        "Eh, Eh/bohr, Eh/bohr²",
        _ACCURACY_SCHEMA,
        _ALL,
        _EGH,
        "none",
        "Accuracy targets (energy, gradient, hessian), turned into the SCF accuracy factor by axtb's "
        "error model and recorded under record.accuracy; GFN-FF meets every target (no thresholds).",
    ),
    Key(
        "accuracy_scale",
        False,
        "expert",
        "convergence",
        None,
        _POS,
        _SCF,
        _EGH,
        "1",
        "Multiplier of every SCF threshold, and of GFN2-xTB's integral screening (xtb's --acc); a floor under "
        "the accuracy targets.",
        "accuracy",
    ),
    Key(
        "conventions",
        False,
        "common",
        "model",
        None,
        {"type": "string"},
        _ALL,
        _EGH,
        "axtb",
        "Whose numbers to reproduce: axtb (CODATA 2022, fixes), or a reference program: GFN2-xTB "
        "xtb-6.8.0, xtb-6.7.1, tblite-0.8.0; g-xTB gxtb-1.1.0, gxtb-2.0.1; GFN-FF xtb-6.8.0, xtb-6.7.1 "
        "(aliases xtb, tblite, gxtb).",
        "conventions",
    ),
    Key(
        "scf.max_iter",
        True,
        "expert",
        "convergence",
        None,
        _INT1,
        _SCF,
        _EGH,
        "GFN2-xTB, g-xTB 2.0.1 250; v1.1 1000",
        "Iteration limit of one SCF.",
        "max_iter",
    ),
    Key(
        "scf.etemp",
        False,
        "common",
        "model",
        "K",
        _NONNEG,
        _SCF,
        _EGH,
        "GFN2-xTB 300; g-xTB 0; the composite each part's (g-xTB 0, GFN2 300; not echoed)",
        "Electronic temperature of the Fermi smearing; given for the composite, it sets every part (and "
        "turns off the g-xTB part's retries at 1000 and 3000 K).",
        "etemp",
    ),
    Key(
        "scf.retries",
        False,
        "expert",
        "convergence",
        None,
        _RETRIES_SCHEMA,
        _SCF,
        _EGH,
        "none",
        "SCF retry ladder: after a calculation whose SCF fails, a fresh SCF at each rung in order until one "
        'converges; a rung is {"etemp": K, "max_iter": n}, either key optional (missing: the calculation\'s). '
        "No rung below a given scf.etemp. A Hessian runs whole at one rung. A result from a rung has the warning "
        "scf_retry (contingent). The composite: replaces its built-in retries, for every part.",
        "scf_retries",
    ),
    Key(
        "scf.mixer",
        False,
        "expert",
        "algorithm",
        None,
        _enum_schema("broyden", "diis", "broyden-then-diis"),
        frozenset({"gfn2", "composite"}),
        _EGH,
        "axtb conventions broyden-then-diis; the programs' broyden",
        "Mixer of GFN2-xTB's SCF (and of the composite's GFN2 parts): broyden-then-diis falls back to DIIS "
        "when Broyden fails (contingent).",
        "mixer",
    ),
    Key(
        "scf.guess",
        True,
        "expert",
        "algorithm",
        None,
        _enum_schema("eeq", "sad", "eeq-bc", "given"),
        frozenset({"gfn2", "gxtb-v1.1", "gxtb-v2.0.1"}),
        _EGH,
        "GFN2-xTB eeq (tblite-0.8.0 sad); g-xTB eeq-bc",
        "Where an SCF without history starts: GFN2-xTB eeq (charges) or sad (neutral atoms), g-xTB eeq-bc; "
        "given: from the input's given wavefunction.",
        "guess",
    ),
    Key(
        "scf.restart",
        False,
        "expert",
        "algorithm",
        None,
        _enum_schema("previous", "off"),
        _SCF,
        _EGH,
        "previous; g-xTB in its program's conventions off",
        "Sessions: whether an SCF starts from its lineage's history (previous) or fresh (off).",
        "restart",
    ),
    Key(
        "scf.e_tol",
        True,
        "expert",
        "convergence",
        "Eh",
        _UNSUPPORTED_SCHEMA,
        frozenset(),
        frozenset(),
        "-",
        "Not supported.",
        unsupported="axtb's SCF thresholds follow accuracy_scale and accuracy, not explicit tolerances",
    ),
    Key(
        "scf.g_tol",
        True,
        "expert",
        "convergence",
        None,
        _UNSUPPORTED_SCHEMA,
        frozenset(),
        frozenset(),
        "-",
        "Not supported.",
        unsupported="axtb's SCF thresholds follow accuracy_scale and accuracy, not explicit tolerances",
    ),
    Key(
        "solvent.name",
        True,
        "common",
        "model",
        None,
        {"type": "string"},
        frozenset({"gfn2", "composite"}),
        _EGH,
        "none",
        "Implicit solvent (ALPB or GBSA parameters of GFN2-xTB). g-xTB has none: with solvent.correction "
        "gfn2 it runs the composite E(g-xTB) + [E(GFN2-xTB + solvent) - E(GFN2-xTB)].",
        "alpb, gbsa",
    ),
    Key(
        "solvent.model",
        True,
        "common",
        "model",
        None,
        _enum_schema("alpb", "gbsa"),
        frozenset({"gfn2", "composite"}),
        _EGH,
        "alpb",
        "The continuum model of solvent.name.",
        "alpb, gbsa",
    ),
    Key(
        "solvent.eps",
        True,
        "expert",
        "model",
        None,
        _UNSUPPORTED_SCHEMA,
        frozenset(),
        frozenset(),
        "-",
        "Not supported.",
        unsupported="axtb's ALPB and GBSA models are parametrized per solvent name "
        "(solvent.name), not by a permittivity",
    ),
    Key(
        "solvent.state",
        False,
        "expert",
        "model",
        None,
        _enum_schema("gsolv", "bar1M", "reference"),
        frozenset({"gfn2", "composite"}),
        _EGH,
        "gsolv",
        "Reference state of the solvation free energy.",
        "solvation_state",
    ),
    Key(
        "solvent.correction",
        False,
        "common",
        "model",
        None,
        _enum_schema("gfn2"),
        frozenset({"composite"}),
        _EGH,
        "none",
        "g-xTB in solvent: the composite with GFN2-xTB's solvation correction.",
        "solvation_correction",
    ),
    Key(
        "solvent.correction_conventions",
        False,
        "expert",
        "model",
        None,
        _enum_schema("axtb", "xtb-6.8.0", "xtb-6.7.1", "tblite-0.8.0"),
        frozenset({"composite"}),
        _EGH,
        "axtb; xtb-6.7.1 with g-xTB's program conventions",
        "Conventions of the composite's GFN2-xTB parts.",
        "correction_conventions",
    ),
    Key(
        "topology.geometry",
        False,
        "expert",
        "model",
        "bohr",
        _GEOMETRY_SCHEMA,
        frozenset({"gfnff"}),
        _EGH,
        "the input's geometry (a session: the opening one; not echoed)",
        "The positions GFN-FF sets its topology up from.",
        "topology_positions",
    ),
    Key(
        "topology.fragment_charges",
        False,
        "expert",
        "model",
        "e",
        {"type": "array", "items": _NUM},
        frozenset({"gfnff"}),
        _EGH,
        "the charge on atom 1's fragment (in record.gfnff; not echoed)",
        "Charges of the molecular fragments, in the order of their first atoms, summing to the charge.",
        "fragment_charges",
    ),
    Key(
        "hessian.method",
        False,
        "common",
        "algorithm",
        None,
        _enum_schema("analytic", "response", "numerical"),
        _ALL,
        _H,
        "analytic; g-xTB v1.1 numerical; the composite each part's",
        "analytic second derivatives, response (central differences of the response gradient) or numerical "
        "(central differences of analytic gradients). g-xTB v1.1: numerical, response (axtb conventions); "
        "GFN-FF: analytic, numerical.",
        "hessian_method",
    ),
    Key(
        "hessian.step",
        False,
        "expert",
        "approximation",
        "bohr",
        _POS,
        _ALL,
        _H,
        "numerical 0.005; response GFN2-xTB 5e-4, g-xTB 2e-4",
        "Finite-difference step of the numerical and response Hessians.",
        "hessian_step",
    ),
    Key(
        "hessian.accuracy_scale",
        False,
        "expert",
        "convergence",
        None,
        _POS,
        _SCF,
        _H,
        "numerical 0.01; analytic and response GFN2-xTB 1e-3, g-xTB 1e-5",
        "SCF accuracy factor of the numerical Hessian's displacements, or of the analytic and response "
        "Hessians' reference; never looser than the calculation's.",
        "hessian_accuracy",
    ),
    Key(
        "freq.analysis",
        False,
        "common",
        "output",
        None,
        _enum_schema("harmonic", "off"),
        _ALL,
        _H,
        "harmonic",
        "The harmonic analysis of the Hessian (vibrational_analysis), and its thermochemistry.",
        "vibrations",
    ),
    Key(
        "freq.thermo",
        False,
        "common",
        "output",
        None,
        _enum_schema("qrrho", "off"),
        _ALL,
        _H,
        "qrrho",
        "Thermochemistry: Grimme's quasi-RRHO (xtb's), or none.",
        "thermo",
    ),
    Key(
        "freq.conventions",
        False,
        "expert",
        "model",
        None,
        _enum_schema("axtb", "xtb-6.8.0", "xtb-6.7.1"),
        _ALL,
        _H,
        "the calculation's (axtb; xtb-6.8.0 or xtb-6.7.1 in those; g-xTB's programs xtb-6.7.1)",
        "Constants, masses and projection of the analysis.",
        "analysis_conventions",
    ),
    Key(
        "freq.masses",
        False,
        "expert",
        "model",
        "u",
        {"type": "array", "items": _POS},
        _ALL,
        _H,
        "the analysis conventions' table (not echoed)",
        "Atomic masses of the analysis.",
        "masses",
    ),
    Key(
        "freq.temperature",
        True,
        "common",
        "model",
        "K",
        {"anyOf": [_POS, {"type": "array", "items": _POS, "minItems": 1}]},
        _ALL,
        _H,
        "298.15",
        "Temperature, or several: one thermochemistry entry each, in order.",
        "temperature",
    ),
    Key(
        "freq.pressure",
        True,
        "common",
        "model",
        "Pa",
        _POS,
        _ALL,
        _H,
        "101325",
        "Pressure of the ideal gas.",
        "pressure",
    ),
    Key(
        "freq.scale",
        True,
        "expert",
        "model",
        None,
        _POS,
        _ALL,
        _H,
        "1 (GFN-FF in the xtb analysis conventions 1.03)",
        "Scale factor of the frequencies in the thermochemistry.",
        "frequency_scale",
    ),
    Key(
        "freq.imaginary_cutoff",
        True,
        "expert",
        "model",
        "cm⁻¹",
        _NONNEG,
        _ALL,
        _H,
        "20",
        "Imaginary modes below this magnitude count as real in the thermochemistry (xtb's ithr).",
        "imaginary_cutoff",
    ),
    Key(
        "freq.rotor_cutoff",
        False,
        "expert",
        "model",
        "cm⁻¹",
        _POS,
        _ALL,
        _H,
        "50",
        "The quasi-RRHO rotor cutoff (xtb's sthr).",
        "rotor_cutoff",
    ),
    Key(
        "freq.symmetry_number",
        False,
        "expert",
        "model",
        None,
        _INT1,
        _ALL,
        _H,
        "from the point group (not echoed)",
        "Rotational symmetry number.",
        "symmetry_number",
    ),
    Key(
        "freq.symmetry_tol",
        False,
        "expert",
        "model",
        "bohr",
        _POS,
        _ALL,
        _H,
        "0.1",
        "Tolerance of the point-group detection (xtb's desy).",
        "symmetry_threshold",
    ),
    Key(
        "output.orbitals",
        True,
        "common",
        "output",
        None,
        _BOOL,
        _SCF,
        _EGH,
        "false",
        "The orbitals' coefficients with their basis set (wavefunction.orbitals).",
        "",
    ),
    Key(
        "output.overlap",
        True,
        "common",
        "output",
        None,
        _BOOL,
        _SCF,
        _EGH,
        "false",
        "The AO overlap matrix (wavefunction.overlap), with the basis set.",
        "overlap",
    ),
    Key(
        "output.hamiltonian",
        True,
        "common",
        "output",
        None,
        _BOOL,
        _SCF,
        _EGH,
        "false",
        "The Hamiltonian whose eigenvectors are the orbitals (wavefunction.hamiltonian_a/_b), with the "
        "overlap and the orbitals; not g-xTB 2.0.1 in its program's conventions.",
        "hamiltonian",
    ),
    Key(
        "output.core_hamiltonian",
        True,
        "expert",
        "output",
        None,
        _UNSUPPORTED_SCHEMA,
        frozenset(),
        frozenset(),
        "-",
        "Not available.",
        unsupported="tight binding's H0 is model data, not a core Hamiltonian",
    ),
    Key(
        "output.bond_orders",
        True,
        "common",
        "output",
        None,
        _BOOL,
        _SCF,
        _EGH,
        "false",
        "Mayer bond orders (bond_orders).",
        "bond_orders",
    ),
    Key(
        "output.molden",
        True,
        "expert",
        "output",
        None,
        _UNSUPPORTED_SCHEMA,
        frozenset(),
        frozenset(),
        "-",
        "Not available.",
        unsupported="request output.orbitals and write the Molden file with qcinf",
    ),
    Key(
        "output.density",
        False,
        "expert",
        "output",
        None,
        _UNSUPPORTED_SCHEMA,
        frozenset(),
        frozenset(),
        "-",
        "Not available.",
        unsupported="the SCF density is Orbitals.density() of output.orbitals",
    ),
    Key(
        "memory_gib",
        True,
        "expert",
        "resource",
        "GiB",
        _UNSUPPORTED_SCHEMA,
        frozenset(),
        frozenset(),
        "-",
        "Not supported.",
        unsupported="axtb has no memory setting",
    ),
    Key(
        "backend",
        True,
        "expert",
        "resource",
        None,
        _UNSUPPORTED_SCHEMA,
        frozenset(),
        frozenset(),
        "-",
        "Not supported.",
        unsupported="axtb runs on the CPU only",
    ),
    Key(
        "device_memory_gib",
        True,
        "expert",
        "resource",
        "GiB",
        _UNSUPPORTED_SCHEMA,
        frozenset(),
        frozenset(),
        "-",
        "Not supported.",
        unsupported="axtb runs on the CPU only",
    ),
)
_BY_PATH = {k.path: k for k in KEYS}
GROUPS = ("scf", "solvent", "topology", "hessian", "freq", "output")

#: The 0.1 names and the other names people try: the setting meant.
RENAMED = {
    "etemp": "scf.etemp",
    "max_iter": "scf.max_iter",
    "maxiter": "scf.max_iter",
    "retries": "scf.retries",
    "scf_retries": "scf.retries",
    "mixer": "scf.mixer",
    "guess": "scf.guess (the start) or scf.restart (previous or off)",
    "restart": "scf.restart",
    "alpb": 'solvent.name (with solvent.model "alpb")',
    "gbsa": 'solvent.name with solvent.model "gbsa"',
    "solvent": "solvent.name",
    "solvation_state": "solvent.state",
    "solvation_correction": "solvent.correction",
    "correction_conventions": "solvent.correction_conventions",
    "fragment_charges": "topology.fragment_charges",
    "topology_positions": "topology.geometry",
    "hessian_method": "hessian.method",
    "hessian_step": "hessian.step",
    "hessian_accuracy": "hessian.accuracy_scale",
    "vibrations": 'freq.analysis ("harmonic" or "off")',
    "thermo": 'freq.thermo ("qrrho" or "off")',
    "analysis_conventions": "freq.conventions",
    "masses": "freq.masses",
    "temperature": "freq.temperature",
    "pressure": "freq.pressure",
    "rotor_cutoff": "freq.rotor_cutoff",
    "imaginary_cutoff": "freq.imaginary_cutoff",
    "frequency_scale": "freq.scale",
    "symmetry_number": "freq.symmetry_number",
    "symmetry_threshold": "freq.symmetry_tol",
    "bond_orders": "output.bond_orders",
    "overlap": "output.overlap",
    "hamiltonian": "output.hamiltonian",
    "orbitals": "output.orbitals",
    "molden": "output.orbitals, and the Molden file from qcinf",
    "molden_file": "output.orbitals, and the Molden file from qcinf",
    "uhf": "the structure's multiplicity (unpaired electrons + 1)",
    "charge": "the structure's charge",
    "multiplicity": "the structure's multiplicity",
    "method": "model.method",
    "numbers": "the structure",
    "positions": "the structure",
    "gradient": 'the calctype "gradient"',
    "hessian": 'the calctype "hessian" (and the group hessian.*)',
    "dispersion": "nothing: switching off dispersion is the CLI's (a validation setting)",
    "scf_trace": "nothing: the SCF trace is the CLI's (programs do not print in process)",
    "analyze_hessian": "qcinf (or the CLI's analyze_hessian)",
    "acc": "accuracy_scale",
    "threads_advice": "threads",
}


def _did_you_mean(path: str) -> str:
    """The hint of an unknown key: the setting it was renamed to, or a close spelling."""
    leaf = path.split(".")[-1]
    for name in (path, leaf):
        if name in RENAMED:
            return f"; did you mean {RENAMED[name]}?"
    close = difflib.get_close_matches(path, list(_BY_PATH) + list(GROUPS), n=1, cutoff=0.75)
    return f"; did you mean {close[0]}?" if close else ""


# --------------------------------------------------------------------------- the plan


@dataclass
class Planned:
    """A checked input: the method, its settings, the driver's request and the plan's parts.

    Attributes:
        method: The canonical method.
        conventions: The canonical conventions.
        composite: Whether it is the g-xTB solvation composite.
        calctypes: The calctypes that may run (``compute``: the input's; a session: its own).
        settings: ``resolved.settings`` without ``threads``.
        reasons: ``resolved.reasons``.
        given: The settings given, by path (canonical values).
        request: The driver's request (flat), of the settings given.
        task: The driver's task fields of a Hessian (the ``hessian.*`` and ``freq.*`` given).
        outputs: The outputs the keywords ask for (``orbitals``, ``overlap``, ``hamiltonian``,
            ``bond_orders``).
        accuracy: The keyword's accuracy target, if any.
        threads: The ``threads`` keyword, if given.
        numbers: Atomic numbers; ``positions`` the input's geometry (bohr).
        charge, multiplicity: The structure's.
        record: The plan's record (method, conventions, calctypes, numerics, contingent).
        restart: Whether an SCF starts from history in a session, per part (one entry, or the
            composite's g-xTB part and its GFN2 parts).
    """

    method: str
    conventions: str
    composite: bool
    calctypes: frozenset[str]
    settings: dict[str, Any]
    reasons: dict[str, str]
    given: dict[str, Any]
    request: dict[str, Any]
    task: dict[str, Any]
    outputs: frozenset[str]
    accuracy: Accuracy | None
    threads: int | None
    numbers: list[int]
    positions: np.ndarray
    charge: int
    multiplicity: int
    record: dict[str, Any]
    restart: tuple[bool, ...] = ()
    contingent: list[str] = field(default_factory=list)
    #: The start of a usable `given` (`given_start`): kind ``density`` (g-xTB) or ``moments``
    #: (GFN2-xTB).
    start: dict[str, Any] | None = None

    @property
    def scf(self) -> bool:
        """Whether the method has an SCF."""
        return self.method != "gfnff"

    @property
    def program_conventions(self) -> bool:
        """Whether the conventions are a reference program's."""
        return self.conventions != "axtb"


def _flatten(keywords: Mapping[str, Any]) -> dict[str, Any]:
    """The keywords by dotted path; ``unknown_key`` for keys outside the table, a group that is
    not an object ``invalid_value``."""
    flat: dict[str, Any] = {}
    for k, v in keywords.items():
        if not isinstance(k, str):
            raise refuse("input", "unknown_key", f"{k!r}: keywords are named by strings")
        if k in GROUPS:
            if not isinstance(v, Mapping):
                hint = f" (the 0.1 key {k!r} is {RENAMED[k]})" if k in RENAMED else ""
                raise refuse(
                    "input",
                    "invalid_value",
                    f'{k}: a group of settings (an object, e.g. {{"{k}": {{...}}}}), not {v!r}{hint}',
                )
            for kk, vv in v.items():
                path = f"{k}.{kk}"
                if path not in _BY_PATH:
                    raise refuse("input", "unknown_key", f"{path}: not a setting of axtb{_did_you_mean(path)}")
                flat[path] = vv
        elif k in _BY_PATH:
            flat[k] = v
        else:
            raise refuse("input", "unknown_key", f"{k}: not a setting of axtb{_did_you_mean(k)}")
    return flat


def _number(value: Any, path: str, *, positive: bool = False, nonnegative: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float) or not math.isfinite(value):
        raise refuse("input", "invalid_value", f"{path}: a finite number, not {value!r}")
    if positive and not value > 0:
        raise refuse("input", "invalid_value", f"{path}: must be positive, not {value!r}")
    if nonnegative and not value >= 0:
        raise refuse("input", "invalid_value", f"{path}: must not be negative, not {value!r}")
    return float(value)


def _integer(value: Any, path: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int | np.integer) or value < 1:
        raise refuse("input", "invalid_value", f"{path}: an integer of at least 1, not {value!r}")
    return int(value)


def _boolean(value: Any, path: str) -> bool:
    if not isinstance(value, bool):
        raise refuse("input", "invalid_value", f"{path}: true or false, not {value!r}")
    return value


def _array(value: Any, path: str, shape: tuple[int, ...]) -> np.ndarray:
    try:
        a = np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        raise refuse("input", "invalid_value", f"{path}: numbers of shape {shape}") from None
    if a.shape != shape:
        raise refuse("input", "invalid_value", f"{path}: shape {shape}, not {a.shape}")
    if not np.all(np.isfinite(a)):
        raise refuse("input", "invalid_value", f"{path}: finite numbers")
    return a


def _hessian_methods(method: str, conventions: str) -> tuple[str, ...]:
    """The Hessian methods offered (the composite: its g-xTB part's, which GFN2-xTB all offers)."""
    if method == "gfnff":
        return ("analytic", "numerical")
    if method == "gxtb-v1.1":
        return ("numerical", "response") if conventions == "axtb" else ("numerical",)
    return ("analytic", "response", "numerical")


def _default_hessian_method(method: str) -> str:
    return "numerical" if method == "gxtb-v1.1" else "analytic"


def _analysis_default(method: str, conventions: str) -> str:
    return _ANALYSIS_OF.get(conventions, "axtb")


def _restart_default(method: str, conventions: str) -> bool:
    """Whether the method's SCFs start from history by default (its restart policy)."""
    if method == "gfn2":
        return True
    return conventions == "axtb"


#: Plans of recent inputs by everything but their geometry (an optimizer's or a list's inputs
#: repeat their settings): planning costs ~35 us, a water gradient ~150 us.
_CACHE: dict[tuple[Any, ...], Planned] = {}
_CACHE_SIZE = 64


def _copy_tree(s: dict[str, Any]) -> dict[str, Any]:
    """A copy of a settings tree (groups, lists and the objects in lists copied; values are
    immutable)."""

    def item(x: Any) -> Any:
        return _copy_tree(x) if isinstance(x, dict) else x

    return {
        k: (_copy_tree(x) if isinstance(x, dict) else [item(e) for e in x] if isinstance(x, list) else x)
        for k, x in s.items()
    }


def plan_input(
    inp: ProgramInput,
    *,
    calctypes: frozenset[str] | None = None,
    session: bool = False,
) -> Planned:
    """Check an input and resolve its settings (see the module documentation); plans of inputs
    that differ only in their geometry are reused."""
    key = None
    if inp.given is None and inp.structure is not None and inp.model is not None:
        try:
            s = inp.structure
            key = (
                inp.model.method,
                inp.model.basis,
                inp.calctype,
                calctypes,
                session,
                json.dumps(inp.keywords, sort_keys=True),
                tuple(s.symbols),
                s.charge,
                s.multiplicity,
                bool(inp.files or inp.cmdline_args or inp.subprograms or inp.structures),
            )
        except (TypeError, ValueError):
            key = None
    if key is not None and key in _CACHE:
        p = _CACHE[key]
        assert inp.structure is not None  # the key includes it
        positions = np.asarray(inp.structure.geometry, dtype=float).reshape(-1, 3)
        if not np.all(np.isfinite(positions)):
            raise refuse("input", "invalid_structure", "structure.geometry: positions must be finite numbers")
        return dataclasses.replace(
            p,
            positions=positions,
            settings=_copy_tree(p.settings),
            reasons=dict(p.reasons),
            request=_copy_tree(p.request),
            task=_copy_tree(p.task),
            record=_copy_tree(p.record),
        )
    p = _plan_input(inp, calctypes=calctypes, session=session)
    if key is not None:
        if len(_CACHE) >= _CACHE_SIZE:
            del _CACHE[next(iter(_CACHE))]
        _CACHE[key] = dataclasses.replace(p, settings=_copy_tree(p.settings), request=_copy_tree(p.request))
    return p


def _plan_input(
    inp: ProgramInput,
    *,
    calctypes: frozenset[str] | None = None,
    session: bool = False,
) -> Planned:
    """Check an input and resolve its settings (see the module documentation).

    Args:
        inp: The input.
        calctypes: The calctypes that may run: a session's (each with its lower derivatives);
            None: the input's.
        session: Whether a session plans it (``scf.restart`` applies, ``hessian.*`` and
            ``freq.*`` are contingent when the session may run Hessians).

    Raises:
        ProgramError: The first refusal, in the order of the module documentation.
    """
    calctype = inp.calctype.value
    allowed = calctypes if calctypes is not None else frozenset({calctype})
    for c in sorted(allowed | {calctype}):
        if c not in CALCTYPES:
            raise refuse(
                "unsupported",
                "unsupported_calctype",
                f"calctype: axtb computes energy, gradient and hessian, not {c}; for an {c} run qcopt (or "
                "another driver) with axtb as its subprogram",
            )
    # 1. unknown keys
    flat = _flatten(inp.keywords)

    # 2. values: the model, the conventions, then every key
    if inp.model is None:
        raise refuse(
            "input", "invalid_value", 'model: axtb needs model.method ("gfn2", "gfnff", "gxtb-v1.1" or "gxtb-v2.0.1")'
        )
    method = canonical_method(inp.model.method)
    if inp.model.basis:
        raise refuse(
            "input",
            "invalid_value",
            f"model.basis: {inp.model.method} has its own basis; leave model.basis out, not {inp.model.basis!r}",
        )
    structure = inp.require_structure()
    nat = len(structure.symbols)
    v: dict[str, Any] = {}
    for path, value in flat.items():
        key = _BY_PATH[path]
        if key.unsupported is not None:
            v[path] = value  # refused in step 7, after every input error
            continue
        v[path] = _value(path, value, method, nat, flat)
    conventions = v.get("conventions", "axtb")
    gxtb = method in _GXTB
    composite = gxtb and "solvent.correction" in v and "solvent.name" in v
    kind = "composite" if composite else method
    if "hessian.method" in v:
        offered = _hessian_methods(method, conventions)
        if v["hessian.method"] not in offered:
            raise refuse(
                "input",
                "invalid_value",
                f"hessian.method: {method} in the {conventions} conventions offers {_quoted(offered)}, not "
                f"{v['hessian.method']!r}",
            )

    # 3. no_effect
    def no_effect(path: str, why: str) -> ProgramError:
        return refuse("input", "no_effect", f"{path}: {why}; leave it out")

    for path, value in v.items():
        key = _BY_PATH[path]
        if key.unsupported is not None:
            continue
        if path.startswith("output.") and value is False and not _applies(key, kind):
            raise no_effect(path, f"{method} gives no such output")
        if method == "gfnff":
            if path == "accuracy_scale" or path.startswith("scf.") or path == "hessian.accuracy_scale":
                raise no_effect(path, "GFN-FF has no SCF")
            continue
        if path.startswith("topology."):
            raise no_effect(path, f"{path} is a GFN-FF setting; {method} has no topology")
        if path == "scf.mixer" and gxtb and not composite:
            raise no_effect(path, f"{method} has one mixer (only GFN2-xTB and the composite's GFN2 parts choose)")
        if path == "scf.guess" and composite:
            raise no_effect(path, "the composite starts each part from its method's own start")
        if method == "gfn2" and path in ("solvent.correction", "solvent.correction_conventions"):
            raise no_effect(path, "GFN2-xTB has its own implicit solvent (solvent.name alone)")
    if v.get("output.hamiltonian") is False and method == "gxtb-v2.0.1" and conventions != "axtb":
        raise no_effect("output.hamiltonian", "g-xTB 2.0.1 in its program's conventions gives no Hamiltonian")
    if v.get("scf.guess") == "sad" and method == "gfn2" and conventions == "tblite-0.8.0":
        pass  # the conventions' own start
    elif v.get("scf.guess") == "eeq" and conventions == "tblite-0.8.0":
        raise no_effect("scf.guess", "tblite-0.8.0 starts every SCF from neutral atoms (sad)")
    hm = v.get("hessian.method", _default_hessian_method(method))
    if "hessian.step" in v and hm == "analytic":
        raise no_effect("hessian.step", "the analytic Hessian takes no step")
    if v.get("freq.analysis") == "off":
        for path in v:
            if path.startswith("freq.") and path != "freq.analysis":
                raise no_effect(path, 'freq.analysis is "off"')
    if v.get("freq.thermo") == "off":
        for path in _THERMO:
            if path in v:
                raise no_effect(path, 'freq.thermo is "off"')
    if "scf.restart" in v and not session:
        raise no_effect(
            "scf.restart", "a calculation of its own has one SCF and no history (scf.restart is a session setting)"
        )
    if "scf.retries" in v and "scf.etemp" in v:
        floor = v["scf.etemp"]
        for i, rung in enumerate(v["scf.retries"]):
            if rung.get("etemp", floor) < floor:
                raise no_effect(
                    f"scf.retries[{i}].etemp",
                    f"{rung['etemp']:g} K is below scf.etemp ({floor:g} K), where the ladder starts and below which "
                    "no rung goes",
                )

    # 4. not_applicable
    def not_applicable(path: str, why: str) -> ProgramError:
        return refuse("input", "not_applicable", f"{path}: {why}")

    if "hessian" not in allowed:
        for path in v:
            if path.startswith(("hessian.", "freq.")):
                what = "this session's calctypes" if session else f"a {calctype} calculation"
                raise not_applicable(path, f"a setting of the hessian calctype, not of {what}")
    for name in ("files", "cmdline_args", "subprograms", "structures"):
        if getattr(inp, name):
            raise not_applicable(name, f"axtb takes no {name}")
    start = None
    if inp.given is not None:
        start = given_start(inp, method, conventions, composite)
        if "scf.guess" in v and v["scf.guess"] != "given":
            raise refuse("input", "no_effect", "scf.guess: the given wavefunction starts the SCF; leave scf.guess out")

    # 5. missing_setting
    def missing(path: str, why: str) -> ProgramError:
        return refuse("input", "missing_setting", f"{path}: {why}")

    if method != "gfnff":
        if "solvent.correction" in v and "solvent.name" not in v:
            raise missing("solvent.correction", 'the composite needs the solvent: solvent.name (e.g. "water")')
        for path in ("solvent.model", "solvent.state"):
            if path in v and "solvent.name" not in v:
                raise missing(path, "a setting of the implicit solvent: give solvent.name")
        if gxtb and "solvent.name" in v and "solvent.correction" not in v:
            hint = (
                f"{method} has no implicit-solvent parameters (the g-xTB program stops with --alpb/--gbsa). For g-xTB "
                "in solvent use the solvation composite E = E(g-xTB) + [E(GFN2-xTB + ALPB/GBSA) - E(GFN2-xTB)]: add "
                'solvent.correction "gfn2"'
            )
            raise missing("solvent.name", hint)
        if gxtb and "solvent.correction_conventions" in v and "solvent.correction" not in v:
            raise missing("solvent.correction_conventions", "a setting of the composite: give solvent.correction")
    if v.get("scf.guess") == "given" and start is None:
        raise missing("scf.guess", "given needs a given wavefunction (ProgramInput.given)")

    # 6. the structure
    numbers = list(structure.atomic_numbers)
    charge, mult = int(structure.charge), int(structure.multiplicity)
    positions = np.asarray(structure.geometry, dtype=float).reshape(-1, 3)
    if not np.all(np.isfinite(positions)):
        raise refuse("input", "invalid_structure", "structure.geometry: positions must be finite numbers")
    supported = _supported(method, composite)
    bad = sorted({z for z in numbers if z not in supported})
    if bad:
        raise refuse(
            "input",
            "unsupported_element",
            f"structure.symbols: {method}{' (the composite)' if composite else ''} has no parameters for Z = "
            + ", ".join(map(str, bad)),
            {"z": bad[0]},
        )
    # the electrons the method's orbitals hold (GFN-FF: all), whose parity the multiplicity must have
    # (the composite: in its g-xTB and its GFN2-xTB parts)
    for m in (method, "gfn2") if composite else (method,):
        electrons = _electrons(m, numbers) - charge
        if electrons < 0 or (electrons - (mult - 1)) % 2 or mult - 1 > electrons:
            ok = [k for k in range(1, 9) if (electrons - (k - 1)) % 2 == 0 and k - 1 <= max(electrons, 0)][:3]
            held = "electrons" if m == "gfnff" else f"electrons in {m}'s orbitals"
            raise refuse(
                "input",
                "invalid_structure",
                f"structure.multiplicity: {mult} with {electrons} {held} is impossible; use "
                + ", ".join(str(k) for k in ok)
                + ", ...",
            )
    if method == "gfnff" and mult > 1:
        raise refuse(
            "input",
            "not_modeled",
            f"structure.multiplicity: GFN-FF has no electronic structure and models closed-shell molecules only, not "
            f"multiplicity {mult}",
        )
    if "topology.geometry" in v and v["topology.geometry"].shape != (nat, 3):
        raise refuse("input", "invalid_value", f"topology.geometry: shape ({nat}, 3)")
    if "freq.masses" in v and len(v["freq.masses"]) != nat:
        raise refuse("input", "invalid_value", f"freq.masses: one mass per atom ({nat})")
    if "topology.fragment_charges" in v and not math.isclose(sum(v["topology.fragment_charges"]), charge, abs_tol=1e-6):
        raise refuse("input", "invalid_value", f"topology.fragment_charges: must sum to the charge ({charge})")

    # 7. what axtb cannot honour
    for path, value in v.items():
        key = _BY_PATH[path]
        if key.unsupported is not None:
            if (path.startswith("output.") and value is False) or value == "off":
                raise refuse("input", "no_effect", f"{path}: {key.unsupported}; leave it out")
            code = "not_applicable" if path == "output.core_hamiltonian" else "unsupported_setting"
            raise refuse("unsupported", code, f"{path}: {key.unsupported}")
        if method == "gfnff" and path.startswith("solvent."):
            raise refuse("unsupported", "unsupported_setting", f"{path}: axtb's GFN-FF has no implicit solvent")
        if path.startswith("output.") and value is True and not _applies(key, kind):
            raise refuse("unsupported", "not_applicable", f"{path}: {method} has no orbitals or SCF matrices")
    if v.get("output.hamiltonian") is True and method == "gxtb-v2.0.1" and conventions != "axtb":
        raise refuse(
            "unsupported",
            "not_applicable",
            "output.hamiltonian: g-xTB 2.0.1 in its program's conventions diagonalizes a mixed matrix last, so no "
            'Hamiltonian has its orbitals as eigenvectors; use conventions "axtb"',
        )
    accuracy = None
    if "accuracy" in v:
        accuracy = v["accuracy"]
    p = _resolve(inp, method, conventions, composite, allowed, session, v, accuracy, numbers, positions, charge, mult)
    if start is not None:
        p.start = start
        if not composite:
            p.settings["scf"]["guess"] = "given"
        p.reasons["scf.guess"] = "given: the input's given " + (
            "shell populations and atomic multipoles" if start["kind"] == "moments" else "wavefunction"
        )
    return p


def given_start(inp: ProgramInput, method: str, conventions: str, composite: bool) -> dict[str, Any]:
    """The SCF start of the input's ``given`` results (decision 19): GFN2-xTB from given shell
    populations and atomic multipoles (shell charges = reference occupation - population); g-xTB
    (the composite: its g-xTB part) from a given wavefunction's orbitals with their basis set,
    P = C n C^T per spin channel. Every other field of ``given`` is left alone. Refused (kind
    ``input``, code ``not_applicable``) when nothing is usable, or when a used field does not fit:
    shells of another layout, another basis set (centers, angular momenta, exponents), other spin
    channels."""
    from qcdata import AOConvention, SphericalOrder

    def not_applicable(path: str, why: str) -> ProgramError:
        return refuse("input", "not_applicable", f"{path}: {why}")

    g = inp.given
    assert g is not None  # called with a given
    structure = inp.require_structure()
    numbers = list(structure.atomic_numbers)
    if method == "gfnff":
        raise not_applicable("given", "GFN-FF has no SCF to start from given results")
    layout = _axtb.shell_layout(method, numbers)
    if layout is None:  # elements without parameters are refused before
        raise not_applicable("given", f"{method} has no shells for these elements")
    if method == "gfn2":
        if g.shell_populations is None or g.atomic_dipoles is None or g.atomic_quadrupoles is None:
            raise not_applicable(
                "given",
                "GFN2-xTB starts from given shell populations and atomic multipoles (shell_populations with "
                "shell_atoms and shell_angular_momenta, atomic_dipoles, atomic_quadrupoles), which this given lacks",
            )
        if list(g.shell_atoms or ()) != [a for a, *_ in layout] or list(g.shell_angular_momenta or ()) != [
            ell for _, ell, *_ in layout
        ]:
            raise not_applicable("given.shell_populations", "indexed by other shells than GFN2-xTB's for these atoms")
        q = [ref - pop for (_, _, ref, _), pop in zip(layout, np.asarray(g.shell_populations), strict=True)]
        qp = np.asarray(g.atomic_quadrupoles, dtype=float)
        packed = np.stack([qp[:, 0, 0], qp[:, 0, 1], qp[:, 1, 1], qp[:, 0, 2], qp[:, 1, 2], qp[:, 2, 2]], axis=1)
        return {"kind": "moments", "moments": (q, np.asarray(g.atomic_dipoles, dtype=float), packed)}
    wf = g.wavefunction
    if wf is None or wf.basis is None or wf.orbitals.coefficients_a is None:
        raise not_applicable(
            "given",
            f"{method} starts from a given wavefunction's orbitals with their basis set (the output orbitals), "
            "which this given lacks",
        )
    restricted = wf.orbitals.spin.value == "restricted"
    if restricted != (structure.multiplicity == 1):
        raise not_applicable(
            "given.wavefunction",
            f"{'restricted' if restricted else 'unrestricted'} orbitals cannot start a "
            f"{'closed' if structure.multiplicity == 1 else 'open'}-shell g-xTB SCF (one channel per spin of the "
            "molecule)",
        )
    ours = AOConvention(spherical_order=SphericalOrder.increasing_m)
    if wf.basis.ao_convention != ours:
        wf = wf.with_convention(ours)
    assert wf.basis is not None  # checked above; with_convention keeps it
    shells = wf.basis.shells
    fits = len(shells) == len(layout) and all(
        sh.center == a and sh.angular_momentum == ell and sh.pure and np.array_equal(np.asarray(sh.exponents), e)
        for sh, (a, ell, _, e) in zip(shells, layout, strict=False)
    )
    if not fits:
        raise not_applicable(
            "given.wavefunction.basis",
            f"not {method}'s basis set (q-vSZP) on these atoms: shells, centers, angular momenta and exponents "
            "must be its own",
        )
    orb = wf.orbitals
    ca = np.asarray(orb.coefficients_a, dtype=float)
    if restricted:
        densities = [(ca * (np.asarray(orb.occupations_a) + np.asarray(orb.occupations_b))) @ ca.T]
    else:
        cb = np.asarray(orb.coefficients_b, dtype=float)
        densities = [(ca * np.asarray(orb.occupations_a)) @ ca.T, (cb * np.asarray(orb.occupations_b)) @ cb.T]
    return {"kind": "density", "densities": densities}


_THERMO = (
    "freq.temperature",
    "freq.pressure",
    "freq.scale",
    "freq.imaginary_cutoff",
    "freq.rotor_cutoff",
    "freq.symmetry_number",
    "freq.symmetry_tol",
)


def _applies(key: Key, kind: str) -> bool:
    return kind in key.methods


_VALENCE: dict[str, list[float | None]] = {}


def _electrons(method: str, numbers: list[int]) -> int:
    """The electrons a neutral molecule brings to `method`'s orbitals (valence-only methods
    leave the cores out: Zn 12 of 30 in GFN2-xTB and g-xTB; GFN-FF: all)."""
    if method not in _VALENCE:
        _VALENCE[method] = _axtb.valence_electrons(method)
    table = _VALENCE[method]
    return round(sum(table[z - 1] or 0.0 for z in numbers))


_SUPPORTED: dict[tuple[str, bool], frozenset[int]] = {}


def _supported(method: str, composite: bool) -> frozenset[int]:
    """The elements `method` (the composite: and GFN2-xTB) has parameters for."""
    k = (method, composite)
    if k not in _SUPPORTED:
        s = set(_axtb.supported_elements(method))
        if composite:
            s &= set(_axtb.supported_elements("gfn2"))
        _SUPPORTED[k] = frozenset(s)
    return _SUPPORTED[k]


def _value(path: str, value: Any, method: str, nat: int, flat: Mapping[str, Any]) -> Any:
    """The canonical value of a setting, else ``invalid_value``."""
    if path == "threads":
        return _integer(value, path)
    if path == "accuracy":
        if isinstance(value, int | float) and not isinstance(value, bool):
            raise refuse(
                "input",
                "invalid_value",
                f"accuracy: an object of targets (energy in Eh, gradient in Eh/bohr, hessian in Eh/bohr²), e.g. "
                f'{{"gradient": 1e-05}}, not the number {value!r}; axtb\'s SCF accuracy factor is accuracy_scale',
            )
        try:
            target = Accuracy.from_keywords({"accuracy": value})
        except ProgramError as e:
            if e.error.code in ("unknown_key", "no_effect"):
                raise
            raise refuse("input", "invalid_value", e.error.message) from None
        return target
    if path in (
        "accuracy_scale",
        "hessian.step",
        "hessian.accuracy_scale",
        "freq.pressure",
        "freq.scale",
        "freq.rotor_cutoff",
        "freq.symmetry_tol",
    ):
        return _number(value, path, positive=True)
    if path in ("scf.etemp", "freq.imaginary_cutoff"):
        return _number(value, path, nonnegative=True)
    if path in ("scf.max_iter", "freq.symmetry_number"):
        return _integer(value, path)
    if path == "scf.retries":
        return _retries(value)
    if path == "conventions":
        return canonical_conventions(value, method)
    if path == "scf.mixer":
        return _enum(
            value, ("broyden", "diis", "broyden-then-diis"), path, {"auto": "broyden-then-diis", "pulay": "diis"}
        )
    if path == "scf.guess":
        allowed = {"gfn2": ("eeq", "sad", "given"), "gfnff": ("eeq", "sad", "eeq-bc", "given")}.get(
            method, ("eeq-bc", "given")
        )
        return _enum(value, allowed, path)
    if path == "scf.restart":
        return _enum(value, ("previous", "off"), path)
    if path == "solvent.model":
        return _enum(value, ("alpb", "gbsa"), path)
    if path == "solvent.name":
        model = flat.get("solvent.model", "alpb")
        model = model if isinstance(model, str) and alias(model) in ("alpb", "gbsa") else "alpb"
        if not isinstance(value, str):
            raise refuse("input", "invalid_value", f"solvent.name: a solvent name, not {value!r}")
        name = _axtb.canonical_solvent(value, alias(model))
        if name is None:
            names = _axtb.solvents(alias(model))
            close = difflib.get_close_matches(value.lower(), names, n=1, cutoff=0.6)
            hint = f"; did you mean {close[0]!r}?" if close else ""
            raise refuse(
                "input", "invalid_value", f"solvent.name: {value!r} has no {alias(model).upper()} parameters{hint}"
            )
        return name
    if path == "solvent.state":
        return _enum(value, ("gsolv", "bar1M", "reference"), path)
    if path == "solvent.correction":
        return _enum(value, ("gfn2",), path, {"gfn2xtb": "gfn2"})
    if path == "solvent.correction_conventions":
        return canonical_conventions(value, "gfn2", path)
    if path == "topology.geometry":
        return _array(value, path, (nat, 3))
    if path == "topology.fragment_charges":
        if not isinstance(value, list | tuple) or not value:
            raise refuse("input", "invalid_value", f"{path}: a list of charges, one per fragment")
        return [_number(x, path) for x in value]
    if path == "hessian.method":
        return _enum(
            value, ("analytic", "response", "numerical"), path, {"analytical": "analytic", "numeric": "numerical"}
        )
    if path == "freq.analysis":
        return _enum(value, ("harmonic", "off"), path)
    if path == "freq.thermo":
        return _enum(value, ("qrrho", "off"), path)
    if path == "freq.conventions":
        return _enum(value, ANALYSIS_CONVENTIONS, path, {"xtb": "xtb-6.8.0"})
    if path == "freq.masses":
        if not isinstance(value, list | tuple | np.ndarray):
            raise refuse("input", "invalid_value", f"{path}: a list of masses (u), one per atom")
        return [_number(x, path, positive=True) for x in value]
    if path == "freq.temperature":
        if isinstance(value, list | tuple):
            if not value:
                raise refuse("input", "invalid_value", f"{path}: at least one temperature")
            return [_number(x, path, positive=True) for x in value]
        return _number(value, path, positive=True)
    if path.startswith("output."):
        return _boolean(value, path)
    raise AssertionError(f"no check for {path}")  # pragma: no cover - every key has one


def _retries(value: Any) -> list[dict[str, Any]]:
    """The canonical rungs of ``scf.retries``: each an object with ``etemp`` (K, a float) and/or
    ``max_iter`` (an int), as given."""
    path = "scf.retries"
    if not isinstance(value, list | tuple) or not value:
        raise refuse(
            "input",
            "invalid_value",
            f'{path}: a list of rungs, e.g. [{{"etemp": 1000, "max_iter": 1000}}, {{"etemp": 3000}}], not {value!r}',
        )
    out = []
    for i, rung in enumerate(value):
        at = f"{path}[{i}]"
        if not isinstance(rung, Mapping):
            raise refuse(
                "input", "invalid_value", f'{at}: a rung is an object {{"etemp": K, "max_iter": n}}, not {rung!r}'
            )
        for k in rung:
            if k not in ("etemp", "max_iter"):
                close = difflib.get_close_matches(str(k), ["etemp", "max_iter"], n=1, cutoff=0.6)
                hint = f"; did you mean {close[0]}?" if close else ""
                raise refuse("input", "unknown_key", f"{at}.{k}: not a setting of a rung (etemp, max_iter){hint}")
        if not rung:
            raise refuse(
                "input",
                "invalid_value",
                f"{at}: a rung sets etemp, max_iter or both (one without either repeats the SCF that failed)",
            )
        c: dict[str, Any] = {}
        if "etemp" in rung:
            c["etemp"] = _number(rung["etemp"], f"{at}.etemp", nonnegative=True)
        if "max_iter" in rung:
            c["max_iter"] = _integer(rung["max_iter"], f"{at}.max_iter")
        out.append(c)
    return out


def _put(tree: dict[str, Any], path: str, value: Any) -> None:
    group, _, leaf = path.rpartition(".")
    (tree.setdefault(group, {}) if group else tree)[leaf] = value


def _resolve(
    inp: ProgramInput,
    method: str,
    conventions: str,
    composite: bool,
    allowed: frozenset[str],
    session: bool,
    v: dict[str, Any],
    accuracy: Accuracy | None,
    numbers: list[int],
    positions: np.ndarray,
    charge: int,
    mult: int,
) -> Planned:
    """The resolved settings, the driver's request and the plan's record of a checked input."""
    s: dict[str, Any] = {}
    reasons: dict[str, str] = {}
    gxtb = method in _GXTB
    scf = method != "gfnff"
    hessian = "hessian" in allowed

    def put(path: str, value: Any) -> None:
        _put(s, path, value)

    if accuracy is not None:
        put("accuracy", accuracy.model_dump(mode="json", exclude_none=True, exclude={"qcdata_version"}))
    put("conventions", conventions)
    if scf:
        put("accuracy_scale", v.get("accuracy_scale", 1.0))
        # iteration limit: one value unless the composite's parts differ (v1.1's 1000, GFN2's 250)
        if "scf.max_iter" in v:
            put("scf.max_iter", v["scf.max_iter"])
        elif not (composite and method == "gxtb-v1.1"):
            put("scf.max_iter", 1000 if method == "gxtb-v1.1" else 250)
        if "scf.etemp" in v:
            put("scf.etemp", v["scf.etemp"])
        elif not composite:
            put("scf.etemp", 300.0 if method == "gfn2" else 0.0)
        else:
            reasons["scf.etemp"] = "the composite: each part its method's (g-xTB 0 K, GFN2-xTB 300 K)"
        if "scf.retries" in v:
            put("scf.retries", v["scf.retries"])
        gfn2_conventions = None
        if method == "gfn2":
            gfn2_conventions = conventions
        elif composite:
            gfn2_conventions = v.get("solvent.correction_conventions", "axtb" if conventions == "axtb" else "xtb-6.7.1")
        if gfn2_conventions is not None:
            mixer = v.get("scf.mixer", "broyden-then-diis" if gfn2_conventions == "axtb" else "broyden")
            put("scf.mixer", mixer)
            if "scf.mixer" not in v:
                reasons["scf.mixer"] = (
                    "conventions axtb: Broyden, then DIIS where Broyden fails"
                    if gfn2_conventions == "axtb"
                    else f"conventions {gfn2_conventions}: Broyden only, as the program"
                )
        if not composite:
            default_guess = "eeq-bc" if gxtb else ("sad" if conventions == "tblite-0.8.0" else "eeq")
            put("scf.guess", v.get("scf.guess", default_guess))
        # the restart policy per part (the composite: g-xTB, then the GFN2 parts)
        restart: tuple[bool, ...]
        if composite:
            restart = (_restart_default(method, conventions), gfn2_conventions == "axtb")
        else:
            restart = (_restart_default(method, conventions),)
        if "scf.restart" in v:
            restart = tuple(v["scf.restart"] == "previous" for _ in restart)
        if session:
            if len(set(restart)) == 1:
                put("scf.restart", "previous" if restart[0] else "off")
            else:
                reasons["scf.restart"] = "the composite: each part as its conventions (g-xTB off, GFN2-xTB previous)"
    else:
        restart = ()
    if "solvent.name" in v:
        put("solvent.name", v["solvent.name"])
        put("solvent.model", v.get("solvent.model", "alpb"))
        put("solvent.state", v.get("solvent.state", "gsolv"))
    if composite:
        put("solvent.correction", v["solvent.correction"])
        put("solvent.correction_conventions", gfn2_conventions)
    if hessian:
        hm_given = v.get("hessian.method")
        hm = hm_given or _default_hessian_method(method)
        if composite and hm_given is None:
            reasons["hessian.method"] = "the composite: each part its method's default"
        else:
            put("hessian.method", hm)
        if hm != "analytic" and not (composite and hm_given is None):
            if method == "gfnff" or hm == "numerical":
                step = 0.005
            else:
                step = 5e-4 if method == "gfn2" else 2e-4
            put("hessian.step", v.get("hessian.step", step))
        elif "hessian.step" in v:
            put("hessian.step", v["hessian.step"])
        if scf and not (composite and hm_given is None):
            acc = 0.01 if hm == "numerical" else (1e-3 if method == "gfn2" else 1e-5)
            put("hessian.accuracy_scale", v.get("hessian.accuracy_scale", acc))
        elif "hessian.accuracy_scale" in v:
            put("hessian.accuracy_scale", v["hessian.accuracy_scale"])
        analysis = v.get("freq.analysis", "harmonic")
        put("freq.analysis", analysis)
        if analysis == "harmonic":
            fc = v.get("freq.conventions", _analysis_default(method, conventions))
            put("freq.conventions", fc)
            thermo = v.get("freq.thermo", "qrrho")
            put("freq.thermo", thermo)
            if thermo == "qrrho":
                put("freq.temperature", v.get("freq.temperature", 298.15))
                put("freq.pressure", v.get("freq.pressure", 101325.0))
                scale = 1.03 if method == "gfnff" and fc != "axtb" else 1.0
                put("freq.scale", v.get("freq.scale", scale))
                put("freq.imaginary_cutoff", v.get("freq.imaginary_cutoff", 20.0))
                put("freq.rotor_cutoff", v.get("freq.rotor_cutoff", 50.0))
                put("freq.symmetry_tol", v.get("freq.symmetry_tol", 0.1))
                if "freq.symmetry_number" in v:
                    put("freq.symmetry_number", v["freq.symmetry_number"])
            if "freq.masses" in v:
                put("freq.masses", v["freq.masses"])
    if "topology.geometry" in v:
        put("topology.geometry", v["topology.geometry"].tolist())
    if "topology.fragment_charges" in v:
        put("topology.fragment_charges", v["topology.fragment_charges"])
    outputs = set()
    if scf:
        for name in ("orbitals", "overlap", "hamiltonian", "bond_orders"):
            if name == "hamiltonian" and method == "gxtb-v2.0.1" and conventions != "axtb":
                continue
            value = v.get(f"output.{name}", False)
            put(f"output.{name}", value)
            if value:
                outputs.add(name)
    # the driver's request: the settings given (the driver applies its defaults itself)
    req: dict[str, Any] = {"numbers": [], "positions": [], "method": method, "charge": float(charge)}
    if scf:
        req["uhf"] = mult - 1
    if "conventions" in v:
        req["conventions"] = CONVENTIONS[method][conventions]
    if "accuracy_scale" in v:
        req["accuracy"] = v["accuracy_scale"]
    if "scf.etemp" in v:
        req["etemp"] = v["scf.etemp"]
    if "scf.max_iter" in v:
        req["max_iter"] = v["scf.max_iter"]
    if "scf.retries" in v:
        req["scf_retries"] = v["scf.retries"]
    if "scf.mixer" in v:
        req["mixer"] = {"broyden-then-diis": "auto"}.get(v["scf.mixer"], v["scf.mixer"])
    if "scf.guess" in v and v["scf.guess"] != "given":
        req["guess"] = v["scf.guess"]
        if session:
            # a named start alone means that start at every geometry: keep the restart policy
            req["restart"] = "previous" if restart[0] else "fresh"
    if "scf.restart" in v:
        req["restart"] = "previous" if v["scf.restart"] == "previous" else "fresh"
    if "solvent.name" in v:
        req["gbsa" if v.get("solvent.model") == "gbsa" else "alpb"] = v["solvent.name"]
        if "solvent.state" in v:
            req["solvation_state"] = v["solvent.state"]
    if composite:
        req["solvation_correction"] = v["solvent.correction"]
        if "solvent.correction_conventions" in v:
            req["correction_conventions"] = CONVENTIONS["gfn2"][v["solvent.correction_conventions"]]
    if "topology.geometry" in v:
        req["topology_positions"] = v["topology.geometry"].tolist()
    if "topology.fragment_charges" in v:
        req["fragment_charges"] = v["topology.fragment_charges"]
    task: dict[str, Any] = {}
    for path, field_name in _TASK_FIELDS.items():
        if path in v:
            value = v[path]
            if path in ("freq.analysis", "freq.thermo"):
                value = value != "off"
            elif path == "freq.conventions":
                value = ANALYSIS_CONVENTIONS[value]
            task[field_name] = value
    given = {p: (x.tolist() if isinstance(x, np.ndarray) else x) for p, x in v.items()}
    contingent = []
    if s.get("scf", {}).get("mixer") == "broyden-then-diis":
        contingent.append("scf.mixer")
    if session and hessian:
        contingent += [p for p in v if p.startswith(("hessian.", "freq."))]
    if session and s.get("scf", {}).get("restart") == "previous":
        contingent.append("scf.restart")
    if "scf.retries" in v:
        contingent.append("scf.retries")
    record: dict[str, Any] = {
        "method": method,
        "conventions": conventions,
        "calctypes": sorted(allowed, key=CALCTYPES.index),
        "numerics": _axtb.numerics(method, composite),
    }
    if composite:
        record["composite"] = {"parts": ["gxtb", "gfn2_solvated", "gfn2_gas"]}
    if contingent:
        record["contingent"] = sorted(set(contingent))
    return Planned(
        method=method,
        conventions=conventions,
        composite=composite,
        calctypes=allowed,
        settings=s,
        reasons=reasons,
        given=given,
        request=req,
        task=task,
        outputs=frozenset(outputs),
        accuracy=accuracy,
        threads=v.get("threads"),
        numbers=numbers,
        positions=positions,
        charge=charge,
        multiplicity=mult,
        record=record,
        restart=restart,
        contingent=sorted(set(contingent)),
    )


#: The driver task fields of the Hessian settings.
_TASK_FIELDS = {
    "hessian.method": "hessian_method",
    "hessian.step": "hessian_step",
    "hessian.accuracy_scale": "hessian_accuracy",
    "freq.analysis": "vibrations",
    "freq.thermo": "thermo",
    "freq.conventions": "analysis_conventions",
    "freq.masses": "masses",
    "freq.temperature": "temperature",
    "freq.pressure": "pressure",
    "freq.scale": "frequency_scale",
    "freq.imaginary_cutoff": "imaginary_cutoff",
    "freq.rotor_cutoff": "rotor_cutoff",
    "freq.symmetry_number": "symmetry_number",
    "freq.symmetry_tol": "symmetry_threshold",
}


# --------------------------------------------------------------------------- threads


def plan_threads(p: Planned, calctype_or_largest: str) -> int:
    """The threads a plan states (`Plan.resources.threads`) for its largest calctype: the
    ``threads`` keyword when given, else the band of the method, calctype and size
    (``axtb_driver::THREAD_BANDS``, `thread_bands`; the composite by its g-xTB method)."""
    if p.threads is not None:
        return p.threads
    return _axtb.planned_threads(p.method, calctype_or_largest == "hessian", len(p.numbers))


def thread_bands() -> str:
    """The plan's thread bands as a Markdown table (``docs/python.md``, generated)."""
    names = {"gfn2": "GFN2-xTB", "gxtb-v2.0.1": "g-xTB 2.0.1", "gxtb-v1.1": "g-xTB v1.1", "gfnff": "GFN-FF"}
    rows = ["| method | calculation | threads by size (atoms) |", "|---|---|---|"]
    for hessian in (False, True):
        for method, h, bands in _axtb.thread_bands():
            if h != hessian:
                continue
            cells = []
            for i, (start, n) in enumerate(bands):
                end = bands[i + 1][0] if i + 1 < len(bands) else None
                span = f"from {start}" if end is None else (f"below {end}" if start == 0 else f"{start}-{end - 1}")
                cells.append(f"{n} ({span})" if len(bands) > 1 else f"{n} (every size)")
            kind = "Hessian" if hessian else "energy, gradient"
            rows.append(f"| {names[method]} | {kind} | {', '.join(cells)} |")
    return "\n".join(rows) + "\n"


def largest(calctypes: Iterable[str]) -> str:
    """The highest derivative among `calctypes`."""
    return max(calctypes, key=CALCTYPES.index)


# --------------------------------------------------------------------------- schema and docs


def _applies_to(key: Key) -> dict[str, Any]:
    methods = sorted(m for m in key.methods if m != "composite")
    out: dict[str, Any] = {"methods": methods, "calctypes": sorted(key.calctypes, key=CALCTYPES.index)}
    if "composite" in key.methods:
        out["composite"] = True
    return out


def settings_schema(calctype: str) -> dict[str, Any]:
    """The JSON Schema (2020-12) of the keywords of `calctype`: every key that applies to it,
    with ``x-unit``, ``x-tier``, ``x-concern`` and ``x-applies-to``."""
    props: dict[str, Any] = {}
    for key in KEYS:
        if key.unsupported is not None or calctype not in key.calctypes:
            continue
        entry = dict(key.schema)
        entry["description"] = key.description + f" Default: {key.default}."
        entry["x-tier"] = key.tier
        entry["x-concern"] = key.concern
        entry["x-shared"] = key.shared
        if key.unit:
            entry["x-unit"] = key.unit
        entry["x-applies-to"] = _applies_to(key)
        group, _, leaf = key.path.rpartition(".")
        if group:
            g = props.setdefault(
                group, {"type": "object", "additionalProperties": False, "properties": {}, "x-tier": "common"}
            )
            g["properties"][leaf] = entry
        else:
            props[leaf] = entry
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "title": f"axtb keywords ({calctype})",
        "type": "object",
        "additionalProperties": False,
        "properties": props,
    }


_METHOD_LABELS = {"gfn2": "G2", "gxtb-v1.1": "X1", "gxtb-v2.0.1": "X2", "gfnff": "FF", "composite": "XC"}


def markdown_table() -> str:
    """The keyword reference as a Markdown table (``docs/python.md``, generated)."""
    rows = [
        "| key | tier | unit | default | methods | calctypes | description |",
        "|---|---|---|---|---|---|---|",
    ]
    for key in KEYS:
        tier = "shared" if key.shared else "axtb"
        if key.unsupported is not None:
            rows.append(f"| `{key.path}` | {tier} | | | refused: kind `unsupported` | | {key.unsupported} |")
            continue
        methods = ", ".join(
            _METHOD_LABELS[m] for m in ("gfn2", "gxtb-v1.1", "gxtb-v2.0.1", "composite", "gfnff") if m in key.methods
        )
        calc = "".join(c[0].upper() for c in CALCTYPES if c in key.calctypes)
        rows.append(
            f"| `{key.path}` | {tier} | {key.unit or ''} | {key.default} | {methods} | {calc} | {key.description} |"
        )
    return "\n".join(rows) + "\n"


def _check_table() -> None:  # pragma: no cover - import-time consistency
    for key in KEYS:
        if key.unsupported is None:
            assert key.methods and key.calctypes, key.path
        assert key.tier in ("common", "expert"), key.path


_check_table()
