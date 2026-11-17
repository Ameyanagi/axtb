# SPDX-License-Identifier: LGPL-3.0-or-later
"""Errors of the Python interface: `qcdata.ProgramError` with axtb's kinds and codes.

Planning refuses with kind ``input`` (codes ``unknown_key``, ``invalid_value``, ``no_effect``,
``not_applicable``, ``missing_setting``, ``invalid_structure``, ``unsupported_element``,
``not_modeled``, ``session_closed``) or ``unsupported`` (``unsupported_setting``,
``unsupported_calctype``, ``not_applicable`` for an output axtb cannot give). A driver error
(``axtb._axtb.Error``) maps by its kind (`from_driver`); its fields go to ``ErrorInfo.details``
where qcdata has them, and always to the record's ``error`` group.
"""

from __future__ import annotations

from typing import Any

from qcdata import ErrorInfo, ErrorKind, ProgramError

__all__ = ["HAS_DETAILS", "error_info", "from_driver", "refuse"]

#: Whether this qcdata's `ErrorInfo` has ``details`` (qcdata 0.24.0a5 and later).
HAS_DETAILS = "details" in ErrorInfo.model_fields


def error_info(kind: str, message: str, code: str | None = None, details: dict[str, Any] | None = None) -> ErrorInfo:
    """An `ErrorInfo` of axtb, with ``details`` where qcdata has the field."""
    fields: dict[str, Any] = {"kind": ErrorKind(kind), "message": message, "code": code, "program": "axtb"}
    if details and HAS_DETAILS:
        fields["details"] = details
    return ErrorInfo(**fields)


def refuse(kind: str, code: str | None, message: str, details: dict[str, Any] | None = None) -> ProgramError:
    """The `ProgramError` of a refused input (raise it)."""
    return ProgramError(error_info(kind, message, code, details))


# the driver's error kinds: (qcdata kind, code)
_DRIVER = {
    "input": ("input", "invalid_input"),
    "unsupported_element": ("input", "unsupported_element"),
    "topology": ("input", "topology_failed"),
    "scf_not_converged": ("convergence", "scf_not_converged"),
    "numerical": ("numerical", "numerical_failure"),
    "not_applicable": ("unsupported", "not_applicable"),
    "structure": ("input", "invalid_structure"),
    "panic": ("internal", "panic"),
    "closed": ("input", "session_closed"),
}


def from_driver(e: Any) -> tuple[ErrorInfo, dict[str, Any]]:
    """The `ErrorInfo` of a driver error (``axtb._axtb.Error``) and its fields (the record's
    ``error`` group: ``z``; ``iterations``, ``energy_change``, ``residual``)."""
    kind, code = _DRIVER.get(e.kind, ("internal", e.code))
    fields = dict(e.fields)
    return error_info(kind, e.message, code, fields), fields
