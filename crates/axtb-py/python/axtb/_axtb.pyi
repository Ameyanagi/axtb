# SPDX-License-Identifier: LGPL-3.0-or-later
"""Type stubs of the private extension module ``axtb._axtb`` (src/lib.rs, src/engine.rs): the
engine under ``axtb.compute`` and ``axtb.Session``. Not an API."""

from collections.abc import Sequence
from typing import Any

import numpy.typing as npt

__version__: str

ArrayLike2 = npt.ArrayLike
"""Any array-like of shape (N, 3) (or the shape the argument names), converted to float64."""

Result = dict[str, Any]
"""A response: the driver response's fields (``energy``, ``gradient``, ``charges``, ``scf``,
``hessian``, ``vibrations``, ``composite``, ``orbitals``, ...), arrays as numpy arrays."""

ANALYSIS_NUMERICS: str

class Error(Exception):
    """A failed engine call."""

    kind: str
    code: str
    message: str
    fields: dict[str, Any]
    attempts: list[dict[str, Any]] | None
    """The tries of an ``scf_retries`` ladder that did not converge, else None."""

class Engine:
    """A session's calculator: the driver calculator of one request and one committed
    wavefunction per lineage."""

    def __init__(self, request: str, numbers: Sequence[int], positions: ArrayLike2) -> None: ...
    @property
    def method(self) -> str: ...
    @property
    def conventions(self) -> str: ...
    @property
    def numerics(self) -> str: ...
    @property
    def threads(self) -> int: ...
    @property
    def closed(self) -> bool: ...
    def compute(
        self,
        positions: ArrayLike2,
        task: str,
        lineage: str,
        seeds: Sequence[str] = (),
        commit: bool = True,
        accuracy: float | None = None,
        orbitals: bool = False,
    ) -> Result: ...
    def set_given(
        self,
        densities: Sequence[ArrayLike2] | None = None,
        moments: tuple[Sequence[float], ArrayLike2, ArrayLike2] | None = None,
    ) -> None: ...
    def reset(self, key: str | None = None) -> None: ...
    def lineages(self) -> list[str]: ...
    def close(self) -> None: ...

def compute_once(
    request: str,
    numbers: Sequence[int],
    positions: ArrayLike2,
    task: str,
    accuracy: float | None = None,
    orbitals: bool = False,
    densities: Sequence[ArrayLike2] | None = None,
    moments: tuple[Sequence[float], ArrayLike2, ArrayLike2] | None = None,
) -> Result: ...
def canonical(kind: str, name: str, method: str | None = None) -> str: ...
def supported_elements(method: str) -> list[int]: ...
def numerics(method: str, composite: bool) -> str: ...
def single_thread_atoms(method: str, hessian: bool) -> int: ...
def valence_electrons(method: str) -> list[float | None]: ...
def shell_layout(method: str, numbers: Sequence[int]) -> list[tuple[int, int, float, list[float]]] | None: ...
def solvents(model: str) -> list[str]: ...
def canonical_solvent(name: str, model: str) -> str | None: ...
def lapack_warning() -> str | None: ...
