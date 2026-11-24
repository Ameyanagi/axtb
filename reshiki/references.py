# SPDX-License-Identifier: LGPL-3.0-or-later
"""Independent, version-pinned reference programs for the Reshiki validation suite.

Structures are ``(symbols, positions_angstrom)`` pairs. All reference calculations
receive a Turbomole ``coord`` file in bohr; returned gradients and Hessians are in
atomic units. tblite is imported only by isolated child interpreters.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Iterable, Sequence
from pathlib import Path

import numpy as np

from axtb_cli import sha256
from cases import ATOMIC_NUMBERS, BOHR_2018, vibrational_mode_count

TIMEOUT = 1800.0
KB_XTB_TBLITE = 3.166808578545117e-6
PROGRAM_VARIABLES = {"xtb": "AXTB_XTB680", "gxtb": "AXTB_GXTB2", "tblite": "RESHIKI_TBLITE_PYTHON"}
EXPECTED_XTB = {"xtb": ("6.8.0", "4b8c987"), "gxtb": ("6.7.1", "26dd68d")}
TBLITE_VERSION = "0.8.0"
FLOAT_PATTERN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
RESULT_PREFIX = "RESHIKI_RESULT "


class ReferenceError(RuntimeError):
    """A reference could not supply a valid result."""


class ReferenceConfigurationError(ReferenceError):
    """A required program is missing, unusable, or has an unapproved version."""


class ReferenceExecutionError(ReferenceError):
    """A launch, signal, ignored option, or incomplete execution failed."""


class ReferenceTimeout(ReferenceExecutionError):
    """A reference exceeded its time limit; this is never a molecular refusal."""


class ReferenceParseError(ReferenceError):
    """A successful reference run produced absent or malformed data."""


class ReferenceRefusal(ReferenceError):
    """The reference itself explicitly refused or failed a calculation.

    ``unsupported_elements`` contains only atomic numbers identified by an actual
    unsupported-element diagnostic. An SCF failure, for example, leaves it empty.
    A caller expecting UF6 refusal must check that 92 occurs in this tuple.
    """

    def __init__(self, message: str, *, unsupported_elements: Iterable[int] = ()) -> None:
        super().__init__(message)
        self.unsupported_elements = tuple(sorted(set(unsupported_elements)))


def _environment(xtb_home: Path | None = None) -> dict[str, str]:
    env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", OMP_STACKSIZE="1G")
    if xtb_home is not None:
        env.pop("XTBPATH", None)
        env["XTBHOME"] = str(xtb_home.resolve())
    return env


def _resolve_program(value: str, variable: str, *, preserve_symlink: bool = False) -> str:
    value = value.strip()
    if not value:
        raise ReferenceConfigurationError(f"{variable} is not set; it must name the required reference executable")
    found = shutil.which(os.path.expanduser(value))
    if found is None:
        raise ReferenceConfigurationError(f"{variable}={value!r} does not name an executable file or command")
    # Executing the resolved target of a venv's Python symlink loses that venv.
    path = Path(found).absolute() if preserve_symlink else Path(found).resolve()
    if not path.is_file() or not os.access(path, os.X_OK):
        raise ReferenceConfigurationError(f"{variable}={value!r} is not an executable file")
    return str(path)


def _preflight(
    command: Sequence[str], label: str, *, directory: Path | None = None
) -> subprocess.CompletedProcess[str]:
    try:
        out = subprocess.run(
            command, capture_output=True, text=True, cwd=directory, env=_environment(directory), timeout=60, check=False
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise ReferenceConfigurationError(f"cannot preflight {label}: {exc}") from exc
    if out.returncode:
        raise ReferenceConfigurationError(f"{label} preflight exited {out.returncode}: {out.stderr.strip()[-1000:]}")
    return out


def reference_programs(required: Iterable[str] = ("xtb", "gxtb", "tblite")) -> tuple[dict[str, str], dict]:
    """Resolve and verify all requested programs, returning paths and provenance.

    Validation calls this with no argument, including when its methods are filtered.
    The geometry preparation tool requests only its required xtb-family binaries.
    There is no PATH fallback for xtb or g-xTB: only the named environment variables
    are honored. tblite's Python defaults to the running interpreter.
    """
    names = tuple(dict.fromkeys(required))
    unknown = set(names) - PROGRAM_VARIABLES.keys()
    if unknown:
        raise ValueError(f"unknown reference programs: {', '.join(sorted(unknown))}")
    programs, metadata = {}, {}
    for name in names:
        variable = PROGRAM_VARIABLES[name]
        value = os.environ.get(variable, sys.executable if name == "tblite" else "")
        binary = _resolve_program(value, variable, preserve_symlink=name == "tblite")
        if name in EXPECTED_XTB:
            with _scratch() as scratch:
                out = _preflight([binary, "--version"], variable, directory=Path(scratch))
            match = re.search(r"\bxtb version\s+(\S+)\s+\(([0-9a-fA-F]+)\)", out.stdout + out.stderr)
            expected_version, expected_commit = EXPECTED_XTB[name]
            if not match or match.group(1) != expected_version or not match.group(2).startswith(expected_commit):
                actual = f"{match.group(1)} ({match.group(2)})" if match else "unrecognized version output"
                raise ReferenceConfigurationError(
                    f"{variable} must be xtb {expected_version} ({expected_commit}); found {actual}"
                )
            version = f"{match.group(1)} ({match.group(2)})"
            info = {"version": version}
        else:
            code = (
                "import importlib.metadata as m, json, sys; import tblite.interface; "
                "print('RESHIKI_RESULT '+json.dumps({'version':m.version('tblite'),"
                "'python_version':sys.version.split()[0]}))"
            )
            out = _preflight([binary, "-c", code], f"tblite ({variable})")
            try:
                info = _child_result(out.stdout)
            except ReferenceParseError as exc:
                raise ReferenceConfigurationError(f"tblite preflight returned invalid version data: {exc}") from exc
            if info.get("version") != TBLITE_VERSION:
                raise ReferenceConfigurationError(
                    f"{variable} must provide tblite {TBLITE_VERSION}; found {info.get('version')!r}"
                )
        programs[name] = binary
        metadata[name] = {**info, "path": binary, "sha256": sha256(binary), "environment_variable": variable}
    return programs, metadata


def _structure(structure: tuple[Sequence[str], np.ndarray]) -> tuple[list[str], np.ndarray, np.ndarray]:
    symbols, positions = structure
    symbols = [str(symbol).capitalize() for symbol in symbols]
    try:
        numbers = np.array([ATOMIC_NUMBERS[symbol] for symbol in symbols], dtype=int)
    except KeyError as exc:
        raise ValueError(f"unknown element symbol {exc.args[0]!r}") from exc
    positions = np.asarray(positions, dtype=float)
    if not symbols or positions.shape != (len(symbols), 3) or not np.all(np.isfinite(positions)):
        raise ValueError("a structure needs nonempty symbols and finite positions with shape (natoms, 3)")
    return symbols, positions / BOHR_2018, numbers


def _write_coord(path: Path, symbols: Sequence[str], positions: np.ndarray) -> None:
    lines = ["$coord"]
    lines.extend(f"{x:.17g} {y:.17g} {z:.17g} {symbol.lower()}" for symbol, (x, y, z) in zip(symbols, positions))
    path.write_text("\n".join([*lines, "$end", ""]), encoding="utf-8")


def _scratch() -> tempfile.TemporaryDirectory:
    return tempfile.TemporaryDirectory(prefix="reshiki-reference-", dir=os.environ.get("AXTB_SCRATCH"))


def _unsupported_elements(text: str, numbers: Iterable[int]) -> tuple[int, ...]:
    elements = {int(z) for z in re.findall(r"no basis found for atom\s+\d+\s+Z\s*=\s*(\d+)", text, re.IGNORECASE)}
    for pattern in (
        r"(?:no parameters|no parametrization|unsupported element)[^\n]*?(?:Z\s*=\s*|element\s+)(\d+)",
        r"no support for elements with Z\s*>\s*(\d+)",
    ):
        for match in re.finditer(pattern, text, re.IGNORECASE):
            if "no support for elements" in match.group(0).lower():
                elements.update(int(z) for z in numbers if int(z) > int(match.group(1)))
            else:
                elements.add(int(match.group(1)))
    return tuple(sorted(elements))


def _execute(command: list[str], directory: Path, numbers: np.ndarray, timeout: float) -> str:
    try:
        out = subprocess.run(
            command,
            cwd=directory,
            capture_output=True,
            text=True,
            env=_environment(directory),
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise ReferenceTimeout(f"{Path(command[0]).name} timed out after {timeout:g} s") from exc
    except OSError as exc:
        raise ReferenceExecutionError(f"could not execute {command[0]}: {exc}") from exc
    text = out.stdout + "\n" + out.stderr
    message = f"{Path(command[0]).name} exited {out.returncode}: {text.strip()[-3000:]}"
    if re.search(r"unknown option|unrecognized (?:option|argument)|unknown flag", text, re.IGNORECASE):
        raise ReferenceExecutionError(message)
    if out.returncode < 0:
        raise ReferenceExecutionError(message)
    refused = re.search(r"abnormal termination|\[ERROR\]|SCF not converged|SCC.*not converged", text, re.IGNORECASE)
    unsupported = _unsupported_elements(text, numbers)
    if refused or unsupported:
        raise ReferenceRefusal(message, unsupported_elements=unsupported)
    if out.returncode or not re.search(r"\bnormal termination\b", text, re.IGNORECASE):
        raise ReferenceExecutionError(message)
    return out.stdout


def _read(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        raise ReferenceParseError(f"cannot read reference {path.name}: {exc}") from exc


def _float(value: str) -> float:
    # ifort omits the exponent letter for three-digit exponents in fixed fields.
    value = re.sub(r"(?<=[\d.])([+-]\d{3})$", r"E\1", value.replace("D", "E").replace("d", "e"))
    result = float(value)
    if not np.isfinite(result):
        raise ValueError("non-finite reference value")
    return result


def _array(value, shape: tuple[int, ...], name: str) -> np.ndarray:
    try:
        result = np.asarray(value, dtype=float)
    except (ValueError, TypeError) as exc:
        raise ReferenceParseError(f"invalid {name}: {exc}") from exc
    if result.shape != shape or not np.all(np.isfinite(result)):
        raise ReferenceParseError(f"{name} has shape {result.shape}, expected finite {shape}")
    return result


def _energy(stdout: str, method: str) -> float:
    # Optimizations and Hessians may print both forms repeatedly; use the last one
    # in output order, including the final boxed energy's full printed precision.
    matches = re.findall(rf"(?:\bTOTAL ENERGY|::\s*total energy)\s+({FLOAT_PATTERN})\s+Eh", stdout)
    if not matches:
        raise ReferenceParseError("reference output has no total energy")
    try:
        return _float(matches[-1])
    except ValueError as exc:
        raise ReferenceParseError(f"invalid reference total energy: {exc}") from exc


def _gradient(path: Path, nat: int) -> np.ndarray:
    text = _read(path)
    blocks = re.findall(r"\$grad\b(.*?)(?=\$end|\Z)", text, re.DOTALL | re.IGNORECASE)
    if not blocks:
        raise ReferenceParseError("gradient file has no $grad section")
    # An optimization can append cycles. A cycle holds one header, N coordinates,
    # and N gradient rows; read the final cycle, never the coordinate rows.
    cycles = re.split(r"\bcycle\s*=", blocks[-1], flags=re.IGNORECASE)
    lines = [line.split() for line in cycles[-1].splitlines() if line.strip()]
    rows = lines[1:]
    if len(rows) != 2 * nat or any(len(row) != 3 for row in rows[nat:]):
        raise ReferenceParseError("gradient file has an invalid final cycle or atom count")
    try:
        return _array([[_float(value) for value in row] for row in rows[nat:]], (nat, 3), "gradient")
    except ValueError as exc:
        raise ReferenceParseError(f"invalid reference gradient: {exc}") from exc


def _iterations(stdout: str, method: str) -> int | None:
    if method == "gfnff":
        return None
    if method == "gfn2":
        matches = re.findall(r"convergence criteria satisfied after\s+(\d+)\s+iterations", stdout)
        if matches:
            return int(matches[-1])
    else:
        tables = re.findall(
            r"cycle\s+total energy\s+energy error\s+density error\s*\n-+\s*\n(.*?)\n-+", stdout, re.DOTALL
        )
        if tables:
            cycles = re.findall(
                rf"^\s*(\d+)\s+{FLOAT_PATTERN}\s+{FLOAT_PATTERN}\s+{FLOAT_PATTERN}\s*$", tables[-1], re.MULTILINE
            )
            if cycles and [int(cycle) for cycle in cycles] == list(range(1, len(cycles) + 1)):
                return int(cycles[-1])
    return None


def _values(directory: Path, stdout: str, method: str, nat: int) -> dict:
    result = {
        "energy": _energy(stdout, method),
        "gradient": _gradient(directory / "gradient", nat),
        "iterations": _iterations(stdout, method),
    }
    if method == "gfnff":
        # Empty list entries can be printed as ****** by the Fortran JSON writer.
        try:
            lists = json.loads(re.sub(r"\*{3,}", "null", _read(directory / "gfnff_lists.json")))
            result["charges"] = _array([row[0] for row in lists["eeq"]], (nat,), "EEQ charges")
        except (ValueError, TypeError, KeyError, IndexError) as exc:
            raise ReferenceParseError(f"invalid GFN-FF EEQ charge data: {exc}") from exc
    elif method == "gfn2" or (directory / "charges").exists():
        try:
            values = [_float(value) for value in _read(directory / "charges").split()]
            result["charges"] = _array(values, (nat,), "charges")
        except ValueError as exc:
            raise ReferenceParseError(f"invalid reference charges: {exc}") from exc
    else:
        # g-xTB prints the atomic charges in its population table (make_reference_v2.py).
        table = re.search(r"#\s+Z\s+q\s+p\(s\).*?\n-+\n(.*?)\n-+", stdout, re.DOTALL)
        if table:
            charges = [_float(line.split()[3]) for line in table[1].splitlines() if line.strip()]
            result["charges"] = _array(charges, (nat,), "g-xTB charges")
    if method == "gfn2":
        data = json.loads(_read(directory / "xtbout.json"))
        result["dipole"] = _array(data["dipole / a.u."], (3,), "dipole")
    return result


def parse_molden(text: str) -> dict[str, np.ndarray]:
    """Read orbital energies by the written Alpha/Beta spin channel in ``[MO]``.

    Channel order and orbital order are preserved. Missing channels are not invented
    and channels are never combined: degeneracies across spins are unrelated.
    """
    channels: dict[str, list[float]] = {}
    section, energy, spin = "", None, None

    def finish() -> None:
        nonlocal energy, spin
        if energy is None and spin is None:
            return
        if energy is None or spin is None:
            raise ReferenceParseError("Molden orbital needs both Ene= and Spin=")
        channels.setdefault(spin, []).append(energy)
        energy, spin = None, None

    for raw in text.splitlines():
        line = raw.strip()
        if line.startswith("["):
            if section == "mo":
                finish()
            match = re.match(r"\[([^]]+)\]", line)
            section = match.group(1).strip().lower() if match else ""
            continue
        if section != "mo" or "=" not in line:
            continue
        key, value = (part.strip() for part in line.split("=", 1))
        key = key.lower()
        if key == "sym":
            finish()
        elif key == "ene":
            if energy is not None:
                finish()
            try:
                energy = _float(value)
            except ValueError as exc:
                raise ReferenceParseError(f"invalid Molden orbital energy {value!r}") from exc
        elif key == "spin":
            if spin is not None:
                raise ReferenceParseError("duplicate Molden Spin= before an orbital energy")
            spin = value.lower()
            if spin not in ("alpha", "beta"):
                raise ReferenceParseError(f"unrecognized Molden spin {value!r}")
    if section == "mo":
        finish()
    if not channels:
        raise ReferenceParseError("Molden file has no orbital energies in [MO]")
    return {channel: np.array(energies) for channel, energies in channels.items()}


def read_molden(text: str | Path) -> dict[str, np.ndarray]:
    """Parse Molden text, or read a ``Path`` explicitly supplied by the caller."""
    return parse_molden(_read(text) if isinstance(text, Path) else text)


def _rotmol(positions: np.ndarray) -> np.ndarray:
    """The g-xTB xtb 6.7.1 Hessian's small rotation, in its actual order."""
    result = positions.copy()
    for (i, j), angle in (((1, 2), 1e-4), ((0, 2), 2e-4), ((0, 1), 3e-4)):
        radians = np.deg2rad(angle)
        x, y = result[:, i].copy(), result[:, j].copy()
        result[:, i] = x * np.cos(radians) + y * np.sin(radians)
        result[:, j] = -x * np.sin(radians) + y * np.cos(radians)
    return result


def _hessian(directory: Path, stdout: str, method: str, positions: np.ndarray) -> dict:
    n3 = 3 * len(positions)
    text = _read(directory / "hessian")
    match = re.search(r"\$hessian\b(.*?)(?=\$end|\Z)", text, re.DOTALL | re.IGNORECASE)
    if not match:
        raise ReferenceParseError("reference Hessian has no $hessian section")
    try:
        values = [_float(value) for value in match.group(1).split()]
        if len(values) != n3 * n3:
            raise ValueError(f"has {len(values)} values, expected {n3 * n3}")
        hessian = _array(np.array(values).reshape(n3, n3), (n3, n3), "Hessian")
    except ValueError as exc:
        raise ReferenceParseError(f"invalid reference Hessian: {exc}") from exc
    if method == "gfnff":
        # vibspectrum prints 0.01 cm^-1, including projected zero modes first.
        frequencies = []
        for line in _read(directory / "vibspectrum").splitlines():
            row = line.split()
            if not row or row[0].startswith(("$", "#")):
                continue
            if not row[0].isdigit() or len(row) < 4:
                raise ReferenceParseError("invalid GFN-FF vibspectrum row")
            try:
                frequencies.append(_float(row[-3]))
            except ValueError as exc:
                raise ReferenceParseError("invalid GFN-FF vibrational frequency") from exc
        if len(frequencies) != n3:
            raise ReferenceParseError(f"vibspectrum has {len(frequencies)} modes, expected {n3}")
        frequencies = frequencies[n3 - vibrational_mode_count(positions) :]
    else:
        frequencies = []
        for line in _read(directory / "g98.out").splitlines():
            match = re.match(r"\s*Frequencies\s+--\s+(.*)", line)
            if match:
                try:
                    frequencies.extend(_float(value) for value in match.group(1).split())
                except ValueError as exc:
                    raise ReferenceParseError("invalid g98.out vibrational frequency") from exc
    return {
        "energy": _energy(stdout, method),
        "hessian": hessian,
        "frequencies": np.array(frequencies),
        "positions": _rotmol(positions) if method == "gxtb" else positions.copy(),
        "iterations": _iterations(stdout, method),
    }


def _arguments(method: str, check: str, charge: int, uhf: int | None) -> list[str]:
    if method not in ("gfn2", "gfnff", "gxtb"):
        raise ValueError(f"unknown reference method {method!r}")
    if check not in ("values", "molden", "hessian", "optimize"):
        raise ValueError(f"unknown reference check {check!r}")
    if method == "gfnff" and check == "molden":
        raise ValueError("GFN-FF has no molecular orbitals")
    args = ["--gfn", "2"] if method == "gfn2" else ["--gfnff" if method == "gfnff" else "--gxtb"]
    if check == "values":
        args += ["--grad", "--json"]
        if method == "gxtb":
            args += ["--norestart"]
    elif check == "hessian":
        args += ["--hess"]
    elif check == "molden":
        args += ["--molden", "--norestart"]
    else:
        args += ["--opt", "tight"]
    if method != "gfnff" and check != "optimize":
        accuracy = "0.001" if method == "gfn2" and check == "hessian" else "1e-4"
        if method == "gxtb" and check == "molden":
            # Upstream Molden comparisons use the g-xTB program's default SCF.
            accuracy = "1.0"
        args += ["--acc", accuracy]
    args += ["--chrg", str(charge)]
    if uhf is not None and method != "gfnff":
        args += ["--uhf", str(uhf)]
    return args


def run_xtb(
    binary: str,
    method: str,
    structure: tuple[Sequence[str], np.ndarray],
    charge: int,
    uhf: int | None,
    check: str,
    *,
    timeout: float = TIMEOUT,
) -> dict:
    """Run one reference value, Molden, or Hessian calculation in fresh scratch."""
    symbols, positions, numbers = _structure(structure)
    args = _arguments(method, check, charge, uhf)
    if check == "optimize":
        raise ValueError("use optimize_structure() for geometry optimization")
    with _scratch() as scratch:
        directory = Path(scratch)
        _write_coord(directory / "coord", symbols, positions)
        stdout = _execute([binary, "coord", *args], directory, numbers, timeout)
        if check == "values":
            return _values(directory, stdout, method, len(symbols))
        if check == "molden":
            return {"orbitals": read_molden(directory / "molden.input"), "iterations": _iterations(stdout, method)}
        return _hessian(directory, stdout, method, positions)


_TBLITE_CHILD = r"""
import json
import sys
import numpy as np
from tblite.interface import Calculator

request = json.load(sys.stdin)
try:
    calculator = Calculator("GFN2-xTB", np.array(request["numbers"], dtype=int),
                            np.array(request["positions"], dtype=float),
                            charge=float(request["charge"]), uhf=request["uhf"])
    calculator.set("verbosity", 0)
    calculator.set("accuracy", 1e-6)
    calculator.set("temperature", 300.0 * 3.166808578545117e-6)
    result = calculator.singlepoint()
    response = {"ok": True, "energy": float(result.get("energy"))}
    for key in ("gradient", "charges", "dipole"):
        response[key] = np.asarray(result.get(key)).tolist()
except Exception as error:
    response = {"ok": False, "error": f"{type(error).__name__}: {error}"}
print("RESHIKI_RESULT " + json.dumps(response, allow_nan=False), flush=True)
"""


def _child_result(stdout: str) -> dict:
    lines = [line[len(RESULT_PREFIX) :] for line in stdout.splitlines() if line.startswith(RESULT_PREFIX)]
    if len(lines) != 1:
        raise ReferenceParseError(f"child returned {len(lines)} result records, expected one")
    try:
        result = json.loads(lines[0])
    except ValueError as exc:
        raise ReferenceParseError(f"invalid child JSON: {exc}") from exc
    if not isinstance(result, dict):
        raise ReferenceParseError("child result is not an object")
    return result


def run_tblite(
    python: str,
    structure: tuple[Sequence[str], np.ndarray],
    charge: int,
    uhf: int | None,
    *,
    timeout: float = TIMEOUT,
) -> dict:
    """Run tblite GFN2-xTB in a child, so native failures cannot kill the suite."""
    symbols, positions, numbers = _structure(structure)
    request = {"numbers": numbers.tolist(), "positions": positions.tolist(), "charge": charge, "uhf": uhf}
    with _scratch() as scratch:
        directory = Path(scratch)
        try:
            out = subprocess.run(
                [python, "-c", _TBLITE_CHILD],
                input=json.dumps(request),
                capture_output=True,
                text=True,
                cwd=directory,
                env=_environment(),
                timeout=timeout,
                check=False,
            )
        except subprocess.TimeoutExpired as exc:
            raise ReferenceTimeout(f"tblite timed out after {timeout:g} s") from exc
        except OSError as exc:
            raise ReferenceExecutionError(f"could not launch tblite Python: {exc}") from exc
    if out.returncode:
        raise ReferenceExecutionError(f"tblite child exited {out.returncode}: {out.stderr.strip()[-2000:]}")
    result = _child_result(out.stdout)
    if result.get("ok") is False:
        error = str(result.get("error", "tblite refused without a diagnostic"))
        raise ReferenceRefusal(error, unsupported_elements=_unsupported_elements(error, numbers))
    if result.get("ok") is not True:
        raise ReferenceParseError("tblite child result has no boolean success status")
    try:
        energy = float(result["energy"])
        if not np.isfinite(energy):
            raise ValueError("non-finite energy")
        return {
            "energy": energy,
            "gradient": _array(result["gradient"], (len(symbols), 3), "tblite gradient"),
            "charges": _array(result["charges"], (len(symbols),), "tblite charges"),
            "dipole": _array(result["dipole"], (3,), "tblite dipole"),
        }
    except (KeyError, ValueError, TypeError) as exc:
        raise ReferenceParseError(f"invalid tblite result: {exc}") from exc


def optimize_structure(
    binary: str,
    method: str,
    structure: tuple[Sequence[str], np.ndarray],
    charge: int,
    uhf: int | None,
    *,
    timeout: float = TIMEOUT,
) -> dict:
    """Optimize with the reference's ``--opt tight``; return geometry in Å."""
    symbols, positions, numbers = _structure(structure)
    args = _arguments(method, "optimize", charge, uhf)
    with _scratch() as scratch:
        directory = Path(scratch)
        _write_coord(directory / "coord", symbols, positions)
        stdout = _execute([binary, "coord", *args], directory, numbers, timeout)
        if (directory / "NOT_CONVERGED").exists() or not re.search(
            r"GEOMETRY OPTIMIZATION CONVERGED", stdout, re.IGNORECASE
        ):
            raise ReferenceRefusal("reference geometry optimization did not converge")
        text = _read(directory / "xtbopt.coord")
        match = re.search(r"\$coord\s*\n(.*?)(?=\$|\Z)", text, re.DOTALL | re.IGNORECASE)
        if not match:
            raise ReferenceParseError("optimized coordinate file has no bohr $coord block")
        rows = [line.split() for line in match.group(1).splitlines() if line.strip()]
        if len(rows) != len(symbols) or [row[-1].capitalize() for row in rows] != symbols:
            raise ReferenceParseError("optimized coordinate file changed atom count or atom order")
        try:
            values = [[_float(value) for value in row[:3]] for row in rows]
            optimized = _array(values, positions.shape, "optimized positions")
        except ValueError as exc:
            raise ReferenceParseError(f"invalid optimized coordinates: {exc}") from exc
        return {
            "symbols": symbols,
            "positions": optimized * BOHR_2018,
            "energy": _energy(stdout, method),
        }
