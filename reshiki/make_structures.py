# SPDX-License-Identifier: LGPL-3.0-or-later
"""Generate deterministic, nonstationary single-point geometries in Angstrom.

Run with ``uv run --no-project --with rdkit --with numpy python reshiki/make_structures.py``.
The generated XYZ files are inputs for single-point calculations.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

Array = NDArray[np.float64]
SEED = 42
DISPLACEMENT = np.array([0.018, -0.012, 0.009])
AXES = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], dtype=float)
TETRAHEDRON = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float) / np.sqrt(3)


@dataclass(frozen=True)
class Structure:
    stem: str
    symbols: tuple[str, ...]
    coordinates: Array
    charge: int
    uhf: int
    description: str
    method: str

    def xyz(self) -> str:
        lines = [str(len(self.symbols)), f"charge={self.charge} uhf={self.uhf} {self.description}"]
        for symbol, point in zip(self.symbols, self.coordinates, strict=True):
            clean = np.where(np.abs(point) < 0.5e-10, 0.0, point)
            lines.append(f"{symbol:2s} {clean[0]: .10f} {clean[1]: .10f} {clean[2]: .10f}")
        return "\n".join(lines) + "\n"


def organic(stem: str, smiles: str, description: str, charge: int = 0, uhf: int = 0) -> Structure:
    """Embed one conformer, relax with MMFF94s, then apply a small fixed displacement."""
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid SMILES for {stem}: {smiles}")
    molecule = Chem.AddHs(molecule)
    if Chem.GetFormalCharge(molecule) != charge:
        raise ValueError(f"SMILES charge disagrees with metadata for {stem}")
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = SEED
    parameters.numThreads = 1
    if AllChem.EmbedMolecule(molecule, parameters) != 0:
        raise RuntimeError(f"ETKDG embedding failed for {stem}")
    if not AllChem.MMFFHasAllMoleculeParams(molecule):
        raise RuntimeError(f"MMFF94s parameters unavailable for {stem}")
    status = AllChem.MMFFOptimizeMolecule(molecule, mmffVariant="MMFF94s", maxIters=2000)
    if status != 0:
        raise RuntimeError(f"MMFF94s did not converge for {stem}: status {status}")
    coordinates = np.array(molecule.GetConformer().GetPositions(), dtype=float)
    coordinates[0] += DISPLACEMENT
    coordinates -= coordinates.mean(axis=0)
    return Structure(
        stem,
        tuple(atom.GetSymbol() for atom in molecule.GetAtoms()),
        coordinates,
        charge,
        uhf,
        description,
        "ETKDGv3 seed 42, MMFF94s, fixed 0.0234 A displacement",
    )


def explicit(
    stem: str,
    symbols: list[str],
    coordinates: list[list[float]] | Array,
    description: str,
    method: str,
    charge: int = 0,
    uhf: int = 0,
) -> Structure:
    return Structure(stem, tuple(symbols), np.array(coordinates, dtype=float), charge, uhf, description, method)


def transverse_basis(axis: Array) -> tuple[Array, Array]:
    reference = np.array([0.0, 0.0, 1.0]) if abs(axis[2]) < 0.9 else np.array([0.0, 1.0, 0.0])
    first = np.cross(axis, reference)
    first /= np.linalg.norm(first)
    return first, np.cross(axis, first)


def ammonia_hydrogens(nitrogen: Array, outward: Array, phase: float = 0.0) -> list[Array]:
    first, second = transverse_basis(outward)
    return [
        nitrogen
        + 1.02 * (outward / 3 + np.sqrt(8 / 9) * (np.cos(angle + phase) * first + np.sin(angle + phase) * second))
        for angle in np.arange(3) * (2 * np.pi / 3)
    ]


def build_structures() -> list[Structure]:
    structures = [
        organic("aspirin", "CC(=O)Oc1ccccc1C(=O)O", "Aspirin, neutral acetylsalicylic acid"),
        organic("ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1", "Ibuprofen, neutral carboxylic acid"),
        organic("paracetamol", "CC(=O)Nc1ccc(O)cc1", "Paracetamol, neutral acetaminophen"),
        organic(
            "l-alanine-zwitterion",
            "C[C@H]([NH3+])C(=O)[O-]",
            "L-alanine zwitterion, S stereocenter, net neutral",
        ),
    ]
    water_angle = np.deg2rad(104.5)
    acceptor = np.array([2.90, 0.05, 0.0])
    structures.append(
        explicit(
            "water-dimer",
            ["O", "H", "H", "O", "H", "H"],
            np.array(
                [
                    [0.0, 0.0, 0.0],
                    [0.97, 0.0, 0.0],
                    [0.97 * np.cos(water_angle), 0.97 * np.sin(water_angle), 0.0],
                    acceptor,
                    acceptor + 0.97 * np.array([np.cos(water_angle / 2), 0.0, np.sin(water_angle / 2)]),
                    acceptor + 0.97 * np.array([np.cos(water_angle / 2), 0.0, -np.sin(water_angle / 2)]),
                ]
            ),
            "Hydrogen-bonded water dimer",
            "Explicit, O-O 2.9004 A, O-H 0.97 A, H-O-H 104.5 deg",
        )
    )
    structures.extend(
        [
            organic("bromobenzene", "Brc1ccccc1", "Bromobenzene, aromatic bromide"),
            organic("iodomethane", "CI", "Iodomethane, methyl iodide"),
            explicit(
                "sf6",
                ["S"] + ["F"] * 6,
                np.vstack([np.zeros(3), AXES * 1.59]),
                "Sulfur hexafluoride, octahedral",
                "Explicit Oh, S-F 1.59 A",
            ),
            organic("acetate", "CC(=O)[O-]", "Acetate anion", charge=-1),
            organic("guanidinium", "NC(=[NH2+])N", "Guanidinium cation", charge=1),
            organic("tetramethylammonium", "C[N+](C)(C)C", "Tetramethylammonium cation", charge=1),
            explicit(
                "sulfate",
                ["S"] + ["O"] * 4,
                np.vstack([np.zeros(3), TETRAHEDRON * 1.48]),
                "Sulfate dianion, tetrahedral",
                "Explicit Td, four equal S-O bonds of 1.48 A",
                charge=-2,
            ),
        ]
    )
    methyl_angles = np.arange(3) * (2 * np.pi / 3)
    methyl_h = 1.085 * np.column_stack([np.cos(methyl_angles), np.sin(methyl_angles), np.zeros(3)])
    no2_half_angle = np.deg2rad(67)
    structures.extend(
        [
            explicit(
                "methyl-radical",
                ["C", "H", "H", "H"],
                np.vstack([np.zeros(3), methyl_h]),
                "Methyl radical, planar doublet",
                "Explicit D3h, C-H 1.085 A",
                uhf=1,
            ),
            explicit(
                "no2",
                ["N", "O", "O"],
                [
                    [0.0, 0.0, 0.0],
                    [1.20 * np.cos(no2_half_angle), 1.20 * np.sin(no2_half_angle), 0.0],
                    [1.20 * np.cos(no2_half_angle), -1.20 * np.sin(no2_half_angle), 0.0],
                ],
                "Nitrogen dioxide radical, bent doublet",
                "Explicit, N-O 1.20 A, O-N-O 134 deg",
                uhf=1,
            ),
            organic("tempo", "CC1(C)CCCC(C)(C)N1[O]", "TEMPO nitroxide radical, doublet", uhf=1),
            explicit(
                "o2-triplet",
                ["O", "O"],
                [[-0.615, 0.0, 0.0], [0.615, 0.0, 0.0]],
                "Dioxygen, triplet",
                "Explicit, O-O 1.23 A",
                uhf=2,
            ),
            explicit(
                "co2",
                ["O", "C", "O"],
                [[-1.18, 0.0, 0.0], [0.0, 0.0, 0.0], [1.18, 0.0, 0.0]],
                "Carbon dioxide, linear",
                "Explicit, C-O 1.18 A",
            ),
            explicit(
                "acetylene",
                ["H", "C", "C", "H"],
                [[-1.68, 0.0, 0.0], [-0.61, 0.0, 0.0], [0.61, 0.0, 0.0], [1.68, 0.0, 0.0]],
                "Acetylene, linear",
                "Explicit, C-C 1.22 A, C-H 1.07 A",
            ),
            explicit(
                "hcn",
                ["H", "C", "N"],
                [[-1.07, 0.0, 0.0], [0.0, 0.0, 0.0], [1.17, 0.0, 0.0]],
                "Hydrogen cyanide, linear",
                "Explicit, H-C 1.07 A, C-N 1.17 A",
            ),
        ]
    )
    nickel_coordinates = [np.zeros(3)]
    for direction in TETRAHEDRON:
        nickel_coordinates.extend([direction * 1.85, direction * (1.85 + 1.14)])
    structures.append(
        explicit(
            "nickel-tetracarbonyl",
            ["Ni"] + ["C", "O"] * 4,
            np.array(nickel_coordinates),
            "Nickel tetracarbonyl, tetrahedral Ni(CO)4",
            "Explicit Td, Ni-C 1.85 A, C-O 1.14 A, linear Ni-C-O",
        )
    )
    platinum_coordinates = [np.zeros(3), np.array([2.30, 0.0, 0.0]), np.array([0.0, 2.30, 0.0])]
    platinum_symbols = ["Pt", "Cl", "Cl"]
    for direction in (AXES[1], AXES[3]):
        nitrogen = 2.05 * direction
        platinum_coordinates.extend([nitrogen, *ammonia_hydrogens(nitrogen, direction)])
        platinum_symbols.extend(["N", "H", "H", "H"])
    structures.append(
        explicit(
            "cisplatin",
            platinum_symbols,
            np.array(platinum_coordinates),
            "Cisplatin, cis square-planar Pt(NH3)2Cl2",
            "Explicit square-planar Pt ligands, Pt-Cl 2.30 A, Pt-N 2.05 A, N-H 1.02 A",
        )
    )
    cobalt_coordinates = [np.zeros(3)]
    cobalt_symbols = ["Co"]
    for index, direction in enumerate(AXES):
        nitrogen = 1.97 * direction
        cobalt_coordinates.extend([nitrogen, *ammonia_hydrogens(nitrogen, direction, index * np.pi / 6)])
        cobalt_symbols.extend(["N", "H", "H", "H"])
    structures.append(
        explicit(
            "cobalt-hexammine",
            cobalt_symbols,
            np.array(cobalt_coordinates),
            "Hexamminecobalt(III), low-spin singlet",
            "Explicit octahedral CoN6, Co-N 1.97 A, N-H 1.02 A, tetrahedral NH3 donors",
            charge=3,
        )
    )
    iron_coordinates = [np.zeros(3)]
    iron_symbols = ["Fe"]
    for index, direction in enumerate(AXES):
        oxygen = 2.15 * direction
        first, second = transverse_basis(direction)
        transverse = np.cos(index * np.pi / 6) * first + np.sin(index * np.pi / 6) * second
        bisector = np.cos(water_angle / 2) * direction
        spread = np.sin(water_angle / 2) * transverse
        iron_coordinates.extend([oxygen, oxygen + 0.98 * (bisector + spread), oxygen + 0.98 * (bisector - spread)])
        iron_symbols.extend(["O", "H", "H"])
    structures.append(
        explicit(
            "iron-hexaaqua",
            iron_symbols,
            np.array(iron_coordinates),
            "Hexaaquairon(II), high-spin quintet",
            "Explicit octahedral FeO6, Fe-O 2.15 A, O-H 0.98 A, H-O-H 104.5 deg",
            charge=2,
            uhf=4,
        )
    )
    ethylene_h_x = 0.67 + 1.09 / 2
    ethylene_h_perpendicular = 1.09 * np.sqrt(3) / 2
    structures.extend(
        [
            explicit(
                "twisted-ethylene",
                ["C", "C", "H", "H", "H", "H"],
                [
                    [-0.67, 0.0, 0.0],
                    [0.67, 0.0, 0.0],
                    [-ethylene_h_x, ethylene_h_perpendicular, 0.0],
                    [-ethylene_h_x, -ethylene_h_perpendicular, 0.0],
                    [ethylene_h_x, 0.0, ethylene_h_perpendicular],
                    [ethylene_h_x, 0.0, -ethylene_h_perpendicular],
                ],
                "Ethylene, 90-degree twist, singlet input",
                "Explicit perpendicular CH2 planes, C-C 1.34 A, C-H 1.09 A, H-C-H 120 deg",
            ),
            explicit(
                "uf6",
                ["U"] + ["F"] * 6,
                np.vstack([np.zeros(3), AXES * 1.99]),
                "Uranium hexafluoride, octahedral",
                "Explicit Oh, six U-F bonds of 1.99 A",
            ),
        ]
    )
    return structures


def readme(structures: list[Structure]) -> str:
    rows = [
        "<!-- SPDX-License-Identifier: LGPL-3.0-or-later -->",
        "# Single-point structure fixtures",
        "",
        "These 25 XYZ inputs use **Angstrom** coordinates. The second line is always",
        "`charge=Q uhf=U <description>`, where `uhf` is the number of unpaired electrons",
        "(0 for a singlet, 1 for a doublet, 2 for a triplet, and 4 for a quintet).",
        "The spin values specify the requested calculation inputs.",
        "",
        "Run this command from the repository root to regenerate the files.",
        "",
        "```sh",
        "uv run --no-project --with rdkit --with numpy python reshiki/make_structures.py",
        "```",
        "",
        "Add `--check` to compare regenerated contents without writing, or use",
        "`--output-dir /path/to/fixtures` to generate a separate copy.",
        "The generator requires Python 3.10 or newer. Identical toolkit versions and seed",
        "give deterministic output. RDKit or force-field changes can alter regenerated coordinates.",
        f"This copy was generated with RDKit {rdBase.rdkitVersion} and NumPy {np.__version__}.",
        "",
        "Organics have explicit hydrogens, ETKDGv3 embedding with seed 42 and one thread,",
        "and MMFF94s relaxation with up to 2000 iterations. A fixed",
        "(0.018, -0.012, 0.009) A displacement of the first atom follows relaxation, then",
        "the coordinates are centered. This keeps the geometries close to sensible molecular",
        "shapes while deliberately moving away from force-field stationary points.",
        "Small molecules and metal complexes use the explicit bond lengths and symmetries below.",
        "The geometries are intended for nonstationary xTB single points. The 90-degree",
        "ethylene twist is intentionally strained.",
        "",
        "| File | Atoms | Charge | uhf | What | How built |",
        "| --- | ---: | ---: | ---: | --- | --- |",
    ]
    rows.extend(
        f"| [{item.stem}.xyz]({item.stem}.xyz) | {len(item.symbols)} | {item.charge} | {item.uhf} | "
        f"{item.description} | {item.method} |"
        for item in structures
    )
    rows.append("")
    return "\n".join(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().with_name("structures"))
    parser.add_argument("--check", action="store_true", help="Compare generated files without writing")
    arguments = parser.parse_args()
    structures = build_structures()
    outputs = {f"{structure.stem}.xyz": structure.xyz() for structure in structures}
    outputs["README.md"] = readme(structures)
    if arguments.check:
        different = [
            name
            for name, content in outputs.items()
            if not (arguments.output_dir / name).is_file()
            or (arguments.output_dir / name).read_text(encoding="utf-8") != content
        ]
        if different:
            print("Missing or changed files: " + ", ".join(different))
            return 1
        print(f"Verified generated contents of {len(structures)} XYZ files and README.md.")
        return 0
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    for name, content in outputs.items():
        (arguments.output_dir / name).write_text(content, encoding="utf-8")
    print(f"Generated {len(structures)} XYZ files and README.md in {arguments.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
