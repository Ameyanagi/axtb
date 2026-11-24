<!-- SPDX-License-Identifier: LGPL-3.0-or-later -->
# Single-point structure fixtures

These 25 XYZ inputs use **Angstrom** coordinates. The second line is always
`charge=Q uhf=U <description>`, where `uhf` is the number of unpaired electrons
(0 for a singlet, 1 for a doublet, 2 for a triplet, and 4 for a quintet).
The spin values specify the requested calculation inputs.

Run this command from the repository root to regenerate the files.

```sh
uv run --no-project --with rdkit --with numpy python reshiki/make_structures.py
```

Add `--check` to compare regenerated contents without writing, or use
`--output-dir /path/to/fixtures` to generate a separate copy.
The generator requires Python 3.10 or newer. Identical toolkit versions and seed
give deterministic output. RDKit or force-field changes can alter regenerated coordinates.
This copy was generated with RDKit 2026.09.1 and NumPy 2.5.5.

Organics have explicit hydrogens, ETKDGv3 embedding with seed 42 and one thread,
and MMFF94s relaxation with up to 2000 iterations. A fixed
(0.018, -0.012, 0.009) A displacement of the first atom follows relaxation, then
the coordinates are centered. This keeps the geometries close to sensible molecular
shapes while deliberately moving away from force-field stationary points.
Small molecules and metal complexes use the explicit bond lengths and symmetries below.
The geometries are intended for nonstationary xTB single points. The 90-degree
ethylene twist is intentionally strained.

| File | Atoms | Charge | uhf | What | How built |
| --- | ---: | ---: | ---: | --- | --- |
| [aspirin.xyz](aspirin.xyz) | 21 | 0 | 0 | Aspirin, neutral acetylsalicylic acid | ETKDGv3 seed 42, MMFF94s, fixed 0.0234 A displacement |
| [ibuprofen.xyz](ibuprofen.xyz) | 33 | 0 | 0 | Ibuprofen, neutral carboxylic acid | ETKDGv3 seed 42, MMFF94s, fixed 0.0234 A displacement |
| [paracetamol.xyz](paracetamol.xyz) | 20 | 0 | 0 | Paracetamol, neutral acetaminophen | ETKDGv3 seed 42, MMFF94s, fixed 0.0234 A displacement |
| [l-alanine-zwitterion.xyz](l-alanine-zwitterion.xyz) | 13 | 0 | 0 | L-alanine zwitterion, S stereocenter, net neutral | ETKDGv3 seed 42, MMFF94s, fixed 0.0234 A displacement |
| [water-dimer.xyz](water-dimer.xyz) | 6 | 0 | 0 | Hydrogen-bonded water dimer | Explicit, O-O 2.9004 A, O-H 0.97 A, H-O-H 104.5 deg |
| [bromobenzene.xyz](bromobenzene.xyz) | 12 | 0 | 0 | Bromobenzene, aromatic bromide | ETKDGv3 seed 42, MMFF94s, fixed 0.0234 A displacement |
| [iodomethane.xyz](iodomethane.xyz) | 5 | 0 | 0 | Iodomethane, methyl iodide | ETKDGv3 seed 42, MMFF94s, fixed 0.0234 A displacement |
| [sf6.xyz](sf6.xyz) | 7 | 0 | 0 | Sulfur hexafluoride, octahedral | Explicit Oh, S-F 1.59 A |
| [acetate.xyz](acetate.xyz) | 7 | -1 | 0 | Acetate anion | ETKDGv3 seed 42, MMFF94s, fixed 0.0234 A displacement |
| [guanidinium.xyz](guanidinium.xyz) | 10 | 1 | 0 | Guanidinium cation | ETKDGv3 seed 42, MMFF94s, fixed 0.0234 A displacement |
| [tetramethylammonium.xyz](tetramethylammonium.xyz) | 17 | 1 | 0 | Tetramethylammonium cation | ETKDGv3 seed 42, MMFF94s, fixed 0.0234 A displacement |
| [sulfate.xyz](sulfate.xyz) | 5 | -2 | 0 | Sulfate dianion, tetrahedral | Explicit Td, four equal S-O bonds of 1.48 A |
| [methyl-radical.xyz](methyl-radical.xyz) | 4 | 0 | 1 | Methyl radical, planar doublet | Explicit D3h, C-H 1.085 A |
| [no2.xyz](no2.xyz) | 3 | 0 | 1 | Nitrogen dioxide radical, bent doublet | Explicit, N-O 1.20 A, O-N-O 134 deg |
| [tempo.xyz](tempo.xyz) | 29 | 0 | 1 | TEMPO nitroxide radical, doublet | ETKDGv3 seed 42, MMFF94s, fixed 0.0234 A displacement |
| [o2-triplet.xyz](o2-triplet.xyz) | 2 | 0 | 2 | Dioxygen, triplet | Explicit, O-O 1.23 A |
| [co2.xyz](co2.xyz) | 3 | 0 | 0 | Carbon dioxide, linear | Explicit, C-O 1.18 A |
| [acetylene.xyz](acetylene.xyz) | 4 | 0 | 0 | Acetylene, linear | Explicit, C-C 1.22 A, C-H 1.07 A |
| [hcn.xyz](hcn.xyz) | 3 | 0 | 0 | Hydrogen cyanide, linear | Explicit, H-C 1.07 A, C-N 1.17 A |
| [nickel-tetracarbonyl.xyz](nickel-tetracarbonyl.xyz) | 9 | 0 | 0 | Nickel tetracarbonyl, tetrahedral Ni(CO)4 | Explicit Td, Ni-C 1.85 A, C-O 1.14 A, linear Ni-C-O |
| [cisplatin.xyz](cisplatin.xyz) | 11 | 0 | 0 | Cisplatin, cis square-planar Pt(NH3)2Cl2 | Explicit square-planar Pt ligands, Pt-Cl 2.30 A, Pt-N 2.05 A, N-H 1.02 A |
| [cobalt-hexammine.xyz](cobalt-hexammine.xyz) | 25 | 3 | 0 | Hexamminecobalt(III), low-spin singlet | Explicit octahedral CoN6, Co-N 1.97 A, N-H 1.02 A, tetrahedral NH3 donors |
| [iron-hexaaqua.xyz](iron-hexaaqua.xyz) | 19 | 2 | 4 | Hexaaquairon(II), high-spin quintet | Explicit octahedral FeO6, Fe-O 2.15 A, O-H 0.98 A, H-O-H 104.5 deg |
| [twisted-ethylene.xyz](twisted-ethylene.xyz) | 6 | 0 | 0 | Ethylene, 90-degree twist, singlet input | Explicit perpendicular CH2 planes, C-C 1.34 A, C-H 1.09 A, H-C-H 120 deg |
| [uf6.xyz](uf6.xyz) | 7 | 0 | 0 | Uranium hexafluoride, octahedral | Explicit Oh, six U-F bonds of 1.99 A |
