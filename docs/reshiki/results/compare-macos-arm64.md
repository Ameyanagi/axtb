# axtb platform comparison — FAIL

BASE (Linux x64 reference): `linux-x64.json` — Linux-7.0.14-arch1-1-x86_64-with-glibc2.43 / x86_64 / AMD Ryzen 9 7940HS w/ Radeon 780M Graphics
OTHER: `macos-arm64.json` — macOS-27.0.1-arm64-arm-64bit-Mach-O / arm64 / Apple M1

Limits are 1e-8 Eh for |dE|, 1e-7 Eh/bohr for max|dG|, and 0.1 cm⁻¹ for max sorted |dν|. Charges and iteration counts are measured without gates. Input fingerprints must match, selected runs must be present, and at least one numeric pair is required. Only declared expected error pairs are accepted.

| Method / run | Runs | Numeric pairs | Matching errors | Max ΔE (Eh) | Max ΔG (Eh/bohr) | Max Δν (cm⁻¹) | Max Δq (e, ungated) | SCF iterations BASE / OTHER | Iteration changes | Failed runs |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| gfn2 / gradient | 40 | 39 | 1 | 1.264e-04 | 9.440e-04 | — | 9.111e-04 | 570 / 586 | 1 | 1 |
| gfn2 / hessian-analytic | 9 | 9 | 0 | 2.132e-14 | 1.043e-14 | 2.519e-10 | 6.205e-13 | 18 / 18 | 0 | 0 |
| gfnff / gradient | 40 | 40 | 0 | 2.132e-14 | 7.286e-15 | — | 6.513e-14 | — | 0 | 0 |
| gfnff / hessian-analytic | 9 | 9 | 0 | 8.882e-16 | 5.996e-16 | 2.354e-11 | 9.614e-15 | — | 0 | 0 |
| gxtb / gradient | 40 | 40 | 0 | 1.839e-03 | 1.115e-01 | — | 4.651e-02 | 557 / 504 | 2 | 2 |
| gxtb / hessian-analytic | 9 | 9 | 0 | 1.137e-13 | 2.467e-15 | 6.048e-09 | 3.708e-14 | 89 / 89 | 0 | 0 |
| **Total** | 147 | 146 | 1 | 1.839e-03 | 1.115e-01 | 6.048e-09 | 4.651e-02 | 1234 / 1197 | 3 | 3 |

Violations **3**. Changed diagnostics for expected error pairs: 0 (ungated).

**Violations**

- gh106/gfn2/gradient — energy delta 1.264e-04 exceeds 1e-08, gradient delta 9.440e-04 exceeds 1e-07
- gh106/gxtb/gradient — energy delta 1.839e-03 exceeds 1e-08, gradient delta 1.842e-02 exceeds 1e-07
- twisted-ethylene/gxtb/gradient — gradient delta 1.115e-01 exceeds 1e-07

**Iteration changes (ungated)**

- gh106/gfn2/gradient: 195 -> 211
- gh106/gxtb/gradient: 177 -> 123
- twisted-ethylene/gxtb/gradient: 9 -> 10

**Every run**

| Case / method / run | ΔE (Eh) | max ΔG (Eh/bohr) | max Δν (cm⁻¹) | Status |
|---|---:|---:|---:|---|
| acetate/gfn2/gradient | 0.000e+00 | 3.155e-16 | — | pass |
| acetate/gfn2/hessian-analytic | 5.329e-15 | 2.522e-16 | 2.437e-12 | pass |
| acetate/gfnff/gradient | 4.441e-16 | 1.839e-16 | — | pass |
| acetate/gfnff/hessian-analytic | 8.882e-16 | 3.459e-16 | 2.354e-11 | pass |
| acetate/gxtb/gradient | 0.000e+00 | 3.920e-16 | — | pass |
| acetate/gxtb/hessian-analytic | 0.000e+00 | 1.266e-15 | 1.511e-09 | pass |
| acetylene/gfn2/gradient | 0.000e+00 | 2.671e-15 | — | pass |
| acetylene/gfnff/gradient | 2.220e-16 | 2.776e-17 | — | pass |
| acetylene/gxtb/gradient | 0.000e+00 | 6.523e-16 | — | pass |
| anthracene2/gfn2/gradient | 0.000e+00 | 8.092e-16 | — | pass |
| anthracene2/gfnff/gradient | 0.000e+00 | 6.184e-16 | — | pass |
| anthracene2/gxtb/gradient | 0.000e+00 | 1.865e-15 | — | pass |
| aspirin/gfn2/gradient | 0.000e+00 | 5.464e-16 | — | pass |
| aspirin/gfn2/hessian-analytic | 2.132e-14 | 2.588e-15 | 5.002e-12 | pass |
| aspirin/gfnff/gradient | 1.776e-15 | 4.250e-16 | — | pass |
| aspirin/gfnff/hessian-analytic | 0.000e+00 | 5.996e-16 | 1.806e-11 | pass |
| aspirin/gxtb/gradient | 0.000e+00 | 1.153e-15 | — | pass |
| aspirin/gxtb/hessian-analytic | 1.137e-13 | 1.031e-15 | 2.071e-09 | pass |
| beni-acetone/gfn2/gradient | 7.105e-15 | 6.713e-16 | — | pass |
| beni-acetone/gfnff/gradient | 4.441e-16 | 1.639e-16 | — | pass |
| beni-acetone/gxtb/gradient | 1.137e-13 | 1.138e-15 | — | pass |
| bromobenzene/gfn2/gradient | 3.553e-15 | 3.834e-16 | — | pass |
| bromobenzene/gfnff/gradient | 1.332e-15 | 1.847e-16 | — | pass |
| bromobenzene/gxtb/gradient | 0.000e+00 | 1.047e-15 | — | pass |
| caffeine/gfn2/gradient | 0.000e+00 | 8.603e-16 | — | pass |
| caffeine/gfn2/hessian-analytic | 0.000e+00 | 8.162e-16 | 1.312e-11 | pass |
| caffeine/gfnff/gradient | 1.776e-15 | 5.135e-16 | — | pass |
| caffeine/gfnff/hessian-analytic | 8.882e-16 | 3.570e-16 | 1.017e-11 | pass |
| caffeine/gxtb/gradient | 0.000e+00 | 1.339e-15 | — | pass |
| caffeine/gxtb/hessian-analytic | 0.000e+00 | 1.561e-15 | 6.048e-09 | pass |
| cf3br-benzene/gfn2/gradient | 0.000e+00 | 3.868e-16 | — | pass |
| cf3br-benzene/gfnff/gradient | 1.776e-15 | 9.541e-17 | — | pass |
| cf3br-benzene/gxtb/gradient | 0.000e+00 | 8.734e-16 | — | pass |
| cisplatin/gfn2/gradient | 3.553e-15 | 8.674e-16 | — | pass |
| cisplatin/gfnff/gradient | 6.661e-16 | 1.006e-16 | — | pass |
| cisplatin/gxtb/gradient | 0.000e+00 | 8.535e-16 | — | pass |
| co2/gfn2/gradient | 0.000e+00 | 1.665e-15 | — | pass |
| co2/gfn2/hessian-analytic | 8.882e-15 | 2.998e-15 | 9.095e-13 | pass |
| co2/gfnff/gradient | 1.110e-16 | 1.665e-16 | — | pass |
| co2/gfnff/hessian-analytic | 1.110e-16 | 3.887e-16 | 1.364e-12 | pass |
| co2/gxtb/gradient | 0.000e+00 | 1.867e-15 | — | pass |
| co2/gxtb/hessian-analytic | 0.000e+00 | 2.467e-15 | 2.015e-10 | pass |
| cobalt-hexammine/gfn2/gradient | 3.553e-15 | 7.065e-16 | — | pass |
| cobalt-hexammine/gfnff/gradient | 2.442e-15 | 1.376e-15 | — | pass |
| cobalt-hexammine/gxtb/gradient | 0.000e+00 | 1.341e-15 | — | pass |
| cu-nh3-4/gfn2/gradient | 3.553e-15 | 6.332e-16 | — | pass |
| cu-nh3-4/gfnff/gradient | 8.882e-16 | 1.190e-15 | — | pass |
| cu-nh3-4/gxtb/gradient | 0.000e+00 | 9.506e-16 | — | pass |
| ferrocene/gfn2/gradient | 7.105e-15 | 1.301e-14 | — | pass |
| ferrocene/gfnff/gradient | 0.000e+00 | 1.717e-16 | — | pass |
| ferrocene/gxtb/gradient | 0.000e+00 | 1.230e-15 | — | pass |
| gap-p-benzyne/gfn2/gradient | 5.329e-15 | 2.934e-09 | — | pass |
| gap-p-benzyne/gfnff/gradient | 2.220e-16 | 2.082e-16 | — | pass |
| gap-p-benzyne/gxtb/gradient | 0.000e+00 | 7.987e-15 | — | pass |
| gh106/gfn2/gradient | 1.264e-04 | 9.440e-04 | — | fail |
| gh106/gfnff/gradient | 0.000e+00 | 1.084e-19 | — | pass |
| gh106/gxtb/gradient | 1.839e-03 | 1.842e-02 | — | fail |
| guanidinium/gfn2/gradient | 3.553e-15 | 2.165e-15 | — | pass |
| guanidinium/gfnff/gradient | 2.220e-16 | 6.193e-16 | — | pass |
| guanidinium/gxtb/gradient | 0.000e+00 | 1.044e-15 | — | pass |
| hcn/gfn2/gradient | 4.441e-15 | 1.066e-16 | — | pass |
| hcn/gfn2/hessian-analytic | 3.553e-15 | 1.615e-16 | 1.364e-12 | pass |
| hcn/gfnff/gradient | 1.110e-16 | 2.776e-17 | — | pass |
| hcn/gfnff/hessian-analytic | 0.000e+00 | 0.000e+00 | 1.364e-12 | pass |
| hcn/gxtb/gradient | 0.000e+00 | 2.207e-15 | — | pass |
| hcn/gxtb/hessian-analytic | 0.000e+00 | 2.015e-15 | 2.101e-10 | pass |
| ibuprofen/gfn2/gradient | 0.000e+00 | 7.676e-16 | — | pass |
| ibuprofen/gfnff/gradient | 1.776e-15 | 5.287e-16 | — | pass |
| ibuprofen/gxtb/gradient | 0.000e+00 | 1.560e-15 | — | pass |
| iodomethane/gfn2/gradient | 2.665e-15 | 1.353e-15 | — | pass |
| iodomethane/gfnff/gradient | 0.000e+00 | 1.388e-17 | — | pass |
| iodomethane/gxtb/gradient | 0.000e+00 | 5.829e-16 | — | pass |
| iron-hexaaqua/gfn2/gradient | 2.132e-14 | 3.116e-15 | — | pass |
| iron-hexaaqua/gfnff/gradient | 0.000e+00 | 2.300e-15 | — | pass |
| iron-hexaaqua/gxtb/gradient | 0.000e+00 | 4.675e-16 | — | pass |
| l-alanine-zwitterion/gfn2/gradient | 3.553e-15 | 1.795e-15 | — | pass |
| l-alanine-zwitterion/gfnff/gradient | 4.441e-16 | 4.163e-16 | — | pass |
| l-alanine-zwitterion/gxtb/gradient | 0.000e+00 | 4.406e-16 | — | pass |
| mb33/gfn2/gradient | 0.000e+00 | 1.481e-15 | — | pass |
| mb33/gfnff/gradient | 6.661e-16 | 1.509e-16 | — | pass |
| mb33/gxtb/gradient | 0.000e+00 | 9.263e-16 | — | pass |
| methyl-radical/gfn2/gradient | 1.776e-15 | 8.335e-17 | — | pass |
| methyl-radical/gfn2/hessian-analytic | 8.882e-16 | 7.102e-17 | 6.764e-12 | pass |
| methyl-radical/gfnff/gradient | 1.110e-16 | 2.114e-17 | — | pass |
| methyl-radical/gfnff/hessian-analytic | 2.220e-16 | 1.408e-17 | 2.274e-12 | pass |
| methyl-radical/gxtb/gradient | 0.000e+00 | 3.757e-15 | — | pass |
| methyl-radical/gxtb/hessian-analytic | 0.000e+00 | 1.011e-15 | 2.615e-10 | pass |
| nickel-tetracarbonyl/gfn2/gradient | 0.000e+00 | 2.663e-15 | — | pass |
| nickel-tetracarbonyl/gfn2/hessian-analytic | 7.105e-15 | 1.043e-14 | 2.519e-10 | pass |
| nickel-tetracarbonyl/gfnff/gradient | 6.661e-16 | 4.944e-16 | — | pass |
| nickel-tetracarbonyl/gfnff/hessian-analytic | 2.220e-16 | 1.385e-16 | 2.871e-12 | pass |
| nickel-tetracarbonyl/gxtb/gradient | 0.000e+00 | 9.298e-16 | — | pass |
| nickel-tetracarbonyl/gxtb/hessian-analytic | 0.000e+00 | 1.501e-15 | 1.026e-09 | pass |
| no2/gfn2/gradient | 1.776e-15 | 9.327e-15 | — | pass |
| no2/gfnff/gradient | 0.000e+00 | 0.000e+00 | — | pass |
| no2/gxtb/gradient | 0.000e+00 | 3.198e-15 | — | pass |
| o2-triplet/gfn2/gradient | 2.665e-15 | 3.261e-16 | — | pass |
| o2-triplet/gfnff/gradient | 0.000e+00 | 0.000e+00 | — | pass |
| o2-triplet/gxtb/gradient | 0.000e+00 | 3.608e-16 | — | pass |
| paracetamol/gfn2/gradient | 1.421e-14 | 5.504e-16 | — | pass |
| paracetamol/gfn2/hessian-analytic | 0.000e+00 | 1.788e-15 | 5.912e-12 | pass |
| paracetamol/gfnff/gradient | 4.441e-16 | 5.378e-16 | — | pass |
| paracetamol/gfnff/hessian-analytic | 8.882e-16 | 3.194e-16 | 1.462e-11 | pass |
| paracetamol/gxtb/gradient | 0.000e+00 | 1.554e-15 | — | pass |
| paracetamol/gxtb/hessian-analytic | 0.000e+00 | 7.201e-16 | 3.937e-09 | pass |
| s229-cl/gfn2/gradient | 3.553e-15 | 5.343e-16 | — | pass |
| s229-cl/gfnff/gradient | 6.661e-15 | 6.412e-15 | — | pass |
| s229-cl/gxtb/gradient | 0.000e+00 | 2.576e-15 | — | pass |
| sf6/gfn2/gradient | 3.553e-15 | 3.916e-16 | — | pass |
| sf6/gfnff/gradient | 4.441e-16 | 8.622e-16 | — | pass |
| sf6/gxtb/gradient | 0.000e+00 | 1.477e-15 | — | pass |
| sulfate/gfn2/gradient | 3.553e-15 | 9.593e-16 | — | pass |
| sulfate/gfnff/gradient | 4.441e-16 | 1.388e-16 | — | pass |
| sulfate/gxtb/gradient | 0.000e+00 | 1.120e-15 | — | pass |
| taxol/gfn2/gradient | 0.000e+00 | 1.455e-15 | — | pass |
| taxol/gfnff/gradient | 2.132e-14 | 7.286e-15 | — | pass |
| taxol/gxtb/gradient | 0.000e+00 | 2.759e-15 | — | pass |
| tempo/gfn2/gradient | 7.105e-15 | 2.512e-15 | — | pass |
| tempo/gfnff/gradient | 2.665e-15 | 1.332e-15 | — | pass |
| tempo/gxtb/gradient | 0.000e+00 | 7.754e-16 | — | pass |
| tetramethylammonium/gfn2/gradient | 7.105e-15 | 4.558e-16 | — | pass |
| tetramethylammonium/gfnff/gradient | 4.441e-16 | 1.023e-16 | — | pass |
| tetramethylammonium/gxtb/gradient | 0.000e+00 | 7.216e-16 | — | pass |
| tm-MnO4-/gfn2/gradient | 0.000e+00 | 1.090e-13 | — | pass |
| tm-MnO4-/gfnff/gradient | 4.441e-16 | 4.814e-17 | — | pass |
| tm-MnO4-/gxtb/gradient | 0.000e+00 | 8.923e-15 | — | pass |
| twisted-ethylene/gfn2/gradient | 3.689e-12 | 6.951e-08 | — | pass |
| twisted-ethylene/gfnff/gradient | 1.110e-16 | 5.898e-17 | — | pass |
| twisted-ethylene/gxtb/gradient | 4.724e-09 | 1.115e-01 | — | fail |
| tyr-xab/gfn2/gradient | 7.105e-15 | 1.731e-15 | — | pass |
| tyr-xab/gfnff/gradient | 4.441e-15 | 7.772e-16 | — | pass |
| tyr-xab/gxtb/gradient | 0.000e+00 | 1.197e-15 | — | pass |
| uf6/gfn2/gradient | — | — | — | expected error |
| uf6/gfnff/gradient | 4.441e-16 | 6.377e-16 | — | pass |
| uf6/gxtb/gradient | 0.000e+00 | 2.773e-15 | — | pass |
| water/gfn2/gradient | 8.882e-16 | 4.441e-16 | — | pass |
| water/gfnff/gradient | 0.000e+00 | 1.735e-18 | — | pass |
| water/gxtb/gradient | 0.000e+00 | 2.984e-16 | — | pass |
| water-dimer/gfn2/gradient | 1.776e-15 | 2.533e-16 | — | pass |
| water-dimer/gfn2/hessian-analytic | 1.776e-15 | 2.020e-16 | 2.248e-11 | pass |
| water-dimer/gfnff/gradient | 4.441e-16 | 8.240e-17 | — | pass |
| water-dimer/gfnff/hessian-analytic | 1.110e-16 | 3.955e-16 | 1.600e-11 | pass |
| water-dimer/gxtb/gradient | 0.000e+00 | 6.488e-16 | — | pass |
| water-dimer/gxtb/hessian-analytic | 0.000e+00 | 3.970e-16 | 2.663e-09 | pass |
| zncl4/gfn2/gradient | 0.000e+00 | 1.776e-16 | — | pass |
| zncl4/gfnff/gradient | 0.000e+00 | 1.735e-18 | — | pass |
| zncl4/gxtb/gradient | 0.000e+00 | 3.945e-16 | — | pass |
