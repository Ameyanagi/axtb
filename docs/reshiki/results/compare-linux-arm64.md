# axtb platform comparison — FAIL

BASE (Linux x64 reference): `linux-x64.json` — Linux-7.0.14-arch1-1-x86_64-with-glibc2.43 / x86_64 / AMD Ryzen 9 7940HS w/ Radeon 780M Graphics
OTHER: `linux-arm64.json` — Linux-6.17.0-1014-nvidia-aarch64-with-glibc2.39 / aarch64 / aarch64

Limits are 1e-8 Eh for |dE|, 1e-7 Eh/bohr for max|dG|, and 0.1 cm⁻¹ for max sorted |dν|. Charges and iteration counts are measured without gates. Input fingerprints must match, selected runs must be present, and at least one numeric pair is required. Only declared expected error pairs are accepted.

| Method / run | Runs | Numeric pairs | Matching errors | Max ΔE (Eh) | Max ΔG (Eh/bohr) | Max Δν (cm⁻¹) | Max Δq (e, ungated) | SCF iterations BASE / OTHER | Iteration changes | Failed runs |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| gfn2 / gradient | 40 | 39 | 1 | 1.264e-04 | 9.436e-04 | — | 8.362e-04 | 570 / 575 | 1 | 1 |
| gfn2 / hessian-analytic | 9 | 9 | 0 | 1.421e-14 | 1.665e-14 | 1.055e-10 | 8.080e-13 | 18 / 18 | 0 | 0 |
| gfnff / gradient | 40 | 40 | 0 | 3.997e-15 | 8.143e-15 | — | 7.006e-14 | — | 0 | 0 |
| gfnff / hessian-analytic | 9 | 9 | 0 | 3.553e-15 | 1.435e-15 | 6.662e-11 | 1.979e-14 | — | 0 | 0 |
| gxtb / gradient | 40 | 40 | 0 | 7.442e-06 | 1.115e-01 | — | 1.369e-03 | 557 / 552 | 2 | 2 |
| gxtb / hessian-analytic | 9 | 9 | 0 | 1.137e-13 | 4.018e-15 | 6.882e-09 | 4.419e-14 | 89 / 89 | 0 | 0 |
| **Total** | 147 | 146 | 1 | 1.264e-04 | 1.115e-01 | 6.882e-09 | 1.369e-03 | 1234 / 1234 | 3 | 3 |

Violations **3**. Changed diagnostics for expected error pairs: 0 (ungated).

**Violations**

- gh106/gfn2/gradient — energy delta 1.264e-04 exceeds 1e-08, gradient delta 9.436e-04 exceeds 1e-07
- gh106/gxtb/gradient — energy delta 7.442e-06 exceeds 1e-08, gradient delta 5.203e-04 exceeds 1e-07
- twisted-ethylene/gxtb/gradient — gradient delta 1.115e-01 exceeds 1e-07

**Iteration changes (ungated)**

- gh106/gfn2/gradient: 195 -> 200
- gh106/gxtb/gradient: 177 -> 167
- twisted-ethylene/gxtb/gradient: 9 -> 14

**Every run**

| Case / method / run | ΔE (Eh) | max ΔG (Eh/bohr) | max Δν (cm⁻¹) | Status |
|---|---:|---:|---:|---|
| acetate/gfn2/gradient | 0.000e+00 | 1.226e-15 | — | pass |
| acetate/gfn2/hessian-analytic | 3.553e-15 | 4.237e-16 | 1.793e-11 | pass |
| acetate/gfnff/gradient | 0.000e+00 | 0.000e+00 | — | pass |
| acetate/gfnff/hessian-analytic | 0.000e+00 | 0.000e+00 | 2.897e-11 | pass |
| acetate/gxtb/gradient | 2.842e-14 | 5.516e-16 | — | pass |
| acetate/gxtb/hessian-analytic | 0.000e+00 | 1.184e-15 | 4.740e-09 | pass |
| acetylene/gfn2/gradient | 0.000e+00 | 2.671e-15 | — | pass |
| acetylene/gfnff/gradient | 0.000e+00 | 0.000e+00 | — | pass |
| acetylene/gxtb/gradient | 0.000e+00 | 1.080e-15 | — | pass |
| anthracene2/gfn2/gradient | 2.842e-14 | 9.211e-16 | — | pass |
| anthracene2/gfnff/gradient | 1.776e-15 | 7.225e-16 | — | pass |
| anthracene2/gxtb/gradient | 0.000e+00 | 1.415e-15 | — | pass |
| aspirin/gfn2/gradient | 7.105e-15 | 6.557e-16 | — | pass |
| aspirin/gfn2/hessian-analytic | 1.421e-14 | 1.967e-15 | 6.636e-12 | pass |
| aspirin/gfnff/gradient | 1.776e-15 | 4.111e-16 | — | pass |
| aspirin/gfnff/hessian-analytic | 0.000e+00 | 1.435e-15 | 3.250e-11 | pass |
| aspirin/gxtb/gradient | 0.000e+00 | 2.051e-15 | — | pass |
| aspirin/gxtb/hessian-analytic | 1.137e-13 | 9.800e-16 | 1.165e-09 | pass |
| beni-acetone/gfn2/gradient | 7.105e-15 | 6.124e-16 | — | pass |
| beni-acetone/gfnff/gradient | 0.000e+00 | 1.292e-16 | — | pass |
| beni-acetone/gxtb/gradient | 1.137e-13 | 1.363e-15 | — | pass |
| bromobenzene/gfn2/gradient | 0.000e+00 | 5.690e-16 | — | pass |
| bromobenzene/gfnff/gradient | 0.000e+00 | 1.735e-16 | — | pass |
| bromobenzene/gxtb/gradient | 0.000e+00 | 1.685e-15 | — | pass |
| caffeine/gfn2/gradient | 1.421e-14 | 1.043e-15 | — | pass |
| caffeine/gfn2/hessian-analytic | 0.000e+00 | 8.303e-16 | 3.056e-11 | pass |
| caffeine/gfnff/gradient | 2.665e-15 | 6.800e-16 | — | pass |
| caffeine/gfnff/hessian-analytic | 3.553e-15 | 8.542e-16 | 6.662e-11 | pass |
| caffeine/gxtb/gradient | 0.000e+00 | 1.346e-15 | — | pass |
| caffeine/gxtb/hessian-analytic | 0.000e+00 | 1.981e-15 | 6.882e-09 | pass |
| cf3br-benzene/gfn2/gradient | 7.105e-15 | 4.580e-16 | — | pass |
| cf3br-benzene/gfnff/gradient | 1.776e-15 | 5.117e-17 | — | pass |
| cf3br-benzene/gxtb/gradient | 0.000e+00 | 8.528e-16 | — | pass |
| cisplatin/gfn2/gradient | 3.553e-15 | 1.516e-15 | — | pass |
| cisplatin/gfnff/gradient | 2.220e-16 | 4.510e-17 | — | pass |
| cisplatin/gxtb/gradient | 0.000e+00 | 7.806e-16 | — | pass |
| co2/gfn2/gradient | 0.000e+00 | 1.665e-15 | — | pass |
| co2/gfn2/hessian-analytic | 8.882e-15 | 2.998e-15 | 9.095e-13 | pass |
| co2/gfnff/gradient | 0.000e+00 | 0.000e+00 | — | pass |
| co2/gfnff/hessian-analytic | 0.000e+00 | 0.000e+00 | 9.095e-13 | pass |
| co2/gxtb/gradient | 0.000e+00 | 5.135e-16 | — | pass |
| co2/gxtb/hessian-analytic | 2.842e-14 | 4.018e-15 | 3.165e-10 | pass |
| cobalt-hexammine/gfn2/gradient | 1.421e-14 | 6.965e-16 | — | pass |
| cobalt-hexammine/gfnff/gradient | 6.661e-16 | 1.020e-15 | — | pass |
| cobalt-hexammine/gxtb/gradient | 0.000e+00 | 1.107e-15 | — | pass |
| cu-nh3-4/gfn2/gradient | 3.553e-15 | 5.753e-16 | — | pass |
| cu-nh3-4/gfnff/gradient | 8.882e-16 | 1.405e-15 | — | pass |
| cu-nh3-4/gxtb/gradient | 0.000e+00 | 7.173e-16 | — | pass |
| ferrocene/gfn2/gradient | 7.105e-15 | 6.386e-15 | — | pass |
| ferrocene/gfnff/gradient | 0.000e+00 | 1.613e-16 | — | pass |
| ferrocene/gxtb/gradient | 0.000e+00 | 1.233e-15 | — | pass |
| gap-p-benzyne/gfn2/gradient | 0.000e+00 | 1.629e-10 | — | pass |
| gap-p-benzyne/gfnff/gradient | 0.000e+00 | 1.804e-16 | — | pass |
| gap-p-benzyne/gxtb/gradient | 0.000e+00 | 1.413e-14 | — | pass |
| gh106/gfn2/gradient | 1.264e-04 | 9.436e-04 | — | fail |
| gh106/gfnff/gradient | 0.000e+00 | 0.000e+00 | — | pass |
| gh106/gxtb/gradient | 7.442e-06 | 5.203e-04 | — | fail |
| guanidinium/gfn2/gradient | 1.776e-15 | 6.835e-16 | — | pass |
| guanidinium/gfnff/gradient | 2.220e-16 | 4.120e-16 | — | pass |
| guanidinium/gxtb/gradient | 2.842e-14 | 7.954e-16 | — | pass |
| hcn/gfn2/gradient | 4.441e-15 | 1.066e-16 | — | pass |
| hcn/gfn2/hessian-analytic | 3.553e-15 | 1.615e-16 | 1.364e-12 | pass |
| hcn/gfnff/gradient | 0.000e+00 | 0.000e+00 | — | pass |
| hcn/gfnff/hessian-analytic | 0.000e+00 | 0.000e+00 | 1.364e-12 | pass |
| hcn/gxtb/gradient | 0.000e+00 | 1.586e-15 | — | pass |
| hcn/gxtb/hessian-analytic | 0.000e+00 | 1.036e-15 | 1.573e-10 | pass |
| ibuprofen/gfn2/gradient | 7.105e-15 | 6.414e-16 | — | pass |
| ibuprofen/gfnff/gradient | 8.882e-16 | 1.928e-15 | — | pass |
| ibuprofen/gxtb/gradient | 0.000e+00 | 1.510e-15 | — | pass |
| iodomethane/gfn2/gradient | 8.882e-16 | 3.851e-16 | — | pass |
| iodomethane/gfnff/gradient | 0.000e+00 | 0.000e+00 | — | pass |
| iodomethane/gxtb/gradient | 0.000e+00 | 1.069e-15 | — | pass |
| iron-hexaaqua/gfn2/gradient | 2.132e-14 | 1.735e-15 | — | pass |
| iron-hexaaqua/gfnff/gradient | 8.882e-16 | 1.492e-15 | — | pass |
| iron-hexaaqua/gxtb/gradient | 0.000e+00 | 3.834e-16 | — | pass |
| l-alanine-zwitterion/gfn2/gradient | 0.000e+00 | 1.119e-15 | — | pass |
| l-alanine-zwitterion/gfnff/gradient | 0.000e+00 | 1.180e-16 | — | pass |
| l-alanine-zwitterion/gxtb/gradient | 0.000e+00 | 1.080e-15 | — | pass |
| mb33/gfn2/gradient | 0.000e+00 | 7.078e-16 | — | pass |
| mb33/gfnff/gradient | 8.882e-16 | 1.388e-16 | — | pass |
| mb33/gxtb/gradient | 0.000e+00 | 6.627e-16 | — | pass |
| methyl-radical/gfn2/gradient | 8.882e-16 | 1.032e-16 | — | pass |
| methyl-radical/gfn2/hessian-analytic | 8.882e-16 | 7.102e-17 | 6.764e-12 | pass |
| methyl-radical/gfnff/gradient | 0.000e+00 | 0.000e+00 | — | pass |
| methyl-radical/gfnff/hessian-analytic | 0.000e+00 | 0.000e+00 | 5.201e-12 | pass |
| methyl-radical/gxtb/gradient | 0.000e+00 | 1.053e-15 | — | pass |
| methyl-radical/gxtb/hessian-analytic | 0.000e+00 | 1.619e-15 | 2.519e-10 | pass |
| nickel-tetracarbonyl/gfn2/gradient | 3.553e-15 | 6.608e-15 | — | pass |
| nickel-tetracarbonyl/gfn2/hessian-analytic | 3.553e-15 | 1.665e-14 | 1.055e-10 | pass |
| nickel-tetracarbonyl/gfnff/gradient | 0.000e+00 | 0.000e+00 | — | pass |
| nickel-tetracarbonyl/gfnff/hessian-analytic | 0.000e+00 | 0.000e+00 | 4.576e-12 | pass |
| nickel-tetracarbonyl/gxtb/gradient | 0.000e+00 | 1.561e-15 | — | pass |
| nickel-tetracarbonyl/gxtb/hessian-analytic | 0.000e+00 | 1.737e-15 | 1.597e-09 | pass |
| no2/gfn2/gradient | 1.776e-15 | 9.327e-15 | — | pass |
| no2/gfnff/gradient | 0.000e+00 | 0.000e+00 | — | pass |
| no2/gxtb/gradient | 0.000e+00 | 3.517e-16 | — | pass |
| o2-triplet/gfn2/gradient | 2.665e-15 | 3.261e-16 | — | pass |
| o2-triplet/gfnff/gradient | 0.000e+00 | 0.000e+00 | — | pass |
| o2-triplet/gxtb/gradient | 0.000e+00 | 3.608e-16 | — | pass |
| paracetamol/gfn2/gradient | 1.421e-14 | 9.784e-16 | — | pass |
| paracetamol/gfn2/hessian-analytic | 7.105e-15 | 1.829e-15 | 7.276e-12 | pass |
| paracetamol/gfnff/gradient | 4.441e-16 | 6.180e-16 | — | pass |
| paracetamol/gfnff/hessian-analytic | 1.776e-15 | 7.456e-16 | 1.467e-11 | pass |
| paracetamol/gxtb/gradient | 0.000e+00 | 1.280e-15 | — | pass |
| paracetamol/gxtb/hessian-analytic | 0.000e+00 | 1.057e-15 | 2.889e-09 | pass |
| s229-cl/gfn2/gradient | 3.553e-15 | 2.591e-15 | — | pass |
| s229-cl/gfnff/gradient | 3.997e-15 | 4.496e-15 | — | pass |
| s229-cl/gxtb/gradient | 0.000e+00 | 1.707e-15 | — | pass |
| sf6/gfn2/gradient | 3.553e-15 | 7.818e-16 | — | pass |
| sf6/gfnff/gradient | 0.000e+00 | 0.000e+00 | — | pass |
| sf6/gxtb/gradient | 0.000e+00 | 1.056e-15 | — | pass |
| sulfate/gfn2/gradient | 3.553e-15 | 9.593e-16 | — | pass |
| sulfate/gfnff/gradient | 0.000e+00 | 0.000e+00 | — | pass |
| sulfate/gxtb/gradient | 0.000e+00 | 1.003e-15 | — | pass |
| taxol/gfn2/gradient | 2.842e-14 | 2.137e-15 | — | pass |
| taxol/gfnff/gradient | 3.553e-15 | 8.143e-15 | — | pass |
| taxol/gxtb/gradient | 0.000e+00 | 2.054e-15 | — | pass |
| tempo/gfn2/gradient | 7.105e-15 | 8.882e-16 | — | pass |
| tempo/gfnff/gradient | 2.665e-15 | 1.790e-15 | — | pass |
| tempo/gxtb/gradient | 0.000e+00 | 6.167e-16 | — | pass |
| tetramethylammonium/gfn2/gradient | 1.066e-14 | 4.341e-16 | — | pass |
| tetramethylammonium/gfnff/gradient | 6.661e-16 | 1.353e-16 | — | pass |
| tetramethylammonium/gxtb/gradient | 0.000e+00 | 6.314e-16 | — | pass |
| tm-MnO4-/gfn2/gradient | 0.000e+00 | 1.357e-13 | — | pass |
| tm-MnO4-/gfnff/gradient | 0.000e+00 | 0.000e+00 | — | pass |
| tm-MnO4-/gxtb/gradient | 0.000e+00 | 7.803e-15 | — | pass |
| twisted-ethylene/gfn2/gradient | 8.251e-13 | 1.145e-08 | — | pass |
| twisted-ethylene/gfnff/gradient | 0.000e+00 | 0.000e+00 | — | pass |
| twisted-ethylene/gxtb/gradient | 1.046e-09 | 1.115e-01 | — | fail |
| tyr-xab/gfn2/gradient | 1.421e-14 | 1.908e-15 | — | pass |
| tyr-xab/gfnff/gradient | 8.882e-16 | 7.563e-16 | — | pass |
| tyr-xab/gxtb/gradient | 0.000e+00 | 1.397e-15 | — | pass |
| uf6/gfn2/gradient | — | — | — | expected error |
| uf6/gfnff/gradient | 0.000e+00 | 0.000e+00 | — | pass |
| uf6/gxtb/gradient | 0.000e+00 | 1.606e-15 | — | pass |
| water/gfn2/gradient | 8.882e-16 | 4.441e-16 | — | pass |
| water/gfnff/gradient | 0.000e+00 | 0.000e+00 | — | pass |
| water/gxtb/gradient | 0.000e+00 | 5.829e-16 | — | pass |
| water-dimer/gfn2/gradient | 0.000e+00 | 9.975e-17 | — | pass |
| water-dimer/gfn2/hessian-analytic | 1.776e-15 | 2.741e-16 | 4.206e-12 | pass |
| water-dimer/gfnff/gradient | 0.000e+00 | 0.000e+00 | — | pass |
| water-dimer/gfnff/hessian-analytic | 0.000e+00 | 0.000e+00 | 3.058e-11 | pass |
| water-dimer/gxtb/gradient | 0.000e+00 | 2.359e-16 | — | pass |
| water-dimer/gxtb/hessian-analytic | 0.000e+00 | 6.757e-16 | 5.752e-09 | pass |
| zncl4/gfn2/gradient | 3.553e-15 | 4.267e-16 | — | pass |
| zncl4/gfnff/gradient | 0.000e+00 | 0.000e+00 | — | pass |
| zncl4/gxtb/gradient | 0.000e+00 | 2.983e-16 | — | pass |
