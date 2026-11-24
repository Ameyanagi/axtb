# ReShiki reference validation

Run date (UTC) `2026-11-24T20:40:45.434341+00:00`
Host `Linux-7.0.14-arch1-1-x86_64-with-glibc2.43` / `x86_64` / `AMD Ryzen 9 7940HS w/ Radeon 780M Graphics`

axtb checkout revision `449483a8fa7979e3568ef25d34863c5e9f3a8f99`

| Program | Version | Binary SHA256 |
|---|---|---|
| axtb | axtb 0.2.0 | `c6e121205c8aafb132a5f2c9769d8e17aa65ed6c372f20b29c885adc6aab0fb1` |
| xtb | 6.8.0 (4b8c987066c7e9bd0c2b727ce026fb9a581c273a) | `95d82a623580bb1ef7e1f9aec6a91d3773b79e5aabf7bb8fada287debf65c93b` |
| gxtb | 6.7.1 (26dd68d) | `1b4e30b68ed4e88b4075f60d97f4756ee440fe92d3294cc20826d53f8121cd26` |
| tblite | 0.8.0 | `0a97d28d3635c7dfd143de6d00b09c0cb382249c020b28ee33fb4002db1efb9c` |

Settings

- **Selection** methods gfn2, gfnff, gxtb, case substring None, jobs 8
- **Threads** one per calculation, axtb SCF eigensolver faer (AXTB_LAPACK and AXTB_LAPACK_BLAS unset)
- **xtb and g-xTB environment** XTBPATH unset and XTBHOME set to each run's scratch directory
- **GFN2** values axtb and tblite 1e-6, xtb 1e-4. Molden both 1e-4. Analytic Hessian 1e-6, xtb 0.001.
- **g-xTB** values both 1e-4 and 250 iterations maximum. Molden default accuracy. Reference Hessian 1e-4.
- **GFN-FF** xtb 6.8.0 conventions in axtb. Reference gradient and Hessian use xtb defaults.
- **Long molecules** cutoff energy and gradient differences compared at their published precision
- **Hessian geometries** method-specific optimized XYZ when present, otherwise original XYZ
- **GFN2 analytic frequencies** matrix comparison uses modes with reference |ν| > 50 cm⁻¹
- **Units** energy and orbital energies Eh, gradient Eh/bohr, charges e, dipole e·bohr, frequency cm⁻¹

## Checks

| Method | Check | Cases | Pass | Fail | Other | Worst deviation by quantity | Tolerance by quantity |
|---|---|---:|---:|---:|---:|---|---|
| gfn2 | xtb_values | 40 | 39 | 0 | 1 | energy 7.38964e-13, gradient 1.02703e-09, charges 9.29496e-09, dipole 1.00977e-07 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 |
| gfn2 | tblite_values | 39 | 38 | 0 | 1 | energy 4.9738e-14, gradient 1.95788e-12, charges 8.25253e-11, dipole 4.25087e-11 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 |
| gfn2 | default_vs_xtb | 39 | 37 | 0 | 2 | energy 9.359e-06, gradient 3.14946e-06, charges 3.30401e-09, dipole 9.17242e-08 | energy 1e-07 / 9.4e-06, gradient 1e-06 / 3.1e-06, charges 1e-05, dipole 0.0001 |
| gfn2 | molden | 39 | 39 | 0 | 0 | orbital_energies 1.56958e-10 | orbital_energies 7e-10 |
| gfn2 | hessian_analytic | 9 | 9 | 0 | 0 | frequencies 0.176754 | frequencies 6.6 |
| gfn2 | hessian_same | 9 | 9 | 0 | 0 | frequencies 5.00486e-05 | frequencies 5e-05 |
| gfnff | xtb_values | 40 | 40 | 0 | 0 | energy 4.99156e-13, gradient 1.17961e-14, charges 7.37188e-14 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 |
| gfnff | molden | 40 | 40 | 0 | 0 | — | — |
| gfnff | hessian_analytic | 9 | 6 | 3 | 0 | frequencies 0.438898 | frequencies 0.27 |
| gfnff | hessian_same | 9 | 9 | 0 | 0 | frequencies 0.00499417 | frequencies 0.005 |
| gfnff | hessian_numerical | 9 | 9 | 0 | 0 | frequencies 0.00499408 | frequencies 0.005 |
| gxtb | gxtb_values | 40 | 36 | 3 | 1 | energy 2.72848e-12, gradient 4.7488e-09, charges 8.93069e-09 | energy 2.1e-12, gradient 6e-11, charges not gated |
| gxtb | molden | 40 | 36 | 4 | 0 | orbital_energies 0.0045183 | orbital_energies 3.3e-12 |
| gxtb | hessian_analytic | 9 | 9 | 0 | 0 | frequencies 0.348004 | frequencies 0.35 |
| gxtb | hessian_same | 9 | 8 | 1 | 0 | frequencies 5.05589e-05 | frequencies 5e-05 |

## Non-pass cases

| Case | Method | Check | Status | Deviation | Details |
|---|---|---|---|---|---|
| acetate | gfnff | hessian_analytic | fail | frequencies 0.41753 | Tolerance exceeded for frequencies |
| caffeine | gfnff | hessian_analytic | fail | frequencies 0.323836 | Tolerance exceeded for frequencies |
| caffeine | gxtb | hessian_same | fail | frequencies 5.05589e-05 | Tolerance exceeded for frequencies |
| cu-nh3-4 | gxtb | molden | fail | orbital_energies 3.61e-12 | Tolerance exceeded for orbital_energies |
| gh106 | gfn2 | tblite_values | both_failed | — | Both SCFs failed to converge (tests/validate_tblite.py KNOWN_DIFFERENCES["gh106"]: tblite 0.8.0 does not converge) / axtb scf_not_converged: SCF not converged in 250 iterations (last energy change -1.074e-3 Eh, residual 1.261e-2) / reference ReferenceRefusal: TBLiteRuntimeError: SCF not converged in 250 cycles |
| gh106 | gfn2 | default_vs_xtb | expected_difference | energy 2.63812e-06, gradient 3.53987e-07, charges 3.30401e-09, dipole 9.17242e-08 | Documented real-space cutoff difference at the source's two-significant-figure precision |
| gh106 | gxtb | gxtb_values | both_failed | — | Both SCFs failed to converge (tests/structures/README.md: gh106 is a pathological SCF case) / axtb scf_not_converged: SCF not converged in 250 iterations (last energy change -6.151e-5 Eh, residual 2.197e-4) / reference ReferenceRefusal: xtb exited 128: 2  -1.3262781E-05   1.6872702E-05     208     -998.8272447028   2.8197478E-05   1.5451514E-05     209     -998.8272467835  -2.0807542E-06   1.5408501E-05     210     -998.8272663678  -1.9584256E-05   1.5279172E-05     211 |
| gh106 | gxtb | molden | fail | orbital_energies 0.0045183 | Tolerance exceeded for orbital_energies |
| taxol | gfn2 | default_vs_xtb | expected_difference | energy 9.359e-06, gradient 3.14946e-06, charges 1.02391e-10, dipole 2.03466e-10 | Documented real-space cutoff difference at the source's two-significant-figure precision |
| taxol | gxtb | gxtb_values | fail | energy 2.72848e-12, gradient 2.6576e-15, charges 4.98381e-09 | Tolerance exceeded for energy |
| tm-MnO4- | gxtb | gxtb_values | fail | energy 2.27374e-13, gradient 1.85479e-10, charges 4.75014e-09 | Tolerance exceeded for gradient |
| tm-MnO4- | gxtb | molden | fail | orbital_energies 3.9299e-11 | Tolerance exceeded for orbital_energies |
| twisted-ethylene | gxtb | gxtb_values | fail | energy 4.40536e-13, gradient 4.7488e-09, charges 8.93069e-09 | Tolerance exceeded for gradient |
| twisted-ethylene | gxtb | molden | fail | orbital_energies 5.19169e-07 | Tolerance exceeded for orbital_energies |
| uf6 | gfn2 | xtb_values | expected_error | — | Both programs refused uranium (Z=92) / axtb unsupported_element: GFN2-xTB has no parameters for element Z=92 / reference ReferenceRefusal: xtb exited 128: 63/1.4959605        for FOD:    * S. Grimme and A. Hansen, Angew. Chem. Int. Ed., 2015, 54, 12308-12313.      DOI: 10.1002/anie.201501887    * C. A. Bauer, A. Hansen, S. Grimme, Chem. Eur. J., 2017, 23, 6150-6164.      DOI: |
| water-dimer | gfnff | hessian_analytic | fail | frequencies 0.438898 | Tolerance exceeded for frequencies |

## Every check

Deviation is max |axtb − reference| per quantity; tolerance is the documented limit (— not gated).

| Case | Method | Check | Reference | Deviation | Tolerance | Iterations axtb / ref | Status |
|---|---|---|---|---|---|---|---|
| acetate | gfn2 | default_vs_xtb | axtb[xtb] | energy 4.27072e-11, gradient 2.98262e-11, charges 8.90288e-12, dipole 1.0117e-11 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 17 / 17 | pass |
| acetate | gfn2 | hessian_analytic | xtb | frequencies 0.0276418 | frequencies 6.6 | 2 / 14 | pass |
| acetate | gfn2 | hessian_same | xtb | frequencies 4.49108e-05 | frequencies 5e-05 | — / 14 | pass |
| acetate | gfn2 | molden | xtb | orbital_energies 8.99281e-15 | orbital_energies 7e-10 | 15 / 15 | pass |
| acetate | gfn2 | tblite_values | tblite | energy 1.77636e-15, gradient 1.35324e-12, charges 2.0804e-11, dipole 1.90066e-11 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 17 / — | pass |
| acetate | gfn2 | xtb_values | xtb | energy 1.63425e-13, gradient 4.08296e-11, charges 4.93387e-09, dipole 2.26716e-09 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 17 / 15 | pass |
| acetate | gfnff | hessian_analytic | xtb | frequencies 0.41753 | frequencies 0.27 | — / — | fail |
| acetate | gfnff | hessian_numerical | xtb | frequencies 0.00499408 | frequencies 0.005 | — / — | pass |
| acetate | gfnff | hessian_same | xtb | frequencies 0.00499417 | frequencies 0.005 | — / — | pass |
| acetate | gfnff | molden | — | — | — | — / — | pass |
| acetate | gfnff | xtb_values | xtb | energy 4.30322e-13, gradient 5.13478e-16, charges 1.22125e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| acetate | gxtb | gxtb_values | gxtb | energy 1.98952e-13, gradient 1.97932e-15, charges 3.61629e-09 | energy 2.1e-12, gradient 6e-11, charges — | 14 / 14 | pass |
| acetate | gxtb | hessian_analytic | gxtb | frequencies 0.348004 | frequencies 0.35 | 9 / 13 | pass |
| acetate | gxtb | hessian_same | gxtb | frequencies 4.96552e-05 | frequencies 5e-05 | — / 13 | pass |
| acetate | gxtb | molden | gxtb | orbital_energies 9.99201e-15 | orbital_energies 3.3e-12 | 9 / 9 | pass |
| acetylene | gfn2 | default_vs_xtb | axtb[xtb] | energy 1.52118e-11, gradient 1.67185e-12, charges 7.99527e-12, dipole 1.50121e-11 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 11 / 11 | pass |
| acetylene | gfn2 | molden | xtb | orbital_energies 9.99201e-15 | orbital_energies 7e-10 | 10 / 10 | pass |
| acetylene | gfn2 | tblite_values | tblite | energy 0, gradient 1.20407e-13, charges 2.17676e-11, dipole 4.25087e-11 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 10 / — | pass |
| acetylene | gfn2 | xtb_values | xtb | energy 9.76996e-15, gradient 6.85398e-12, charges 3.88425e-09, dipole 7.24913e-11 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 11 / 10 | pass |
| acetylene | gfnff | molden | — | — | — | — / — | pass |
| acetylene | gfnff | xtb_values | xtb | energy 3.70259e-13, gradient 4.51028e-16, charges 1.04777e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| acetylene | gxtb | gxtb_values | gxtb | energy 7.10543e-14, gradient 1.05332e-14, charges 7.77359e-11 | energy 2.1e-12, gradient 6e-11, charges — | 10 / 10 | pass |
| acetylene | gxtb | molden | gxtb | orbital_energies 2.10942e-15 | orbital_energies 3.3e-12 | 7 / 7 | pass |
| anthracene2 | gfn2 | default_vs_xtb | axtb[xtb] | energy 3.65446e-10, gradient 9.2868e-12, charges 3.50728e-11, dipole 4.93111e-12 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 16 / 16 | pass |
| anthracene2 | gfn2 | molden | xtb | orbital_energies 9.99201e-15 | orbital_energies 7e-10 | 14 / 14 | pass |
| anthracene2 | gfn2 | tblite_values | tblite | energy 4.26326e-14, gradient 1.6851e-14, charges 6.66689e-13, dipole 2.18044e-12 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 16 / — | pass |
| anthracene2 | gfn2 | xtb_values | xtb | energy 9.9476e-14, gradient 6.77628e-11, charges 5.65214e-09, dipole 6.12776e-12 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 16 / 14 | pass |
| anthracene2 | gfnff | molden | — | — | — | — / — | pass |
| anthracene2 | gfnff | xtb_values | xtb | energy 1.66978e-13, gradient 6.08888e-16, charges 2.52975e-14 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| anthracene2 | gxtb | gxtb_values | gxtb | energy 2.27374e-13, gradient 1.35655e-15, charges 4.93581e-09 | energy 2.1e-12, gradient 6e-11, charges — | 14 / 14 | pass |
| anthracene2 | gxtb | molden | gxtb | orbital_energies 6.99441e-15 | orbital_energies 3.3e-12 | 8 / 8 | pass |
| aspirin | gfn2 | default_vs_xtb | axtb[xtb] | energy 1.56149e-10, gradient 3.17277e-11, charges 7.40952e-11, dipole 1.26033e-10 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 24 / 24 | pass |
| aspirin | gfn2 | hessian_analytic | xtb | frequencies 0.124089 | frequencies 6.6 | 2 / 18 | pass |
| aspirin | gfn2 | hessian_same | xtb | frequencies 5.00486e-05 | frequencies 5e-05 | — / 18 | pass |
| aspirin | gfn2 | molden | xtb | orbital_energies 6.99441e-15 | orbital_energies 7e-10 | 20 / 20 | pass |
| aspirin | gfn2 | tblite_values | tblite | energy 2.84217e-14, gradient 2.17968e-15, charges 1.78746e-14, dipole 4.04121e-14 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 25 / — | pass |
| aspirin | gfn2 | xtb_values | xtb | energy 2.41585e-13, gradient 1.02703e-09, charges 6.05333e-09, dipole 9.14224e-09 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 24 / 20 | pass |
| aspirin | gfnff | hessian_analytic | xtb | frequencies 0.255261 | frequencies 0.27 | — / — | pass |
| aspirin | gfnff | hessian_numerical | xtb | frequencies 0.00493679 | frequencies 0.005 | — / — | pass |
| aspirin | gfnff | hessian_same | xtb | frequencies 0.00493684 | frequencies 0.005 | — / — | pass |
| aspirin | gfnff | molden | — | — | — | — / — | pass |
| aspirin | gfnff | xtb_values | xtb | energy 2.4869e-13, gradient 1.23512e-15, charges 6.32827e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| aspirin | gxtb | gxtb_values | gxtb | energy 2.27374e-13, gradient 3.31506e-15, charges 4.9602e-09 | energy 2.1e-12, gradient 6e-11, charges — | 17 / 17 | pass |
| aspirin | gxtb | hessian_analytic | gxtb | frequencies 0.174144 | frequencies 0.35 | 9 / 17 | pass |
| aspirin | gxtb | hessian_same | gxtb | frequencies 4.97489e-05 | frequencies 5e-05 | — / 17 | pass |
| aspirin | gxtb | molden | gxtb | orbital_energies 9.99201e-15 | orbital_energies 3.3e-12 | 9 / 9 | pass |
| beni-acetone | gfn2 | default_vs_xtb | axtb[xtb] | energy 1.18781e-10, gradient 2.70709e-11, charges 5.96305e-11, dipole 9.54267e-11 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 26 / 26 | pass |
| beni-acetone | gfn2 | molden | xtb | orbital_energies 6.99441e-15 | orbital_energies 7e-10 | 22 / 22 | pass |
| beni-acetone | gfn2 | tblite_values | tblite | energy 1.42109e-14, gradient 1.06946e-15, charges 1.82077e-14, dipole 1.7597e-13 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 26 / — | pass |
| beni-acetone | gfn2 | xtb_values | xtb | energy 3.05533e-13, gradient 4.75782e-10, charges 5.01462e-09, dipole 9.37917e-09 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 26 / 22 | pass |
| beni-acetone | gfnff | molden | — | — | — | — / — | pass |
| beni-acetone | gfnff | xtb_values | xtb | energy 4.74731e-13, gradient 3.01842e-16, charges 2.95597e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| beni-acetone | gxtb | gxtb_values | gxtb | energy 4.54747e-13, gradient 1.72778e-15, charges 4.61836e-09 | energy 2.1e-12, gradient 6e-11, charges — | 15 / 15 | pass |
| beni-acetone | gxtb | molden | gxtb | orbital_energies 9.99201e-15 | orbital_energies 3.3e-12 | 9 / 9 | pass |
| bromobenzene | gfn2 | default_vs_xtb | axtb[xtb] | energy 7.4408e-11, gradient 1.53938e-11, charges 5.27602e-11, dipole 2.47788e-10 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 24 / 24 | pass |
| bromobenzene | gfn2 | molden | xtb | orbital_energies 2.9976e-15 | orbital_energies 7e-10 | 20 / 20 | pass |
| bromobenzene | gfn2 | tblite_values | tblite | energy 0, gradient 5.88071e-16, charges 2.4647e-14, dipole 5.06262e-14 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 24 / — | pass |
| bromobenzene | gfn2 | xtb_values | xtb | energy 3.69482e-13, gradient 9.34908e-11, charges 4.59214e-09, dipole 1.13375e-08 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 24 / 20 | pass |
| bromobenzene | gfnff | molden | — | — | — | — / — | pass |
| bromobenzene | gfnff | xtb_values | xtb | energy 2.70894e-14, gradient 1.50921e-16, charges 1.94289e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| bromobenzene | gxtb | gxtb_values | gxtb | energy 0, gradient 2.27075e-15, charges 4.89609e-09 | energy 2.1e-12, gradient 6e-11, charges — | 15 / 15 | pass |
| bromobenzene | gxtb | molden | gxtb | orbital_energies 2.9976e-15 | orbital_energies 3.3e-12 | 7 / 7 | pass |
| caffeine | gfn2 | default_vs_xtb | axtb[xtb] | energy 2.14577e-10, gradient 3.35124e-11, charges 5.86653e-11, dipole 3.58329e-11 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 22 / 22 | pass |
| caffeine | gfn2 | hessian_analytic | xtb | frequencies 0.176754 | frequencies 6.6 | 2 / 17 | pass |
| caffeine | gfn2 | hessian_same | xtb | frequencies 4.97259e-05 | frequencies 5e-05 | — / 17 | pass |
| caffeine | gfn2 | molden | xtb | orbital_energies 5.9952e-15 | orbital_energies 7e-10 | 19 / 19 | pass |
| caffeine | gfn2 | tblite_values | tblite | energy 2.84217e-14, gradient 2.19012e-15, charges 1.22125e-14, dipole 5.14033e-14 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 23 / — | pass |
| caffeine | gfn2 | xtb_values | xtb | energy 2.34479e-13, gradient 5.17978e-10, charges 6.26843e-09, dipole 2.73657e-09 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 22 / 19 | pass |
| caffeine | gfnff | hessian_analytic | xtb | frequencies 0.323836 | frequencies 0.27 | — / — | fail |
| caffeine | gfnff | hessian_numerical | xtb | frequencies 0.00497484 | frequencies 0.005 | — / — | pass |
| caffeine | gfnff | hessian_same | xtb | frequencies 0.00497498 | frequencies 0.005 | — / — | pass |
| caffeine | gfnff | molden | — | — | — | — / — | pass |
| caffeine | gfnff | xtb_values | xtb | energy 1.74083e-13, gradient 5.23886e-16, charges 6.55032e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| caffeine | gxtb | gxtb_values | gxtb | energy 3.41061e-13, gradient 1.00902e-14, charges 4.37517e-09 | energy 2.1e-12, gradient 6e-11, charges — | 17 / 17 | pass |
| caffeine | gxtb | hessian_analytic | gxtb | frequencies 0.302762 | frequencies 0.35 | 9 / 17 | pass |
| caffeine | gxtb | hessian_same | gxtb | frequencies 5.05589e-05 | frequencies 5e-05 | — / 17 | fail |
| caffeine | gxtb | molden | gxtb | orbital_energies 9.99201e-15 | orbital_energies 3.3e-12 | 9 / 9 | pass |
| cf3br-benzene | gfn2 | default_vs_xtb | axtb[xtb] | energy 1.06773e-10, gradient 3.67675e-11, charges 1.34623e-11, dipole 6.39147e-10 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 23 / 23 | pass |
| cf3br-benzene | gfn2 | molden | xtb | orbital_energies 1.09912e-14 | orbital_energies 7e-10 | 20 / 20 | pass |
| cf3br-benzene | gfn2 | tblite_values | tblite | energy 2.84217e-14, gradient 1.69656e-15, charges 1.58762e-14, dipole 1.76359e-13 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 24 / — | pass |
| cf3br-benzene | gfn2 | xtb_values | xtb | energy 3.76588e-13, gradient 4.92234e-10, charges 5.40515e-09, dipole 1.60877e-08 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 23 / 20 | pass |
| cf3br-benzene | gfnff | molden | — | — | — | — / — | pass |
| cf3br-benzene | gfnff | xtb_values | xtb | energy 2.42029e-13, gradient 4.33681e-16, charges 1.8735e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| cf3br-benzene | gxtb | gxtb_values | gxtb | energy 4.54747e-13, gradient 1.14405e-15, charges 4.94392e-09 | energy 2.1e-12, gradient 6e-11, charges — | 14 / 14 | pass |
| cf3br-benzene | gxtb | molden | gxtb | orbital_energies 1.9984e-15 | orbital_energies 3.3e-12 | 9 / 9 | pass |
| cisplatin | gfn2 | default_vs_xtb | axtb[xtb] | energy 5.48646e-11, gradient 7.99813e-11, charges 2.53753e-11, dipole 2.79844e-10 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 24 / 24 | pass |
| cisplatin | gfn2 | molden | xtb | orbital_energies 9.99201e-15 | orbital_energies 7e-10 | 23 / 23 | pass |
| cisplatin | gfn2 | tblite_values | tblite | energy 7.10543e-15, gradient 3.4469e-15, charges 2.39808e-14, dipole 1.22569e-13 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 25 / — | pass |
| cisplatin | gfn2 | xtb_values | xtb | energy 2.94875e-13, gradient 3.43847e-11, charges 4.52069e-09, dipole 3.8549e-09 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 24 / 23 | pass |
| cisplatin | gfnff | molden | — | — | — | — / — | pass |
| cisplatin | gfnff | xtb_values | xtb | energy 4.99156e-13, gradient 5.3256e-16, charges 2.9976e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| cisplatin | gxtb | gxtb_values | gxtb | energy 4.54747e-13, gradient 1.00267e-15, charges 4.92442e-09 | energy 2.1e-12, gradient 6e-11, charges — | 16 / 16 | pass |
| cisplatin | gxtb | molden | gxtb | orbital_energies 9.99201e-16 | orbital_energies 3.3e-12 | 9 / 9 | pass |
| co2 | gfn2 | default_vs_xtb | axtb[xtb] | energy 4.71161e-11, gradient 1.32575e-11, charges 8.90121e-12, dipole 5.83977e-14 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 12 / 12 | pass |
| co2 | gfn2 | hessian_analytic | xtb | frequencies 0.0638748 | frequencies 6.6 | 2 / 10 | pass |
| co2 | gfn2 | hessian_same | xtb | frequencies 3.73716e-05 | frequencies 5e-05 | — / 10 | pass |
| co2 | gfn2 | molden | xtb | orbital_energies 1.9984e-15 | orbital_energies 7e-10 | 10 / 10 | pass |
| co2 | gfn2 | tblite_values | tblite | energy 0, gradient 4.85723e-16, charges 2.66454e-15, dipole 1.12688e-14 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 13 / — | pass |
| co2 | gfn2 | xtb_values | xtb | energy 1.5099e-13, gradient 1.68736e-10, charges 4.69351e-09, dipole 2.68674e-14 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 12 / 10 | pass |
| co2 | gfnff | hessian_analytic | xtb | frequencies 0.0580212 | frequencies 0.27 | — / — | pass |
| co2 | gfnff | hessian_numerical | xtb | frequencies 0.00326724 | frequencies 0.005 | — / — | pass |
| co2 | gfnff | hessian_same | xtb | frequencies 0.00326723 | frequencies 0.005 | — / — | pass |
| co2 | gfnff | molden | — | — | — | — / — | pass |
| co2 | gfnff | xtb_values | xtb | energy 1.91291e-13, gradient 1.9869e-16, charges 1.05471e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| co2 | gxtb | gxtb_values | gxtb | energy 5.11591e-13, gradient 3.07797e-15, charges 4.20842e-09 | energy 2.1e-12, gradient 6e-11, charges — | 12 / 12 | pass |
| co2 | gxtb | hessian_analytic | gxtb | frequencies 0.0379133 | frequencies 0.35 | 8 / 12 | pass |
| co2 | gxtb | hessian_same | gxtb | frequencies 3.26523e-05 | frequencies 5e-05 | — / 12 | pass |
| co2 | gxtb | molden | gxtb | orbital_energies 2.9976e-15 | orbital_energies 3.3e-12 | 8 / 8 | pass |
| cobalt-hexammine | gfn2 | default_vs_xtb | axtb[xtb] | energy 2.36174e-10, gradient 2.84135e-11, charges 4.48969e-10, dipole 1.34937e-12 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 19 / 19 | pass |
| cobalt-hexammine | gfn2 | molden | xtb | orbital_energies 9.99201e-15 | orbital_energies 7e-10 | 15 / 15 | pass |
| cobalt-hexammine | gfn2 | tblite_values | tblite | energy 7.10543e-15, gradient 6.95631e-13, charges 8.25253e-11, dipole 1.30232e-11 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 19 / — | pass |
| cobalt-hexammine | gfn2 | xtb_values | xtb | energy 4.90274e-13, gradient 8.74223e-11, charges 4.43811e-09, dipole 4.6063e-09 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 19 / 15 | pass |
| cobalt-hexammine | gfnff | molden | — | — | — | — / — | pass |
| cobalt-hexammine | gfnff | xtb_values | xtb | energy 1.88738e-14, gradient 1.69136e-15, charges 2.81997e-14 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| cobalt-hexammine | gxtb | gxtb_values | gxtb | energy 4.54747e-13, gradient 1.22753e-15, charges 4.91044e-09 | energy 2.1e-12, gradient 6e-11, charges — | 15 / 15 | pass |
| cobalt-hexammine | gxtb | molden | gxtb | orbital_energies 1.9984e-15 | orbital_energies 3.3e-12 | 9 / 9 | pass |
| cu-nh3-4 | gfn2 | default_vs_xtb | axtb[xtb] | energy 7.33635e-11, gradient 7.4117e-12, charges 2.67306e-10, dipole 1.30451e-14 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 16 / 16 | pass |
| cu-nh3-4 | gfn2 | molden | xtb | orbital_energies 2.9976e-15 | orbital_energies 7e-10 | 14 / 14 | pass |
| cu-nh3-4 | gfn2 | tblite_values | tblite | energy 1.77636e-14, gradient 1.2218e-13, charges 1.37467e-11, dipole 3.32512e-14 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 18 / — | pass |
| cu-nh3-4 | gfn2 | xtb_values | xtb | energy 2.73559e-13, gradient 8.04975e-11, charges 2.97293e-09, dipole 2.76446e-14 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 16 / 14 | pass |
| cu-nh3-4 | gfnff | molden | — | — | — | — / — | pass |
| cu-nh3-4 | gfnff | xtb_values | xtb | energy 3.35287e-13, gradient 9.88792e-16, charges 1.26565e-14 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| cu-nh3-4 | gxtb | gxtb_values | gxtb | energy 4.54747e-13, gradient 3.78949e-12, charges 4.77982e-09 | energy 2.1e-12, gradient 6e-11, charges — | 15 / 15 | pass |
| cu-nh3-4 | gxtb | molden | gxtb | orbital_energies 3.61e-12 | orbital_energies 3.3e-12 | 9 / 9 | fail |
| ferrocene | gfn2 | default_vs_xtb | axtb[xtb] | energy 2.56364e-11, gradient 1.66759e-11, charges 1.22226e-09, dipole 6.88477e-14 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 18 / 18 | pass |
| ferrocene | gfn2 | molden | xtb | orbital_energies 6.10623e-15 | orbital_energies 7e-10 | 15 / 15 | pass |
| ferrocene | gfn2 | tblite_values | tblite | energy 1.42109e-14, gradient 1.06165e-15, charges 1.29896e-14, dipole 5.3933e-14 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 18 / — | pass |
| ferrocene | gfn2 | xtb_values | xtb | energy 1.52767e-13, gradient 6.16515e-11, charges 6.87011e-09, dipole 6.15197e-09 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 18 / 15 | pass |
| ferrocene | gfnff | molden | — | — | — | — / — | pass |
| ferrocene | gfnff | xtb_values | xtb | energy 4.6807e-13, gradient 5.34295e-16, charges 5.62397e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| ferrocene | gxtb | gxtb_values | gxtb | energy 1.59162e-12, gradient 1.47105e-15, charges 4.3063e-09 | energy 2.1e-12, gradient 6e-11, charges — | 14 / 14 | pass |
| ferrocene | gxtb | molden | gxtb | orbital_energies 9.99201e-15 | orbital_energies 3.3e-12 | 9 / 9 | pass |
| gap-p-benzyne | gfn2 | default_vs_xtb | axtb[xtb] | energy 6.09237e-11, gradient 9.57635e-12, charges 9.50873e-12, dipole 9.4918e-12 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 15 / 15 | pass |
| gap-p-benzyne | gfn2 | molden | xtb | orbital_energies 8.99281e-15 | orbital_energies 7e-10 | 13 / 13 | pass |
| gap-p-benzyne | gfn2 | tblite_values | tblite | energy 3.55271e-15, gradient 5.64132e-15, charges 1.75304e-13, dipole 4.31308e-13 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 15 / — | pass |
| gap-p-benzyne | gfn2 | xtb_values | xtb | energy 3.85469e-13, gradient 3.40754e-11, charges 4.77658e-09, dipole 1.4681e-11 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 15 / 13 | pass |
| gap-p-benzyne | gfnff | molden | — | — | — | — / — | pass |
| gap-p-benzyne | gfnff | xtb_values | xtb | energy 9.65894e-14, gradient 4.82253e-16, charges 1.70697e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| gap-p-benzyne | gxtb | gxtb_values | gxtb | energy 1.13687e-13, gradient 4.0952e-13, charges 3.60609e-09 | energy 2.1e-12, gradient 6e-11, charges — | 12 / 12 | pass |
| gap-p-benzyne | gxtb | molden | gxtb | orbital_energies 2.9976e-15 | orbital_energies 3.3e-12 | 7 / 7 | pass |
| gh106 | gfn2 | default_vs_xtb | axtb[xtb] | energy 2.63812e-06, gradient 3.53987e-07, charges 3.30401e-09, dipole 9.17242e-08 | energy 9.4e-06, gradient 3.1e-06, charges 1e-05, dipole 0.0001 | 240 / 200 | expected_difference |
| gh106 | gfn2 | molden | xtb | orbital_energies 1.56958e-10 | orbital_energies 7e-10 | 187 / 193 | pass |
| gh106 | gfn2 | tblite_values | tblite | — | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | — / — | both_failed |
| gh106 | gfn2 | xtb_values | xtb | energy 6.39488e-14, gradient 1.03856e-10, charges 9.29496e-09, dipole 1.00977e-07 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 200 / 193 | pass |
| gh106 | gfnff | molden | — | — | — | — / — | pass |
| gh106 | gfnff | xtb_values | xtb | energy 2.68174e-13, gradient 5.27356e-16, charges 1e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| gh106 | gxtb | gxtb_values | gxtb | — | energy 2.1e-12, gradient 6e-11, charges — | — / — | both_failed |
| gh106 | gxtb | molden | gxtb | orbital_energies 0.0045183 | orbital_energies 3.3e-12 | 140 / 96 | fail |
| guanidinium | gfn2 | default_vs_xtb | axtb[xtb] | energy 1.00636e-10, gradient 5.15413e-12, charges 2.28493e-11, dipole 1.51838e-12 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 17 / 17 | pass |
| guanidinium | gfn2 | molden | xtb | orbital_energies 6.99441e-15 | orbital_energies 7e-10 | 14 / 14 | pass |
| guanidinium | gfn2 | tblite_values | tblite | energy 1.42109e-14, gradient 3.97609e-13, charges 8.41749e-12, dipole 4.83374e-12 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 17 / — | pass |
| guanidinium | gfn2 | xtb_values | xtb | energy 1.38556e-13, gradient 8.49023e-11, charges 5.33594e-09, dipole 1.18638e-09 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 17 / 14 | pass |
| guanidinium | gfnff | molden | — | — | — | — / — | pass |
| guanidinium | gfnff | xtb_values | xtb | energy 7.08322e-14, gradient 4.38885e-16, charges 3.55271e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| guanidinium | gxtb | gxtb_values | gxtb | energy 8.52651e-14, gradient 1.55084e-15, charges 4.77976e-09 | energy 2.1e-12, gradient 6e-11, charges — | 14 / 14 | pass |
| guanidinium | gxtb | molden | gxtb | orbital_energies 9.99201e-15 | orbital_energies 3.3e-12 | 8 / 8 | pass |
| hcn | gfn2 | default_vs_xtb | axtb[xtb] | energy 2.23563e-11, gradient 1.096e-11, charges 2.35623e-12, dipole 1.34763e-11 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 15 / 15 | pass |
| hcn | gfn2 | hessian_analytic | xtb | frequencies 0.081048 | frequencies 6.6 | 2 / 12 | pass |
| hcn | gfn2 | hessian_same | xtb | frequencies 2.32401e-05 | frequencies 5e-05 | — / 12 | pass |
| hcn | gfn2 | molden | xtb | orbital_energies 4.996e-15 | orbital_energies 7e-10 | 13 / 13 | pass |
| hcn | gfn2 | tblite_values | tblite | energy 4.44089e-15, gradient 3.747e-16, charges 6.10623e-15, dipole 4.44089e-15 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 15 / — | pass |
| hcn | gfn2 | xtb_values | xtb | energy 2.45137e-13, gradient 1.80648e-11, charges 2.90617e-09, dipole 3.02196e-09 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 15 / 13 | pass |
| hcn | gfnff | hessian_analytic | xtb | frequencies 0.0805272 | frequencies 0.27 | — / — | pass |
| hcn | gfnff | hessian_numerical | xtb | frequencies 0.00304755 | frequencies 0.005 | — / — | pass |
| hcn | gfnff | hessian_same | xtb | frequencies 0.0030476 | frequencies 0.005 | — / — | pass |
| hcn | gfnff | molden | — | — | — | — / — | pass |
| hcn | gfnff | xtb_values | xtb | energy 4.27658e-13, gradient 4.44089e-16, charges 1.38778e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| hcn | gxtb | gxtb_values | gxtb | energy 3.41061e-13, gradient 1.02696e-15, charges 2.41044e-09 | energy 2.1e-12, gradient 6e-11, charges — | 13 / 13 | pass |
| hcn | gxtb | hessian_analytic | gxtb | frequencies 0.0726175 | frequencies 0.35 | 9 / 13 | pass |
| hcn | gxtb | hessian_same | gxtb | frequencies 3.94421e-05 | frequencies 5e-05 | — / 13 | pass |
| hcn | gxtb | molden | gxtb | orbital_energies 5.9952e-15 | orbital_energies 3.3e-12 | 9 / 9 | pass |
| ibuprofen | gfn2 | default_vs_xtb | axtb[xtb] | energy 2.0578e-10, gradient 2.6244e-11, charges 8.00432e-11, dipole 2.54571e-11 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 25 / 25 | pass |
| ibuprofen | gfn2 | molden | xtb | orbital_energies 7.99361e-15 | orbital_energies 7e-10 | 20 / 20 | pass |
| ibuprofen | gfn2 | tblite_values | tblite | energy 2.13163e-14, gradient 1.65579e-15, charges 1.4877e-14, dipole 9.67004e-14 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 24 / — | pass |
| ibuprofen | gfn2 | xtb_values | xtb | energy 2.13163e-14, gradient 4.0122e-10, charges 7.61131e-09, dipole 2.97904e-09 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 25 / 20 | pass |
| ibuprofen | gfnff | molden | — | — | — | — / — | pass |
| ibuprofen | gfnff | xtb_values | xtb | energy 1.15463e-14, gradient 1.80498e-15, charges 1.5099e-14 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| ibuprofen | gxtb | gxtb_values | gxtb | energy 3.41061e-13, gradient 2.38351e-15, charges 4.83952e-09 | energy 2.1e-12, gradient 6e-11, charges — | 15 / 15 | pass |
| ibuprofen | gxtb | molden | gxtb | orbital_energies 9.99201e-15 | orbital_energies 3.3e-12 | 8 / 8 | pass |
| iodomethane | gfn2 | default_vs_xtb | axtb[xtb] | energy 1.25793e-11, gradient 1.01507e-11, charges 2.12097e-11, dipole 3.68143e-10 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 17 / 17 | pass |
| iodomethane | gfn2 | molden | xtb | orbital_energies 9.99201e-16 | orbital_energies 7e-10 | 15 / 15 | pass |
| iodomethane | gfn2 | tblite_values | tblite | energy 4.44089e-15, gradient 3.00107e-16, charges 9.65894e-15, dipole 7.77156e-15 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 17 / — | pass |
| iodomethane | gfn2 | xtb_values | xtb | energy 4.45866e-13, gradient 1.36706e-10, charges 3.41107e-09, dipole 3.94459e-09 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 17 / 15 | pass |
| iodomethane | gfnff | molden | — | — | — | — / — | pass |
| iodomethane | gfnff | xtb_values | xtb | energy 2.78888e-13, gradient 4.23273e-16, charges 6.38378e-16 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| iodomethane | gxtb | gxtb_values | gxtb | energy 3.97904e-13, gradient 7.84095e-16, charges 4.25363e-09 | energy 2.1e-12, gradient 6e-11, charges — | 12 / 12 | pass |
| iodomethane | gxtb | molden | gxtb | orbital_energies 7.99361e-15 | orbital_energies 3.3e-12 | 8 / 8 | pass |
| iron-hexaaqua | gfn2 | default_vs_xtb | axtb[xtb] | energy 1.3193e-09, gradient 4.3982e-10, charges 1.70334e-09, dipole 2.45006e-09 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 21 / 21 | pass |
| iron-hexaaqua | gfn2 | molden | xtb | orbital_energies 1.79856e-14 | orbital_energies 7e-10 | 18 / 18 | pass |
| iron-hexaaqua | gfn2 | tblite_values | tblite | energy 2.13163e-14, gradient 1.49414e-12, charges 2.179e-11, dipole 6.79445e-12 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 23 / — | pass |
| iron-hexaaqua | gfn2 | xtb_values | xtb | energy 4.76064e-13, gradient 2.95301e-10, charges 4.58467e-09, dipole 4.28121e-09 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 21 / 18 | pass |
| iron-hexaaqua | gfnff | molden | — | — | — | — / — | pass |
| iron-hexaaqua | gfnff | xtb_values | xtb | energy 2.55018e-13, gradient 1.2941e-15, charges 7.10543e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| iron-hexaaqua | gxtb | gxtb_values | gxtb | energy 6.82121e-13, gradient 2.03631e-12, charges 4.53932e-09 | energy 2.1e-12, gradient 6e-11, charges — | 30 / 30 | pass |
| iron-hexaaqua | gxtb | molden | gxtb | orbital_energies 9.99201e-15 | orbital_energies 3.3e-12 | 9 / 9 | pass |
| l-alanine-zwitterion | gfn2 | default_vs_xtb | axtb[xtb] | energy 1.09829e-10, gradient 3.06072e-11, charges 8.44771e-11, dipole 9.43143e-11 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 20 / 20 | pass |
| l-alanine-zwitterion | gfn2 | molden | xtb | orbital_energies 6.99441e-15 | orbital_energies 7e-10 | 17 / 17 | pass |
| l-alanine-zwitterion | gfn2 | tblite_values | tblite | energy 1.06581e-14, gradient 1.08073e-15, charges 1.39888e-14, dipole 1.02141e-14 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 20 / — | pass |
| l-alanine-zwitterion | gfn2 | xtb_values | xtb | energy 1.3145e-13, gradient 2.22673e-10, charges 5.56015e-09, dipole 4.37969e-09 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 20 / 17 | pass |
| l-alanine-zwitterion | gfnff | molden | — | — | — | — / — | pass |
| l-alanine-zwitterion | gfnff | xtb_values | xtb | energy 2.66454e-14, gradient 5.22152e-16, charges 7.21645e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| l-alanine-zwitterion | gxtb | gxtb_values | gxtb | energy 1.7053e-13, gradient 9.47159e-16, charges 4.83698e-09 | energy 2.1e-12, gradient 6e-11, charges — | 15 / 15 | pass |
| l-alanine-zwitterion | gxtb | molden | gxtb | orbital_energies 5.9952e-15 | orbital_energies 3.3e-12 | 9 / 9 | pass |
| mb33 | gfn2 | default_vs_xtb | axtb[xtb] | energy 5.90958e-11, gradient 1.60794e-10, charges 1.24515e-10, dipole 5.00023e-10 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 29 / 29 | pass |
| mb33 | gfn2 | molden | xtb | orbital_energies 4.996e-15 | orbital_energies 7e-10 | 25 / 25 | pass |
| mb33 | gfn2 | tblite_values | tblite | energy 7.10543e-15, gradient 9.47159e-16, charges 9.65894e-15, dipole 4.01901e-14 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 29 / — | pass |
| mb33 | gfn2 | xtb_values | xtb | energy 3.55271e-13, gradient 2.64898e-10, charges 6.82847e-09, dipole 2.5055e-08 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 29 / 25 | pass |
| mb33 | gfnff | molden | — | — | — | — / — | pass |
| mb33 | gfnff | xtb_values | xtb | energy 2.82885e-13, gradient 5.16948e-16, charges 2.48412e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| mb33 | gxtb | gxtb_values | gxtb | energy 2.27374e-13, gradient 1.61286e-15, charges 4.80877e-09 | energy 2.1e-12, gradient 6e-11, charges — | 20 / 20 | pass |
| mb33 | gxtb | molden | gxtb | orbital_energies 9.99201e-15 | orbital_energies 3.3e-12 | 10 / 10 | pass |
| methyl-radical | gfn2 | default_vs_xtb | axtb[xtb] | energy 1.02176e-11, gradient 1.29423e-12, charges 6.40155e-13, dipole 2.33147e-15 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 10 / 10 | pass |
| methyl-radical | gfn2 | hessian_analytic | xtb | frequencies 0.0601081 | frequencies 6.6 | 2 / 9 | pass |
| methyl-radical | gfn2 | hessian_same | xtb | frequencies 4.47182e-05 | frequencies 5e-05 | — / 9 | pass |
| methyl-radical | gfn2 | molden | xtb | orbital_energies 9.99201e-16 | orbital_energies 7e-10 | 9 / 9 | pass |
| methyl-radical | gfn2 | tblite_values | tblite | energy 8.88178e-16, gradient 2.08844e-16, charges 1.55431e-15, dipole 1.81799e-15 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 9 / — | pass |
| methyl-radical | gfn2 | xtb_values | xtb | energy 1.62537e-13, gradient 1.23896e-11, charges 4.47399e-09, dipole 1.5245e-12 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 10 / 9 | pass |
| methyl-radical | gfnff | hessian_analytic | xtb | frequencies 0.131012 | frequencies 0.27 | — / — | pass |
| methyl-radical | gfnff | hessian_numerical | xtb | frequencies 0.00446319 | frequencies 0.005 | — / — | pass |
| methyl-radical | gfnff | hessian_same | xtb | frequencies 0.00446301 | frequencies 0.005 | — / — | pass |
| methyl-radical | gfnff | molden | — | — | — | — / — | pass |
| methyl-radical | gfnff | xtb_values | xtb | energy 4.40092e-13, gradient 5.16948e-16, charges 1.00614e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| methyl-radical | gxtb | gxtb_values | gxtb | energy 3.76588e-13, gradient 4.52478e-15, charges 2.00113e-09 | energy 2.1e-12, gradient 6e-11, charges — | 12 / 12 | pass |
| methyl-radical | gxtb | hessian_analytic | gxtb | frequencies 0.0592252 | frequencies 0.35 | 9 / 12 | pass |
| methyl-radical | gxtb | hessian_same | gxtb | frequencies 4.8115e-05 | frequencies 5e-05 | — / 12 | pass |
| methyl-radical | gxtb | molden | gxtb | orbital_energies 4.79616e-14 | orbital_energies 3.3e-12 | 9 / 9 | pass |
| nickel-tetracarbonyl | gfn2 | default_vs_xtb | axtb[xtb] | energy 9.47509e-12, gradient 9.02384e-12, charges 2.00088e-10, dipole 4.74065e-14 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 18 / 18 | pass |
| nickel-tetracarbonyl | gfn2 | hessian_analytic | xtb | frequencies 0.0776858 | frequencies 6.6 | 2 / 14 | pass |
| nickel-tetracarbonyl | gfn2 | hessian_same | xtb | frequencies 4.1115e-05 | frequencies 5e-05 | — / 14 | pass |
| nickel-tetracarbonyl | gfn2 | molden | xtb | orbital_energies 3.39728e-14 | orbital_energies 7e-10 | 15 / 15 | pass |
| nickel-tetracarbonyl | gfn2 | tblite_values | tblite | energy 4.9738e-14, gradient 3.67316e-14, charges 1.94023e-12, dipole 4.40203e-13 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 16 / — | pass |
| nickel-tetracarbonyl | gfn2 | xtb_values | xtb | energy 2.70006e-13, gradient 6.60736e-11, charges 1.36782e-09, dipole 2.02061e-14 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 18 / 15 | pass |
| nickel-tetracarbonyl | gfnff | hessian_analytic | xtb | frequencies 0.0304733 | frequencies 0.27 | — / — | pass |
| nickel-tetracarbonyl | gfnff | hessian_numerical | xtb | frequencies 0.00413706 | frequencies 0.005 | — / — | pass |
| nickel-tetracarbonyl | gfnff | hessian_same | xtb | frequencies 0.00413709 | frequencies 0.005 | — / — | pass |
| nickel-tetracarbonyl | gfnff | molden | — | — | — | — / — | pass |
| nickel-tetracarbonyl | gfnff | xtb_values | xtb | energy 4.08784e-13, gradient 1.48319e-16, charges 1.9984e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| nickel-tetracarbonyl | gxtb | gxtb_values | gxtb | energy 4.54747e-13, gradient 6.42686e-12, charges 4.67566e-09 | energy 2.1e-12, gradient 6e-11, charges — | 14 / 14 | pass |
| nickel-tetracarbonyl | gxtb | hessian_analytic | gxtb | frequencies 0.0825005 | frequencies 0.35 | 9 / 14 | pass |
| nickel-tetracarbonyl | gxtb | hessian_same | gxtb | frequencies 4.62145e-05 | frequencies 5e-05 | — / 14 | pass |
| nickel-tetracarbonyl | gxtb | molden | gxtb | orbital_energies 9.2204e-13 | orbital_energies 3.3e-12 | 9 / 9 | pass |
| no2 | gfn2 | default_vs_xtb | axtb[xtb] | energy 5.48273e-11, gradient 1.9015e-11, charges 5.93636e-12, dipole 3.23749e-11 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 14 / 14 | pass |
| no2 | gfn2 | molden | xtb | orbital_energies 2.9976e-15 | orbital_energies 7e-10 | 12 / 12 | pass |
| no2 | gfn2 | tblite_values | tblite | energy 1.06581e-14, gradient 1.91877e-14, charges 6.19504e-14, dipole 2.27485e-13 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 14 / — | pass |
| no2 | gfn2 | xtb_values | xtb | energy 1.75859e-13, gradient 2.03805e-10, charges 4.48798e-09, dipole 7.7999e-09 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 14 / 12 | pass |
| no2 | gfnff | molden | — | — | — | — / — | pass |
| no2 | gfnff | xtb_values | xtb | energy 3.1658e-13, gradient 4.02456e-16, charges 1.22125e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| no2 | gxtb | gxtb_values | gxtb | energy 3.12639e-13, gradient 4.12864e-16, charges 4.1296e-09 | energy 2.1e-12, gradient 6e-11, charges — | 15 / 15 | pass |
| no2 | gxtb | molden | gxtb | orbital_energies 9.99201e-15 | orbital_energies 3.3e-12 | 10 / 10 | pass |
| o2-triplet | gfn2 | default_vs_xtb | axtb[xtb] | energy 2.9007e-11, gradient 2.16422e-11, charges 1.13243e-14, dipole 1.83464e-14 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 8 / 8 | pass |
| o2-triplet | gfn2 | molden | xtb | orbital_energies 9.99201e-16 | orbital_energies 7e-10 | 7 / 7 | pass |
| o2-triplet | gfn2 | tblite_values | tblite | energy 4.44089e-15, gradient 5.0307e-16, charges 9.10383e-15, dipole 1.35447e-14 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 8 / — | pass |
| o2-triplet | gfn2 | xtb_values | xtb | energy 4.88498e-14, gradient 5.9566e-12, charges 3.77476e-15, dipole 5.44009e-15 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 8 / 7 | pass |
| o2-triplet | gfnff | molden | — | — | — | — / — | pass |
| o2-triplet | gfnff | xtb_values | xtb | energy 3.41949e-14, gradient 3.91354e-15, charges 0 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| o2-triplet | gxtb | gxtb_values | gxtb | energy 1.42109e-13, gradient 8.71733e-14, charges 4.21885e-15 | energy 2.1e-12, gradient 6e-11, charges — | 11 / 11 | pass |
| o2-triplet | gxtb | molden | gxtb | orbital_energies 9.99201e-16 | orbital_energies 3.3e-12 | 7 / 7 | pass |
| paracetamol | gfn2 | default_vs_xtb | axtb[xtb] | energy 1.53946e-10, gradient 1.89541e-11, charges 4.85918e-11, dipole 3.24836e-11 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 23 / 23 | pass |
| paracetamol | gfn2 | hessian_analytic | xtb | frequencies 0.158652 | frequencies 6.6 | 2 / 18 | pass |
| paracetamol | gfn2 | hessian_same | xtb | frequencies 4.73414e-05 | frequencies 5e-05 | — / 18 | pass |
| paracetamol | gfn2 | molden | xtb | orbital_energies 7.99361e-15 | orbital_energies 7e-10 | 20 / 20 | pass |
| paracetamol | gfn2 | tblite_values | tblite | energy 1.42109e-14, gradient 3.01495e-15, charges 1.82077e-14, dipole 7.43849e-14 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 24 / — | pass |
| paracetamol | gfn2 | xtb_values | xtb | energy 4.26326e-13, gradient 1.41786e-10, charges 5.00334e-09, dipole 3.88485e-09 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 23 / 20 | pass |
| paracetamol | gfnff | hessian_analytic | xtb | frequencies 0.151992 | frequencies 0.27 | — / — | pass |
| paracetamol | gfnff | hessian_numerical | xtb | frequencies 0.00492316 | frequencies 0.005 | — / — | pass |
| paracetamol | gfnff | hessian_same | xtb | frequencies 0.00492304 | frequencies 0.005 | — / — | pass |
| paracetamol | gfnff | molden | — | — | — | — / — | pass |
| paracetamol | gfnff | xtb_values | xtb | energy 2.70006e-13, gradient 5.10876e-16, charges 4.32987e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| paracetamol | gxtb | gxtb_values | gxtb | energy 5.68434e-13, gradient 1.29584e-15, charges 4.71088e-09 | energy 2.1e-12, gradient 6e-11, charges — | 16 / 16 | pass |
| paracetamol | gxtb | hessian_analytic | gxtb | frequencies 0.215186 | frequencies 0.35 | 9 / 16 | pass |
| paracetamol | gxtb | hessian_same | gxtb | frequencies 4.9772e-05 | frequencies 5e-05 | — / 16 | pass |
| paracetamol | gxtb | molden | gxtb | orbital_energies 9.99201e-15 | orbital_energies 3.3e-12 | 9 / 9 | pass |
| s229-cl | gfn2 | default_vs_xtb | axtb[xtb] | energy 8.71978e-11, gradient 1.94829e-10, charges 2.0649e-10, dipole 9.93952e-11 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 21 / 21 | pass |
| s229-cl | gfn2 | molden | xtb | orbital_energies 4.996e-15 | orbital_energies 7e-10 | 18 / 18 | pass |
| s229-cl | gfn2 | tblite_values | tblite | energy 0, gradient 2.34882e-15, charges 1.04361e-14, dipole 4.84057e-14 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 23 / — | pass |
| s229-cl | gfn2 | xtb_values | xtb | energy 1.49214e-13, gradient 4.32371e-10, charges 6.81202e-09, dipole 4.23889e-09 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 21 / 18 | pass |
| s229-cl | gfnff | molden | — | — | — | — / — | pass |
| s229-cl | gfnff | xtb_values | xtb | energy 2.33147e-13, gradient 1.17961e-14, charges 7.37188e-14 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| s229-cl | gxtb | gxtb_values | gxtb | energy 2.27374e-13, gradient 6.48787e-15, charges 4.93083e-09 | energy 2.1e-12, gradient 6e-11, charges — | 16 / 16 | pass |
| s229-cl | gxtb | molden | gxtb | orbital_energies 9.99201e-15 | orbital_energies 3.3e-12 | 9 / 9 | pass |
| sf6 | gfn2 | default_vs_xtb | axtb[xtb] | energy 2.57039e-11, gradient 1.55046e-11, charges 2.64677e-11, dipole 8.42957e-15 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 14 / 14 | pass |
| sf6 | gfn2 | molden | xtb | orbital_energies 9.99201e-15 | orbital_energies 7e-10 | 13 / 13 | pass |
| sf6 | gfn2 | tblite_values | tblite | energy 1.06581e-14, gradient 1.12334e-15, charges 7.99361e-15, dipole 1.3517e-14 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 14 / — | pass |
| sf6 | gfn2 | xtb_values | xtb | energy 5.08038e-13, gradient 1.59038e-11, charges 4.22181e-09, dipole 6.43929e-15 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 14 / 13 | pass |
| sf6 | gfnff | molden | — | — | — | — / — | pass |
| sf6 | gfnff | xtb_values | xtb | energy 1.55542e-13, gradient 3.61239e-16, charges 2.38698e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| sf6 | gxtb | gxtb_values | gxtb | energy 3.41061e-13, gradient 4.44804e-12, charges 3.97461e-09 | energy 2.1e-12, gradient 6e-11, charges — | 13 / 13 | pass |
| sf6 | gxtb | molden | gxtb | orbital_energies 1.91003e-12 | orbital_energies 3.3e-12 | 8 / 8 | pass |
| sulfate | gfn2 | default_vs_xtb | axtb[xtb] | energy 4.79012e-11, gradient 1.93686e-12, charges 1.56342e-10, dipole 1.46549e-14 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 13 / 13 | pass |
| sulfate | gfn2 | molden | xtb | orbital_energies 2.9976e-15 | orbital_energies 7e-10 | 11 / 11 | pass |
| sulfate | gfn2 | tblite_values | tblite | energy 0, gradient 1.28201e-13, charges 4.87455e-12, dipole 1.09912e-14 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 14 / — | pass |
| sulfate | gfn2 | xtb_values | xtb | energy 1.98952e-13, gradient 2.14524e-11, charges 3.09963e-09, dipole 1.22125e-14 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 13 / 11 | pass |
| sulfate | gfnff | molden | — | — | — | — / — | pass |
| sulfate | gfnff | xtb_values | xtb | energy 4.2899e-13, gradient 2.25514e-16, charges 3.33067e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| sulfate | gxtb | gxtb_values | gxtb | energy 6.82121e-13, gradient 6.01268e-12, charges 4.30996e-09 | energy 2.1e-12, gradient 6e-11, charges — | 13 / 13 | pass |
| sulfate | gxtb | molden | gxtb | orbital_energies 1.95606e-12 | orbital_energies 3.3e-12 | 9 / 9 | pass |
| taxol | gfn2 | default_vs_xtb | axtb[xtb] | energy 9.359e-06, gradient 3.14946e-06, charges 1.02391e-10, dipole 2.03466e-10 | energy 9.4e-06, gradient 3.1e-06, charges 1e-05, dipole 0.0001 | 27 / 27 | expected_difference |
| taxol | gfn2 | molden | xtb | orbital_energies 5.9952e-15 | orbital_energies 7e-10 | 22 / 22 | pass |
| taxol | gfn2 | tblite_values | tblite | energy 0, gradient 2.78944e-15, charges 1.11022e-14, dipole 2.00284e-13 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 27 / — | pass |
| taxol | gfn2 | xtb_values | xtb | energy 7.38964e-13, gradient 5.53656e-10, charges 8.09545e-09, dipole 3.0605e-08 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 27 / 22 | pass |
| taxol | gfnff | molden | — | — | — | — / — | pass |
| taxol | gfnff | xtb_values | xtb | energy 4.1922e-13, gradient 3.90833e-15, charges 5.02931e-14 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| taxol | gxtb | gxtb_values | gxtb | energy 2.72848e-12, gradient 2.6576e-15, charges 4.98381e-09 | energy 2.1e-12, gradient 6e-11, charges — | 17 / 17 | fail |
| taxol | gxtb | molden | gxtb | orbital_energies 2.02061e-14 | orbital_energies 3.3e-12 | 9 / 9 | pass |
| tempo | gfn2 | default_vs_xtb | axtb[xtb] | energy 1.67063e-10, gradient 7.00065e-11, charges 1.88431e-10, dipole 7.93761e-11 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 22 / 22 | pass |
| tempo | gfn2 | molden | xtb | orbital_energies 4.996e-15 | orbital_energies 7e-10 | 18 / 18 | pass |
| tempo | gfn2 | tblite_values | tblite | energy 2.84217e-14, gradient 5.82867e-16, charges 8.99281e-15, dipole 4.45616e-14 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 22 / — | pass |
| tempo | gfn2 | xtb_values | xtb | energy 1.35003e-13, gradient 3.41265e-10, charges 6.12726e-09, dipole 1.04026e-08 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 22 / 18 | pass |
| tempo | gfnff | molden | — | — | — | — / — | pass |
| tempo | gfnff | xtb_values | xtb | energy 1.5099e-14, gradient 1.52656e-15, charges 1.08524e-14 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| tempo | gxtb | gxtb_values | gxtb | energy 3.97904e-13, gradient 3.53884e-15, charges 4.82204e-09 | energy 2.1e-12, gradient 6e-11, charges — | 21 / 21 | pass |
| tempo | gxtb | molden | gxtb | orbital_energies 3.05311e-15 | orbital_energies 3.3e-12 | 11 / 11 | pass |
| tetramethylammonium | gfn2 | default_vs_xtb | axtb[xtb] | energy 9.37881e-11, gradient 2.30866e-11, charges 4.93539e-11, dipole 3.33222e-12 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 16 / 16 | pass |
| tetramethylammonium | gfn2 | molden | xtb | orbital_energies 1.9984e-15 | orbital_energies 7e-10 | 14 / 14 | pass |
| tetramethylammonium | gfn2 | tblite_values | tblite | energy 7.10543e-15, gradient 1.95788e-12, charges 6.40015e-11, dipole 6.69667e-12 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 16 / — | pass |
| tetramethylammonium | gfn2 | xtb_values | xtb | energy 1.3145e-13, gradient 7.30786e-11, charges 5.05278e-09, dipole 2.58624e-09 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 16 / 14 | pass |
| tetramethylammonium | gfnff | molden | — | — | — | — / — | pass |
| tetramethylammonium | gfnff | xtb_values | xtb | energy 2.10942e-13, gradient 4.57967e-16, charges 7.28584e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| tetramethylammonium | gxtb | gxtb_values | gxtb | energy 3.41061e-13, gradient 1.10328e-15, charges 4.82508e-09 | energy 2.1e-12, gradient 6e-11, charges — | 13 / 13 | pass |
| tetramethylammonium | gxtb | molden | gxtb | orbital_energies 3.9968e-15 | orbital_energies 3.3e-12 | 9 / 9 | pass |
| tm-MnO4- | gfn2 | default_vs_xtb | axtb[xtb] | energy 1.92227e-10, gradient 3.26292e-11, charges 2.19189e-11, dipole 8.12491e-11 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 13 / 13 | pass |
| tm-MnO4- | gfn2 | molden | xtb | orbital_energies 7.87037e-13 | orbital_energies 7e-10 | 12 / 12 | pass |
| tm-MnO4- | gfn2 | tblite_values | tblite | energy 1.42109e-14, gradient 5.63626e-13, charges 4.58418e-11, dipole 1.71582e-12 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 17 / — | pass |
| tm-MnO4- | gfn2 | xtb_values | xtb | energy 1.77636e-14, gradient 7.08902e-12, charges 4.27251e-09, dipole 5.38857e-11 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 13 / 12 | pass |
| tm-MnO4- | gfnff | molden | — | — | — | — / — | pass |
| tm-MnO4- | gfnff | xtb_values | xtb | energy 5.72875e-14, gradient 1.30737e-16, charges 1.16573e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| tm-MnO4- | gxtb | gxtb_values | gxtb | energy 2.27374e-13, gradient 1.85479e-10, charges 4.75014e-09 | energy 2.1e-12, gradient 6e-11, charges — | 15 / 15 | fail |
| tm-MnO4- | gxtb | molden | gxtb | orbital_energies 3.9299e-11 | orbital_energies 3.3e-12 | 11 / 11 | fail |
| twisted-ethylene | gfn2 | default_vs_xtb | axtb[xtb] | energy 2.50053e-09, gradient 2.795e-12, charges 2.88458e-12, dipole 3.22881e-13 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 14 / 14 | pass |
| twisted-ethylene | gfn2 | molden | xtb | orbital_energies 9.28591e-12 | orbital_energies 7e-10 | 11 / 11 | pass |
| twisted-ethylene | gfn2 | tblite_values | tblite | energy 4.44089e-15, gradient 3.80654e-13, charges 1.4443e-11, dipole 2.21438e-11 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 12 / — | pass |
| twisted-ethylene | gfn2 | xtb_values | xtb | energy 1.57208e-13, gradient 5.33784e-11, charges 1.32887e-09, dipole 1.52743e-11 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 14 / 11 | pass |
| twisted-ethylene | gfnff | molden | — | — | — | — / — | pass |
| twisted-ethylene | gfnff | xtb_values | xtb | energy 4.55191e-15, gradient 4.96131e-16, charges 1.58901e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| twisted-ethylene | gxtb | gxtb_values | gxtb | energy 4.40536e-13, gradient 4.7488e-09, charges 8.93069e-09 | energy 2.1e-12, gradient 6e-11, charges — | 15 / 16 | fail |
| twisted-ethylene | gxtb | molden | gxtb | orbital_energies 5.19169e-07 | orbital_energies 3.3e-12 | 11 / 11 | fail |
| tyr-xab | gfn2 | default_vs_xtb | axtb[xtb] | energy 2.41137e-10, gradient 2.01944e-11, charges 5.74595e-11, dipole 7.89251e-11 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 25 / 25 | pass |
| tyr-xab | gfn2 | molden | xtb | orbital_energies 1.19904e-14 | orbital_energies 7e-10 | 21 / 21 | pass |
| tyr-xab | gfn2 | tblite_values | tblite | energy 2.84217e-14, gradient 2.30718e-15, charges 1.19904e-14, dipole 4.76286e-14 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 26 / — | pass |
| tyr-xab | gfn2 | xtb_values | xtb | energy 2.55795e-13, gradient 1.85729e-10, charges 7.24041e-09, dipole 3.68824e-09 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 25 / 21 | pass |
| tyr-xab | gfnff | molden | — | — | — | — / — | pass |
| tyr-xab | gfnff | xtb_values | xtb | energy 2.20268e-13, gradient 8.96852e-16, charges 1.38223e-14 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| tyr-xab | gxtb | gxtb_values | gxtb | energy 1.13687e-13, gradient 1.61069e-15, charges 4.91911e-09 | energy 2.1e-12, gradient 6e-11, charges — | 16 / 16 | pass |
| tyr-xab | gxtb | molden | gxtb | orbital_energies 9.99201e-15 | orbital_energies 3.3e-12 | 9 / 9 | pass |
| uf6 | gfn2 | xtb_values | xtb | — | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | — / — | expected_error |
| uf6 | gfnff | molden | — | — | — | — / — | pass |
| uf6 | gfnff | xtb_values | xtb | energy 3.91243e-13, gradient 4.80518e-16, charges 2.9976e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| uf6 | gxtb | gxtb_values | gxtb | energy 2.27374e-13, gradient 2.21898e-12, charges 1.78652e-09 | energy 2.1e-12, gradient 6e-11, charges — | 14 / 14 | pass |
| uf6 | gxtb | molden | gxtb | orbital_energies 3.12017e-12 | orbital_energies 3.3e-12 | 9 / 9 | pass |
| water | gfn2 | default_vs_xtb | axtb[xtb] | energy 1.89671e-11, gradient 2.41336e-11, charges 3.75506e-11, dipole 3.92731e-11 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 12 / 12 | pass |
| water | gfn2 | molden | xtb | orbital_energies 2.9976e-15 | orbital_energies 7e-10 | 10 / 10 | pass |
| water | gfn2 | tblite_values | tblite | energy 8.88178e-16, gradient 1.82146e-16, charges 1.77636e-15, dipole 9.99201e-16 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 13 / — | pass |
| water | gfn2 | xtb_values | xtb | energy 1.7053e-13, gradient 3.24809e-11, charges 4.08446e-09, dipole 1.33606e-09 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 12 / 10 | pass |
| water | gfnff | molden | — | — | — | — / — | pass |
| water | gfnff | xtb_values | xtb | energy 3.94129e-15, gradient 7.76289e-17, charges 4.996e-16 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| water | gxtb | gxtb_values | gxtb | energy 1.42109e-14, gradient 5.6205e-16, charges 3.45216e-09 | energy 2.1e-12, gradient 6e-11, charges — | 11 / 11 | pass |
| water | gxtb | molden | gxtb | orbital_energies 2.79776e-14 | orbital_energies 3.3e-12 | 9 / 9 | pass |
| water-dimer | gfn2 | default_vs_xtb | axtb[xtb] | energy 3.75486e-11, gradient 3.26145e-11, charges 4.3092e-11, dipole 9.40292e-11 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 15 / 15 | pass |
| water-dimer | gfn2 | hessian_analytic | xtb | frequencies 0.102789 | frequencies 6.6 | 2 / 12 | pass |
| water-dimer | gfn2 | hessian_same | xtb | frequencies 4.8457e-05 | frequencies 5e-05 | — / 12 | pass |
| water-dimer | gfn2 | molden | xtb | orbital_energies 5.05151e-15 | orbital_energies 7e-10 | 13 / 13 | pass |
| water-dimer | gfn2 | tblite_values | tblite | energy 0, gradient 3.07046e-16, charges 4.88498e-15, dipole 1.9984e-15 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 16 / — | pass |
| water-dimer | gfn2 | xtb_values | xtb | energy 2.07834e-13, gradient 1.77715e-10, charges 4.48701e-09, dipole 2.3138e-09 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 15 / 13 | pass |
| water-dimer | gfnff | hessian_analytic | xtb | frequencies 0.438898 | frequencies 0.27 | — / — | fail |
| water-dimer | gfnff | hessian_numerical | xtb | frequencies 0.00387746 | frequencies 0.005 | — / — | pass |
| water-dimer | gfnff | hessian_same | xtb | frequencies 0.00387696 | frequencies 0.005 | — / — | pass |
| water-dimer | gfnff | molden | — | — | — | — / — | pass |
| water-dimer | gfnff | xtb_values | xtb | energy 3.22742e-13, gradient 1.00614e-16, charges 5.55112e-16 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| water-dimer | gxtb | gxtb_values | gxtb | energy 1.42109e-13, gradient 8.6129e-16, charges 4.23842e-09 | energy 2.1e-12, gradient 6e-11, charges — | 12 / 12 | pass |
| water-dimer | gxtb | hessian_analytic | gxtb | frequencies 0.146197 | frequencies 0.35 | 9 / 12 | pass |
| water-dimer | gxtb | hessian_same | gxtb | frequencies 4.77495e-05 | frequencies 5e-05 | — / 12 | pass |
| water-dimer | gxtb | molden | gxtb | orbital_energies 3.9968e-15 | orbital_energies 3.3e-12 | 8 / 8 | pass |
| zncl4 | gfn2 | default_vs_xtb | axtb[xtb] | energy 1.95897e-11, gradient 2.10226e-11, charges 2.72093e-11, dipole 1.57652e-14 | energy 1e-07, gradient 1e-06, charges 1e-05, dipole 0.0001 | 15 / 15 | pass |
| zncl4 | gfn2 | molden | xtb | orbital_energies 7.99361e-15 | orbital_energies 7e-10 | 14 / 14 | pass |
| zncl4 | gfn2 | tblite_values | tblite | energy 1.77636e-14, gradient 6.80175e-13, charges 2.68503e-11, dipole 1.12688e-14 | energy 1e-09, gradient 1e-08, charges 1e-07, dipole 1e-06 | 17 / — | pass |
| zncl4 | gfn2 | xtb_values | xtb | energy 3.09086e-13, gradient 5.66825e-11, charges 2.56989e-09, dipole 1.37668e-14 | energy 1e-07, gradient 1e-07, charges 1e-06, dipole 1e-05 | 15 / 14 | pass |
| zncl4 | gfnff | molden | — | — | — | — / — | pass |
| zncl4 | gfnff | xtb_values | xtb | energy 2.91323e-13, gradient 4.68375e-17, charges 1.88738e-15 | energy 7e-13, gradient 4.3e-12, charges 2.3e-12 | — / — | pass |
| zncl4 | gxtb | gxtb_values | gxtb | energy 9.09495e-13, gradient 5.77901e-13, charges 3.09349e-09 | energy 2.1e-12, gradient 6e-11, charges — | 12 / 12 | pass |
| zncl4 | gxtb | molden | gxtb | orbital_energies 9.58011e-13 | orbital_energies 3.3e-12 | 9 / 9 | pass |

40 cases, 380 checks: 364 pass, 11 fail, 1 expected_error, 2 expected_difference, 2 both_failed, 0 error
