# axtb for ReShiki: validation

This fork (branch `reshiki`) is qualified as ReShiki's xTB engine. This page records what was
compared, against which programs, with which tolerances, and every result. The scripts are in
`reshiki/` (see `reshiki/README.md`); the raw results are in `docs/reshiki/results/`.

- Engine source validated: `449483a8fa7979e3568ef25d34863c5e9f3a8f99` (upstream axtb 0.2.0,
  `65d0f64`, plus the optional-g-xTB feature; `CHANGES-reshiki.md`). The commit that adds this
  page changes no engine code.
- Build on every platform: `cargo +1.99.1 build --release --locked -p axtb-cli`, default
  features, no `RUSTFLAGS` or target-CPU flags, `AXTB_LAPACK` and `AXTB_LAPACK_BLAS` unset
  (every response reports the faer eigensolver), one thread per calculation.
- Date: 2026-11-25.

## Summary

| Area | Result |
|---|---|
| Reference comparison, Linux x64 | 40 structures, 380 checks: **364 pass, 11 fail, 1 expected error, 2 expected differences, 2 both failed, 0 errors** |
| GFN2-xTB vs xtb 6.8.0 (energy, gradient, charges, dipole, Molden orbital energies) | 39 of 39 pass (UF6 refused by both, as expected); worst ΔE 7.4e-13 Eh, ΔG 1.0e-9 Eh/bohr |
| GFN2-xTB vs tblite 0.8.0 | 38 of 38 pass; gh106 does not converge in either (documented) |
| GFN2-xTB Hessians vs `xtb --hess` | analytic 9/9 (worst 0.18 cm⁻¹), same Hessian 9/9 (5.0e-5 cm⁻¹) |
| GFN-FF vs xtb 6.8.0 | values 40/40 (ΔE 5.0e-13 Eh), same Hessian 9/9, axtb numerical Hessian 9/9 (≤ 0.005 cm⁻¹); analytic Hessian 6/9 (3 at 0.32–0.44 cm⁻¹ vs the documented 0.27) |
| g-xTB 2.0.1 vs the g-xTB 2.0.1 program | values 36/39 + 1 both failed; Molden 36/40; Hessians analytic 9/9 (0.348 cm⁻¹ vs 0.35), same Hessian 8/9 |
| Cross-platform (axtb vs axtb on Linux x64) | Linux arm64 and macOS arm64: 144 of 147 runs within the gates (worst ΔE 3.7e-12 Eh, ΔG 7.0e-8 Eh/bohr, Δν 6.9e-9 cm⁻¹); 3 runs differ: SCFs with several solutions (below) |
| Windows x64 | `cargo check --target x86_64-pc-windows-msvc -p axtb-driver -p axtb-cli` and `-p axtb-driver --no-default-features`: pass (Rust 1.99.1, MSVC); no run-time comparison |
| `gxtb` cargo feature | default on; `--no-default-features` removes `axtb-gxtb` from the graph (`cargo tree`: 0 matches); tests pass both ways |
| Run-time LAPACK loader | Linux only (`cfg(target_os = "linux")`; elsewhere `open` returns an error), and used only when `AXTB_LAPACK` is set |

No failure is an axtb implementation error. Every failure is one of: the reference program's
own finite-difference or print precision, a documented limit measured on other molecules and
exceeded at round-off level, or an SCF with more than one solution. Details below.

## Reference programs (official release binaries, Linux x64)

| Program | Release asset | Asset SHA-256 | Executable SHA-256 |
|---|---|---|---|
| xtb 6.8.0 (4b8c987) | https://github.com/grimme-lab/xtb/releases/download/v6.8.0/xtb-6.8.0-linux-x86_64-ifx-mkl.tar.xz | `db090333361749ae0d9ccf0ce7c0035c0000f2ebad0847f6d77191d1954f4d68` | `bin/xtb` `95d82a623580bb1ef7e1f9aec6a91d3773b79e5aabf7bb8fada287debf65c93b` |
| g-xTB 2.0.1 (xtb 6.7.1 26dd68d, `xtb --gxtb`) | https://github.com/grimme-lab/g-xtb/releases/download/v2.0.1/xtb-6.7.1-gxtb-140526-linux-x86_64.tar.xz | `cdfb81164f8d3300fc451dbbf2414528b533e71f7a30ec824b564cd30d6d9e5c` | `bin/xtb` `1b4e30b68ed4e88b4075f60d97f4756ee440fe92d3294cc20826d53f8121cd26` |
| tblite 0.8.0 (no CLI in the release; the release wheel, Python 3.12) | https://github.com/tblite/tblite/releases/download/v0.8.0/tblite-0.8.0-cp312-cp312-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl | `e82fece754c5423eca7d49e2c9a328f66df0386b075893590007a4a51a744f73` | — |

The asset hashes equal the digests GitHub lists for the releases. `validate.py` refuses to run
when a program is missing or reports another version. The xtb programs run with `XTBPATH`
unset and `XTBHOME` set to their own scratch directory, `OMP_NUM_THREADS=1`. The host CPU is an
AMD Ryzen 9 7940HS: the xtb 6.8.0 binary takes its baseline (non-AVX2) code path there.

## Method and tolerances

`reshiki/validate.py` runs each of the 40 structures (`reshiki/structures/`: 25 ReShiki
molecules: drug-like organics, charged, open-shell, linear, metal complexes, a hard SCF and
UF6; plus axtb's `tests/structures/` and two of its `tests/scf_hard/` cases) through axtb and
the reference program with the same bohr coordinates, charge and unpaired electrons, in the
conventions that reproduce that program. Hessians use geometries optimized by the reference
program (`reshiki/optimize_structures.py`, `--opt tight`; `reshiki/structures/opt/`).

Every tolerance is a number that axtb's own documentation or gates state; the `source` of each
is in `reshiki/cases.py`. As in axtb's `tests/gate.py`, a deviation passes when it is at most
its limit at two significant digits.

| Method | Check | Settings | Tolerance (source) |
|---|---|---|---|
| GFN2-xTB | values vs xtb 6.8.0, conventions `xtb` | axtb accuracy 1e-6, xtb `--acc 1e-4` | E 1e-7 Eh, G 1e-7 Eh/bohr, q 1e-6 e, μ 1e-5 e·bohr (`tests/reference_matrix/run.py` TOL["xtb"]) |
| GFN2-xTB | values vs tblite 0.8.0, conventions `tblite` | accuracy 1e-6 | 1e-9 / 1e-8 / 1e-7 / 1e-6 (TOL["tblite"]) |
| GFN2-xTB | default conventions vs `xtb` conventions (axtb only) | same bohr input | 1e-7 / 1e-6 / 1e-5 / 1e-4 (TOL["axtb_vs_xtbmode"]); molecules longer than 25 bohr: 9.4e-6 Eh, 3.1e-6 Eh/bohr (documented cutoff difference) |
| GFN2-xTB | Molden orbital energies | both `--acc 1e-4` | 7.0e-10 Eh (docs/validation.md, Molden) |
| GFN2-xTB | analytic Hessian vs `xtb --hess --acc 0.001` | modes above 50 cm⁻¹ | 6.6 cm⁻¹ (reference matrix, analytic_xtb) |
| GFN2-xTB | same Hessian: axtb's analysis of xtb's `hessian` file | `analysis_conventions` `xtb` | 5.0e-5 cm⁻¹ (`tests/vibrational_check.py`) |
| GFN-FF | values vs `xtb --gfnff` | conventions `xtb` | E 7.0e-13 Eh, G 4.3e-12 Eh/bohr, EEQ q 2.3e-12 e (docs/validation.md, GFN-FF) |
| GFN-FF | analytic Hessian / numerical Hessian (xtb's step 0.005) / same Hessian | | 0.27 / 0.005 / 0.005 cm⁻¹ (docs/validation.md; `tests/gfnff_reproduction.py`) |
| g-xTB 2.0.1 | values vs `xtb --gxtb` | conventions `gxtb`, both `--acc 1e-4` | E 2.1e-12 Eh, G 6e-11 Eh/bohr (`tests/gxtb2_reproduction.py`); charges reported, not gated |
| g-xTB 2.0.1 | Molden orbital energies | program default accuracy | 3.3e-12 Eh (docs/validation.md, Molden) |
| g-xTB 2.0.1 | analytic Hessian vs `xtb --gxtb --hess --acc 1e-4` / same Hessian | | 0.35 cm⁻¹ / 5.0e-5 cm⁻¹ (docs/validation.md; `tests/vibrational_check.py`) |

## The failures, explained

| Case | Check | Measured | Limit | Why |
|---|---|---|---|---|
| acetate, caffeine, water dimer | GFN-FF analytic Hessian vs `xtb --hess` | 0.42, 0.32, 0.44 cm⁻¹ | 0.27 | xtb's finite-difference Hessian error. axtb's numerical Hessian with xtb's step matches xtb to ≤ 0.005 cm⁻¹ on all 9 cases, and the same-Hessian check passes, so the analytic Hessian is exact and xtb's is not. The 0.27 limit was measured on other molecules. |
| caffeine | g-xTB same Hessian | 5.06e-5 cm⁻¹ | 5.0e-5 | xtb prints frequencies with 4 decimals; the limit is half the last digit. 6e-7 cm⁻¹ over, at a rounding boundary. |
| [Cu(NH3)4]²⁺ | g-xTB Molden orbital energies | 3.6e-12 Eh | 3.3e-12 | Round-off (the documented 3.3e-12 includes this molecule, measured on another CPU). |
| taxol | g-xTB energy | 2.7e-12 Eh | 2.1e-12 | 9e-16 relative (|E| = 2930 Eh); axtb's gate gives large molecules their own limits (C50H102 2.05e-12, (H2O)100 9.1e-12). |
| [MnO4]⁻, twisted ethylene | g-xTB gradient, Molden | 1.9e-10 / 4.7e-9 Eh/bohr; 3.9e-11 / 5.2e-7 Eh | 6e-11; 3.3e-12 | Hard SCFs (small gap; twisted ethylene has a degenerate π pair filled by aufbau at g-xTB's 0 K). Energies agree to 2e-13 / 4e-13 Eh. Same behaviour as the documented near-degenerate cases (CeCl3 4.0e-9 Eh/bohr). |
| gh106 | g-xTB Molden (program default accuracy 1) | 4.5e-3 Eh | 3.3e-12 | gh106 is axtb's documented pathological SCF (several solutions, tests/structures/README.md); at accuracy 1 the two programs land on different solutions; at `--acc 1e-4` neither converges in 250 iterations (both_failed). |

The two expected differences (taxol, gh106: default vs `xtb` conventions, 9.4e-6 and 2.6e-6
Eh) are the documented real-space cutoff change of xtb 6.8.0. UF6 with GFN2-xTB is refused by
both programs (no parameters for Z = 92); GFN-FF and g-xTB compute it, and both match.

## Cross-platform

`reshiki/platform_dump.py` runs the axtb CLI alone (default conventions, the three methods,
energy, gradient, charges; analytic Hessians for the 9 Hessian cases) and
`reshiki/platform_compare.py` compares each platform with Linux x64. Gates (chosen for ReShiki,
not upstream numbers): |ΔE| ≤ 1e-8 Eh, max |ΔG| ≤ 1e-7 Eh/bohr, max |Δν| ≤ 0.1 cm⁻¹, identical
inputs (fingerprints), identical errors only where declared (UF6 / GFN2).

| Platform | Binary SHA-256 | Runs | Within gates | Worst of those (ΔE / ΔG / Δν) | Outside gates |
|---|---|---:|---:|---|---|
| Linux x64 (Arch, AMD Ryzen 9 7940HS) | `c6e121205c8aafb132a5f2c9769d8e17aa65ed6c372f20b29c885adc6aab0fb1` | 147 | baseline | — | — |
| Linux arm64 (Ubuntu, NVIDIA GB10, Cortex-X925) | `f82157f507d5d66fd00663215fd7d277bf9fa0dede65d27ec20999fbd9225375` | 147 | 144 | 8.3e-13 Eh / 1.1e-8 Eh/bohr / 6.9e-9 cm⁻¹ | gh106 GFN2, gh106 g-xTB, twisted ethylene g-xTB |
| macOS arm64 (Apple M1) | `6d513300578ee4eb0a3dd2ec9fb4883effd7ca5a38ca0df0caf33f6239716099` | 147 | 144 | 3.7e-12 Eh / 7.0e-8 Eh/bohr / 6.0e-9 cm⁻¹ | the same three |
| Windows x64 | — (cargo check only) | — | — | — | not run |

The three runs outside the gates are SCFs with several solutions, where round-off picks the
solution: gh106 (ΔE up to 1.8e-3 Eh, different iteration counts) and twisted ethylene with
g-xTB (energy equal to 1e-9 Eh, gradients differ by 0.11 Eh/bohr: aufbau at 0 K occupies one
of two degenerate π orbitals, and which one depends on the eigensolver's round-off; at 300 K
the platforms still reach different fractional-occupation solutions, −78.3231 vs −78.3168 Eh).
No tolerance was widened for them. Every other run, including all Hessians, agrees to
round-off.

## Notes for ReShiki

- **Multiple SCF solutions are platform-dependent.** For molecules with degenerate frontier
  orbitals (twisted or square π systems, radicals at symmetric geometries, metal clusters)
  g-xTB at its default 0 K, and sometimes GFN2-xTB, can give different gradients (and so
  different optimized structures) on different machines. ReShiki should report the HOMO–LUMO
  gap and fractional occupations with each result and warn when the gap is near zero.
- **`[patch.crates-io]` does not carry over.** axtb's workspace patches `gemm-common` and
  `private-gemm-x86` (`third_party/`, start-up cost on many-CPU Linux machines). Cargo ignores
  patches of a dependency's workspace; ReShiki's root `Cargo.toml` must repeat them (git source
  of this fork at the pinned revision) to get the same dependency graph that was validated, or
  accept the upstream crates.
- **Clear `AXTB_LAPACK` and `AXTB_LAPACK_BLAS`** in the environment of the calculation worker,
  so a user's environment cannot load a LAPACK library into ReShiki (Linux).
- **Default conventions** (`axtb`) differ from xtb 6.8.0 by ≤ 1e-7 Eh except for molecules
  longer than 25 bohr (≤ 9.4e-6 Eh, the cutoffs). Use `conventions: "xtb"` when results must
  match the xtb program.

## Upstream-worthy findings (not reported upstream)

- xtb 6.8.0 leaves orbitals above 2 Eh out of its Molden file (ferrocene: 58 of 59); axtb's
  own `tests/molden_check.py` already compares the common ones.
- The documented GFN-FF "analytic Hessian ≤ 0.27 cm⁻¹" bound is exceeded by xtb's own
  finite-difference error on other molecules (up to 0.44 cm⁻¹ here).

## Reproduce

```sh
# Linux x64, with the three reference programs
export AXTB_XTB680=.../xtb680/stage/bin/xtb AXTB_GXTB2=.../xtb-6.7.1/bin/xtb
export RESHIKI_TBLITE_PYTHON=.../venv/bin/python AXTB_BIN=target/release/axtb
$RESHIKI_TBLITE_PYTHON reshiki/validate.py --jobs 8 --out out
# any platform
uv run --no-project --with numpy python reshiki/platform_dump.py --out this-platform.json
uv run --no-project --with numpy python reshiki/platform_compare.py linux-x64.json this-platform.json --md compare.md
```

The full per-check and per-run tables follow (generated: `reference-linux-x64.md`,
`compare-linux-arm64.md`, `compare-macos-arm64.md` in `docs/reshiki/results/`).


## Appendix A: reference comparison, every check (Linux x64)


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


## Appendix: platform comparison, linux-arm64 vs Linux x64, every run


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


## Appendix: platform comparison, macos-arm64 vs Linux x64, every run


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
