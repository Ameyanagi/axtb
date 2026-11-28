# Vendored, patched crates

| crate | version | upstream | license | patch |
|---|---|---|---|---|
| gemm-common | 0.18.2 | https://github.com/sarah-quinones/gemm | MIT | `src/cache.rs`: read only cpu0's cache files |
| private-gemm-x86 | 0.1.22 | https://github.com/sarah-quinones/gemm-x64-v2 | MIT | `src/cache.rs`: read only cpu0's cache files; `src/lib.rs` and `src/packing.rs`: bound packing buffers by matrix panels |

Both crates (dependencies of faer 0.22) detect cache sizes by opening every file under
`/sys/devices/system/cpu/cpu*/cache/index*/` the first time a matrix product runs: about 5,800
opens per process on a 128-thread machine. When many processes start together the kernel
serializes these, and each process spends ~0.4 s in start-up (measured on a 64-core / 128-thread
AMD node: 1024 axtb CLI runs on 128 workers took 5.87 s; with this patch 0.52 s; see
bench/cluster/cli_throughput.py and bench/results/cluster-9985wx/cli_throughput.json).

The upstream loop overwrites its result with each CPU it visits, so the answer is the last CPU
in directory order; on machines whose cores share one cache layout, reading cpu0 gives the
same answer. The `gemm-common` copy can be removed once upstream releases its cache-detection fix. The
`private-gemm-x86` copy also carries the packing-storage fix below; both patches must be
present upstream before removing that copy and its `[patch.crates-io]` entry.

Reported upstream as [sarah-quinones/gemm#53](https://github.com/sarah-quinones/gemm/issues/53)
(gemm-common) and [sarah-quinones/gemm-x64-v2#8](https://github.com/sarah-quinones/gemm-x64-v2/issues/8)
(private-gemm-x86); both are open, and upstream's `cache.rs` is unpatched.

The patch applies only to builds of this workspace: the CLI, the Python wheels and anything else
built inside this repository, and builds from the Python sdist, which carries this directory
(`crates/axtb-py/third_party` is a symbolic link to it, for maturin; see `pyproject.toml`). Cargo ignores the `[patch]` sections of dependencies, so a project
that depends on the axtb crates from crates.io or from this repository's URL builds the
unpatched gemm crates and pays the start-up cost on many-CPU machines, unless it adds the same
`[patch.crates-io]` entries to its own workspace.

## Matrix-packing storage

The x86 kernel previously reserved two whole detected L3 caches in thread-local storage
before deciding whether a matrix product needed packing. With a 260 MiB detected cache, even
the six-basis-function water density product requested 520 MiB and failed a 256 MiB heap
budget. The copy here bounds each existing cache-sized buffer by the padded whole matrix
panel times the packed depth. Whole dimensions account for recursive outer offsets; LHS
padding is at most 16 elements and RHS storage rounds to the microkernel width. TLS starts empty and
grows to the required storage. Cache detection, tiling and arithmetic are unchanged.

Run the focused regression on x86-64 with AVX2 and FMA:

```sh
python3 scripts/check_packing_budget.py
```

It builds the real vendored kernel in a temporary source copy with a synthetic 260 MiB L3
cache. A 256 MiB allocator budget covers 40 exact-result cases with one and four threads,
including odd dimensions, depth tails beyond 512 and a tall matrix that forces both packed buffers
and recursive offsets. The original kernel exits with status
75 on the first water-sized product. This test uses Rust dependencies only and avoids the
vendored crate's BLAS-dependent development tests. CI runs it separately because
`third_party/` is excluded from the workspace. `--kernel-source PATH` selects an older kernel
copy for a before/after check; the synthetic cache override never changes production source.

Portable sizing and padding tests can also run without building x86 assembly:

```sh
rustc --edition=2021 --test third_party/private-gemm-x86/src/packing.rs -o /tmp/axtb-packing-tests
/tmp/axtb-packing-tests
```
