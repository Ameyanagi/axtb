// SPDX-License-Identifier: MIT
//! Standalone regression executable driven by scripts/check_packing_budget.py.

use private_gemm_x86::{Accum, DType, DstKind, IType, InstrSet, gemm};
use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering};

const HEAP_BYTES: usize = 256 * 1024 * 1024;
static USED: AtomicUsize = AtomicUsize::new(0);
static LARGEST: AtomicUsize = AtomicUsize::new(0);
struct LimitedHeap;

#[global_allocator]
static ALLOCATOR: LimitedHeap = LimitedHeap;

unsafe impl GlobalAlloc for LimitedHeap {
	// fetch_update supports this vendored crate's Rust 1.84 minimum.
	#[allow(deprecated)]
	unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
		let bytes = layout.size();
		LARGEST.fetch_max(bytes, Ordering::Relaxed);
		if let Err(used) = USED.fetch_update(Ordering::Relaxed, Ordering::Relaxed, |used| {
			used.checked_add(bytes).filter(|&total| total <= HEAP_BYTES)
		}) {
			eprintln!("PACKING_HEAP_LIMIT {HEAP_BYTES} {used} {bytes}");
			std::process::exit(75);
		}
		let data = unsafe { System.alloc(layout) };
		if data.is_null() {
			USED.fetch_sub(bytes, Ordering::Relaxed);
		}
		data
	}

	unsafe fn dealloc(&self, data: *mut u8, layout: Layout) {
		unsafe { System.dealloc(data, layout) };
		USED.fetch_sub(layout.size(), Ordering::Relaxed);
	}
}

fn check(m: usize, n: usize, k: usize, kind: DstKind, threads: usize) {
	let a: Vec<f64> = (0..m * k).map(|i| ((i % m + 3 * (i / m)) % 5) as f64 - 2.0).collect();
	let b: Vec<f64> = (0..k * n).map(|i| ((2 * (i % k) + i / k) % 7) as f64 - 3.0).collect();
	let mut c = vec![123.0f64; m * n];
	let alpha = 1.0f64;
	unsafe {
		gemm(
			DType::F64,
			IType::U32,
			InstrSet::Avx256,
			m,
			n,
			k,
			c.as_mut_ptr().cast(),
			1,
			m as isize,
			std::ptr::null(),
			std::ptr::null(),
			kind,
			Accum::Replace,
			a.as_ptr().cast(),
			1,
			m as isize,
			false,
			std::ptr::null(),
			0,
			b.as_ptr().cast(),
			1,
			k as isize,
			false,
			(&alpha as *const f64).cast(),
			threads,
		);
	}
	for col in 0..n {
		for row in 0..m {
			let active = match kind {
				DstKind::Full => true,
				DstKind::Lower => row >= col,
				DstKind::Upper => row <= col,
			};
			let expected = if active {
				(0..k).map(|z| a[row + z * m] * b[z + col * k]).sum()
			} else {
				123.0
			};
			assert_eq!(c[row + col * m], expected, "{m}x{n} k{k} ({row},{col}) {kind:?} threads={threads}");
		}
	}
}

fn suite(threads: usize) {
	// The first product is the six-basis-function water density that exposed the failure.
	check(6, 6, 5, DstKind::Lower, threads);
	for (m, n, k) in [(1, 1, 1), (7, 9, 17), (17, 31, 65), (49, 67, 513), (131, 29, 529), (51, 137, 129)] {
		for kind in [DstKind::Full, DstKind::Lower, DstKind::Upper] {
			check(m, n, k, kind, threads);
		}
	}
	// With the synthetic 260 MiB cache and k=1, the tall-kernel cutoff is 532480
	// rows. Crossing it forces real RHS packing and recursive outer row offsets.
	let tall_rows = (260 * 1024 * 1024 / std::mem::size_of::<f64>() / 64) + 1;
	check(tall_rows, 3, 1, DstKind::Full, threads);
}

fn main() {
	assert!(
		std::arch::is_x86_feature_detected!("avx2") && std::arch::is_x86_feature_detected!("fma"),
		"This regression requires x86 AVX2 and FMA"
	);
	suite(1);
	rayon::ThreadPoolBuilder::new().num_threads(4).build().unwrap().install(|| suite(4));
	println!(
		"40 exact-result packing cases passed; largest allocation {} bytes",
		LARGEST.load(Ordering::Relaxed)
	);
}
