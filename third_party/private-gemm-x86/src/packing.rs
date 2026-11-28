// SPDX-License-Identifier: MIT
//! Packing storage is bounded by both the cache tile and the actual matrix panel.

#[repr(align(4096))]
pub(crate) struct Page([u8; 4096]);

pub(crate) fn buffer_pages(cache_elements: usize, extent: usize, padding: usize, depth: usize, element_size: usize) -> usize {
	// Recursive packed offsets advance by logical matrix coordinates, so the
	// whole padded dimension bounds both outer offsets and the final leaf pack.
	// LHS packs pad to at most 16 elements; RHS storage rounds to microkernel width.
	let cache_elements = cache_elements.next_multiple_of(padding);
	let panel_elements = extent.next_multiple_of(padding).saturating_mul(depth);
	(cache_elements.min(panel_elements) * element_size).div_ceil(core::mem::size_of::<Page>())
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn water_panel_does_not_reserve_the_hosts_large_cache() {
		let l3_elements = 260 * 1024 * 1024 / 8;
		assert_eq!(buffer_pages(l3_elements, 6, 16, 5, 8), 1);
		assert_eq!(buffer_pages(l3_elements, 6, 4, 5, 8), 1);
	}

	#[test]
	fn panel_and_cache_each_bound_the_allocation() {
		let l3_elements = 8 * 1024 * 1024 / 8;
		assert_eq!(buffer_pages(l3_elements, 1024, 16, 512, 8), 1024);
		assert_eq!(buffer_pages(l3_elements, 4096, 16, 512, 8), 2048);
	}

	#[test]
	fn complex_and_partial_micro_panels_include_padding() {
		let l3_elements = 260 * 1024 * 1024 / 16;
		assert_eq!(buffer_pages(l3_elements, 17, 16, 33, 16), 5);
		assert_eq!(buffer_pages(l3_elements, 17, 4, 33, 16), 3);
	}

	#[test]
	fn partial_simd_panels_pad_before_multiplying_by_depth() {
		let l3 = 260 * 1024 * 1024;
		assert_eq!(buffer_pages(l3 / 4, 17, 16, 512, 4), 16);
		assert_eq!(buffer_pages(l3 / 8, 15, 16, 512, 8), 16);
	}

	#[test]
	fn cache_bound_avoids_overflow_from_a_large_panel() {
		assert_eq!(buffer_pages(1024, usize::MAX / 32, 16, 512, 8), 2);
	}

	#[test]
	fn whole_dimensions_cover_recursive_outer_offsets() {
		let l3_elements = 2 * 1024 * 1024 / 8;
		assert_eq!(buffer_pages(l3_elements, 10001, 16, 512, 8), 512);
		assert_eq!(buffer_pages(l3_elements, 777, 4, 512, 8), 512);
	}

	#[test]
	fn empty_panels_need_no_storage() {
		assert_eq!(buffer_pages(1024, 0, 16, 5, 8), 0);
		assert_eq!(buffer_pages(1024, 6, 16, 0, 8), 0);
	}
}
