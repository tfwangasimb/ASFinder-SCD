// @ts-check


/**
 * @param {string} ss AA sequence
 * @param {number} min_num_AA number of Q each polyQ region default=4
 * @param {string} AA
 * @param {number} max_num_not_AA pos(Q2) - pos(Q1) < max_num_not_AA
 */
function any_polyAA_impure(ss, min_num_AA, AA, max_num_not_AA = 2) {
	// if (max_num_not_AA != 3) {
	// 	throw new Error();
	// }
	if (!AA) {
		throw new Error();
	}
	const mm_arr = [...ss.matchAll(RegExp(`[^${AA}]{${max_num_not_AA},}`, "g"))].map(mm => { // search not QQ
		/** @type {Number} */
		const i = /** @type {any} */(mm.index); // dynamic_cast
		return {
			idx_start: i,
			idx_end: i + mm[0].length,
		};
	});
	if (mm_arr.length) {
		/** @type {{ idx_start: number; idx_end: number; }[]} */
		const region = [];
		if (mm_arr[0].idx_start != 0) {// first fragment
			region.push({
				idx_start: 0,
				idx_end: mm_arr[0].idx_start - 1,
			});
		}
		if (mm_arr[mm_arr.length - 1].idx_end != ss.length) {// last fragment
			region.push({
				idx_start: mm_arr[mm_arr.length - 1].idx_end,
				idx_end: ss.length,
			});
		}
		// 111 xxx 111 xxx 111 xxx 111
		// 111 xxx
		//     xxx 111 xxx
		//             xxx 111 xxx
		//                     xxx 111
		mm_arr.slice(0, -1).forEach((m, idx) => {
			const next = mm_arr[idx + 1];
			const idx_start = m.idx_end;
			const idx_end = next.idx_start;
			region.push({
				idx_start,
				idx_end,
			});
		});

		// console.log(...region);

		return region.map(a => {
			const frag = ss.slice(a.idx_start, a.idx_end);
			const num_AA = calc_num_AA(frag, AA);
			if (num_AA >= min_num_AA) {
				const len = frag.length;
				const is_pure = num_AA == len;
				return {
					index: a.idx_start,
					end_idx: a.idx_end,
					num_AA: num_AA,
					length: len,
					type: is_pure ? num_AA.toString() : `${num_AA}/${len}`,
					frag: frag,
					is_pure: is_pure,
				};
			}
		}).filter(a => a);
	}
	return [];// fill empty
}
/**
 * @param {string} ss
 * @param {string} AA
 */
function calc_num_AA(ss, AA) {
	let num = 0;
	for (let i = 0; i < ss.length; ++i) {
		if (ss[i] == AA) {
			++num;
		}
	}
	return num;
}

/**
 * @param {string} ss
 * @param {number} min_num_AA default=4
 * @param {string[]} AA_list
 * @param {number} max_num_not_AA
 */
function any_polyXX_impure(ss, min_num_AA, AA_list, max_num_not_AA = 2) {
	if (!Array.isArray(AA_list)) {
		throw new Error();
	}
	// QNQN^^QN -> 6/8
	const mm_arr = [...ss.matchAll(RegExp(`[^${AA_list.join("")}]{${max_num_not_AA},}`, "g"))].map(mm => {
		const i = mm.index;
		if (typeof i === "undefined") {
			throw new TypeError();
		}
		return {
			idx_start: i,
			idx_end: i + mm[0].length,
		};
	});
	if (mm_arr.length) {
		/** @type {{ idx_start: number; idx_end: number; }[]} */
		const region = [];
		if (mm_arr[0].idx_start != 0) {// first fragment
			region.push({
				idx_start: 0,
				idx_end: mm_arr[0].idx_start - 1,
			});
		}
		if (mm_arr[mm_arr.length - 1].idx_end != ss.length) {// last fragment
			region.push({
				idx_start: mm_arr[mm_arr.length - 1].idx_end,
				idx_end: ss.length,
			});
		}
		// 111 xxx 111 xxx 111 xxx 111
		// 111 xxx
		//     xxx 111 xxx
		//             xxx 111 xxx
		//                     xxx 111
		mm_arr.slice(0, -1).forEach((m, idx) => {
			const next = mm_arr[idx + 1];
			const idx_start = m.idx_end;
			const idx_end = next.idx_start;
			region.push({
				idx_start,
				idx_end,
			});
		});

		// console.log(...region);

		return region.map(a => {
			const frag = ss.slice(a.idx_start, a.idx_end);
			const num_AA = calc_num_XX(frag, AA_list);
			if (num_AA >= min_num_AA) {
				const len = frag.length;
				const is_pure = num_AA == len;
				return {
					index: a.idx_start,
					end_idx: a.idx_end,
					num_AA: num_AA,
					length: frag.length,
					type: is_pure ? num_AA : `${num_AA}/${len}`,
					frag: frag,
					is_pure: is_pure,
				};
			}
		}).filter(a => a);
	}
	return [];// fill empty
}
/**
 * @param {string} ss
 * @param {string[]} AA_list
 */
function calc_num_XX(ss, AA_list) {
	let num = 0;
	for (let i = 0; i < ss.length; ++i) {
		if (AA_list.includes(ss[i])) {
			++num;
		}
	}
	return num;
}

/**
 * @param {string} input_seq
 * @param {number} min_num_AA
 * @param {string|string[]} amino_acid
 * @param {number} max_num_not_AA
 */
function findPolyAminoAcid(input_seq, min_num_AA, amino_acid, max_num_not_AA = 2) {
	if (Array.isArray(amino_acid)) {
		return any_polyXX_impure(input_seq, min_num_AA, amino_acid, max_num_not_AA);
	}
	else {
		return any_polyAA_impure(input_seq, min_num_AA, amino_acid, max_num_not_AA);
	}
}

module.exports.any_polyAA_impure = any_polyAA_impure;
module.exports.any_polyXX_impure = any_polyXX_impure;
module.exports.findPolyAminoAcid = findPolyAminoAcid;
