// @ts-check

const MAX_DISTANCE = 100;

/**
 * @param {string} sequence
 * @param {number} max_distance
 */
function find_SCD_v4(sequence, max_distance) {
	const motif_len = "SQ".length;
	const sq_tq = [...sequence.matchAll(/[ST]Q/g)];
	const group = sq_tq.map((a, head_idx) => {
		const list = sq_tq.slice(head_idx).filter(b => {
			const end_pos = /** @type {number} */(b.index) + motif_len - 1;
			const len = end_pos - /** @type {number} */(a.index) + 1;
			return len <= max_distance;
		}).map(a => /** @type {number} */(a.index));
		if (list.length >= 3) {
			const end = list[list.length - 1] + motif_len -1
			const seq = sequence.slice(list[0], end + 1);
			// console.log(a.index + 1, end + 1, seq.length, seq);
			return {
				start: list[0] + 1,
				end : end + 1,
				list: list,
				seq: seq,
			};
		}
	}).filter(a => a);
	return group;
}

module.exports.find_SCD_v4 = find_SCD_v4;
module.exports.MAX_DISTANCE = MAX_DISTANCE;
