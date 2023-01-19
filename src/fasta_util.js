//@ts-check

/**
 * 2020/12/15: add reverseComplementEx
 * 2021/11/25: output NCBI format
 */

const fs = require("fs");
const Path = require("path");

const VERBOSE = process.argv.indexOf("--verbose") >= 0;

if (typeof TextDecoder != "function") {
	if (require) {
		global.TextDecoder = require("util").TextDecoder;
	}
}

class FastaFAI {
	/**
	 * @param {string|Partial<FastaFAI>} seq_name
	 * @param {number} len
	 * @param {number} offset
	 */
	constructor(seq_name, len, offset) {
		if (typeof seq_name == "string") {
			/** @type {string} */
			this.seq_name = String(seq_name);

			// TODO: Number.isSafeInteger(len) && Number.isSafeInteger(offset)
			/** @type {number} */

			this.len = Number(len | 0); // to Safe Integer

			/** @type {number} */
			this.offset = Number(offset | 0); // to Safe Integer
		}
		else {
			const oo = seq_name;
			if (oo && typeof oo == "object" && oo.seq_name && oo.len && oo.offset) {
				/** @type {string} */
				this.seq_name = String(seq_name.seq_name ?? "");
				
				/** @type {number} */
				this.len = (seq_name.len ?? 0) | 0;

				/** @type {number} */
				this.offset = (seq_name.offset ?? 0) | 0;
			}
			else {
				throw new TypeError();
			}
		}
	}
}

/**
 * @param {string} fname
 */
function readFastaFAI(fname) {
	const tab = fs.readFileSync(fname).toString().split("\n").map(a => a.split("\t").map(b => b.trim()));
	return Object.fromEntries(tab.map(v => {
		const seq_name = v[0];
		const len = v[1];
		const offset = v[2];
		return [
			seq_name,
			{
				seq_name,
				len,
				offset,
			},
		];
	}))
}

class FastaAttribute {
	/**
	 * @param {Partial<FastaAttribute>} [data]
	 */
	constructor(data = null) {
		/** @type {null|"chromosome"|"mitochondrion"} */
		this.location = undefined;

		/**
		 * I, II, III, IV, V, VI, VII...
		 * @type {string}
		 * @see {@link romanize}
		 */
		this.chromosome = undefined;
		
		/** @type {null|"circular"} */
		this.top = undefined;

		if (data != null) {
 			Object.assign(this, data);
		}
	}

	toString() {
		return this._toString();
	}

	_toString() {
		return Object.entries(this).map(a => {
			return a[1] != null ? `[${a[0]}=${a[1]}]` : null;
		}).filter(a => a != null).join(" ");
	}
}

/**
 * @typedef FastaCollection
 * @type {{ [name: string]: FastaData }}
 */
class FastaData {
	/**
	 * @param {Partial<FastaData>} [data]
	 */
	constructor(data) {
		this.name = "";
		
		this.header = "";
		
		this.seq = "";
		
		Object.assign(this, data);
		
		/** @type {Partial<FastaAttribute>} */
		this.attr = new FastaAttribute();

		if (data != null) {
			for (let key of Object.keys(data.attr)) {
				this.attr[key] = data.attr[key];
			}
		}
	}
}

if (__filename && process.argv[1] && Path.basename(process.argv[1]) == Path.basename(__filename)) {
	const run_cmd = process.argv[2];

	if (run_cmd == "--sum") {// faSize
		const fl = process.argv.slice(3);
		fl.forEach(ff => {
			const fa = readFasta(ff);
			const fa_size = Object.keys(fa).map(k => {
				// console.log(ff, k, fa[k].length);
				return fa[k].length;
			}).reduce((acc, v) => acc + v, 0);
			console.log(ff, fa_size);
		});
	}
	else if (run_cmd == "--get-seq-by-gff") {
		// --get-seq-by-gff --list=list.txt --fasta=QM6a.fa --gff=QM6a.gff3
		cmd_getSeqByGFF();
	}
	else if (run_cmd == "--get-all-isoform") {
		// node fasta_util.js --get-isoform <geneId> [<geneId>.fa]

		const geneId = process.argv[3];
		const genome = readFasta_Ex(process.argv[4]);
		// const fa = genome[geneId];
		Object.values(genome).filter(fa => fa.name.startsWith(geneId) || fa.header.indexOf(geneId) >= 0).forEach(fa => {
			saveFasta(`${fa.name}.fa`, {
				[fa.header]: fa.seq,
			});
		});
	}
	else if (run_cmd == "--getAll") {
		// node fasta_util.js --getAll input.fa

		const genome = readFasta_Ex(process.argv[3]);
		Object.keys(genome).forEach(name => {
			const fa = genome[name];
			saveFasta(`${name}.fa`, {
				[fa.header]: fa.seq,
			});
		});
	}
	else if (run_cmd == "--getByIndex") {
		// node fasta_util.js --getByIndex input.fa <index>

		const genome = readFasta_Ex(process.argv[3]);
		const keys = Object.keys(genome);
		const name = keys[0];
		const fa = genome[name];
		saveFasta(`${name}.fa`, {
			[fa.header]: fa.seq,
		});
	}
	else if (run_cmd == "--getByExcludeList") {
		// node fasta_util.js --getByExcludeList input.fa list.txt output.fa

		const genome = readFasta_Ex(process.argv[3], false);

		// only use column 1
		const _exclude_list = fs.readFileSync(process.argv[4]).toString().trim().split("\n").map(line => line.split("\t")[0]);
		const exclude_set = new Set(_exclude_list);

		const output_file_name = process.argv[5];

		/** @type {{ [seqName: string]: string; }} */
		const output_fa = {
		};

		const keys = Object.keys(genome);

		console.log(_exclude_list[0]);
		console.log(keys[0]);

		keys.forEach(name => {
			if (!exclude_set.has(name)) {
				const fa = genome[name];
				output_fa[fa.header] = fa.seq;
			}
		});

		saveFasta(output_file_name, output_fa);
	}
	else if (run_cmd == "--getByList") {
		(async function () {
			// 0    1             2           3        4        5
			// node fasta_util.js --getByList input.fa list.txt output.fa

			const input_fa = process.argv[3];
			const input_list = process.argv[4];
			const output_file_name = process.argv[5];

			// console.log({
			// 	input_fa,
			// 	input_list,
			// 	output_file_name,
			// });

			const genome = readFasta_Ex(input_fa, false);

			// only use column 1
			const list = await (async function () {
				if (input_list && input_list != "-") {
					return fs.readFileSync(input_list).toString().trim().split("\n").map(line => line.split("\t")[0]);
				}
				else {
					const readline = require("readline");
					const rl = readline.createInterface({
						input: process.stdin,
						crlfDelay: Infinity,
					});
					const lines = [];
					for await (const line of rl) {
						lines.push(line.trim());
					}
					// console.log(lines);
					return lines;
				}
			})();

			/** @type {{ [seqName: string]: string; }} */
			const output_fa = {
			};

			if (VERBOSE) {
				console.log(list[0]);
				console.log(Object.keys(genome)[0]);
			}

			list.forEach(name => {
				const fa = genome[name];
				if (fa) {
					output_fa[fa.header] = fa.seq;
				}
			});

			if (output_file_name && output_file_name != "-") {
				saveFasta(output_file_name, output_fa);
			}
			else {
				const fa_text = _saveFasta_encode(output_fa);
				console.log(fa_text);
			}
		})();
	}
	else if (run_cmd == "--splice" || run_cmd == "--split") {
		// node fasta_util.js --splice input.fa <splice_len>
		// node fasta_util.js --split input.fa <split_len>
		
		const in_fa = readFasta_Ex(process.argv[3], true);
		const keys = Object.keys(in_fa);

		const splice_len = (() => {
			if (run_cmd == "--split") {
				const s_len = Number(process.argv[4]) | 0;
				if (!Number.isNaN(s_len) &&
					Number.isSafeInteger(s_len) &&
					s_len > 0
				) {
					return Math.trunc(keys.length / s_len) + 1;
				}
			}
			else {
				return Number(process.argv[4]) | 0;
			}
		})();

		console.log({
			splice_len,
			"keys.length": keys.length,
		});
		
		const basename = Path.basename(process.argv[3], Path.extname(process.argv[3]));

		if (!Number.isNaN(splice_len) &&
			Number.isSafeInteger(splice_len) &&
			splice_len > 0
		) {
			const pad_len = Math.trunc(Math.log10(Math.ceil(keys.length / splice_len)) + 1);

			for (let splice_id = 1; keys.length; ++splice_id) {
				const list = keys.splice(0, splice_len);

				// console.log(keys.length);
				
				const id = String(splice_id).padStart(pad_len, "0");

				/** @type {{ [seqName: string]: string; }} */
				const out_fa = {};

				list.forEach(name => out_fa[name] = in_fa[name].seq);

				saveFasta(`${basename}.${id}.fa`, out_fa);

				// console.log(`${basename}.${id}.fa`);
			}

			// console.log({ keys });
		}
		// else {
		// 	console.error({
		// 		splice_len
		// 	});
		// }
	}
}

/**
 * ASCII only, uint8
 * @param {FastaCollectionSimple} fasta
 */
function _saveFasta_encode(fasta, line_length = 50) {
	const decoder = new TextDecoder("utf-8");
	const newline = "\n".charCodeAt(0);
	
	// const encoder = new TextEncoder();
	// const newline = encoder.encode("\n");

	let str = "";
	Object.keys(fasta).forEach(key => {
		if (fasta[key] && fasta[key].length) {
			str += ">";
			str += key;
			str += "\n";

			// str += fasta[key];
			const seq = fasta[key];

			const len = seq.length + Math.ceil(seq.length / line_length) - 1;
			const ui8a = new Uint8Array(len);
			for (let i = 0, ii = 0; i < seq.length; i += line_length, ii += (line_length + 1)) {
				for (let j = 0; j < line_length; ++j) {
					ui8a[ii + j] = seq.charCodeAt(i + j);
				}
				ui8a[ii + line_length] = newline;
			}
			str += decoder.decode(ui8a);

			str += "\n";
		}
		else {
			if (VERBOSE) {
				console.warn("saveFasta:", key, "=>", fasta[key]);
			}
		}
	});
	return str;
 }

/**
 * ASCII only, uint8
 * @param {string} filename
 * @param {FastaCollectionSimple} fasta
 */
function saveFasta(filename, fasta, line_length = 50) {
	if (filename == null) {
		throw new Error("filename is null");
	}
	try {
		fs.writeFileSync(filename, _saveFasta_encode(fasta));
	}
	catch (ex) {
		console.error("write file error:", filename);
		throw ex;
	}
}

/**
 * ASCII only, uint8
 * @param {string} filename
 * @param {FastaCollectionSimple} fasta
 */
async function saveFastaAsync(filename, fasta, line_length = 50) {
	if (filename == null) {
		throw new Error("filename");
	}
	try {
		await fs.promises.writeFile(filename, _saveFasta_encode(fasta));
	}
	catch (ex) {
		console.error("write file error:", filename);
		throw ex;
	}
}

/**
 * @param {string} filename
 * @returns {{[chr:string]:string}}
 */
function readFasta(filename) {
	try {
		let text = fs.readFileSync(filename).toString();
		return _parseFasta(text);
	}
	catch (ex) {
		console.error("read file error:", filename, process.cwd());
		throw ex;
	}
}

/**
 * @param {string} in_seq
 * @returns {{[chr:string]:string}}
 */
function _parseFasta(in_seq) {
	let all = in_seq.split("\n").map(a => a.trim()).map(a => {
		if (a.startsWith(">")) {
			return a.split(" ")[0];
		}
		else {
			return a;
		}
	}).join("\n").split(">");
	
	/** @type {{[chr:string]:string}} */
	const results = {};

	if (all.length) {
		all.map(a => a.trim()).filter(a => a.length).map(sub_fa => {
			// const lines = sub_fa.split("\n");
			// const head = lines.shift();
			// const name = head.split(" ")[0];
			// const seq = lines.map(ss => ss.trim().replace(/ /g, "")).join("\n");
			// results[name] = seq;
	
			let li = sub_fa.indexOf("\n");
			try {
				let out_name = sub_fa.slice(0, li).trim().replace(/:/g, "");
				let out_file_name = out_name.match(/^([^ ]+)/)[1];
				let out_seq = sub_fa.slice(li).trim();
	
				out_seq = out_seq.replace(/\n/g, "");
				out_seq = out_seq.replace(/\r/g, "");
				out_seq = out_seq.replace(/ /g, "");
	
				results[out_file_name] = out_seq.toUpperCase();
			}
			catch (ex) {
				console.error(ex);
				console.error("err seq:", sub_fa.slice(0, li));
			}
		});
	}
	else {
		throw new Error("not fasta file");
	}

	return results;
}

/**
 * @param {string} file_path
 * @param {(name: string, header: string) => string} fn_name_transform
 */
function readFA(file_path, fn_name_transform) {
	const text = fs.readFileSync(file_path).toString();
	try {
		return _parseFA(text, fn_name_transform);
	}
	catch (m) {
		const { cc, line, idx } = m;
		const err = "Unexpected char:" + line ? line.charCodeAt(0) : "";
		console.error(`${err}.\t${m}\t${file_path}: ${idx}`);
		console.error("\t" + JSON.stringify(line));
	}
}

/**
 * @param {string} text
 * @param {(name: string, header: string) => string} fn_name_transform
 */
function _parseFA(text, fn_name_transform) {
	/** @type {"OUT_READ" | "IN_SEQ" | "SEQ_END"} */
	let state = "OUT_READ";

	const lines = text.split("\n");
	lines.push("\0");//EOF

	/** @type {string} */
	let head = null;
	/** @type {string[]} */
	let seq_frag_list = [];

	/**
	 * old type {{ [seqID: string]: string; }}
	 * @type {FastaCollection}
	 */
	const map = {
	};

	for (let idx = 0; idx < lines.length;) {
		const line = lines[idx];
		switch (state) {// OUT_READ, IN_NAME, IN_SEQ, SEQ_END, IN_QUAL
			case "OUT_READ":
				switch (line[0]) {
					case ">":
						head = line;
						state = "IN_SEQ";
						++idx;
						break;
					case null:
					case undefined:
					case "":
					case "\n":
						++idx;
						break;
					case "\0":// eof
						++idx;
						break;
					default:
						const err = "Unexpected char:" + line.charCodeAt(0);
						// console.error(`${err}.\t${file_path}: ${idx}`);
						// console.error("\t" + JSON.stringify(line));
						throw {
							error: new Error(err),
							cc: line.charCodeAt(0),
							line: line,
							idx: idx,
						};
				}
				break;
			case "IN_SEQ":
				switch (line[0]) {
					case "\0":// eof
					case ">":
						state = "SEQ_END";
						break;
					case null:
					case undefined:
					case "":
					case "\n":
						++idx;
						break;
					default:
						seq_frag_list.push(line);
						++idx;
				}
				break;
			case "SEQ_END":
				const header = head.slice(1);
				const head_cols = header.split(" ");
				const seqID = fn_name_transform ? fn_name_transform(head_cols[0], header) : head_cols[0];
				map[seqID] = {
					name: seqID,
					header: header,
					seq: seq_frag_list.join("").replace(/\s/g, ""),
					attr: null,
				};
				head = null;
				seq_frag_list = [];
				state = "OUT_READ";
				break;
		}
	}
	return map;
}

/**
 * @param {string} filename
 * @param {boolean} [use_origin_header]
 * @returns {FastaCollection}
 */
function readFasta_Ex(filename, use_origin_header = false) {
	try {
		let text = fs.readFileSync(filename).toString();
		return _parseFasta_Ex(text, use_origin_header, undefined);
	}
	catch (ex) {
		console.error("read file error:", filename);
		throw ex;
	}
}

/**
 * @param {string} in_seq
 * @param {boolean} [use_origin_header]
 * @param {"toUpperCase"|"toLowerCase"} [to_case]
 * @returns {FastaCollection}
 */
function _parseFasta_Ex(in_seq, use_origin_header = undefined, to_case = undefined) {
	let all = in_seq.split("\n").map(a => a.trim()).map(a => {
		// if (a.startsWith(">")) {
		// 	return a.split(" ")[0];
		// }
		// else {
			return a;
		// }
	}).join("\n").split(">");

	/** @type {FastaCollection} */
	let results = {};

	if (all.length) {
		all.filter(a => a.length).map(sub_fa => {
			const li = sub_fa.indexOf("\n");
			const header = sub_fa.slice(0, li).trim();

			const out_name = header.split(" ")[0].replace(/:/g, "");
			const out_file_name = out_name.match(/^([^ ]+)/)[1];
			
			let out_seq = sub_fa.slice(li).trim().replace(/\n|\r| /g, "");
			
			if (to_case == "toUpperCase") {
				out_seq = out_seq.toUpperCase();
			}
			else if (to_case == "toLowerCase") {
				out_seq = out_seq.toLowerCase();
			}
			
			results[use_origin_header ? header : out_file_name] = new FastaData({
				name: out_file_name,
				header: header,
				seq: out_seq,
				attr: new FastaAttribute({
					location: getFastaAttrVal(header, "location"),   // chromosome, mitochondrion
					chromosome: getFastaAttrVal(header, "chromosome"), // I, II, III...
					top: getFastaAttrVal(header, "top"),        // circular
				}),
			});
		});
	}
	else {
		throw new Error("not fasta file");
	}

	return results;
}

/**
 * @template T
 * @param {string} header
 * @param {string} attr
 * @returns {T}
 * @example
getFastaAttrVal("[location=mitochondrion] [chromosome=III] [top=circular]", "top")
getFastaAttrVal("[location=mitochondrion] [chromosome=III] [top=circular]", "chromosome")
getFastaAttrVal("[location=mitochondrion] [chromosome=III] [top=circular]", "location")
 */
function getFastaAttrVal(header, attr) {
	const m = header.match(new RegExp(`\\[${attr}=([^\\]]+)\\]`));
	if (m) {
		// console.log(header, attr, m[1]);
		/** @type {any} */
		const v = m[1];
		return v;
	}
}

function save_fasta_async(path, fasta) {
	return new Promise(function (resolve, reject) {
		try {
			let fp = fs.createWriteStream(path);
			Object.keys(fasta).forEach(key => {
				fp.write("> ");
				fp.write(key);
				fp.write("\n");

				fp.write(fasta[key]);
				fp.write("\n");
			});
			fp.end();
			fp.on("close", resolve);
		}
		catch (ex) {
			reject(ex);
		}
	});
}


/**
 * @param {{[chr:string]:string}[]} fasta_list
 * @returns {{[chr:string]:string}}
 */
function joinFastaSeq(fasta_list) {
	/** @type {{[chr:string]:string}} */
	let out = {};
	Object.keys(fasta_list[0]).forEach(k => {
		if (!out[k]) {
			out[k] = "";
		}
		fasta_list.forEach(seq => {
			out[k] += seq[k];
		});
	});
	return out;
}


/**
 * index: 0 ~ length - 1
 * @param {string} ma_seq
 * @returns {number[]}
 */
function chrPos_to_multialign_posMap_old(ma_seq) {
	let nPos = 1;
	let posmap = [];
	for (let index = 0; index < ma_seq.length; index++) {
		const element = ma_seq[index];
		if (element != "-") {
			posmap[nPos] = index + 1;
			++nPos;
		}
	}
	return posmap;
}

/**
 * index: 0 ~ length - 1
 * @param {string} ma_seq
 * @returns {Uint32Array}
 * @version 20220627 webView
 */
function chrPos_to_multialign_posMap(ma_seq) {
	let nPos = 1;
	let posmap = new Uint32Array(ma_seq.length);
	// for (let index = 0; index < ma_seq.length; ++index) {
	// 	const element = ma_seq[index];
	// 	// if (element != "-") {
	// 	// 	posmap[nPos] = index + 1;
	// 	// 	++nPos;
	// 	// }
	// }
	// 20220627
	for (let pos = 0; pos < ma_seq.length; ++pos) {
		const element = ma_seq[pos];
		(function (ref_pos) {
			if (!posmap[ref_pos]) {
				posmap[ref_pos] = pos + 1;
			}
		})(Math.max(1, element != "-" ? nPos : (nPos - 1)));
		
		if (element != "-") {
			++nPos;
		}
		// if (nPos == 10) {
		// 	console.log(posmap.slice(0, 10));
		// 	break;
		// }
	}
	return posmap;
}

/**
 * index: 0 ~ length - 1
 * @param {string} ma_seq
 * @returns {number[]}
 */
function multialign_to_chrPos_posMap_old(ma_seq) {
	let nPos = 0;
	let posmap = [];

	posmap[0] = 0;

	for (let index = 0; index < ma_seq.length; index++) {
		const element = ma_seq[index];
		
		if (element != "-") {
			++nPos;
		}
		
		posmap[index + 1] = nPos;
	}

	return posmap;
}

/**
 * index: 0 ~ length - 1
 * posmap[index + 1] = ++nPos;
 * @param {string} ma_seq
 * @returns {Uint32Array}
 * @version 20220627 webView
 */
function multialign_to_chrPos_posMap(ma_seq) {
	let nPos = 1;
	let posmap = new Uint32Array(ma_seq.length);

	posmap[0] = 0;

	for (let pos = 0; pos < ma_seq.length; ++pos) {
		const element = ma_seq[pos];

		// if (element != "-") {
		// 	++nPos;
		// }
		// // posmap[pos + 1] = nPos;
		// posmap[pos] = nPos - 1;// idx to idx

		(function (ref_pos) {
			posmap[pos] = ref_pos;
		})(Math.max(1, element != "-" ? nPos : (nPos - 1)));

		if (element != "-") {
			++nPos;
		}
	}
	
	// const a = posmap;
	// const b = pos_ref1_uint32array;
	// for (let i = 0; i < a.length; ++i) {
	// 	if (a[i] != b[i]) {
	// 		console.log(i, a[i], b[i]);
	// 		break;
	// 	}
	// }

	return posmap;
}

const complement_map = {};
complement_map["A"] = "T";
complement_map["T"] = "A";
complement_map["C"] = "G";
complement_map["G"] = "C";
complement_map["N"] = "N";
complement_map["a"] = "t";
complement_map["t"] = "a";
complement_map["c"] = "g";
complement_map["g"] = "c";
complement_map["n"] = "n";

/**
 * @param {string} seq
 */
function reverseComplement(seq) {
	/** @type {Map<string, number>} */
	const hasUnknow = new Map();
	
	const result = [...seq].reverse().map(a => {
		if (complement_map[a] == null || complement_map[a] == "N" || complement_map[a] == "n") {
			hasUnknow.set(a, (hasUnknow.get(a) || 0) + 1);
		}
		return complement_map[a];
	}).join("");

	if (result.length != seq.length || hasUnknow.size) {
		for (let ent of hasUnknow.entries()) {
			console.warn(`found ${ent[0]}: ${ent[1]}, ${(ent[1] * 100 / seq.length).toFixed(2)}%`);
		}
	}

	return result;
}

/**
 * @param {string} seq
 */
function reverseComplementEx(seq) {
	return [...seq].reverse().map(a => complement_map[a] || a).join("");
}

/**
 * @param {FastaCollectionSimple} fasta
 * @param {string[]} ordered_list
 * @returns {FastaCollectionSimple}
 */
function sortFasta(fasta, ordered_list) {
	// @ts-ignore
	return Object.fromEntries(ordered_list.map(k => [k, fasta[k]]));
}

async function cmd_getSeqByGFF() {
	const { ArgvParser } = require("./tools/ArgvParser.js");

	const arg_p = new ArgvParser(process.argv);

	const arg_ls = arg_p.get(/--list=(.+)/, (arg, args) => args[0]);
	const path_to_fasta = arg_p.get(/--fasta=(.+)/, (arg, args) => args[0]);
	const path_to_gff = arg_p.get(/--gff=(.+)/, (arg, args) => args[0]);
	const path_to_output = arg_p.get(/--output=(.+)/, (arg, args) => args[0]);
	const path_to_options = arg_p.get(/--options=(.+)/, (arg, args) => args[0]);

	if (path_to_fasta && path_to_gff) {
		const ls = await (async () => {
			if (arg_ls && arg_ls != "-") {
				return arg_ls.split(",");
			}
			else {
				const readline = require("readline");
				const rl = readline.createInterface({
					input: process.stdin,
					crlfDelay: Infinity,
				});
				const lines = [];
				for await (const line of rl) {
					lines.push(line.trim());
				}
				// console.log(lines);
				return lines;
			}
		})();

		const get_seq_opt = path_to_options ? require(path_to_options) : null;
		const text = await _getSeqByGFF(ls, path_to_fasta, path_to_gff, get_seq_opt);

		if (path_to_output) {
			fs.writeFileSync(path_to_output, text);
		}
		else {
			console.log(text);
		}
	}
}

/**
 * @param {string[]} ls
 * @param {string} path_to_fasta
 * @param {string} path_to_gff
 */
async function _getSeqByGFF(ls, path_to_fasta, path_to_gff, get_seq_opt) {
	const { GFF_ROW, parseGFF } = require("./gff.js");

	const fa = readFasta(path_to_fasta);
	const gff = Object.values(parseGFF(fs.readFileSync(path_to_gff).toString())).flat(1);

	// feature
	// auto_RC
	// fwd: s-500, s+100
	// rev: e-100, e+500

	let vv = await Promise.all(ls.map(id => getGeneSeqById(id, get_seq_opt)));
	return vv.filter(arr => arr).map(arr => {
		if (arr) {
			if (arr[0].gene.strand < 0) {
				return [...arr].sort((a, b) => b.gene.end - a.gene.end)[0];
			}
			else if (arr[0].gene.strand > 0) {
				return arr[0];
			}
			else {
				throw new Error("arr[0].strand");
			}
		}
	}).map(a => a?.fa_header + "\n" + a?.seq).join("\n");

	/**
	 * @param {string} gene_id
	 * @param {null|{ feature?: string, auto_RC?: boolean; ext_func?: (start: number, end: number, is_rev: boolean) => { start: number; end: number; relative_5_prime: number; relative_3_prime: number; }; }} [__options]
	 */
	async function getGeneSeqById(gene_id, __options = {}) {
		/**
		 * @type {{ feature?: string, auto_RC?: boolean; ext_func?: (start: number, end: number, is_rev: boolean) => { start: number; end: number; relative_5_prime: number; relative_3_prime: number; }; }}
		 */
		const options = {
			feature: "gene",// mRNA
			...(__options ?? {}),
		};

		const feature = options.feature ?? "gene";// mRNA

		const predicate = row => row.type == feature && row.attributes?.ID?.startsWith(gene_id);
		const gg = gff.filter(predicate);

		if (gg) {
			// const b_out_fa = options.b_out_fa;
			
			// const node_parent = await findGeneByGffAll(refName, nChr, row => row.type == "gene" && row.attributes?.ID?.startsWith(gene_id));
			// const node_parent_strand = node_parent[0].strand;
			// return (() => {
			// 	if (node_parent_strand > 0) {
			// 		return gg;
			// 	}
			// 	else {
			// 		return [...gg].sort((a, b) => b.start - a.start);
			// 	}
			// })().map(gene => {
			// 	const has_rc = options?.auto_RC && (gene.strand ?? node_parent_strand) < 0;
			return gg.map(gene => {
				const chr_seq = fa[gene.seqid];
				const has_rc = !!options?.auto_RC && gene.strand < 0;
			
				// gff_data_map['QM6a']['ChVII_QM6a'].find(a => a.type == "gene" && a.attributes.Name)
				
				// const ext_ATG_fwd = options.ext_ATG_fwd ?? 0;
				// const ext_ATG_rev = options.ext_ATG_rev ?? 0;
				// const pos_start = gene.start;// has_rc ? (gene.start - ext_ATG_fwd) : (gene.start - ext_ATG_rev);
				// const pos_end = gene.end;// has_rc ? (gene.end + ext_ATG_rev) : (gene.end + ext_ATG_fwd);

				// const pos_start = has_rc ? (gene.start - 100) : (gene.start - 500);
				// const pos_end = has_rc ? (gene.start + 500) : (gene.start + 100);

				const {
					start: pos_start,
					end: pos_end,
					relative_5_prime,
					relative_3_prime,
				} = (() => {
					if (options.ext_func) {
						return options.ext_func(gene.start, gene.end, has_rc);
					}
					else {
						return {
							relative_5_prime: null,
							relative_3_prime: null,
							start: gene.start,
							end: gene.end,
						};
					}
				})();

				const seq = (s => has_rc ? reverseComplement(s) : s)(chr_seq.slice(pos_start - 1, pos_end));
				const strand = gene.strand > 0 ? "+" : (gene.strand < 0 ? "-" : "");
				const fa_header = ">" + [
					gene.attributes.ID,
					gene.attributes.Name,
					[
						gene.geneID,
						`${gene.start}-${gene.end}`,
						strand
					].join(":"),
					has_rc ? "[reverse] [complement]" : null,
					
					`${gene.seqid}_${pos_start}_${pos_end}`,

					`[relative 5' ${relative_5_prime}]-[relative 3' ${relative_3_prime}`,
				].filter(v => v).join(" ");

				return {
					gene,
					is_gene_seq: pos_start == gene.start && pos_end == gene.end ? true : false,
					pos_start,
					pos_end,
					has_reverse: has_rc,
					has_complement: has_rc,
					fa_header,
					seq,
				};
			});// for ... of
		}
	}
}

module.exports.FastaFAI = FastaFAI;
module.exports.readFastaFAI = readFastaFAI;

module.exports.readFA = readFA;
module.exports.readFasta = readFasta;
module.exports.readFasta_Ex = readFasta_Ex;
module.exports._parseFasta = _parseFasta;
module.exports._parseFA = _parseFA;

module.exports.saveFasta = saveFasta;
module.exports.saveFastaAsync = saveFastaAsync;
module.exports._saveFasta_encode = _saveFasta_encode;
module.exports.sortFasta = sortFasta;
module.exports.reverseComplement = reverseComplement;
module.exports.reverseComplementEx = reverseComplementEx;

module.exports.joinFastaSeq = joinFastaSeq;

module.exports.chrPos_to_multialign_posMap = chrPos_to_multialign_posMap;
module.exports.multialign_to_chrPos_posMap = multialign_to_chrPos_posMap;


/**
 * @typedef FastaCollectionSimple
 * @type {{ [name: string]: string }}
 */
