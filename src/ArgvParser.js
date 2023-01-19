// @ts-check

class ArgvParser {
	/**
	 * @private
	 * @type {string[]}
	 */
	_argv = [];
	
	/**
	 * @param {string[]} [argv]
	 */
	constructor(argv = process.argv) {
		this._argv = [...argv];
	}

	/**
	 * @template T
	 * @param {RegExp} regexp
	 * @param {(arg: string, args: string[]) => T|null|undefined} trx
	 * @example
	 * get(/--aaa=(\d+)/, (arg, args) => Number(args[0]))
	 */
	get(regexp, trx) {
		const aa = this._argv.map(s => s.match(regexp)).filter(a => a)[0];
		if (aa) {
			const i = this._argv.indexOf(aa[0]);
			if (i >= 0) {
				const args = this._argv[i];

				const r = trx?.(aa[0], aa.slice(1));

				this._argv.splice(i, 1);

				return r;
			}
		}
	}
}

module.exports.ArgvParser = ArgvParser;

