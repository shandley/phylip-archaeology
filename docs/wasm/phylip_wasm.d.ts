/* tslint:disable */
/* eslint-disable */

/**
 * Bootstrap NJ analysis with majority-rule consensus.
 *
 * `nreps`: number of bootstrap replicates
 * `seed`: random seed (0 = use default 42)
 * Returns JSON: { original_newick, consensus_newick, nreps, split_support }
 */
export function bootstrap_nj(fasta: string, nreps: number, seed: number): string;

/**
 * Evaluate tree likelihood for an alignment under a substitution model.
 *
 * `model`: "jc69" or "f84"
 * Returns JSON: { lnl, ntaxa, nsites, model, newick }
 */
export function compute_likelihood(fasta: string, newick_str: string, model: string): string;

/**
 * Compute a distance matrix and Neighbor-Joining tree from FASTA input.
 *
 * `model`: "jc69" or "k2p"
 * Returns JSON: { names, matrix, newick, ntaxa, nsites, model }
 */
export function compute_nj(fasta: string, model: string): string;

/**
 * Run maximum parsimony tree search on FASTA input.
 *
 * `seed`: random seed (0 = use default seed 42)
 * Returns JSON: { score, newick, ntaxa, nsites }
 */
export function compute_parsimony(fasta: string, seed: number): string;

/**
 * Simulate the Felsenstein Zone: demonstrate long-branch attraction.
 *
 * Parameters: branch lengths, number of sites, replicates, and seed.
 * Returns JSON with per-replicate results and accuracy summary.
 */
export function felsenstein_zone(long_branch: number, short_branch: number, internal: number, nsites: number, nreps: number, seed: number): string;

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly bootstrap_nj: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly compute_likelihood: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number, number, number];
    readonly compute_nj: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly compute_parsimony: (a: number, b: number, c: number) => [number, number, number, number];
    readonly felsenstein_zone: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number, number, number];
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
    readonly __externref_table_dealloc: (a: number) => void;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
