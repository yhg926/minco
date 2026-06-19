#ifndef SKETCH_WRAPPER_H
#define SKETCH_WRAPPER_H

#include "global_basic.h"
#include "pairwise_graph.h"
#include <stdbool.h>
#include <argp.h>

typedef pairwise_metric_expr_t sketch_dedup_metric_t;

#ifndef MINCO_DEFAULT_SKETCH_SIZE
#ifdef MINCO_SKETCH_SIZE
#define MINCO_DEFAULT_SKETCH_SIZE ((uint32_t)(MINCO_SKETCH_SIZE))
#else
#define MINCO_DEFAULT_SKETCH_SIZE 10000u
#endif
#endif

typedef enum minco_ctxmeta_mode
{
	MINCO_CTXMETA_NONE = 0,
	MINCO_CTXMETA_PRECONFLICT = 1,
	MINCO_CTXMETA_POSTCONFLICT = 2,
	MINCO_CTXMETA_BOTH = 3
} minco_ctxmeta_mode_t;

typedef struct sketch_opt
{
	// sketch stats
	int hclen;	// half context length, 1..16
	int holen;	// half outer object length, 0..8
	int iolen;	// half outer object length, 0..8
		int compat_filter_shift; // hidden source hash prefilter shift; new stat files store 0
	uint32_t sketch_size; // final bottom-k context hashes per sample
	int kmerocrs;
	double npercentile; // dynamic k-mer occurrence threshold percentile, 0 disables
	int ncap; // maximum effective lower k-mer occurrence threshold, 0 disables
	bool reads_qc; // infer a read k-mer count range from the frequency distribution
	int p; // threads counts
	bool abundance;
	bool asone;	// treat input genomes as parts of final genome.
	bool conflict; // keep conflict context-object or not 
	bool anno; // write FASTA/FASTQ header annotations to minco.anno
	bool compute_meta; // write per-input metadata to minco.infilemeta
	minco_ctxmeta_mode_t ctxmeta_mode; // write minco context-cardinality metadata
	bool density_threshold_enabled; // internal: keep all hashed contexts up to density_threshold
	uint64_t density_threshold; // internal: reference-derived context hash threshold
	uint32_t qc_hash_mult; // readsQC hash-sampling multiplier relative to final sketch size
	uint32_t qc_hash_target; // readsQC hash-sampling target; 0 uses qc_hash_mult
	bool position; // write per-context-object sequence positions to minco.ctxobj64.position
	bool drop_position; // drop minco.ctxobj64.position during filtering maintenance modes
	bool merge_minco_mode;
	bool append_minco_mode;
	bool remove_minco_mode;
	bool keep_minco_mode;
	bool dedup_minco_mode;
	bool dedup_raw_build_from_inputs;
	bool append_copy_mode;
	bool remove_copy_mode;
	bool keep_copy_mode;
	bool dedup_copy_mode;
	bool sketch_qc;
	double dedup_cutoff;
	double dedup_max_afcut;
	uint32_t dedup_ctxcut;
	sketch_dedup_metric_t dedup_metric;
	pairwise_dedup_strategy_t dedup_strategy;
	bool dedup_index;
	uint32_t dedup_index_max_ctx_freq;
	uint32_t dedup_index_min_votes;
	uint32_t dedup_index_sample_step;
	int print_mode; // 1 samples, 2 sketch entries, 3 sorted index, 4 positions, 5 ctxmeta, 6 ctxsetmeta
	bool split_mfa;
	bool coden_ctxobj_pattern; 	 
	char index[PATHLEN];
	char *remove_list;
	char *remove_source;
	char *keep_list;
	char *keep_source;
	char *dedup_source;
	char *fpath;
	char *outdir;  // results dir
	char *pipecmd; // pipe command
	int num_remaining_args;
	char **remaining_args;
} sketch_opt_t;

int cmd_sketch(struct argp_state *state);
/*
extern void compute_sketch(sketch_opt_t*, infile_tab_t*);
extern void gen_inverted_index_for_minco(const char* sketchdir);
extern int merge_minco_sketches (sketch_opt_t * sketch_opt_val);
extern uint32_t get_sketching_id(uint32_t hclen, uint32_t holen,uint32_t iolen,uint32_t compat_filter_shift,uint32_t FILTER);
*/
#endif
