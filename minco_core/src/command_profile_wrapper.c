#include "command_profile_wrapper.h"
#include "command_ani.h"
#include "global_basic.h"

#include <argp.h>
#include <ctype.h>
#include <err.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

typedef struct profile_args
{
	char ref[PATHLEN];
	char reads[PATHLEN];
	char outf[PATHLEN];
	char cami_profile[PATHLEN];
	char cami_taxmap[PATHLEN];
	char cami_sample_id[PATHLEN];
	char track[PATHLEN];
	char track_summary[PATHLEN];
	char edge_out[PATHLEN];
	char gtdb_taxmap[PATHLEN];
	char ncbi_taxmap[PATHLEN];
	char pipecmd[PATHLEN];
	int p;
	int ntop;
	uint64_t edge_max;
	uint32_t edge_max_candidates;
	uint32_t density_block_ctx;
	bool density_block_set;
	bool report_all;
	bool dual_evidence;
	bool edge_ambiguous_only;
	bool edge_selected_only;
	ani_readwise_assign_mode_t assign_mode;
	ani_readwise_ani_model_t ani_model;
	ani_readwise_ctx_filter_model_t ctx_filter_model;
	ani_readwise_taxonomy_mode_t taxonomy_mode;
	double fake_threshold;
} profile_args_t;

enum
{
	PROFILE_CAMI_PROFILE = 1100,
	PROFILE_CAMI_TAXMAP,
	PROFILE_CAMI_SAMPLE_ID,
	PROFILE_TRACK,
	PROFILE_TRACK_SUMMARY,
	PROFILE_TAXONOMY,
	PROFILE_EDGE_OUT,
	PROFILE_EDGE_MAX,
	PROFILE_EDGE_MAX_CANDIDATES,
	PROFILE_EDGE_ALL,
	PROFILE_EDGE_SELECTED_ONLY,
	PROFILE_GTDB_TAXMAP,
	PROFILE_NCBI_TAXMAP,
	PROFILE_PIPECMD,
	PROFILE_DENSITY_BLOCK_CTX,
	PROFILE_REPORT_ALL,
	PROFILE_DUAL_EVIDENCE,
	PROFILE_ANI_MODEL,
	PROFILE_ASSIGN,
	PROFILE_CTX_FILTER,
	PROFILE_FAKE_THRESHOLD,
	PROFILE_PRESET
};

static struct argp_option opt_profile[] =
{
	{0, 0, 0, 0, "Inputs:", 1},
	{"ref", 'r', "<SKETCH>", 0, "Reference sketch directory.", 1},
	{"reads", 'q', "<FASTA|FASTQ[.gz]|->", 0, "Raw reads/query sequence input.", 1},
	{"qraw", 'x', "<FASTA|FASTQ[.gz]|->", 0, "Alias for --reads.", 1},

	{0, 0, 0, 0, "Output:", 2},
	{"outfile", 'o', "<FILE>", 0, "Write the profile table. [stdout]", 2},
	{"top", 'N', "<INT>", 0, "Report at most top N references.", 2},
	{"cami-profile", PROFILE_CAMI_PROFILE, "<FILE>", 0, "Write a CAMI taxonomic profile from called rows.", 2},
	{"cami-taxmap", PROFILE_CAMI_TAXMAP, "<TSV>", 0, "CAMI taxonomy map for --cami-profile.", 2},
	{"cami-sample-id", PROFILE_CAMI_SAMPLE_ID, "<ID>", 0, "Sample ID for --cami-profile.", 2},

	{0, 0, 0, 0, "Read tracking:", 3},
	{"track", PROFILE_TRACK, "<FILE>", 0, "Write Kraken-like per-read hit tracking TSV.", 3},
	{"track-summary", PROFILE_TRACK_SUMMARY, "<FILE>", 0, "Write read tracking summary TSV.", 3},
	{"taxonomy", PROFILE_TAXONOMY, "<none|gtdb|ncbi|both>", 0, "Taxonomy namespace for --track LCA labels. [none]", 3},
	{"edge-out", PROFILE_EDGE_OUT, "<FILE>", 0, "Experimental: write context-level readwise ambiguity edges for EM/debugging. Forces exact per-read density units.", 3},
	{"edge-max", PROFILE_EDGE_MAX, "<INT>", 0, "Maximum ambiguity edge rows to write; 0 disables the cap. [10000000]", 3},
	{"edge-max-candidates", PROFILE_EDGE_MAX_CANDIDATES, "<INT>", 0, "Skip context groups with more candidate refs than this; 0 disables the cap. [64]", 3},
	{"edge-all", PROFILE_EDGE_ALL, 0, 0, "With --edge-out, include unambiguous context groups too. Default writes ambiguous groups only.", 3},
	{"edge-selected-only", PROFILE_EDGE_SELECTED_ONLY, 0, 0, "With --edge-out, write only candidate edges selected by the assignment mode.", 3},
	{"gtdb-taxmap", PROFILE_GTDB_TAXMAP, "<TSV>", 0, "GTDB ref taxonomy map for --track.", 3},
	{"ncbi-taxmap", PROFILE_NCBI_TAXMAP, "<TSV>", 0, "NCBI ref taxonomy map for --track.", 3},

	{0, 0, 0, 0, "Execution and tuning:", 4},
	{"threads", 'p', "<INT>", 0, "Number of threads. [1]", 4},
	{"pipecmd", PROFILE_PIPECMD, "<CMD>", 0, "Stream raw input through CMD; '{}' is replaced by the input path.", 4},
	{"density-block-ctx", PROFILE_DENSITY_BLOCK_CTX, "<INT>", 0, "Read block size in retained density contexts; 0 keeps exact per-read mode. [" MINCO_DEFAULT_DENSITY_BLOCK_CTX_STR "]", 4},
	{"preset", PROFILE_PRESET, "<universal|debug-all>", 0, "Profile preset. [universal]", 4},
	{"report-all", PROFILE_REPORT_ALL, 0, 0, "Write every ref comparison instead of default called rows.", 4},
	{"dual-evidence", PROFILE_DUAL_EVIDENCE, 0, 0, "Append Marker_* columns from contexts unique to one reference in a full index.", 4},
	{"ani-model", PROFILE_ANI_MODEL, "<naive|zip-aaf>", 0, "Selected readwise ANI model. [naive]", 4},
	{"assign", PROFILE_ASSIGN, "<best-diff-split|best-diff-unique|best-diff|all>", 0, "Readwise shared-context assignment. [best-diff-split]", 4},
	{"ctx-filter", PROFILE_CTX_FILTER, "<none|product-topfrac-median|product1-topfrac-median>", 0, "Context defake/reliability filter. [product-topfrac-median]", 4},
	{"fake-threshold", PROFILE_FAKE_THRESHOLD, "<FLOAT>", 0, "Threshold for --ctx-filter; 0.25 means top 25 percent for topfrac modes. [0.25]", 4},
	{0}
};

static char doc_profile[] =
	"\n"
	"Profile raw reads against a minco reference sketch with conservative direct readwise depth defaults."
	"\v"
	"`profile` is the conservative direct front end for metagenome/gene-panel profiling.\n"
	"It runs the direct readwise path with reference-density extraction,\n"
	"depth abundance, profile-only memory bounds, best-diff-split assignment,\n"
	"naive readwise ANI, and product-topfrac-median context defaking.\n"
	"It is not the calibrated universal-auto-exact species wrapper; use\n"
	"`scripts/minco_profile_default.py` for the current F1-priority species\n"
	"default candidate when species taxmap and training features are available.\n"
	"\n"
	"Examples:\n"
	"  minco profile -r gtdb.s2000.ctxmarker.minco reads.fq.gz -o profile.tsv\n"
	"  minco profile -p16 -r ref.minco --reads sample.fq.gz \\\n"
	"    --cami-taxmap ref.taxmap.tsv --cami-profile sample.profile -o profile.tsv\n"
	"  minco profile -r ref.minco reads.fq.gz --track reads.track.tsv \\\n"
	"    --taxonomy gtdb --gtdb-taxmap ref.gtdb.tsv -o profile.tsv\n"
	"\n"
	"Use `ani` directly for assembled-genome ANI, matrices, and experimental\n"
	"benchmark sweeps that need every low-level flag.";

static int profile_parse_int_range(struct argp_state *state, const char *option_name,
								   const char *arg, int min_value, int max_value)
{
	char *end = NULL;
	errno = 0;
	long value = strtol(arg, &end, 10);
	if (errno != 0 || end == arg || *end != '\0' ||
		value < min_value || value > max_value)
		argp_error(state, "%s requires an integer in range %d..%d",
				   option_name, min_value, max_value);
	return (int)value;
}

static uint64_t profile_parse_uint64(struct argp_state *state,
									 const char *option_name,
									 const char *arg,
									 uint64_t max_value)
{
	char *end = NULL;
	errno = 0;
	unsigned long long value = strtoull(arg, &end, 10);
	if (errno != 0 || end == arg || *end != '\0' ||
		(uint64_t)value > max_value)
		argp_error(state, "%s requires an integer in range 0..%llu",
				   option_name, (unsigned long long)max_value);
	return (uint64_t)value;
}

static double profile_parse_nonnegative_double(struct argp_state *state,
											   const char *option_name,
											   const char *arg)
{
	char *end = NULL;
	errno = 0;
	double value = strtod(arg, &end);
	if (errno != 0 || end == arg || *end != '\0' || !isfinite(value) ||
		value < 0.0)
		argp_error(state, "%s requires a non-negative number", option_name);
	return value;
}

static void profile_copy_path(struct argp_state *state, const char *option_name,
							  char *dest, size_t dest_size, const char *arg)
{
	if (strlen(arg) >= dest_size)
		argp_error(state, "%s path is too long; maximum supported length is %zu bytes",
				   option_name, dest_size - 1);
	snprintf(dest, dest_size, "%s", arg);
}

static ani_readwise_taxonomy_mode_t profile_parse_taxonomy(struct argp_state *state,
														  const char *arg)
{
	if (strcasecmp(arg, "none") == 0 || strcasecmp(arg, "off") == 0 ||
		strcmp(arg, "0") == 0)
		return ANI_READWISE_TAXONOMY_NONE;
	if (strcasecmp(arg, "gtdb") == 0 || strcmp(arg, "1") == 0)
		return ANI_READWISE_TAXONOMY_GTDB;
	if (strcasecmp(arg, "ncbi") == 0 || strcmp(arg, "2") == 0)
		return ANI_READWISE_TAXONOMY_NCBI;
	if (strcasecmp(arg, "both") == 0 || strcasecmp(arg, "all") == 0 ||
		strcmp(arg, "3") == 0)
		return ANI_READWISE_TAXONOMY_BOTH;
	argp_error(state, "--taxonomy must be one of none, gtdb, ncbi, or both");
	return ANI_READWISE_TAXONOMY_NONE;
}

static ani_readwise_ani_model_t profile_parse_ani_model(struct argp_state *state,
														const char *arg)
{
	if (strcasecmp(arg, "naive") == 0 || strcasecmp(arg, "raw") == 0 ||
		strcmp(arg, "0") == 0)
		return ANI_READWISE_ANI_NAIVE;
	if (strcasecmp(arg, "zip-aaf") == 0 || strcasecmp(arg, "zip") == 0 ||
		strcmp(arg, "1") == 0)
		return ANI_READWISE_ANI_ZIP_AAF;
	argp_error(state, "--ani-model must be one of naive or zip-aaf");
	return ANI_READWISE_ANI_NAIVE;
}

static ani_readwise_assign_mode_t profile_parse_assign(struct argp_state *state,
													   const char *arg)
{
	if (strcasecmp(arg, "best-diff-split") == 0 ||
		strcasecmp(arg, "best-split") == 0 ||
		strcasecmp(arg, "split") == 0 || strcmp(arg, "2") == 0)
		return ANI_READWISE_ASSIGN_BEST_DIFF_SPLIT;
	if (strcasecmp(arg, "best-diff-unique") == 0 ||
		strcasecmp(arg, "best-unique") == 0 ||
		strcasecmp(arg, "unique") == 0 || strcmp(arg, "3") == 0)
		return ANI_READWISE_ASSIGN_BEST_DIFF_UNIQUE;
	if (strcasecmp(arg, "best-diff") == 0 || strcasecmp(arg, "best") == 0 ||
		strcmp(arg, "1") == 0)
		return ANI_READWISE_ASSIGN_BEST_DIFF;
	if (strcasecmp(arg, "all") == 0 || strcasecmp(arg, "none") == 0 ||
		strcmp(arg, "0") == 0)
		return ANI_READWISE_ASSIGN_ALL;
	argp_error(state, "--assign must be one of best-diff-split, best-diff-unique, best-diff, or all");
	return ANI_READWISE_ASSIGN_BEST_DIFF_SPLIT;
}

static ani_readwise_ctx_filter_model_t profile_parse_ctx_filter(struct argp_state *state,
																const char *arg)
{
	if (strcasecmp(arg, "none") == 0 || strcasecmp(arg, "off") == 0 ||
		strcmp(arg, "0") == 0)
		return ANI_READWISE_CTX_FILTER_NONE;
	if (strcasecmp(arg, "product-topfrac-median") == 0 ||
		strcasecmp(arg, "topfrac-median") == 0 ||
		strcasecmp(arg, "product-median") == 0 || strcmp(arg, "7") == 0)
		return ANI_READWISE_CTX_FILTER_PRODUCT_TOPFRAC_MEDIAN;
	if (strcasecmp(arg, "product1-topfrac-median") == 0 ||
		strcasecmp(arg, "product-plus1-topfrac-median") == 0 ||
		strcasecmp(arg, "product1-median") == 0 || strcmp(arg, "8") == 0)
		return ANI_READWISE_CTX_FILTER_PRODUCT1_TOPFRAC_MEDIAN;
	if (strcasecmp(arg, "poisson-product") == 0 ||
		strcasecmp(arg, "product") == 0 || strcmp(arg, "4") == 0)
		return ANI_READWISE_CTX_FILTER_POISSON_PRODUCT;
	if (strcasecmp(arg, "product-nb") == 0 ||
		strcasecmp(arg, "nb-product") == 0 || strcmp(arg, "5") == 0)
		return ANI_READWISE_CTX_FILTER_PRODUCT_NB;
	argp_error(state, "--ctx-filter must be one of none, product-topfrac-median, product1-topfrac-median, poisson-product, or product-nb");
	return ANI_READWISE_CTX_FILTER_PRODUCT_TOPFRAC_MEDIAN;
}

static void profile_apply_preset(struct argp_state *state, profile_args_t *args,
								 const char *arg)
{
	if (strcasecmp(arg, "universal") == 0 || strcasecmp(arg, "default") == 0)
		return;
	if (strcasecmp(arg, "debug-all") == 0 || strcasecmp(arg, "benchmark") == 0)
	{
		args->report_all = true;
		return;
	}
	argp_error(state, "--preset must be one of universal or debug-all");
}

static error_t parse_profile(int key, char *arg, struct argp_state *state)
{
	profile_args_t *args = state->input;

	switch (key)
	{
	case 'r':
		profile_copy_path(state, "-r/--ref", args->ref, sizeof(args->ref), arg);
		break;
	case 'q':
	case 'x':
		profile_copy_path(state, "--reads", args->reads, sizeof(args->reads), arg);
		break;
	case 'o':
		profile_copy_path(state, "-o/--outfile", args->outf, sizeof(args->outf), arg);
		break;
	case 'p':
		args->p = profile_parse_int_range(state, "-p/--threads", arg, 1, 65536);
		break;
	case 'N':
		args->ntop = profile_parse_int_range(state, "-N/--top", arg, 1, INT_MAX);
		break;
	case PROFILE_CAMI_PROFILE:
		profile_copy_path(state, "--cami-profile", args->cami_profile,
						  sizeof(args->cami_profile), arg);
		break;
	case PROFILE_CAMI_TAXMAP:
		profile_copy_path(state, "--cami-taxmap", args->cami_taxmap,
						  sizeof(args->cami_taxmap), arg);
		break;
	case PROFILE_CAMI_SAMPLE_ID:
		profile_copy_path(state, "--cami-sample-id", args->cami_sample_id,
						  sizeof(args->cami_sample_id), arg);
		break;
	case PROFILE_TRACK:
		profile_copy_path(state, "--track", args->track, sizeof(args->track), arg);
		break;
	case PROFILE_TRACK_SUMMARY:
		profile_copy_path(state, "--track-summary", args->track_summary,
						  sizeof(args->track_summary), arg);
		break;
	case PROFILE_TAXONOMY:
		args->taxonomy_mode = profile_parse_taxonomy(state, arg);
		break;
	case PROFILE_EDGE_OUT:
		profile_copy_path(state, "--edge-out", args->edge_out,
						  sizeof(args->edge_out), arg);
		break;
	case PROFILE_EDGE_MAX:
		args->edge_max = profile_parse_uint64(state, "--edge-max", arg,
											  UINT64_MAX);
		break;
	case PROFILE_EDGE_MAX_CANDIDATES:
		args->edge_max_candidates =
			(uint32_t)profile_parse_uint64(state, "--edge-max-candidates",
										   arg, UINT32_MAX);
		break;
	case PROFILE_EDGE_ALL:
		args->edge_ambiguous_only = false;
		break;
	case PROFILE_EDGE_SELECTED_ONLY:
		args->edge_selected_only = true;
		break;
	case PROFILE_GTDB_TAXMAP:
		profile_copy_path(state, "--gtdb-taxmap", args->gtdb_taxmap,
						  sizeof(args->gtdb_taxmap), arg);
		break;
	case PROFILE_NCBI_TAXMAP:
		profile_copy_path(state, "--ncbi-taxmap", args->ncbi_taxmap,
						  sizeof(args->ncbi_taxmap), arg);
		break;
	case PROFILE_PIPECMD:
		profile_copy_path(state, "--pipecmd", args->pipecmd,
						  sizeof(args->pipecmd), arg);
		break;
	case PROFILE_DENSITY_BLOCK_CTX:
		args->density_block_ctx =
			(uint32_t)profile_parse_int_range(state, "--density-block-ctx",
											  arg, 0, INT32_MAX);
		args->density_block_set = true;
		break;
	case PROFILE_PRESET:
		profile_apply_preset(state, args, arg);
		break;
	case PROFILE_REPORT_ALL:
		args->report_all = true;
		break;
	case PROFILE_DUAL_EVIDENCE:
		args->dual_evidence = true;
		break;
	case PROFILE_ANI_MODEL:
		args->ani_model = profile_parse_ani_model(state, arg);
		break;
	case PROFILE_ASSIGN:
		args->assign_mode = profile_parse_assign(state, arg);
		break;
	case PROFILE_CTX_FILTER:
		args->ctx_filter_model = profile_parse_ctx_filter(state, arg);
		break;
	case PROFILE_FAKE_THRESHOLD:
		args->fake_threshold =
			profile_parse_nonnegative_double(state, "--fake-threshold", arg);
		break;
	case ARGP_KEY_ARG:
		if (args->reads[0] != '\0')
			argp_error(state, "only one reads input is supported");
		profile_copy_path(state, "reads", args->reads, sizeof(args->reads), arg);
		break;
	case ARGP_KEY_END:
		if (args->ref[0] == '\0')
			argp_error(state, "-r/--ref is required");
		if (args->reads[0] == '\0')
			argp_error(state, "--reads or one positional reads input is required");
		if (args->track[0] == '\0' && args->track_summary[0] != '\0')
			argp_error(state, "--track-summary requires --track");
		if (args->track[0] == '\0' &&
			args->taxonomy_mode != ANI_READWISE_TAXONOMY_NONE)
			argp_error(state, "--taxonomy requires --track");
		if (args->track[0] != '\0' &&
			(args->taxonomy_mode == ANI_READWISE_TAXONOMY_GTDB ||
			 args->taxonomy_mode == ANI_READWISE_TAXONOMY_BOTH) &&
			args->gtdb_taxmap[0] == '\0')
			argp_error(state, "--taxonomy gtdb/both requires --gtdb-taxmap");
		if (args->track[0] != '\0' &&
			(args->taxonomy_mode == ANI_READWISE_TAXONOMY_NCBI ||
			 args->taxonomy_mode == ANI_READWISE_TAXONOMY_BOTH) &&
			args->ncbi_taxmap[0] == '\0')
			argp_error(state, "--taxonomy ncbi/both requires --ncbi-taxmap");
		if (args->cami_profile[0] != '\0' && args->cami_taxmap[0] == '\0')
			argp_error(state, "--cami-profile requires --cami-taxmap");
		if (strcmp(args->cami_profile, "-") == 0 && args->outf[0] == '\0')
			argp_error(state, "--cami-profile - cannot be combined with profile table on stdout; use -o or write the CAMI profile to a file");
		if (strcmp(args->track, "-") == 0 && args->outf[0] == '\0')
			argp_error(state, "--track - cannot be combined with profile table on stdout; use -o or write tracking to a file");
		if (strcmp(args->edge_out, "-") == 0 && args->outf[0] == '\0')
			argp_error(state, "--edge-out - cannot be combined with profile table on stdout; use -o or write edges to a file");
		break;
	default:
		return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

static struct argp argp_profile =
{
	opt_profile,
	parse_profile,
	"-r REF_SKETCH [--reads] READS [options]",
	doc_profile,
	0,
	0,
	0
};

static void apply_profile_args_to_ani(const profile_args_t *args)
{
	ani_opt.fmt = 0;
	ani_opt.p = args->p;
	ani_opt.v = true;
	ani_opt.unassembled = true;
	ani_opt.readwise_profile_only = true;
	ani_opt.readwise_dual_evidence = args->dual_evidence;
	ani_opt.query_density_model = ANI_QUERY_DENSITY_REF;
	ani_opt.abundance_model = ANI_ABUNDANCE_DEPTH;
	ani_opt.readwise_assign_mode = args->assign_mode;
	ani_opt.readwise_ani_model = args->ani_model;
	ani_opt.readwise_ctx_filter_model = args->ctx_filter_model;
	ani_opt.readwise_fake_ctx_threshold = args->fake_threshold;
	ani_opt.ntop = args->ntop;

	if (args->density_block_set)
		ani_opt.density_block_ctx = args->density_block_ctx;
	if (args->track[0] != '\0' && ani_opt.density_block_ctx != 0)
	{
		warnx("--track forces --density-block-ctx 0 for exact read IDs and context offsets");
		ani_opt.density_block_ctx = 0;
	}
	if (args->edge_out[0] != '\0' && ani_opt.density_block_ctx != 0)
	{
		warnx("--edge-out forces --density-block-ctx 0 for exact read/context ambiguity groups");
		ani_opt.density_block_ctx = 0;
	}

	snprintf(ani_opt.refdir, sizeof(ani_opt.refdir), "%s", args->ref);
	snprintf(ani_opt.qrydir, sizeof(ani_opt.qrydir), "%s", args->reads);
	snprintf(ani_opt.outf, sizeof(ani_opt.outf), "%s", args->outf);
	snprintf(ani_opt.cami_profile, sizeof(ani_opt.cami_profile), "%s",
			 args->cami_profile);
	snprintf(ani_opt.cami_taxmap, sizeof(ani_opt.cami_taxmap), "%s",
			 args->cami_taxmap);
	snprintf(ani_opt.cami_sample_id, sizeof(ani_opt.cami_sample_id), "%s",
			 args->cami_sample_id);
	snprintf(ani_opt.readwise_track, sizeof(ani_opt.readwise_track), "%s",
			 args->track);
	snprintf(ani_opt.readwise_track_summary,
			 sizeof(ani_opt.readwise_track_summary), "%s", args->track_summary);
	snprintf(ani_opt.readwise_edge_out, sizeof(ani_opt.readwise_edge_out), "%s",
			 args->edge_out);
	ani_opt.readwise_edge_max = args->edge_max;
	ani_opt.readwise_edge_max_candidates = args->edge_max_candidates;
	ani_opt.readwise_edge_ambiguous_only = args->edge_ambiguous_only;
	ani_opt.readwise_edge_selected_only = args->edge_selected_only;
	snprintf(ani_opt.gtdb_taxmap, sizeof(ani_opt.gtdb_taxmap), "%s",
			 args->gtdb_taxmap);
	snprintf(ani_opt.ncbi_taxmap, sizeof(ani_opt.ncbi_taxmap), "%s",
			 args->ncbi_taxmap);
	snprintf(ani_opt.sketch_pipecmd, sizeof(ani_opt.sketch_pipecmd), "%s",
			 args->pipecmd);
	ani_opt.readwise_taxonomy_mode = args->taxonomy_mode;

	if (args->report_all)
	{
		ani_opt.ctxcut = 0;
		ani_opt.ctxcut_set = true;
		ani_opt.afcut = 0.0f;
		ani_opt.afcut_set = true;
		ani_opt.anicut = 0.0f;
		ani_opt.anicut_set = true;
	}
}

int cmd_profile(struct argp_state *state)
{
	profile_args_t args = {
		.p = 1,
		.ntop = -1,
		.edge_max = 10000000ULL,
		.edge_max_candidates = 64u,
		.density_block_ctx = MINCO_DEFAULT_DENSITY_BLOCK_CTX,
		.assign_mode = ANI_READWISE_ASSIGN_BEST_DIFF_SPLIT,
		.ani_model = ANI_READWISE_ANI_NAIVE,
		.ctx_filter_model = ANI_READWISE_CTX_FILTER_PRODUCT_TOPFRAC_MEDIAN,
		.taxonomy_mode = ANI_READWISE_TAXONOMY_NONE,
		.fake_threshold = 0.25,
		.edge_ambiguous_only = true,
		.edge_selected_only = false,
	};

	int argc = state->argc - state->next + 1;
	char **argv = &state->argv[state->next - 1];
	argp_parse(&argp_profile, argc, argv, ARGP_IN_ORDER, &argc, &args);
	state->next += argc - 1;

	apply_profile_args_to_ani(&args);
	return run_ani_configured(&ani_opt);
}
