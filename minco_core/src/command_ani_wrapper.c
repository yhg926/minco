#include "command_ani.h"
#include "command_sketch_wrapper.h"
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include <argp.h>
#include <argz.h>
#include <err.h>
#include <errno.h>
#include <math.h>
#include <libgen.h>
#include <dirent.h>
#include <omp.h>
#include <stdatomic.h>
#include <ctype.h>
#include <limits.h>
#include <sys/stat.h>
#include <unistd.h>

/*** argp wrapper ***/
struct arg_ani
{
	struct arg_global *global;

	char *name;
	int manual_pattern_seen;
	int coden_pattern_seen;
};

enum
{
	ANI_SKETCH_USE_CODEN = 781,
	ANI_SKETCH_CTXLEN,
	ANI_SKETCH_OUTEROBJLEN,
	ANI_SKETCH_INNEROBJLEN,
	ANI_SKETCH_FILTER_SHIFT,
	ANI_SKETCH_SIZE,
	ANI_SKETCH_LST_KMER_OCRS,
	ANI_SKETCH_NPERCENTILE,
	ANI_SKETCH_NCAP,
	ANI_SKETCH_READS_QC,
	ANI_SKETCH_ABUNDANCE,
	ANI_SKETCH_ANNO,
	ANI_SKETCH_CONFLICT,
	ANI_SKETCH_ASONE,
	ANI_SKETCH_SPLITMFA,
	ANI_SKETCH_PIPECMD,
	ANI_RAW_OUTPUT,
	ANI_UNIFIED_METRIC,
	ANI_QUERY_DENSITY,
	ANI_DENSITY_BLOCK_CTX,
	ANI_SAVE_QUERY_SKETCH,
	ANI_ABUNDANCE_EST
};

enum
{
	ANI_GROUP_INPUT = 1,
	ANI_GROUP_MODE,
	ANI_GROUP_FILTER,
	ANI_GROUP_REPORT,
	ANI_GROUP_EXECUTION,
	ANI_GROUP_AUTO_MODE,
	ANI_GROUP_AUTOSKETCH
};

static struct argp_option opt_ani[] =
	{
		{0, 0, 0, 0, "Inputs:", ANI_GROUP_INPUT},
		{"ref", 'r', "<SKETCH>", 0, "Reference sketch directory.", ANI_GROUP_INPUT},
		{"query", 'q', "<SKETCH|FASTA|FASTQ[.gz]>", 0, "Query sketch directory or sequence input; sequence inputs are auto-sketched for normal assembled-genome ANI.", ANI_GROUP_INPUT},
		{"qraw", 'x', "<SKETCH|FASTA|FASTQ[.gz]>", 0, "Raw-read query sketch or sequence input; keeps query conflicts and uses read/query-aware ANI.", ANI_GROUP_INPUT},
		{"reflist", 779, "<FILE>", 0, "Reference path list for auto ANI; entries may be sketches, FASTA, or FASTQ.", ANI_GROUP_INPUT},
		{"qrylist", 780, "<FILE>", 0, "Query path list for auto ANI; entries may be sketches, FASTA, or FASTQ.", ANI_GROUP_INPUT},
		{"index", 'i', "<FILE>", OPTION_HIDDEN, "Reserved; currently unused.", 2},
		//		{"model", 'M', "<FILE>", 0, "specify the trained ani model.", 2},

		{0, 0, 0, 0, "ANI mode and conflict handling:", ANI_GROUP_MODE},
		{"naive", 'v', 0, 0, "Use naive distance/ANI calculation.", ANI_GROUP_MODE},
		{"ignoreconflict", 778, 0, 0, "Ignore reference-side contexts that contain conflicting objects.", ANI_GROUP_MODE},

		{0, 0, 0, 0, "Filtering and metrics:", ANI_GROUP_FILTER},
		{"afcut", 'f', "<FLOAT>", 0, "Skip reports with max effective aligned fraction below this value; uses density-estimated Real_* AF when ctxmeta is available. [0.5; 0.2 for --qraw]", ANI_GROUP_FILTER},
		{"anicut", 'n', "<FLOAT>", 0, "Skip reports with selected ANI below this value. [0.95]", ANI_GROUP_FILTER},
		{"control", 'c', "<FLOAT>", 0, "Skip duplicated samples with distance below this value. [0]", ANI_GROUP_FILTER},
		{"ctxcut", 't', "<INT>", 0, "Skip reports below this context-support cutoff. In readwise FASTQ density mode this is unique context overlap. [3]", ANI_GROUP_FILTER},
		{"slmetrics", 's', "<+-1..8>", 0, "Selected metric: Best(1), Recalibrated(2), CtxMoE(3), Naive(4), MashD(5), AafD(6), MashD_if_far(7), AafD_if_far(8). Positive values report distance in matrix/triangle formats; negative values report ANI. Detail output prints both. [1]", ANI_GROUP_FILTER},
		{"unified-metric", ANI_UNIFIED_METRIC, 0, 0, "Expert: in --qraw mode, honor -s instead of forcing 1..4 to Naive; Best/Recalibrated still fall back when unavailable.", ANI_GROUP_FILTER},

		{0, 0, 0, 0, "Reporting and output:", ANI_GROUP_REPORT},
		{"diagonal", 'd', 0, 0, "Set diagonal values.", ANI_GROUP_REPORT},
		{"exception", 'e', "<INT>", 0, "Distance value to use when skipped. [1]", ANI_GROUP_REPORT},
		{"glist", 'g', "<FILE>", 0, "Write sample grouping/report file.", ANI_GROUP_REPORT},
		{"outfmt", 'm', "<0/1/2>", 0, "Output format: detail(0), full matrix(1), or lower triangle(2). [0]", ANI_GROUP_REPORT},
		{"outfile", 'o', "<FILE>", 0, "Output file. [STDOUT]", ANI_GROUP_REPORT},
		{"raw-output", ANI_RAW_OUTPUT, 0, 0, "Skip calibrated/best ANI computation; selected calibrated metrics fall back to raw distances when unavailable.", ANI_GROUP_REPORT},
		{"abundance-est", ANI_ABUNDANCE_EST, "<none|depth>", 0, "Experimental: append readwise breadth/depth abundance estimates for direct --qraw FASTQ --query-density ref; default output reports major and low-abundance calls unless filters are set. [none]", ANI_GROUP_REPORT},
		{"top", 'N', "<INT>", 0, "Report at most top N references per query. [all]", ANI_GROUP_REPORT},

		{0, 0, 0, 0, "Execution:", ANI_GROUP_EXECUTION},
		{"threads", 'p', "<INT>", 0, "Number of threads to use.", ANI_GROUP_EXECUTION},

		{0, 0, 0, 0, "Direct FASTA/FASTQ/sketch input mode:", ANI_GROUP_AUTO_MODE},
		{"pair", 777, 0, 0, "Treat positional inputs as one reference followed by one or more queries; optional when -r/-q/--qraw are omitted.", ANI_GROUP_AUTO_MODE},
		{"pipecmd", ANI_SKETCH_PIPECMD, "<cmd>", 0, "Stream each raw input through command; '{}' is replaced by the input path, otherwise the path is appended.", ANI_GROUP_AUTO_MODE},
		{"save-query-sketch", ANI_SAVE_QUERY_SKETCH, "<DIR>", 0, "Debug: save the generated query sketch used for ANI. DIR must not already exist.", ANI_GROUP_AUTO_MODE},
		{"density-block-ctx", ANI_DENSITY_BLOCK_CTX, "<INT>", 0, "Experimental: in direct readwise FASTQ density mode, accumulate consecutive reads until at least this many retained density contexts before lookup; 0 keeps per-read mode. [" MINCO_DEFAULT_DENSITY_BLOCK_CTX_STR "]", ANI_GROUP_AUTO_MODE},

		{0, 0, 0, 0, "Auto-sketch options for FASTA/FASTQ inputs:", ANI_GROUP_AUTOSKETCH},
		{"use_coden_ctxobj", ANI_SKETCH_USE_CODEN, 0, 0, "Use the default coden context-object structure pattern. [default]", ANI_GROUP_AUTOSKETCH},
		{"hash-filter-shift", ANI_SKETCH_FILTER_SHIFT, "<INT>", OPTION_HIDDEN, "Internal source hash prefilter shift for auto-sketching; hidden from public help.", ANI_GROUP_AUTOSKETCH},
		{"sketch-size", 'S', "<INT>", 0, "Number of bottom context hashes to keep per auto-sketched sample. [10000]", ANI_GROUP_AUTOSKETCH},
		{"query-density", ANI_QUERY_DENSITY, "<fixed|ref>", 0, "Sequence query extraction for ANI: fixed bottom-k or reference-density threshold. Temporary sketches are not kept. [fixed]", ANI_GROUP_AUTOSKETCH},
		{"ctxlen", ANI_SKETCH_CTXLEN, "<INT>", 0, "Manual half context length; disables default coden mode.", ANI_GROUP_AUTOSKETCH},
		{"outerobjlen", ANI_SKETCH_OUTEROBJLEN, "<INT>", 0, "Manual half outer object length; disables default coden mode.", ANI_GROUP_AUTOSKETCH},
		{"innerobjlen", ANI_SKETCH_INNEROBJLEN, "<INT>", 0, "Manual inner object length; disables default coden mode.", ANI_GROUP_AUTOSKETCH},
		{"LstKmerOcrs", ANI_SKETCH_LST_KMER_OCRS, "<INT>", 0, "Minimum k-mer/context-object count to retain before final bottom-k. [1]", ANI_GROUP_AUTOSKETCH},
		{"npercentile", ANI_SKETCH_NPERCENTILE, "<0..1>", 0, "Weighted lower-tail count percentile among entries passing --LstKmerOcrs.", ANI_GROUP_AUTOSKETCH},
		{"ncap", ANI_SKETCH_NCAP, "<INT>", OPTION_HIDDEN, "Cap for effective lower k-mer occurrence cutoff; 0 disables. [0]", ANI_GROUP_AUTOSKETCH},
		{"readsQC", ANI_SKETCH_READS_QC, 0, 0, "For noisy reads, infer a count range from a larger hash sample, filter by that range, then take final bottom-k.", ANI_GROUP_AUTOSKETCH},
		{"abundance", ANI_SKETCH_ABUNDANCE, 0, 0, "Keep per-entry counts in temporary sketches.", ANI_GROUP_AUTOSKETCH},
		{"anno", ANI_SKETCH_ANNO, 0, 0, "Write FASTA/FASTQ header annotations to minco.anno.", ANI_GROUP_AUTOSKETCH},
		{"conflict", ANI_SKETCH_CONFLICT, 0, 0, "Keep conflicting context-objects.", ANI_GROUP_AUTOSKETCH},
		{"asone", ANI_SKETCH_ASONE, 0, 0, "Treat each raw input path as one final genome.", ANI_GROUP_AUTOSKETCH},
		{"splitmfa", ANI_SKETCH_SPLITMFA, 0, 0, "Treat each multi-FASTA record as one genome.", ANI_GROUP_AUTOSKETCH},
		{0}};

static char doc_ani[] =
	"\n"
	"Estimate average nucleotide identity (ANI) from minco sketches or direct sequence inputs."
	"\v"
	"Input modes:\n"
	"  -r REF -q QRY\n"
	"      Normal assembled-genome ANI; QRY may be sketch, FASTA, or FASTQ.\n"
	"  -r REF --qraw QRY\n"
	"      Raw-read/query-aware ANI; query conflicts are kept for sequence inputs.\n"
	"  REF QRY...\n"
	"      Positional shortcut: first input is reference, rest are queries.\n"
	"  --reflist/--qrylist\n"
	"      Read reference and query paths from files.\n"
	"\n"
	"Direct FASTA/FASTQ inputs are auto-sketched into temporary directories.\n"
	"If any input is already a sketch, sequence inputs reuse that sketch's\n"
	"minco.stat parameters. Otherwise they use public defaults: coden\n"
	"context/object pattern and --sketch-size 10000.\n"
	"`--query-density ref` extracts raw sequence queries at the reference\n"
	"density inside ANI and deletes the temporary sketch: one-sample references\n"
	"use that sample threshold; combined references use the largest sample\n"
	"density stored in minco.stat, with binary minco.ctxmeta as fallback.\n"
	"Direct --qraw FASTQ with --query-density ref streams reads in density\n"
	"blocks by default, tracks coverage-corrected AF, and appends read-support\n"
	"columns. By default --density-block-ctx 100 batches consecutive reads into\n"
	"pseudo-read blocks with at least 100 retained density contexts before lookup;\n"
	"use --density-block-ctx 0 for exact per-read lookup.\n"
	"OpenMP builds with -p > 1 process this direct FASTQ path in parallel\n"
	"read batches and merge batch state while streaming to keep memory bounded.\n"
	"It reports reads processed, rate, and file percent when available to stderr.\n"
	"Use --save-query-sketch DIR to keep the generated query sketch for debug.\n"
	"For sequence files, -q keeps conflicts only with --conflict; --qraw keeps\n"
	"query conflicts by default. Use --readsQC for noisy long-read FASTQ.\n"
	"\n"
	"Default filters are -n 0.95, -f 0.5, and -t 3. In direct readwise\n"
	"--qraw FASTQ density mode with --abundance-est depth and no custom\n"
	"-f/-n/-t/--top, minco applies the default abundance report: ctxcut=S/100,\n"
	"anicut=0.95, no hidden AF cutoff, and only major or low_abundance calls are printed.\n"
	"In readwise mode, -t is unique context overlap. Use -f0 -n0 -t0 when every\n"
	"comparison must be reported. In --qraw mode without this abundance report,\n"
	"-f defaults to 0.2.\n"
	"When minco.ctxmeta is present, -f uses density-estimated Real_* aligned\n"
	"fractions; otherwise it uses the fixed-sketch aligned fractions.\n"
	"Build an index with `minco sketch -i DIR` before large reference or\n"
	"all-vs-all runs.\n"
	"\n"
	"Output:\n"
	"  -m0 detail: Qry, Ref, ANI, Distance, Confidence, Selected_metric,\n"
	"      diagnostics, Ref_annotation, Real_*_align_fraction, AF_source.\n"
	"  -m1 full matrix and -m2 lower triangle; one sketch gives self output.\n"
		"  --abundance-est depth appends experimental breadth/depth abundance\n"
		"  columns plus Default_call for direct readwise FASTQ density ANI.\n"
		"  With no custom filters it prints default major/low_abundance calls;\n"
		"  normalized abundance sums to 1 over printed rows. It is incompatible\n"
		"  with --readsQC, --abundance, and --save-query-sketch.\n"
	"  Positive -s values print distance in matrix/triangle formats.\n"
	"  Negative -s values print ANI in matrix/triangle formats.\n"
	"  --raw-output skips calibrated/best ANI fields when unavailable.\n"
	"  CtxMoE is a minco-trained Nayfach ANIm MoE calibration for assembled\n"
	"  genomes. Best/Recalibrated apply a matching HGB layer on top of it.\n"
	"  Best may use a density-aware exact-context AAF fallback and report\n"
	"  guarded_low_confidence in a narrow low-AF complete-assembly case.\n"
	"  Use --raw-output -s3 to inspect raw CtxMoE.\n"
	"\n"
	"Use '-' as one raw input to read FASTA/FASTQ from stdin. Use --pipecmd\n"
	"CMD to stream each raw input through a command. '{}' is replaced by the\n"
	"input path; otherwise the path is appended.\n"
	"\n"
	"Examples:\n"
	"  minco sketch -i ref\n"
	"  minco ani -p8 -r ref -q qry -m0 -o ani.tsv\n"
	"  minco ani -p8 -r ref -q query.fna.gz -m0 -o ani.tsv\n"
	"  minco ani -p8 -r ref -q query.fna.gz --query-density ref -m0 -o ani.tsv\n"
	"  minco ani -r ref -q query.fna.gz --query-density ref \\\n"
	"    --save-query-sketch query.debug.minco -o ani.tsv\n"
	"  minco ani -S 20000 -p8 ref.fna query.fna -o ani.20k.tsv\n"
	"  minco ani -p8 -r ref --qraw reads -m0 -o raw.tsv\n"
	"  minco ani -p8 -r ref --qraw reads.fq.gz --readsQC -m0 -o raw.tsv\n"
	"  minco ani -p8 -q genomes -m2 -s -1 -d -o ani.tri.tsv\n"
	"  minco ani -f0 -n0 -t0 -o pair.tsv ref.fna qry.fna\n"
	"  samtools fastq reads.bam | \\\n"
	"    minco ani -r ref --qraw - --readsQC -o raw.tsv\n"
	"  minco ani --pipecmd 'samtools fastq {}' -r ref --qraw reads.bam -o raw.tsv\n"
	"  minco ani --reflist refs.txt --qrylist qrys.txt -o ani.tsv";

ani_opt_t ani_opt = {
	.fmt = 0, // 0:detail, 1:matrix 2: triangle
	.c = 0.0, // control duplicated sample by skip distance < c;
	.p = 1,
	.d = 0,	   // diagonal
	.v = 0,	   // naive model
	.s = 1,	   // select metrics: 1:standard, 2:MashD, 3:AafD, 4:MashD_if_far, 5:AafD_if_far
	.pair = 0, // pairwise compute
	.unassembled = 0,
	.unified_metric = 0,
	.readwise_query = 0,
	.ignoreconflict = 0,
	.raw_output = 0,
	.ctxcut = 3,
	.ctxcut_set = false,
	.afcut = 0.5,
	.afcut_set = false,
	.anicut = 0.95,
	.anicut_set = false,
	.sketch_hclen = 11,
	.sketch_holen = 0,
	.sketch_iolen = 0,
	.sketch_filter_shift = 8,
	.sketch_size = MINCO_DEFAULT_SKETCH_SIZE,
	.sketch_kmerocrs = 1,
	.sketch_ncap = 0,
	.sketch_npercentile = 0.0,
	.sketch_reads_qc = false,
	.sketch_abundance = false,
	.sketch_asone = false,
	.sketch_conflict = false,
	.sketch_anno = false,
	.sketch_split_mfa = false,
	.sketch_coden_ctxobj_pattern = true,
	.density_block_ctx = MINCO_DEFAULT_DENSITY_BLOCK_CTX,
	.query_density_model = ANI_QUERY_DENSITY_FIXED,
	.abundance_model = ANI_ABUNDANCE_NONE,
	.ntop = -1, //all 
	.e = 1, // abort
	.index[0] = '\0',
	.refdir[0] = '\0',
	.qrydir[0] = '\0',
	.reflist[0] = '\0',
	.qrylist[0] = '\0',
	.sketch_pipecmd[0] = '\0',
	.save_query_sketch[0] = '\0',
	.outf[0] = '\0',
	.gl[0] = '\0',
	//	.model[0] = '\0',
	.num_remaining_args = 0, // int num_remaining_args; no option arguments num.
	.remaining_args = NULL	 // char **remaining_args; no option arguments array.
};

static int parse_int_range(struct argp_state *state, const char *option_name,
						   const char *arg, int min_value, int max_value)
{
	char *end = NULL;
	errno = 0;
	long value = strtol(arg, &end, 10);
	if (errno != 0 || end == arg || *end != '\0' || value < min_value || value > max_value)
		argp_error(state, "%s requires an integer in range %d..%d", option_name, min_value, max_value);
	return (int)value;
}

static double parse_nonnegative_double(struct argp_state *state, const char *option_name,
									   const char *arg)
{
	char *end = NULL;
	errno = 0;
	double value = strtod(arg, &end);
	if (errno != 0 || end == arg || *end != '\0' || !isfinite(value) || value < 0.0)
		argp_error(state, "%s requires a non-negative number", option_name);
	return value;
}

static double parse_double_range(struct argp_state *state, const char *option_name,
								 const char *arg, double min_value, double max_value)
{
	char *end = NULL;
	errno = 0;
	double value = strtod(arg, &end);
	if (errno != 0 || end == arg || *end != '\0' || !isfinite(value) || value < min_value || value > max_value)
		argp_error(state, "%s requires a number in range %.6g..%.6g", option_name, min_value, max_value);
	return value;
}

static ani_query_density_model_t parse_query_density_model(struct argp_state *state, const char *arg)
{
	if (strcasecmp(arg, "fixed") == 0 || strcmp(arg, "0") == 0)
		return ANI_QUERY_DENSITY_FIXED;
	if (strcasecmp(arg, "ref") == 0 || strcasecmp(arg, "reference") == 0 || strcmp(arg, "1") == 0)
		return ANI_QUERY_DENSITY_REF;
	argp_error(state, "--query-density must be one of fixed or ref");
	return ANI_QUERY_DENSITY_FIXED;
}

static ani_abundance_model_t parse_abundance_model(struct argp_state *state, const char *arg)
{
	if (strcasecmp(arg, "none") == 0 || strcmp(arg, "0") == 0)
		return ANI_ABUNDANCE_NONE;
	if (strcasecmp(arg, "depth") == 0 || strcasecmp(arg, "breadth") == 0 ||
		strcasecmp(arg, "depth-breadth") == 0 || strcmp(arg, "1") == 0)
		return ANI_ABUNDANCE_DEPTH;
	argp_error(state, "--abundance-est must be one of none or depth");
	return ANI_ABUNDANCE_NONE;
}

static void copy_path_arg(struct argp_state *state, const char *option_name,
						  char *dest, size_t dest_size, const char *arg)
{
	if (strlen(arg) >= dest_size)
		argp_error(state, "%s path is too long; maximum supported length is %zu bytes", option_name, dest_size - 1);
	snprintf(dest, dest_size, "%s", arg);
}

static error_t parse_ani(int key, char *arg, struct argp_state *state)
{
	struct arg_ani *ani = state->input;

	switch (key)
	{
	case 'm':
	{
		ani_opt.fmt = parse_int_range(state, "-m/--outfmt", arg, 0, 2);
		break;
	}
	case 'c':
	{
		ani_opt.c = parse_nonnegative_double(state, "-c/--control", arg);
		break;
	}
	case 'e':
	{
		ani_opt.e = parse_int_range(state, "-e/--exception", arg, 1, INT_MAX);
		break;
	}
	case 's':
	{
		ani_opt.s = parse_int_range(state, "-s/--slmetrics", arg, -8, 8);
		if (ani_opt.s == 0)
			argp_error(state, "-s/--slmetrics accepts -8..-1 or 1..8, but not 0");
		break;
	}
	case 'p':
	{
		ani_opt.p = parse_int_range(state, "-p/--threads", arg, 1, 65536);
		break;
	}
	case ANI_SKETCH_USE_CODEN:
	{
		ani->coden_pattern_seen = 1;
		ani_opt.sketch_coden_ctxobj_pattern = true;
		break;
	}
	case ANI_SKETCH_CTXLEN:
	{
		ani->manual_pattern_seen = 1;
		ani_opt.sketch_coden_ctxobj_pattern = false;
		ani_opt.sketch_hclen = parse_int_range(state, "--ctxlen", arg, 4, 16);
		break;
	}
	case ANI_SKETCH_OUTEROBJLEN:
	{
		ani->manual_pattern_seen = 1;
		ani_opt.sketch_coden_ctxobj_pattern = false;
		ani_opt.sketch_holen = parse_int_range(state, "--outerobjlen", arg, 0, 8);
		break;
	}
	case ANI_SKETCH_INNEROBJLEN:
	{
		ani->manual_pattern_seen = 1;
		ani_opt.sketch_coden_ctxobj_pattern = false;
		ani_opt.sketch_iolen = parse_int_range(state, "--innerobjlen", arg, 0, 8);
		break;
	}
	case ANI_SKETCH_FILTER_SHIFT:
	{
		ani_opt.sketch_filter_shift = parse_int_range(state, "--hash-filter-shift", arg, 0, 24);
		break;
	}
	case 'S':
	case ANI_SKETCH_SIZE:
	{
		ani_opt.sketch_size = (uint32_t)parse_int_range(state, "-S/--sketch-size", arg, 1, INT32_MAX);
		break;
	}
	case ANI_SKETCH_LST_KMER_OCRS:
	{
		ani_opt.sketch_kmerocrs = parse_int_range(state, "--LstKmerOcrs", arg, 1, 65536);
		break;
	}
	case ANI_SKETCH_NPERCENTILE:
	{
		ani_opt.sketch_npercentile = parse_double_range(state, "--npercentile", arg, 0.0, 1.0);
		break;
	}
	case ANI_SKETCH_NCAP:
	{
		ani_opt.sketch_ncap = parse_int_range(state, "--ncap", arg, 0, 65536);
		break;
	}
	case ANI_SKETCH_READS_QC:
	{
		ani_opt.sketch_reads_qc = true;
		if (ani_opt.sketch_kmerocrs < 2)
			ani_opt.sketch_kmerocrs = 2;
		break;
	}
	case ANI_SKETCH_ABUNDANCE:
	{
		ani_opt.sketch_abundance = true;
		break;
	}
	case ANI_SKETCH_ANNO:
	{
		ani_opt.sketch_anno = true;
		break;
	}
	case ANI_SKETCH_CONFLICT:
	{
		ani_opt.sketch_conflict = true;
		break;
	}
	case ANI_SKETCH_ASONE:
	{
		ani_opt.sketch_asone = true;
		break;
	}
	case ANI_SKETCH_SPLITMFA:
	{
		ani_opt.sketch_split_mfa = true;
		break;
	}
	case ANI_SKETCH_PIPECMD:
	{
		copy_path_arg(state, "--pipecmd", ani_opt.sketch_pipecmd, sizeof(ani_opt.sketch_pipecmd), arg);
		break;
	}
	case ANI_SAVE_QUERY_SKETCH:
	{
		copy_path_arg(state, "--save-query-sketch", ani_opt.save_query_sketch, sizeof(ani_opt.save_query_sketch), arg);
		break;
	}
	case ANI_DENSITY_BLOCK_CTX:
	{
		ani_opt.density_block_ctx = (uint32_t)parse_int_range(state, "--density-block-ctx", arg, 0, INT32_MAX);
		break;
	}
	case ANI_RAW_OUTPUT:
	{
		ani_opt.raw_output = true;
		break;
	}
	case ANI_UNIFIED_METRIC:
	{
		ani_opt.unified_metric = true;
		break;
	}
	case ANI_QUERY_DENSITY:
	{
		ani_opt.query_density_model = parse_query_density_model(state, arg);
		break;
	}
	case ANI_ABUNDANCE_EST:
	{
		ani_opt.abundance_model = parse_abundance_model(state, arg);
		break;
	}
	case 't':
	{
		ani_opt.ctxcut = parse_int_range(state, "-t/--ctxcut", arg, 0, INT_MAX);
		ani_opt.ctxcut_set = true;
		break;
	}
	case 'f':
	{
		ani_opt.afcut = (float)parse_nonnegative_double(state, "-f/--afcut", arg);
		ani_opt.afcut_set = true;
		break;
	}
	case 'n':
	{
		ani_opt.anicut = (float)parse_nonnegative_double(state, "-n/--anicut", arg);
		ani_opt.anicut_set = true;
		break;
	}
	case 'N':
	{
		ani_opt.ntop = parse_int_range(state, "-N/--top", arg, 1, INT_MAX);
		break;
	}
	case 'i':
	{
		copy_path_arg(state, "-i/--index", ani_opt.index, sizeof(ani_opt.index), arg);
		break;
	}
	case 'M':
	{
		copy_path_arg(state, "-M/--model", ani_opt.model, sizeof(ani_opt.model), arg);
		break;
	}
	case 'q':
	{
		copy_path_arg(state, "-q/--query", ani_opt.qrydir, sizeof(ani_opt.qrydir), arg);
		break;
	}
	case 780:
	{
		copy_path_arg(state, "--qrylist", ani_opt.qrylist, sizeof(ani_opt.qrylist), arg);
		break;
	}
	case 'r':
	{
		copy_path_arg(state, "-r/--ref", ani_opt.refdir, sizeof(ani_opt.refdir), arg);
		break;
	}
	case 779:
	{
		copy_path_arg(state, "--reflist", ani_opt.reflist, sizeof(ani_opt.reflist), arg);
		break;
	}
	case 'o':
	{
		copy_path_arg(state, "-o/--outfile", ani_opt.outf, sizeof(ani_opt.outf), arg);
		break;
	}
	case 'g':
	{
		copy_path_arg(state, "-g/--glist", ani_opt.gl, sizeof(ani_opt.gl), arg);
		break;
	}
	case 'd':
	{
		ani_opt.d = 1;
		break;
	}
	case 'v':
	{
		ani_opt.v = 1;
		break;
	}
	case 778:
	{
		ani_opt.ignoreconflict = 1;
		break;
	}
	case 'x':
	{
		copy_path_arg(state, "--qraw", ani_opt.qrydir, sizeof(ani_opt.qrydir), arg);
		ani_opt.unassembled = 1; 
		ani_opt.v = 1;
		break;
	}
	case 777:
	{
		ani_opt.pair = 1;
		break;
	}
	case ARGP_KEY_ARGS:
	{
		ani_opt.num_remaining_args = state->argc - state->next;
		ani_opt.remaining_args = state->argv + state->next;
		break;
	}
	case ARGP_KEY_END:
	{
		const bool has_reflist = ani_opt.reflist[0] != '\0';
		const bool has_qrylist = ani_opt.qrylist[0] != '\0';

		if (ani->coden_pattern_seen && ani->manual_pattern_seen)
			argp_error(state, "--use_coden_ctxobj selects the coden pattern; do not combine it with --ctxlen, --outerobjlen, or --innerobjlen.");
		if (!ani_opt.sketch_coden_ctxobj_pattern)
		{
			int sketch_klen = ani_opt.sketch_iolen + 2 * (ani_opt.sketch_holen + ani_opt.sketch_hclen);
			if (sketch_klen > 32)
				argp_error(state, "auto-sketch k-mer length %d should be at most 32", sketch_klen);
		}
		if (ani_opt.sketch_pipecmd[0] != '\0' && ani_opt.sketch_split_mfa)
			argp_error(state, "--splitmfa does not support --pipecmd streaming inputs");

		if (has_reflist || has_qrylist)
		{
			if (!has_reflist || !has_qrylist)
			{
				printf("\nError: --reflist and --qrylist must be used together.\n\n");
				argp_state_help(state, stdout, ARGP_HELP_STD_HELP);
				argp_usage(state);
			}
			if (ani_opt.refdir[0] != '\0' || ani_opt.qrydir[0] != '\0')
			{
				printf("\nError: --reflist/--qrylist use auto ANI inputs; do not combine them with -r, -q, or --qraw.\n\n");
				argp_state_help(state, stdout, ARGP_HELP_STD_HELP);
				argp_usage(state);
			}
			if (ani_opt.pair || ani_opt.num_remaining_args > 0)
			{
				printf("\nError: --reflist/--qrylist cannot be combined with --pair or positional inputs.\n\n");
				argp_state_help(state, stdout, ARGP_HELP_STD_HELP);
				argp_usage(state);
			}
		}
		else if (ani_opt.pair)
		{
			if (ani_opt.refdir[0] != '\0' || ani_opt.qrydir[0] != '\0')
			{
				printf("\nError: --pair uses positional FASTA/FASTQ/sketch inputs; do not combine it with -r, -q, or --qraw.\n\n");
				argp_state_help(state, stdout, ARGP_HELP_STD_HELP);
				argp_usage(state);
			}
			if (ani_opt.num_remaining_args < 2)
			{
				printf("\nError: --pair requires one reference input followed by at least one query input.\n\n");
				argp_state_help(state, stdout, ARGP_HELP_STD_HELP);
				argp_usage(state);
			}
		}
			else if (ani_opt.refdir[0] == '\0' && ani_opt.qrydir[0] == '\0' &&
					 ani_opt.num_remaining_args >= 2)
			{
				ani_opt.pair = 1;
			}
			else if (ani_opt.refdir[0] == '\0' && ani_opt.qrydir[0] == '\0' &&
					 ani_opt.num_remaining_args == 1 && ani_opt.fmt != 0)
			{
				copy_path_arg(state, "sketch", ani_opt.qrydir, sizeof(ani_opt.qrydir),
							  ani_opt.remaining_args[0]);
			}
			else if (ani_opt.qrydir[0] == '\0')
			{
				printf("\nError: Mandatory options: '-q' or '--qraw' are missing unless positional inputs are used.\n\n");
				argp_state_help(state, stdout, ARGP_HELP_STD_HELP);
				argp_usage(state);
			}
			else if (ani_opt.refdir[0] == '\0')
			{
				if (ani_opt.fmt == 0)
				{
					printf("\nError: -r/--ref is required with -q/--query or --qraw for detail output. Use -m1 or -m2 for one-sketch self matrices.\n\n");
					argp_state_help(state, stdout, ARGP_HELP_STD_HELP);
					argp_usage(state);
				}
			}
		if (ani_opt.s < -8 || ani_opt.s > 8 || ani_opt.s == 0)
			argp_error(state, "-s option should be within range 1..8 or -8..-1");

		break;
		/*
					if(state->argc<2)
					{
				printf("\v");
						argp_state_help(state,stdout,ARGP_HELP_SHORT_USAGE);
						printf("\v");
				argp_state_help(state,stdout,ARGP_HELP_LONG);
				printf("\v");
				return EINVAL;
					}
		*/
	}
	default:
		return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

static struct argp argp_ani =
	{
		opt_ani,
		parse_ani,
			"[-r REF_SKETCH -q QRY | --qraw RAW_QRY | --reflist REF_LIST --qrylist QRY_LIST | [--pair] REF_INPUT QUERY_INPUT... | SELF_SKETCH with -m1/-m2]",
			doc_ani};
extern const char sorted_comb_ctxgid64obj32[];
extern void gen_inverted_index_for_minco(const char *sketchdir);

get_generic_dist_from_features_fn get_generic_dist_from_features = NULL;
int ani_model_compat_filter_shift = ANI_MODEL_REFERENCE_FILTER_SHIFT;
uint32_t ani_model_target_sketch_size = ANI_MODEL_REFERENCE_SKETCH_SIZE;

static int qraw_multi_query_index_threshold(void)
{
	const char *env = getenv("MINCO_QRAW_MULTI_THRESHOLD");
	if (!env || env[0] == '\0')
		return 20;

	char *end = NULL;
	errno = 0;
	long value = strtol(env, &end, 10);
	if (errno != 0 || end == env || *end != '\0' || value < 1 || value > INT_MAX)
		errx(EXIT_FAILURE, "invalid MINCO_QRAW_MULTI_THRESHOLD='%s' (use positive integer)", env);
	return (int)value;
}

static int auto_ref_index_threshold(void)
{
	const char *env = getenv("MINCO_AUTO_REF_INDEX_THRESHOLD");
	if (!env || env[0] == '\0')
		return 100;

	char *end = NULL;
	errno = 0;
	long value = strtol(env, &end, 10);
	if (errno != 0 || end == env || *end != '\0' || value < 0 || value > INT_MAX)
		errx(EXIT_FAILURE, "invalid MINCO_AUTO_REF_INDEX_THRESHOLD='%s' (use non-negative integer; 0 disables)", env);
	return (int)value;
}

static int ani_matrix_direct_threshold(void)
{
	const char *name = "MINCO_ANI_MATRIX_DIRECT_THRESHOLD";
	const char *env = getenv(name);
	if (!env || env[0] == '\0') {
		name = "MINCO_ANI_M1_DIRECT_THRESHOLD";
		env = getenv(name);
	}
	if (!env || env[0] == '\0')
		return 1000;

	char *end = NULL;
	errno = 0;
	long value = strtol(env, &end, 10);
	if (errno != 0 || end == env || *end != '\0' || value < 0 || value > INT_MAX)
		errx(EXIT_FAILURE, "invalid %s='%s' (use non-negative integer)", name, env);
	return (int)value;
}

static bool ani_matrix_should_use_index(int ref_infile_num, int qry_infile_num, int direct_threshold)
{
	if (ref_infile_num < 1 || qry_infile_num < 1)
		return false;
	if (direct_threshold == 0)
		return true;
	return ref_infile_num > direct_threshold || qry_infile_num > direct_threshold;
}

static void ensure_sorted_ref_index_for_ani_m1(const char *refdir)
{
	if (file_exists_in_folder(refdir, sorted_comb_ctxgid64obj32))
		return;
	gen_inverted_index_for_minco(refdir);
	if (!file_exists_in_folder(refdir, sorted_comb_ctxgid64obj32))
		errx(EXIT_FAILURE, "failed to create reference index '%s/%s'",
			 refdir, sorted_comb_ctxgid64obj32);
}

static int sketch_infile_num(const char *sketch_dir)
{
	size_t stat_size = 0;
	char *stat_path = test_get_fullpath(sketch_dir, sketch_stat);
	void *stat = read_from_file(stat_path, &stat_size);
	free(stat_path);
	minco_sketch_stat_t decoded = {0};
	if (!minco_stat_decode_mem(stat, stat_size, &decoded, NULL))
		errx(EINVAL, "%s(): malformed %s/%s", __func__, sketch_dir, sketch_stat);
	const int infile_num = decoded.infile_num;
	free_read_from_file(stat, stat_size);
	return infile_num;
}

static bool force_ref_index_requested(void)
{
	const char *env = getenv("MINCO_FORCE_REF_INDEX");
	return env && env[0] != '\0' && strcmp(env, "0") != 0;
}

typedef struct
{
	char **items;
	int n;
	int cap;
} path_vec_t;

static void path_vec_push(path_vec_t *vec, const char *path)
{
	if (vec->n == vec->cap)
	{
		int new_cap = vec->cap ? vec->cap * 2 : 8;
		char **new_items = realloc(vec->items, (size_t)new_cap * sizeof(new_items[0]));
		if (!new_items)
			err(errno, "%s(): OOM path vector", __func__);
		vec->items = new_items;
		vec->cap = new_cap;
	}
	vec->items[vec->n] = strdup(path);
	if (!vec->items[vec->n])
		err(errno, "%s(): OOM path copy", __func__);
	vec->n++;
}

static void path_vec_free(path_vec_t *vec)
{
	for (int i = 0; i < vec->n; ++i)
		free(vec->items[i]);
	free(vec->items);
	vec->items = NULL;
	vec->n = 0;
	vec->cap = 0;
}

static char *trim_path_line(char *line)
{
	char *start = line;
	while (isspace((unsigned char)*start))
		start++;

	char *end = start + strlen(start);
	while (end > start && isspace((unsigned char)end[-1]))
		end--;
	*end = '\0';
	return start;
}

static void read_path_list_to_vec(const char *list_path, path_vec_t *out, const char *option_name)
{
	FILE *list = fopen(list_path, "r");
	if (!list)
		err(errno, "can't open %s file %s", option_name, list_path);

	char *line = malloc(LMAX);
	if (!line)
		err(errno, "%s(): OOM list line buffer", __func__);

	int line_no = 0;
	while (fgets(line, LMAX, list) != NULL)
	{
		line_no++;
		size_t len = strlen(line);
		if (len == LMAX - 1 && line[len - 1] != '\n' && line[len - 1] != '\r')
			errx(EXIT_FAILURE, "%s file %s line %d exceeds maximum line length %d",
				 option_name, list_path, line_no, LMAX - 1);

		char *path = trim_path_line(line);
		if (path[0] == '\0')
			continue;
		if (strlen(path) >= PATHLEN)
			errx(EXIT_FAILURE, "%s file %s line %d path exceeds maximum length %d",
				 option_name, list_path, line_no, PATHLEN - 1);
		path_vec_push(out, path);
	}
	if (ferror(list))
		err(errno, "failed while reading %s file %s", option_name, list_path);

	free(line);
	fclose(list);

	if (out->n < 1)
		errx(EXIT_FAILURE, "%s file %s has no input paths", option_name, list_path);
}

static bool is_minco_sketch_dir(const char *path)
{
	struct stat st;
	if (stat(path, &st) != 0 || !S_ISDIR(st.st_mode))
		return false;
	return file_exists_in_folder(path, sketch_stat) &&
		   file_exists_in_folder(path, combined_sketch_suffix) &&
		   file_exists_in_folder(path, idx_sketch_suffix);
}

static minco_sketch_stat_t read_minco_sketch_stat_info(const char *sketch_dir,
												   minco_sketch_info_t *info_out)
{
	size_t stat_size = 0;
	char *stat_path = test_get_fullpath(sketch_dir, sketch_stat);
	void *stat = read_from_file(stat_path, &stat_size);
	free(stat_path);
	minco_sketch_stat_t copy = {0};
	minco_sketch_info_t info = {0};
	if (!minco_stat_decode_mem(stat, stat_size, &copy, &info))
		err(EINVAL, "%s(): %s/%s has %zu bytes, expected at least %zu",
			__func__, sketch_dir, sketch_stat, stat_size, sizeof(copy));
	if (info.has_minco_ext)
		copy.hash_id = info.sketch_id;
	if (info_out)
		*info_out = info;
	free_read_from_file(stat, stat_size);
	return copy;
}

static minco_sketch_stat_t read_minco_sketch_stat(const char *sketch_dir)
{
	return read_minco_sketch_stat_info(sketch_dir, NULL);
}

typedef struct ani_ref_density_threshold
{
	bool enabled;
	uint64_t threshold;
	uint32_t hash_bits;
} ani_ref_density_threshold_t;

static const char ani_ctxmeta_legacy_tsv_name[] = "minco.ctxmeta.tsv";

static int ani_split_tsv_fields(char *line, char **fields, int max_fields)
{
	int n = 0;
	char *p = line;
	while (n < max_fields)
	{
		fields[n++] = p;
		char *tab = strchr(p, '\t');
		if (!tab)
			break;
		*tab = '\0';
		p = tab + 1;
	}
	return n;
}

static bool ani_parse_u64_field(const char *s, uint64_t *out)
{
	if (!s || !*s)
		return false;
	errno = 0;
	char *end = NULL;
	unsigned long long value = strtoull(s, &end, 10);
	if (errno != 0 || end == s || *end != '\0')
		return false;
	*out = (uint64_t)value;
	return true;
}

static long double ani_threshold_density(uint64_t threshold, uint32_t hash_bits)
{
	if (hash_bits == 0 || hash_bits > 64)
		return 0.0L;
	const uint64_t max_threshold =
		hash_bits >= 64 ? UINT64_MAX : ((1ULL << hash_bits) - 1ULL);
	if (threshold >= max_threshold)
		return 1.0L;
	const long double density =
		((long double)threshold + 1.0L) / ldexpl(1.0L, (int)hash_bits);
	return density > 1.0L ? 1.0L : density;
}

static uint32_t ani_hash_bits_from_stat(const minco_sketch_stat_t *stat)
{
	const int ctx_bits = stat->coden_len > 0 ? 4 * stat->coden_len : 4 * stat->hclen;
	const int obj_bits = 2 * stat->klen - ctx_bits;
	if (obj_bits < 0 || obj_bits > 64)
		errx(EXIT_FAILURE, "invalid sketch context/object lengths in reference template");
	return (uint32_t)(64 - obj_bits);
}

static bool ani_read_ctxset_largest_threshold(const char *sketch_dir,
											  ani_ref_density_threshold_t *out)
{
	if (!file_exists_in_folder(sketch_dir, minco_ctxsetmeta_legacy_tsv_stat))
		return false;
	char *path = test_get_fullpath(sketch_dir, minco_ctxsetmeta_legacy_tsv_stat);
	FILE *fp = fopen(path, "r");
	if (!fp)
		err(errno, "%s(): cannot open %s", __func__, path);

	bool have_threshold = false;
	bool have_hash_bits = false;
	uint64_t threshold = 0;
	uint64_t hash_bits_u64 = 0;
	char *line = NULL;
	size_t cap = 0;
	while (getline(&line, &cap, fp) >= 0)
	{
		line[strcspn(line, "\r\n")] = '\0';
		if (line[0] == '\0')
			continue;
		char *fields[2] = {0};
		if (ani_split_tsv_fields(line, fields, 2) < 2)
			continue;
		if (strcmp(fields[0], "largest_sample_threshold") == 0)
			have_threshold = ani_parse_u64_field(fields[1], &threshold);
		else if (strcmp(fields[0], "largest_sample_hash_bits") == 0)
			have_hash_bits = ani_parse_u64_field(fields[1], &hash_bits_u64);
	}
	free(line);
	if (fclose(fp) != 0)
		err(errno, "%s(): cannot close %s", __func__, path);
	free(path);

	if (!have_threshold || !have_hash_bits || hash_bits_u64 > UINT32_MAX)
		return false;
	*out = (ani_ref_density_threshold_t){
		.enabled = true,
		.threshold = threshold,
		.hash_bits = (uint32_t)hash_bits_u64,
	};
	return true;
}

static bool ani_read_ctxmeta_threshold(const char *sketch_dir, bool choose_largest_density,
									   ani_ref_density_threshold_t *out)
{
	if (file_exists_in_folder(sketch_dir, minco_ctxmeta_bin_stat))
	{
		minco_sketch_stat_t stat = read_minco_sketch_stat(sketch_dir);
		char *path = test_get_fullpath(sketch_dir, minco_ctxmeta_bin_stat);
		size_t file_size = 0;
		minco_ctxmeta_record_t *records = read_from_file(path, &file_size);
		free(path);
		const size_t expected_size = (size_t)stat.infile_num * sizeof(records[0]);
		if (file_size != expected_size)
			errx(EXIT_FAILURE, "%s(): %s/%s has %zu bytes, expected %zu",
				 __func__, sketch_dir, minco_ctxmeta_bin_stat,
				 file_size, expected_size);
		bool found = false;
		long double best_density = 0.0L;
		for (int i = 0; i < stat.infile_num; ++i)
		{
			if (!records[i].valid)
				continue;
			const long double density =
				ani_threshold_density(records[i].threshold, records[i].hash_bits);
			if (!choose_largest_density || !found || density > best_density)
			{
				*out = (ani_ref_density_threshold_t){
					.enabled = true,
					.threshold = records[i].threshold,
					.hash_bits = records[i].hash_bits,
				};
				best_density = density;
				found = true;
				if (!choose_largest_density)
					break;
			}
		}
		free_read_from_file(records, file_size);
		return found;
	}

	if (!file_exists_in_folder(sketch_dir, ani_ctxmeta_legacy_tsv_name))
		return false;
	char *path = test_get_fullpath(sketch_dir, ani_ctxmeta_legacy_tsv_name);
	FILE *fp = fopen(path, "r");
	if (!fp)
		err(errno, "%s(): cannot open %s", __func__, path);

	char *line = NULL;
	size_t cap = 0;
	ssize_t len = getline(&line, &cap, fp);
	if (len < 0)
	{
		free(line);
		fclose(fp);
		free(path);
		return false;
	}

	bool found = false;
	long double best_density = 0.0L;
	while ((len = getline(&line, &cap, fp)) >= 0)
	{
		line[strcspn(line, "\r\n")] = '\0';
		if (line[0] == '\0')
			continue;
		char *fields[13] = {0};
		if (ani_split_tsv_fields(line, fields, 13) < 13)
			continue;
		uint64_t valid = 0;
		uint64_t hash_bits_u64 = 0;
		uint64_t threshold = 0;
		if (!ani_parse_u64_field(fields[3], &valid) ||
			!ani_parse_u64_field(fields[4], &hash_bits_u64) ||
			!ani_parse_u64_field(fields[5], &threshold))
			continue;
		if (!valid || hash_bits_u64 > UINT32_MAX)
			continue;
		const uint32_t hash_bits = (uint32_t)hash_bits_u64;
		const long double density = ani_threshold_density(threshold, hash_bits);
		if (!choose_largest_density)
		{
			*out = (ani_ref_density_threshold_t){
				.enabled = true,
				.threshold = threshold,
				.hash_bits = hash_bits,
			};
			found = true;
			break;
		}
		if (!found || density > best_density)
		{
			*out = (ani_ref_density_threshold_t){
				.enabled = true,
				.threshold = threshold,
				.hash_bits = hash_bits,
			};
			best_density = density;
			found = true;
		}
	}

	free(line);
	if (fclose(fp) != 0)
		err(errno, "%s(): cannot close %s", __func__, path);
	free(path);
	return found;
}

static bool ani_read_stat_density_threshold(const minco_sketch_info_t *info,
											ani_ref_density_threshold_t *out)
{
	if (!info || !out || !info->has_minco_ext ||
		info->density_valid_sample_count == 0 || info->density_hash_bits == 0)
		return false;
	*out = (ani_ref_density_threshold_t){
		.enabled = true,
		.threshold = info->density_universal_threshold,
		.hash_bits = info->density_hash_bits,
	};
	return true;
}

static ani_ref_density_threshold_t ani_select_ref_density_threshold(const char *refdir,
																	const minco_sketch_stat_t *template_stat)
{
	minco_sketch_info_t ref_info = {0};
	minco_sketch_stat_t ref_stat = read_minco_sketch_stat_info(refdir, &ref_info);
	ani_ref_density_threshold_t selected = {0};
	if (ani_read_stat_density_threshold(&ref_info, &selected))
	{
		/* minco.stat v1 stores one operational threshold: sample density for
		 * single-reference sketches and largest sample density for combined refs. */
	}
	else if (ref_stat.infile_num <= 1)
	{
		if (!ani_read_ctxmeta_threshold(refdir, false, &selected))
			errx(EXIT_FAILURE,
				 "--query-density ref requires density metadata in %s/%s or %s/%s for a one-sample reference sketch",
				 refdir, sketch_stat, refdir, minco_ctxmeta_bin_stat);
	}
	else if (!ani_read_ctxset_largest_threshold(refdir, &selected) &&
			 !ani_read_ctxmeta_threshold(refdir, true, &selected))
	{
		errx(EXIT_FAILURE,
			 "--query-density ref requires density metadata in %s/%s or %s/%s for a combined reference sketch",
			 refdir, sketch_stat, refdir, minco_ctxmeta_bin_stat);
	}

	const uint32_t expected_hash_bits = ani_hash_bits_from_stat(template_stat);
	if (selected.hash_bits != expected_hash_bits)
		errx(EXIT_FAILURE,
			 "--query-density ref metadata hash_bits=%u does not match sketch template hash_bits=%u",
			 selected.hash_bits, expected_hash_bits);
	return selected;
}

static bool same_ani_auto_sketch_stat(const minco_sketch_stat_t *a, const minco_sketch_stat_t *b)
{
	return a->hash_id == b->hash_id &&
		   a->koc == b->koc &&
		   a->conflict == b->conflict &&
		   a->coden_len == b->coden_len &&
		   a->klen == b->klen &&
		   a->hclen == b->hclen &&
		   a->holen == b->holen &&
		   a->compat_filter_shift == b->compat_filter_shift;
}

static void describe_sketch_stat_mismatch(
	const char *path,
	const minco_sketch_stat_t *expected,
	const minco_sketch_stat_t *actual)
{
	errx(EXIT_FAILURE,
		 "%s sketch parameters do not match first auto ANI sketch "
		 "(hash_id %u/%u, koc %d/%d, conflict %d/%d, coden_len %d/%d, "
		 "klen %d/%d, hclen %d/%d, holen %d/%d, compat_filter_shift %d/%d)",
		 path,
		 actual->hash_id, expected->hash_id,
		 actual->koc, expected->koc,
		 actual->conflict, expected->conflict,
		 actual->coden_len, expected->coden_len,
		 actual->klen, expected->klen,
		 actual->hclen, expected->hclen,
		 actual->holen, expected->holen,
		 actual->compat_filter_shift, expected->compat_filter_shift);
}

extern uint32_t FILTER;
extern uint32_t hash_id;
extern minco_sketch_stat_t minco_stat_one;
extern uint32_t minco_target_sketch_size;
extern uint32_t minco_selection_mode;
extern uint64_t minco_density_threshold;
extern uint32_t minco_compiled_stat_flags(void);
extern uint32_t get_sketching_id(uint32_t hclen, uint32_t holen, uint32_t iolen, uint32_t compat_filter_shift, uint32_t FILTER);
extern void compute_sketch(sketch_opt_t *sketch_opt_val, infile_tab_t *infile_stat);
extern int merge_minco_sketches(sketch_opt_t *sketch_opt_val);

static minco_sketch_stat_t default_auto_ani_sketch_stat(const ani_opt_t *opt)
{
	minco_sketch_stat_t stat = {0};
	stat.koc = opt->sketch_abundance;
	stat.conflict = opt->sketch_conflict;
	stat.compat_filter_shift = 0;
	if (opt->sketch_coden_ctxobj_pattern)
	{
		stat.coden_len = NUM_CODENS;
#if NUM_CODENS < 11
		stat.klen = 3 * NUM_CODENS + 1;
#else
		stat.klen = 32;
#endif
		stat.hclen = 0;
		stat.holen = 0;
	}
	else
	{
		stat.coden_len = 0;
		stat.hclen = opt->sketch_hclen;
		stat.holen = opt->sketch_holen;
		stat.klen = opt->sketch_iolen + 2 * (opt->sketch_hclen + opt->sketch_holen);
	}
	stat.hash_id = minco_stat_sketch_id_from_dim(&stat,
												 opt->sketch_size ? opt->sketch_size : MINCO_DEFAULT_SKETCH_SIZE,
												 MINCO_STAT_SELECTION_BOTTOMK,
												 minco_compiled_stat_flags());
	stat.infile_num = 0;
	return stat;
}

static int stat_inner_object_len(const minco_sketch_stat_t *stat)
{
	return stat->klen - 2 * (stat->hclen + stat->holen);
}

static void validate_supported_auto_sketch_stat(const minco_sketch_stat_t *stat)
{
	if (stat->compat_filter_shift < 0 || stat->compat_filter_shift > 32)
		errx(EXIT_FAILURE, "unsupported sketch compat_filter_shift %d in first auto ANI sketch", stat->compat_filter_shift);
	if (stat->klen > 32)
		errx(EXIT_FAILURE, "first auto ANI sketch has unsupported klen=%d; expected at most 32", stat->klen);
	if (stat->coden_len > 0)
	{
		if (stat->coden_len != NUM_CODENS)
			errx(EXIT_FAILURE,
				 "first auto ANI sketch uses coden_len=%d, but this binary was built with NUM_CODENS=%d",
				 stat->coden_len, NUM_CODENS);
#if NUM_CODENS < 11
		const int expected_klen = 3 * NUM_CODENS + 1;
#else
		const int expected_klen = 32;
#endif
		if (stat->klen != expected_klen)
			errx(EXIT_FAILURE,
				 "first auto ANI sketch klen=%d is not compatible with coden_len=%d in this binary",
				 stat->klen, stat->coden_len);
	}
	else if (stat_inner_object_len(stat) < 0)
	{
		errx(EXIT_FAILURE,
			 "first auto ANI sketch has invalid manual context-object lengths: klen=%d hclen=%d holen=%d",
			 stat->klen, stat->hclen, stat->holen);
	}
}

static sketch_opt_t sketch_opt_from_stat(const minco_sketch_stat_t *stat, const char *outdir,
										 const ani_opt_t *ani_opt_val,
										 uint32_t target_sketch_size)
{
	sketch_opt_t opt = {
		.hclen = stat->hclen,
		.holen = stat->holen,
		.iolen = stat_inner_object_len(stat),
		.compat_filter_shift = 0,
		.sketch_size = target_sketch_size ? target_sketch_size : ani_opt_val->sketch_size,
		.kmerocrs = ani_opt_val->sketch_kmerocrs,
		.npercentile = ani_opt_val->sketch_npercentile,
		.ncap = ani_opt_val->sketch_ncap,
		.reads_qc = ani_opt_val->sketch_reads_qc,
		.p = ani_opt_val->p > 0 ? ani_opt_val->p : 1,
		.abundance = stat->koc,
		.asone = ani_opt_val->sketch_asone,
		.conflict = stat->conflict,
		.anno = ani_opt_val->sketch_anno,
		.compute_meta = true,
		.ctxmeta_mode = MINCO_CTXMETA_PRECONFLICT,
		.density_threshold_enabled = false,
		.density_threshold = UINT64_MAX,
		.merge_minco_mode = false,
		.sketch_qc = false,
		.split_mfa = ani_opt_val->sketch_split_mfa,
		.coden_ctxobj_pattern = stat->coden_len > 0,
		.index[0] = '\0',
		.fpath = NULL,
		.outdir = (char *)outdir,
		.pipecmd = ani_opt_val->sketch_pipecmd[0] ? (char *)ani_opt_val->sketch_pipecmd : NULL,
		.num_remaining_args = 0,
		.remaining_args = NULL,
	};
	if (opt.coden_ctxobj_pattern)
		opt.iolen = 0;
	return opt;
}

static void apply_auto_sketch_stat(const minco_sketch_stat_t *stat, int infile_num,
								   uint32_t target_sketch_size,
								   uint32_t selection_mode,
								   uint64_t density_threshold)
{
	FILTER = UINT32_MAX;
	minco_target_sketch_size = target_sketch_size ? target_sketch_size : MINCO_DEFAULT_SKETCH_SIZE;
	minco_selection_mode = selection_mode ? selection_mode : MINCO_STAT_SELECTION_BOTTOMK;
	minco_density_threshold = density_threshold;
	minco_stat_one = *stat;
	minco_stat_one.compat_filter_shift = 0;
	minco_stat_one.infile_num = infile_num;
	hash_id = minco_stat_sketch_id_from_dim(&minco_stat_one,
											minco_target_sketch_size,
											minco_selection_mode,
											minco_compiled_stat_flags());
	minco_stat_one.hash_id = hash_id;
	const_comask_init(&minco_stat_one);
	set_uint64kmer2generic_ctxobj(stat->coden_len > 0);
}

static void free_infile_tab(infile_tab_t *tab)
{
	if (!tab)
		return;
	for (int i = 0; i < tab->infile_num; ++i)
		free(tab->organized_infile_tab[i].fpath);
	free(tab->organized_infile_tab);
	free(tab);
}

static void remove_tree(const char *path)
{
	DIR *dir = opendir(path);
	if (!dir)
		return;
	struct dirent *ent;
	while ((ent = readdir(dir)) != NULL)
	{
		if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0)
			continue;
		char child[PATHLEN * 2];
		snprintf(child, sizeof(child), "%s/%s", path, ent->d_name);
		struct stat st;
		if (lstat(child, &st) != 0)
			continue;
		if (S_ISDIR(st.st_mode))
			remove_tree(child);
		else
			unlink(child);
	}
	closedir(dir);
	rmdir(path);
}

static bool path_exists_any(const char *path)
{
	struct stat st;
	return path && lstat(path, &st) == 0;
}

static void ensure_save_query_sketch_path_available(const char *path)
{
	if (!path || path[0] == '\0')
		return;
	if (path_exists_any(path))
		errx(EXIT_FAILURE, "--save-query-sketch target already exists: %s", path);
}

static void mkdir_parent_dirs_for_path(const char *path)
{
	char parent[PATHLEN * 2];
	if (strlen(path) >= sizeof(parent))
		errx(EXIT_FAILURE, "path is too long: %s", path);
	snprintf(parent, sizeof(parent), "%s", path);
	char *slash = strrchr(parent, '/');
	if (!slash || slash == parent)
		return;
	*slash = '\0';
	mkdir_p(parent);
}

static void copy_regular_file(const char *src, const char *dst)
{
	FILE *in = fopen(src, "rb");
	if (!in)
		err(errno, "%s(): cannot open %s", __func__, src);
	FILE *out = fopen(dst, "wb");
	if (!out)
		err(errno, "%s(): cannot create %s", __func__, dst);
	setvbuf(in, NULL, _IOFBF, 8u << 20);
	setvbuf(out, NULL, _IOFBF, 8u << 20);

	const size_t buf_size = 8u << 20;
	char *buf = (char *)malloc(buf_size);
	if (!buf)
		err(errno, "%s(): OOM copy buffer", __func__);
	for (;;)
	{
		size_t got = fread(buf, 1, buf_size, in);
		if (got > 0 && fwrite(buf, 1, got, out) != got)
			err(errno, "%s(): failed writing %s", __func__, dst);
		if (got < buf_size)
		{
			if (ferror(in))
				err(errno, "%s(): failed reading %s", __func__, src);
			break;
		}
	}
	free(buf);
	if (fclose(in) != 0)
		err(errno, "%s(): failed closing %s", __func__, src);
	if (fclose(out) != 0)
		err(errno, "%s(): failed closing %s", __func__, dst);
}

static void copy_tree_to_new_dir(const char *src, const char *dst)
{
	struct stat st;
	if (lstat(src, &st) != 0)
		err(errno, "%s(): cannot stat %s", __func__, src);
	if (!S_ISDIR(st.st_mode))
		errx(EXIT_FAILURE, "%s(): source is not a directory: %s", __func__, src);
	mkdir_parent_dirs_for_path(dst);
	if (mkdir(dst, 0777) != 0)
		err(errno, "%s(): cannot create %s", __func__, dst);

	DIR *dir = opendir(src);
	if (!dir)
		err(errno, "%s(): cannot open %s", __func__, src);
	struct dirent *ent;
	while ((ent = readdir(dir)) != NULL)
	{
		if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0)
			continue;
		char child_src[PATHLEN * 2];
		char child_dst[PATHLEN * 2];
		int written = snprintf(child_src, sizeof(child_src), "%s/%s", src, ent->d_name);
		if (written < 0 || (size_t)written >= sizeof(child_src))
			errx(EXIT_FAILURE, "%s(): source path too long under %s", __func__, src);
		written = snprintf(child_dst, sizeof(child_dst), "%s/%s", dst, ent->d_name);
		if (written < 0 || (size_t)written >= sizeof(child_dst))
			errx(EXIT_FAILURE, "%s(): destination path too long under %s", __func__, dst);
		if (lstat(child_src, &st) != 0)
			err(errno, "%s(): cannot stat %s", __func__, child_src);
		if (S_ISDIR(st.st_mode))
			copy_tree_to_new_dir(child_src, child_dst);
		else if (S_ISREG(st.st_mode))
			copy_regular_file(child_src, child_dst);
		else
			errx(EXIT_FAILURE, "%s(): unsupported file type in sketch directory: %s", __func__, child_src);
	}
	if (closedir(dir) != 0)
		err(errno, "%s(): cannot close %s", __func__, src);
}

static void save_generated_query_sketch_if_requested(const ani_opt_t *opt, const char *query_sketch_dir)
{
	if (!opt || opt->save_query_sketch[0] == '\0')
		return;
	ensure_save_query_sketch_path_available(opt->save_query_sketch);
	copy_tree_to_new_dir(query_sketch_dir, opt->save_query_sketch);
}

static void sketch_sequence_arg_to_dir(
	const char *input,
	const minco_sketch_stat_t *template_stat,
	const minco_sketch_info_t *template_info,
	const char *outdir,
	const ani_opt_t *ani_opt_val,
	const ani_ref_density_threshold_t *density_threshold)
{
	char *arg = (char *)input;
	int fmt_ck = (ani_opt_val->sketch_pipecmd[0] == '\0' && strcmp(input, "-") != 0) ? 1 : 0;
	infile_tab_t *tab = organize_infile_frm_arg(1, &arg, fmt_ck);
	if (!tab || tab->infile_num < 1)
		errx(EXIT_FAILURE, "%s is not a valid FASTA/FASTQ input for auto ANI", input);

	int infile_num = ani_opt_val->sketch_asone ? (tab->infile_num < 1 ? 0 : 1) : tab->infile_num;
	const uint32_t target_sketch_size =
		(template_info && template_info->target_sketch_size)
			? template_info->target_sketch_size
			: (ani_opt_val->sketch_size ? ani_opt_val->sketch_size : MINCO_DEFAULT_SKETCH_SIZE);
	uint32_t selection_mode = MINCO_STAT_SELECTION_BOTTOMK;
	uint64_t stat_density_threshold = UINT64_MAX;
	if (density_threshold && density_threshold->enabled)
	{
		selection_mode = MINCO_STAT_SELECTION_DENSITY_THRESHOLD;
		stat_density_threshold = density_threshold->threshold;
	}
	apply_auto_sketch_stat(template_stat, infile_num, target_sketch_size,
						   selection_mode, stat_density_threshold);
	sketch_opt_t opt = sketch_opt_from_stat(template_stat, outdir, ani_opt_val,
											target_sketch_size);
	if (density_threshold && density_threshold->enabled)
	{
		opt.density_threshold_enabled = true;
		opt.density_threshold = density_threshold->threshold;
	}
	mkdir_p(outdir);
	compute_sketch(&opt, tab);
	free_infile_tab(tab);
}

static char *merge_auto_sketches(path_vec_t *sketch_dirs, const char *outdir)
{
	if (sketch_dirs->n == 1)
	{
		char *copy = strdup(sketch_dirs->items[0]);
		if (!copy)
			err(errno, "%s(): OOM sketch path copy", __func__);
		return copy;
	}

	sketch_opt_t merge_opt = {
		.outdir = (char *)outdir,
		.num_remaining_args = sketch_dirs->n,
		.remaining_args = sketch_dirs->items,
	};
	mkdir_p(outdir);
	merge_minco_sketches(&merge_opt);
	char *copy = strdup(outdir);
	if (!copy)
		err(errno, "%s(): OOM merged sketch path copy", __func__);
	return copy;
}

static bool paths_refer_to_same_file(const char *a, const char *b)
{
	struct stat sa;
	struct stat sb;
	return a && b &&
		   stat(a, &sa) == 0 &&
		   stat(b, &sb) == 0 &&
		   sa.st_dev == sb.st_dev &&
		   sa.st_ino == sb.st_ino;
}

static int run_sketch_ani(ani_opt_t *opt)
{
	if (opt->unassembled && !opt->afcut_set)
		opt->afcut = 0.2f;

	size_t file_size;
	char *qry_stat_path = test_get_fullpath(opt->qrydir, sketch_stat);
	void *qry_stat_mem = read_from_file(qry_stat_path, &file_size);
	free(qry_stat_path);
	minco_sketch_stat_t qry_sketch_stat_val = {0};
	minco_sketch_info_t qry_info = {0};
	if (!minco_stat_decode_mem(qry_stat_mem, file_size, &qry_sketch_stat_val, &qry_info))
		errx(EINVAL, "%s(): malformed %s/%s", __func__, opt->qrydir, sketch_stat);
	if (qry_info.has_minco_ext)
		qry_sketch_stat_val.hash_id = qry_info.sketch_id;
	minco_sketch_stat_t *qry_sketch_stat = &qry_sketch_stat_val;
	const int qry_infile_num = qry_sketch_stat->infile_num;

	ani_model_target_sketch_size = qry_info.target_sketch_size
										? qry_info.target_sketch_size
										: ANI_MODEL_REFERENCE_SKETCH_SIZE;
	ani_model_compat_filter_shift = qry_sketch_stat->compat_filter_shift;
	const_comask_init(qry_sketch_stat);
	const bool force_ref_index = force_ref_index_requested();
	const int matrix_direct_n = opt->fmt ? ani_matrix_direct_threshold() : 0;
	if (opt->refdir[0] == '\0')
	{
		if (opt->fmt == 0)
			errx(EXIT_FAILURE, "one-sketch ANI detail output is not supported; use -m1 full matrix or -m2 triangle");
		if (opt->fmt == 1 &&
			(force_ref_index || ani_matrix_should_use_index(qry_infile_num, qry_infile_num, matrix_direct_n)))
		{
			if (force_ref_index && !file_exists_in_folder(opt->qrydir, sorted_comb_ctxgid64obj32))
				errx(EXIT_FAILURE, "MINCO_FORCE_REF_INDEX requires reference index '%s/%s'",
					 opt->qrydir, sorted_comb_ctxgid64obj32);
			else
				ensure_sorted_ref_index_for_ani_m1(opt->qrydir);
			free_read_from_file(qry_stat_mem, file_size);
			if (comb_sortedsketch64_indexed_self_full(opt))
				return 1;
			comb_sortedsketch64_self_matrix(opt);
			return 1;
		}
		if (opt->fmt == 2 &&
			(force_ref_index || ani_matrix_should_use_index(qry_infile_num, qry_infile_num, matrix_direct_n)))
		{
			free_read_from_file(qry_stat_mem, file_size);
			if (comb_sortedsketch64_indexed_self_triangle(opt))
				return 1;
			comb_sortedsketch64_self_matrix(opt);
			return 1;
		}
		free_read_from_file(qry_stat_mem, file_size);
		comb_sortedsketch64_self_matrix(opt);
		return 1;
	}
	if (opt->fmt && paths_refer_to_same_file(opt->refdir, opt->qrydir))
	{
		if (opt->fmt == 1 &&
			(force_ref_index || ani_matrix_should_use_index(qry_infile_num, qry_infile_num, matrix_direct_n)))
		{
			if (force_ref_index && !file_exists_in_folder(opt->qrydir, sorted_comb_ctxgid64obj32))
				errx(EXIT_FAILURE, "MINCO_FORCE_REF_INDEX requires reference index '%s/%s'",
					 opt->qrydir, sorted_comb_ctxgid64obj32);
			else
				ensure_sorted_ref_index_for_ani_m1(opt->qrydir);
			free_read_from_file(qry_stat_mem, file_size);
			if (comb_sortedsketch64_indexed_self_full(opt))
				return 1;
			comb_sortedsketch64_self_matrix(opt);
			return 1;
		}
		if (opt->fmt == 2 &&
			(force_ref_index || ani_matrix_should_use_index(qry_infile_num, qry_infile_num, matrix_direct_n)))
		{
			free_read_from_file(qry_stat_mem, file_size);
			if (comb_sortedsketch64_indexed_self_triangle(opt))
				return 1;
			comb_sortedsketch64_self_matrix(opt);
			return 1;
		}
		free_read_from_file(qry_stat_mem, file_size);
		comb_sortedsketch64_self_matrix(opt);
		return 1;
	}
	if (opt->fmt == 2)
		errx(EXIT_FAILURE, "ANI triangle output requires one sketch or identical -r/-q sketches");
	const int ref_infile_num = sketch_infile_num(opt->refdir);
	const bool use_small_query_stream =
		!force_ref_index && opt->fmt == 0 && ref_infile_num > 0 &&
		qry_infile_num < ref_infile_num;
	if (use_small_query_stream && qry_infile_num == 1)
	{
		free_read_from_file(qry_stat_mem, file_size);
		return stream_ref_sketches_one_qraw_lookup(opt);
	}
	const char *multi_qraw_mode = getenv("MINCO_QRAW_MULTI");
	if (opt->fmt == 0 && qry_infile_num > 1)
	{
		bool use_query_index = use_small_query_stream &&
							   qry_infile_num <= qraw_multi_query_index_threshold();
		if (opt->unassembled && !force_ref_index && multi_qraw_mode && multi_qraw_mode[0] != '\0')
		{
			if (strcmp(multi_qraw_mode, "sortedindex") == 0)
				use_query_index = true;
			else if (strcmp(multi_qraw_mode, "legacy") == 0 ||
					 strcmp(multi_qraw_mode, "refindex") == 0)
				use_query_index = false;
			else
				errx(EXIT_FAILURE, "invalid MINCO_QRAW_MULTI='%s' (use sortedindex, refindex, or legacy)",
					 multi_qraw_mode);
		}
		if (use_query_index)
		{
			free_read_from_file(qry_stat_mem, file_size);
			return stream_ref_sketches_multi_qraw_sortedindex(opt);
		}
	}
	free_read_from_file(qry_stat_mem, file_size);
	const bool ref_index_exists = file_exists_in_folder(opt->refdir, sorted_comb_ctxgid64obj32);
	if (opt->fmt == 1)
	{
		if (force_ref_index)
		{
			if (!ref_index_exists)
				errx(EXIT_FAILURE, "MINCO_FORCE_REF_INDEX requires reference index '%s/%s'",
					 opt->refdir, sorted_comb_ctxgid64obj32);
			return mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(opt);
		}
		if (ani_matrix_should_use_index(ref_infile_num, qry_infile_num, matrix_direct_n))
		{
			ensure_sorted_ref_index_for_ani_m1(opt->refdir);
			return mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(opt);
		}
	}
	else if (ref_index_exists)
		return mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(opt);
	const int auto_index_n = auto_ref_index_threshold();
	if (!opt->v && !opt->unassembled && opt->fmt == 0 &&
		auto_index_n > 0 && ref_infile_num > auto_index_n)
	{
		gen_inverted_index_for_minco(opt->refdir);
		return mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(opt);
	}
	if (force_ref_index)
		errx(EXIT_FAILURE, "MINCO_FORCE_REF_INDEX requires reference index '%s/%s'",
			 opt->refdir, sorted_comb_ctxgid64obj32);

	if (opt->fmt)
	{
		comb_sortedsketch64Xcomb_sortedsketch64(opt);
		return 1;
	}
	comb_manysmall_sortedsketch64Xcomb_fewlarge_sortedsketch64_filter_and_sort_survivors(opt);
	return 1;
}

static minco_sketch_info_t default_auto_ani_sketch_info(const minco_sketch_stat_t *stat,
														const ani_opt_t *opt)
{
	const uint32_t target = opt->sketch_size ? opt->sketch_size : MINCO_DEFAULT_SKETCH_SIZE;
	minco_stat_ext_v1_t ext = minco_stat_make_ext(stat, target,
												  MINCO_STAT_SELECTION_BOTTOMK,
												  minco_compiled_stat_flags(),
												  UINT64_MAX);
	return (minco_sketch_info_t){
		.has_minco_ext = true,
		.stat_version = ext.version,
		.target_sketch_size = ext.target_sketch_size,
		.feature_id = ext.feature_id,
		.sketch_id = ext.sketch_id,
		.hash_function = ext.hash_function,
		.hash_bits = ext.hash_bits,
		.hash_seed = ext.hash_seed,
		.sketch_model = ext.sketch_model,
		.selection_mode = ext.selection_mode,
		.flags = ext.flags,
		.density_threshold = ext.density_threshold,
		.density_min_threshold = ext.density_min_threshold,
		.density_max_threshold = ext.density_max_threshold,
		.density_universal_threshold = ext.density_universal_threshold,
		.density_min_sample_id = ext.density_min_sample_id,
		.density_max_sample_id = ext.density_max_sample_id,
		.density_universal_sample_id = ext.density_universal_sample_id,
		.density_valid_sample_count = ext.density_valid_sample_count,
		.density_hash_bits = ext.density_hash_bits,
		.density_universal_policy = ext.density_universal_policy,
		.density_flags = ext.density_flags,
	};
}

static bool select_auto_template_from_inputs(const path_vec_t *inputs,
											 minco_sketch_stat_t *template_stat,
											 minco_sketch_info_t *template_info)
{
	for (int i = 0; i < inputs->n; ++i)
	{
		const char *path = inputs->items[i];
		if (!is_minco_sketch_dir(path))
			continue;
		*template_stat = read_minco_sketch_stat_info(path, template_info);
		validate_supported_auto_sketch_stat(template_stat);
		return true;
	}
	return false;
}

static minco_sketch_stat_t select_auto_template_stat(const ani_opt_t *opt,
												   const path_vec_t *ref_inputs,
												   const path_vec_t *qry_inputs,
												   minco_sketch_info_t *template_info)
{
	minco_sketch_stat_t template_stat = default_auto_ani_sketch_stat(opt);
	if (template_info)
		*template_info = default_auto_ani_sketch_info(&template_stat, opt);
	bool have_template = select_auto_template_from_inputs(ref_inputs, &template_stat,
														 template_info);
	if (!have_template)
		have_template = select_auto_template_from_inputs(qry_inputs, &template_stat,
														 template_info);
	if (!have_template)
		validate_supported_auto_sketch_stat(&template_stat);
	if (template_info && template_info->target_sketch_size == 0)
		*template_info = default_auto_ani_sketch_info(&template_stat, opt);
	return template_stat;
}

static void validate_auto_sketch_inputs(const path_vec_t *inputs, const minco_sketch_stat_t *template_stat)
{
	for (int i = 0; i < inputs->n; ++i)
	{
		const char *path = inputs->items[i];
		if (!is_minco_sketch_dir(path))
			continue;
		minco_sketch_stat_t stat = read_minco_sketch_stat(path);
		if (!same_ani_auto_sketch_stat(template_stat, &stat))
			describe_sketch_stat_mismatch(path, template_stat, &stat);
	}
}

static int path_vec_stdin_count(const path_vec_t *inputs)
{
	int count = 0;
	for (int i = 0; i < inputs->n; ++i)
		if (strcmp(inputs->items[i], "-") == 0)
			count++;
	return count;
}

static bool auto_input_is_streamed_raw(const ani_opt_t *opt, const char *path)
{
	return strcmp(path, "-") == 0 || opt->sketch_pipecmd[0] != '\0';
}

static bool auto_input_is_fastq_path(const char *path)
{
	if (strcmp(path, "-") == 0 || is_minco_sketch_dir(path))
		return false;
	return isOK_fmt_infile((char *)path, fastq_fmt, FQ_FMT_SZ);
}

static bool auto_query_side_is_unassembled(const ani_opt_t *opt, const path_vec_t *qry_inputs)
{
	for (int i = 0; i < qry_inputs->n; ++i)
	{
		const char *path = qry_inputs->items[i];
		if (is_minco_sketch_dir(path))
			continue;
		if (auto_input_is_fastq_path(path))
			return true;
		if (opt->sketch_conflict && auto_input_is_streamed_raw(opt, path))
			return true;
	}
	return false;
}

static bool path_vec_has_sequence_inputs(const path_vec_t *inputs)
{
	for (int i = 0; i < inputs->n; ++i)
		if (!is_minco_sketch_dir(inputs->items[i]))
			return true;
	return false;
}

static void prepare_auto_input_side(
	const path_vec_t *inputs,
	path_vec_t *sketches,
	const char *side_prefix,
	const minco_sketch_stat_t *template_stat,
	const minco_sketch_info_t *template_info,
	const char *tmp_root,
	const ani_opt_t *ani_opt_val,
	const ani_ref_density_threshold_t *density_threshold)
{
	for (int i = 0; i < inputs->n; ++i)
	{
		const char *input = inputs->items[i];
		if (is_minco_sketch_dir(input))
		{
			path_vec_push(sketches, input);
			continue;
		}

		char outdir[PATHLEN * 2];
		int written = snprintf(outdir, sizeof(outdir), "%s/%sseq_%d", tmp_root, side_prefix, i);
		if (written < 0 || (size_t)written >= sizeof(outdir))
			errx(EXIT_FAILURE, "auto ANI temporary output path is too long for %s", input);
		sketch_sequence_arg_to_dir(input, template_stat, template_info, outdir,
								   ani_opt_val, density_threshold);
		path_vec_push(sketches, outdir);
	}
}

static void copy_runtime_path(char *dest, size_t dest_size, const char *path, const char *label)
{
	if (strlen(path) >= dest_size)
		errx(EXIT_FAILURE, "%s path is too long; maximum supported length is %zu bytes", label, dest_size - 1);
	snprintf(dest, dest_size, "%s", path);
}

static int run_auto_ani_from_sides(ani_opt_t *opt, const path_vec_t *ref_inputs, const path_vec_t *qry_inputs)
{
	if (ref_inputs->n < 1)
		errx(EXIT_FAILURE, "auto ANI requires at least one reference input");
	if (qry_inputs->n < 1)
		errx(EXIT_FAILURE, "auto ANI requires at least one query input");
	if (path_vec_stdin_count(ref_inputs) + path_vec_stdin_count(qry_inputs) > 1)
		errx(EXIT_FAILURE, "stdin input '-' can be used only once in auto ANI");
	const bool qry_has_sequence_inputs = path_vec_has_sequence_inputs(qry_inputs);
	if (opt->save_query_sketch[0] != '\0')
	{
		if (!qry_has_sequence_inputs)
			errx(EXIT_FAILURE, "--save-query-sketch requires at least one FASTA/FASTQ query input");
		ensure_save_query_sketch_path_available(opt->save_query_sketch);
	}

	minco_sketch_info_t template_info = {0};
	minco_sketch_stat_t template_stat = select_auto_template_stat(opt, ref_inputs, qry_inputs,
																&template_info);
	validate_auto_sketch_inputs(ref_inputs, &template_stat);
	validate_auto_sketch_inputs(qry_inputs, &template_stat);

	char tmp_template[] = "/tmp/minco_ani_auto_XXXXXX";
	char *tmp_root = mkdtemp(tmp_template);
	if (!tmp_root)
		err(errno, "%s(): mkdtemp", __func__);

	path_vec_t ref_sketches = {0};
	path_vec_t qry_sketches = {0};
	prepare_auto_input_side(ref_inputs, &ref_sketches, "ref", &template_stat,
							&template_info, tmp_root, opt, NULL);
	if (ref_sketches.n < 1)
		errx(EXIT_FAILURE, "auto ANI failed to prepare reference sketches");

	char ref_merge[PATHLEN * 2];
	char qry_merge[PATHLEN * 2];
	int written = snprintf(ref_merge, sizeof(ref_merge), "%s/ref_merged", tmp_root);
	if (written < 0 || (size_t)written >= sizeof(ref_merge))
		errx(EXIT_FAILURE, "auto ANI reference merge path is too long");
	char *refdir = merge_auto_sketches(&ref_sketches, ref_merge);

	ani_ref_density_threshold_t query_density_threshold = {0};
	if (opt->query_density_model == ANI_QUERY_DENSITY_REF &&
		qry_has_sequence_inputs)
		query_density_threshold = ani_select_ref_density_threshold(refdir, &template_stat);

	prepare_auto_input_side(qry_inputs, &qry_sketches, "qry", &template_stat,
							&template_info, tmp_root, opt,
							query_density_threshold.enabled ? &query_density_threshold : NULL);
	if (qry_sketches.n < 1)
		errx(EXIT_FAILURE, "auto ANI failed to prepare query sketches");

	written = snprintf(qry_merge, sizeof(qry_merge), "%s/qry_merged", tmp_root);
	if (written < 0 || (size_t)written >= sizeof(qry_merge))
		errx(EXIT_FAILURE, "auto ANI query merge path is too long");
	char *qrydir = merge_auto_sketches(&qry_sketches, qry_merge);
	save_generated_query_sketch_if_requested(opt, qrydir);

	char old_refdir[PATHLEN];
	char old_qrydir[PATHLEN];
	const bool old_unassembled = opt->unassembled;
	const bool old_v = opt->v;
	get_generic_dist_from_features_fn old_dist_fn = get_generic_dist_from_features;
	bool auto_unassembled = auto_query_side_is_unassembled(opt, qry_inputs);
	snprintf(old_refdir, sizeof(old_refdir), "%s", opt->refdir);
	snprintf(old_qrydir, sizeof(old_qrydir), "%s", opt->qrydir);
	copy_runtime_path(opt->refdir, sizeof(opt->refdir), refdir, "auto ANI reference sketch");
	copy_runtime_path(opt->qrydir, sizeof(opt->qrydir), qrydir, "auto ANI query sketch");
	opt->unassembled = auto_unassembled;
	if (auto_unassembled)
	{
		opt->v = true;
		get_generic_dist_from_features = get_naive_dist;
	}

	int rc = run_sketch_ani(opt);

	snprintf(opt->refdir, sizeof(opt->refdir), "%s", old_refdir);
	snprintf(opt->qrydir, sizeof(opt->qrydir), "%s", old_qrydir);
	opt->unassembled = old_unassembled;
	opt->v = old_v;
	get_generic_dist_from_features = old_dist_fn;
	free(refdir);
	free(qrydir);
	path_vec_free(&ref_sketches);
	path_vec_free(&qry_sketches);
	remove_tree(tmp_root);
	return rc;
}

static int run_auto_positional_ani(ani_opt_t *opt)
{
	if (opt->num_remaining_args < 2)
		errx(EXIT_FAILURE, "auto ANI requires one reference input followed by at least one query input");

	path_vec_t ref_inputs = {0};
	path_vec_t qry_inputs = {0};
	path_vec_push(&ref_inputs, opt->remaining_args[0]);
	for (int i = 1; i < opt->num_remaining_args; ++i)
		path_vec_push(&qry_inputs, opt->remaining_args[i]);

	int rc = run_auto_ani_from_sides(opt, &ref_inputs, &qry_inputs);
	path_vec_free(&ref_inputs);
	path_vec_free(&qry_inputs);
	return rc;
}

static int run_auto_list_ani(ani_opt_t *opt)
{
	path_vec_t ref_inputs = {0};
	path_vec_t qry_inputs = {0};
	read_path_list_to_vec(opt->reflist, &ref_inputs, "--reflist");
	read_path_list_to_vec(opt->qrylist, &qry_inputs, "--qrylist");

	int rc = run_auto_ani_from_sides(opt, &ref_inputs, &qry_inputs);
	path_vec_free(&ref_inputs);
	path_vec_free(&qry_inputs);
	return rc;
}

static bool path_is_directory(const char *path)
{
	struct stat st;
	return path && stat(path, &st) == 0 && S_ISDIR(st.st_mode);
}

static bool direct_query_input_is_fastq_like(const ani_opt_t *opt, const char *path)
{
	if (!path)
		return false;
	if (strcmp(path, "-") == 0)
		return true;
	if (opt && opt->sketch_pipecmd[0] != '\0')
		return true;
	return isOK_fmt_infile((char *)path, fastq_fmt, FQ_FMT_SZ);
}

static int run_direct_query_sequence_ani(ani_opt_t *opt, bool force_raw_query, const char *option_name)
{
	if (!is_minco_sketch_dir(opt->refdir))
		errx(EXIT_FAILURE,
			 "%s sequence input requires -r/--ref to be a sketch directory; "
			 "for direct reference sequence input use positional auto ANI",
			 option_name);
	if (path_is_directory(opt->qrydir))
		errx(EXIT_FAILURE, "%s input '%s' is a directory but not a valid sketch", option_name, opt->qrydir);
	ensure_save_query_sketch_path_available(opt->save_query_sketch);

	minco_sketch_info_t qry_template_info = {0};
	minco_sketch_stat_t qry_template = read_minco_sketch_stat_info(opt->refdir,
															   &qry_template_info);
	validate_supported_auto_sketch_stat(&qry_template);

	char input[PATHLEN];
	char old_qrydir[PATHLEN];
	snprintf(input, sizeof(input), "%s", opt->qrydir);
	snprintf(old_qrydir, sizeof(old_qrydir), "%s", opt->qrydir);
	const bool raw_query = force_raw_query;
	const bool keep_query_conflicts = force_raw_query || opt->sketch_conflict;
	if (keep_query_conflicts)
		qry_template.conflict = true;

	ani_ref_density_threshold_t query_density_threshold = {0};
	if (opt->query_density_model == ANI_QUERY_DENSITY_REF)
		query_density_threshold = ani_select_ref_density_threshold(opt->refdir, &qry_template);

	const bool use_readwise_density =
		raw_query &&
		query_density_threshold.enabled &&
		opt->save_query_sketch[0] == '\0' &&
		!opt->sketch_reads_qc &&
		!opt->sketch_abundance &&
		direct_query_input_is_fastq_like(opt, input);
	if (use_readwise_density)
	{
		return stream_fastq_query_readwise_density_ani(opt, input, query_density_threshold.threshold);
	}
	if (opt->abundance_model != ANI_ABUNDANCE_NONE)
		errx(EXIT_FAILURE,
			 "--abundance-est currently requires direct --qraw FASTQ --query-density ref "
			 "without --save-query-sketch, --readsQC, or --abundance");

	char tmp_template[] = "/tmp/minco_ani_query_XXXXXX";
	char *tmp_root = mkdtemp(tmp_template);
	if (!tmp_root)
		err(errno, "%s(): mkdtemp", __func__);

	char qry_out[PATHLEN * 2];
	int written = snprintf(qry_out, sizeof(qry_out), "%s/qry_seq", tmp_root);
	if (written < 0 || (size_t)written >= sizeof(qry_out))
		errx(EXIT_FAILURE, "temporary %s sketch output path is too long", option_name);

	sketch_sequence_arg_to_dir(input, &qry_template, &qry_template_info, qry_out, opt,
							   query_density_threshold.enabled ? &query_density_threshold : NULL);
	save_generated_query_sketch_if_requested(opt, qry_out);

	const bool old_unassembled = opt->unassembled;
	const bool old_v = opt->v;
	char old_refdir[PATHLEN];
	char old_outf[PATHLEN];
	snprintf(old_refdir, sizeof(old_refdir), "%s", opt->refdir);
	snprintf(old_outf, sizeof(old_outf), "%s", opt->outf);
	get_generic_dist_from_features_fn old_dist_fn = get_generic_dist_from_features;
	copy_runtime_path(opt->qrydir, sizeof(opt->qrydir), qry_out, "temporary query sketch");
	opt->unassembled = raw_query;
	if (raw_query)
	{
		opt->v = true;
		get_generic_dist_from_features = get_naive_dist;
	}

	int rc = run_sketch_ani(opt);

	snprintf(opt->qrydir, sizeof(opt->qrydir), "%s", old_qrydir);
	snprintf(opt->refdir, sizeof(opt->refdir), "%s", old_refdir);
	snprintf(opt->outf, sizeof(opt->outf), "%s", old_outf);
	opt->unassembled = old_unassembled;
	opt->v = old_v;
	get_generic_dist_from_features = old_dist_fn;
	remove_tree(tmp_root);
	return rc;
}

int cmd_ani(struct argp_state *state)
{

	struct arg_ani ani = {
		0,
	};
	int argc = state->argc - state->next + 1;
	char **argv = &state->argv[state->next - 1];
	ani.global = state->input;
	argp_parse(&argp_ani, argc, argv, ARGP_IN_ORDER, &argc, &ani);
	state->next += argc - 1;

	// instantilize distance fn according selection of naive model or not
	ani_model_compat_filter_shift = ani_opt.sketch_filter_shift;
	ani_model_target_sketch_size = ani_opt.sketch_size ? ani_opt.sketch_size : ANI_MODEL_REFERENCE_SKETCH_SIZE;
	get_generic_dist_from_features = ani_opt.v ? get_naive_dist : lm3ways_dist_from_features;

	if (ani_opt.abundance_model != ANI_ABUNDANCE_NONE &&
		(ani_opt.reflist[0] != '\0' || ani_opt.qrylist[0] != '\0' ||
		 ani_opt.num_remaining_args >= 2))
		errx(EXIT_FAILURE,
			 "--abundance-est currently supports only direct -r REF --qraw FASTQ --query-density ref");

	if (ani_opt.reflist[0] != '\0' || ani_opt.qrylist[0] != '\0')
		return run_auto_list_ani(&ani_opt);
	else if (ani_opt.qrydir[0] != '\0')
	{
		const bool qry_is_sketch = is_minco_sketch_dir(ani_opt.qrydir);
		if (ani_opt.abundance_model != ANI_ABUNDANCE_NONE && qry_is_sketch)
			errx(EXIT_FAILURE,
				 "--abundance-est currently requires a direct FASTQ query, not an existing query sketch");
		if (ani_opt.save_query_sketch[0] != '\0' && qry_is_sketch)
			errx(EXIT_FAILURE, "--save-query-sketch requires a FASTA/FASTQ query input; query sketch %s is already persistent",
				 ani_opt.qrydir);
		if (!qry_is_sketch)
			return run_direct_query_sequence_ani(&ani_opt,
												 ani_opt.unassembled,
												 ani_opt.unassembled ? "--qraw" : "-q/--query");
		return run_sketch_ani(&ani_opt);
	}
	else if (ani_opt.num_remaining_args >= 2)
		return run_auto_positional_ani(&ani_opt);
	return 1;
}
