#include "command_ani.h"
#include "command_matrix.h"
#include "global_basic.h"
#include "minco_sort.h"
#include "pairwise_graph.h"
#include "sketch_rearrange.h"
// #include "command_sketch.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <err.h>
#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <libgen.h>
#include <dirent.h>
#include <omp.h>
#include <stdatomic.h>
#include <ctype.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif
#include "../klib/kstring.h" // from klib
#include "../klib/khash.h"
#include "../klib/kvec.h"      // klib: dynamic arrays
#include "../klib/kseq.h"
// #include "../klib/khashl.h"
KSEQ_INIT(gzFile, gzread)
#define GID_NBITS 20 // 2^20, 1M
#define CONFLICT_OBJ UINT32_MAX
// pulic vars
const char gid_obj_prefix[] = "gidobj", ctx_idx_prefix[] = "ctx.index";
static const char minco_ctxmeta_legacy_tsv_stat[] = "minco.ctxmeta.tsv";
extern const char sorted_comb_ctxgid64obj32[];
extern uint32_t FILTER;
extern void gen_inverted_index_for_minco(const char *sketchdir);
extern double C9O7_98[6], C9O7_96[6];
size_t file_size;

const char unified_detail_header[] = "Qry\tRef\tANI\tDistance\tConfidence\tSelected_metric\tXnY_ctx\tQry_align_fraction\tblastn_Qry_align_fraction\tRef_align_fraction\tblastn_Ref_align_fraction\tN_diff_obj\tN_diff_obj_section\tN_mut2_ctx\tRef_annotation\tReal_Qry_align_fraction\tReal_Ref_align_fraction\tReal_min_align_fraction\tAF_source";
const char readwise_detail_extra_header[] = "Reads_with_ctx_match\tTotal_reads\tRead_match_fraction\tUnique_query_ctx\tUnique_query_ctx_hit\tUnique_ref_ctx_hit\tDensity_block_ctx\tTotal_density_blocks\tBlocks_with_ctx_match\tBlock_match_fraction\tRaw_XnY_ctx\tRejected_ctx\tRejected_diff_ctx\tFake_ctx_fraction\tFake_ctx_prob_mean\tFake_ctx_prob_weighted";
const char readwise_abundance_extra_header[] = "Ref_breadth\tRef_mean_depth\tRef_hit_mean_depth\tRef_depth_variance\tRef_depth_cv\tRef_zero_fraction\tRelative_abundance_depth\tNormalized_abundance_depth\tEffective_abundance_depth\tNormalized_effective_abundance_depth\tRef_zip_af\tRef_zip_aaf_ani\tReliable_Ref_breadth\tReliable_Ref_mean_depth\tReliable_Ref_hit_ctx\tReliable_Ref_hit_mean_depth\tReliable_Ref_hit_median_depth\tReliable_Ref_hit_depth_variance\tReliable_Ref_zip_af\tDefault_call\tDefault_call_rule";
const char readwise_marker_extra_header[] = "Marker_XnY_ctx\tMarker_Raw_XnY_ctx\tMarker_N_diff_obj\tMarker_N_diff_obj_section\tMarker_N_mut2_ctx\tMarker_Ref_ctx_total\tctx_marker_size\tMarker_Ref_breadth\tMarker_Ref_mean_depth\tMarker_Ref_hit_mean_depth\tMarker_Ref_depth_variance\tMarker_Ref_depth_cv\tMarker_Ref_zero_fraction\tMarker_Relative_abundance_depth\tMarker_Effective_abundance_depth\tMarker_Ref_zip_af\tMarker_Ref_zip_aaf_ani\tMarker_Reliable_Ref_breadth\tMarker_Reliable_Ref_mean_depth\tMarker_Reliable_Ref_hit_ctx\tMarker_Reliable_Ref_hit_mean_depth\tMarker_Reliable_Ref_hit_median_depth\tMarker_Reliable_Ref_hit_depth_variance\tMarker_Reliable_Ref_zip_af";

#define MINCO_READWISE_EFFECTIVE_MEDIAN_DEPTH_CUTOFF 20.0
#define MINCO_READWISE_EFFECTIVE_AF_EXPONENT 1.05
#define ANI_SELECTED_METRIC_COUNT 8
const char select_metrics_header[ANI_SELECTED_METRIC_COUNT][20] = {
	"BestDist", "RecalDist", "CtxMoE", "Naive",
	"MashD", "AafD", "MashD_if_far", "AafD_if_far"
};

typedef enum {
	ANI_CONF_LOW = 0,
	ANI_CONF_MEDIUM = 1,
	ANI_CONF_HIGH = 2
} ani_confidence_t;

static inline const char *ani_confidence_label(unsigned char confidence)
{
	switch ((ani_confidence_t)confidence) {
	case ANI_CONF_HIGH:
		return "high";
	case ANI_CONF_MEDIUM:
		return "medium";
	case ANI_CONF_LOW:
	default:
		return "low";
	}
}

static inline unsigned char ani_confidence_from_values(double raw_ani, double calibrated_ani, double ref_af)
{
	const double drop = raw_ani - calibrated_ani;
	if (ref_af < 0.5 || (raw_ani >= 0.95 && drop > 0.02))
		return ANI_CONF_LOW;
	if (fabs(drop) > 0.01)
		return ANI_CONF_MEDIUM;
	return ANI_CONF_HIGH;
}

static inline const infile_meta_t *infile_meta_at(const unify_sketch_t *sketch,
															uint32_t idx)
{
	if (!sketch || !sketch->infile_meta || idx >= (uint32_t)sketch->infile_num)
		return NULL;
	return &sketch->infile_meta[idx];
}

static inline bool infile_meta_complete_like_assembly(const infile_meta_t *s)
{
	if (!s || s->meta_fmt_version != MINCO_INFILE_META_VERSION ||
		s->infile_fmt != MINCO_INFILE_FMT_FASTA ||
		s->total_length_bp == 0 || s->median_length_bp == 0)
		return false;
	return s->total_length_bp >= 3500000ULL &&
		   s->record_count < 5 &&
		   s->median_length_bp > 100000U &&
		   (double)s->total_length_bp / (double)s->median_length_bp < 5.0;
}

static inline bool has_assembly_meta_record(const infile_meta_t *s)
{
	return s && s->meta_fmt_version == MINCO_INFILE_META_VERSION &&
		   s->infile_fmt == MINCO_INFILE_FMT_FASTA &&
		   s->total_length_bp > 0 && s->record_count > 0;
}

static inline double bounded_ani(double value)
{
	if (!isfinite(value))
		return NAN;
	if (value < 0.0)
		return 0.0;
	if (value > 1.0)
		return 1.0;
	return value;
}

static inline double bounded_align_fraction(double value)
{
	if (!isfinite(value))
		return value;
	if (value < 0.0)
		return 0.0;
	if (value > 1.0)
		return 1.0;
	return value;
}

static inline double ctx_k_for_ani(void)
{
	return Bitslen.ctx > 0 ? (double)Bitslen.ctx / 2.0 : 22.0;
}

static inline double mash_ani_from_counts(double overlap, double qry_ctx, double ref_ctx)
{
	if (overlap <= 0.0 || qry_ctx <= 0.0 || ref_ctx <= 0.0)
		return 0.0;
	const double denom = qry_ctx + ref_ctx - overlap;
	if (denom <= 0.0)
		return 0.0;
	const double jaccard = overlap / denom;
	if (jaccard <= 0.0)
		return 0.0;
	return bounded_ani(1.0 + log(2.0 * jaccard / (1.0 + jaccard)) / ctx_k_for_ani());
}

static inline double aaf_ani_from_counts(double overlap, double qry_ctx, double ref_ctx)
{
	if (overlap <= 0.0 || qry_ctx <= 0.0 || ref_ctx <= 0.0)
		return 0.0;
	const double min_ctx = qry_ctx < ref_ctx ? qry_ctx : ref_ctx;
	if (min_ctx <= 0.0)
		return 0.0;
	const double containment = overlap / min_ctx;
	if (containment <= 0.0)
		return 0.0;
	return bounded_ani(1.0 + log(containment) / ctx_k_for_ani());
}

static inline double aaf_ani_from_containment(double containment)
{
	if (containment <= 0.0)
		return 0.0;
	if (containment > 1.0)
		containment = 1.0;
	return bounded_ani(1.0 + log(containment) / ctx_k_for_ani());
}

static inline double ani_zip_observed_breadth(double latent_af, double mean_depth)
{
	if (latent_af <= 0.0 || mean_depth <= 0.0)
		return 0.0;
	const double lambda = mean_depth / latent_af;
	if (lambda > 700.0)
		return latent_af;
	return latent_af * (1.0 - exp(-lambda));
}

static double ani_zip_corrected_af(double observed_breadth, double mean_depth)
{
	if (!isfinite(observed_breadth) || observed_breadth <= 0.0)
		return 0.0;
	if (observed_breadth >= 1.0)
		return 1.0;
	if (!isfinite(mean_depth) || mean_depth <= 0.0)
		return observed_breadth;
	double lo = observed_breadth;
	double hi = 1.0;
	const double f_lo = ani_zip_observed_breadth(lo, mean_depth);
	const double f_hi = ani_zip_observed_breadth(hi, mean_depth);
	if (observed_breadth <= f_lo + 1e-12)
		return lo;
	if (observed_breadth >= f_hi - 1e-12)
		return hi;
	for (int i = 0; i < 48; ++i) {
		const double mid = 0.5 * (lo + hi);
		const double f_mid = ani_zip_observed_breadth(mid, mean_depth);
		if (f_mid < observed_breadth)
			lo = mid;
		else
			hi = mid;
	}
	return 0.5 * (lo + hi);
}

static inline void print_ani_detail_header(FILE *outfp, const ani_opt_t *ani_opt, bool include_selected_metric)
{
	(void)include_selected_metric;
	if (ani_opt && ani_opt->readwise_query)
	{
		fprintf(outfp, "%s\t%s", unified_detail_header, readwise_detail_extra_header);
		if (ani_opt->abundance_model != ANI_ABUNDANCE_NONE)
		{
			fprintf(outfp, "\t%s", readwise_abundance_extra_header);
			if (ani_opt->readwise_dual_evidence)
				fprintf(outfp, "\t%s", readwise_marker_extra_header);
		}
		fputc('\n', outfp);
	}
	else
		fprintf(outfp, "%s\n", unified_detail_header);
}

static inline const char *annotation_or_na(const char *annotation)
{
	return annotation && annotation[0] ? annotation : "NA";
}

static inline const char *annotation_at(char (*annotations)[PATHLEN], uint32_t idx)
{
	return annotation_or_na(annotations ? annotations[idx] : NULL);
}

static inline const char *unify_annotation_at(const unify_sketch_t *sketch, uint32_t idx)
{
	return annotation_at(sketch ? sketch->annotation : NULL, idx);
}

static char (*read_optional_sketch_annotations(const char *sketch_dir, int infile_num))[PATHLEN]
{
	if (infile_num <= 0 || !file_exists_in_folder(sketch_dir, sketch_anno_stat))
		return NULL;
	size_t anno_file_size = 0;
	char *anno_path = test_get_fullpath(sketch_dir, sketch_anno_stat);
	char (*annotations)[PATHLEN] =
		read_from_file(anno_path, &anno_file_size);
	free(anno_path);
	const size_t expected_size = (size_t)infile_num * PATHLEN;
	if (anno_file_size != expected_size)
		err(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
			__func__, sketch_dir, sketch_anno_stat, anno_file_size, expected_size);
	return annotations;
}

static infile_meta_t *read_optional_sketch_infile_meta_stats(const char *sketch_dir, int infile_num)
{
	if (infile_num <= 0 || !file_exists_in_folder(sketch_dir, sketch_infile_meta_stat))
		return NULL;
	size_t meta_file_size = 0;
	char *meta_path = test_get_fullpath(sketch_dir, sketch_infile_meta_stat);
	infile_meta_t *stats = read_from_file(meta_path, &meta_file_size);
	free(meta_path);
	const size_t expected_size = (size_t)infile_num * sizeof(stats[0]);
	if (meta_file_size != expected_size)
		err(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
			__func__, sketch_dir, sketch_infile_meta_stat, meta_file_size, expected_size);
	return stats;
}

static uint8_t *read_optional_domain_profiles(const char *sketch_dir, int infile_num)
{
	if (infile_num <= 0 || !file_exists_in_folder(sketch_dir, sketch_domain_stat))
		return NULL;
	size_t domain_file_size = 0;
	char *domain_path = test_get_fullpath(sketch_dir, sketch_domain_stat);
	uint8_t *profiles = read_from_file(domain_path, &domain_file_size);
	free(domain_path);
	const size_t expected_size = (size_t)infile_num * sizeof(profiles[0]);
	if (domain_file_size != expected_size)
		err(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
			__func__, sketch_dir, sketch_domain_stat, domain_file_size, expected_size);
	for (int i = 0; i < infile_num; ++i)
		if (profiles[i] > MINCO_DOMAIN_PROFILE_MAX)
			errx(EINVAL, "%s(): %s/%s has invalid profile %u at sample %d",
				 __func__, sketch_dir, sketch_domain_stat,
				 (unsigned)profiles[i], i);
	return profiles;
}

static int split_tsv_fields(char *line, char **fields, int max_fields)
{
	int n = 0;
	char *p = line;
	while (n < max_fields) {
		fields[n++] = p;
		char *tab = strchr(p, '\t');
		if (!tab)
			break;
		*tab = '\0';
		p = tab + 1;
	}
	return n;
}

static bool parse_u64_field(const char *s, uint64_t *out)
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

static ani_ctxmeta_rec_t *read_optional_ani_ctxmeta_stats(const char *sketch_dir, int infile_num)
{
	if (infile_num <= 0)
		return NULL;

	if (file_exists_in_folder(sketch_dir, minco_ctxmeta_bin_stat)) {
		char *ctxmeta_path = test_get_fullpath(sketch_dir, minco_ctxmeta_bin_stat);
		size_t ctxmeta_size = 0;
		minco_ctxmeta_record_t *records = read_from_file(ctxmeta_path, &ctxmeta_size);
		free(ctxmeta_path);
		const size_t expected_size = (size_t)infile_num * sizeof(records[0]);
		if (ctxmeta_size != expected_size)
			errx(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
				 __func__, sketch_dir, minco_ctxmeta_bin_stat,
				 ctxmeta_size, expected_size);
		ani_ctxmeta_rec_t *stats = calloc((size_t)infile_num, sizeof(stats[0]));
		if (!stats)
			err(EXIT_FAILURE, "%s(): calloc ctxmeta", __func__);
		for (int i = 0; i < infile_num; ++i) {
			stats[i] = (ani_ctxmeta_rec_t){
				.valid = records[i].valid,
				.hash_bits = records[i].hash_bits,
				.threshold = records[i].threshold,
				.sketch_entries = records[i].sketch_entries,
				.selected_observed_ctx = records[i].selected_observed_ctx,
				.selected_estimated_unique_ctx = records[i].selected_estimated_unique_ctx,
				.preconflict_observed_ctx = records[i].preconflict_observed_ctx,
				.preconflict_estimated_unique_ctx = records[i].preconflict_estimated_unique_ctx,
				.postconflict_observed_ctx = records[i].postconflict_observed_ctx,
				.postconflict_estimated_unique_ctx = records[i].postconflict_estimated_unique_ctx,
			};
		}
		free_read_from_file(records, ctxmeta_size);
		return stats;
	}

	if (!file_exists_in_folder(sketch_dir, minco_ctxmeta_legacy_tsv_stat))
		return NULL;

	char *ctxmeta_path = test_get_fullpath(sketch_dir, minco_ctxmeta_legacy_tsv_stat);
	FILE *fp = fopen(ctxmeta_path, "r");
	if (!fp)
		err(errno, "%s", ctxmeta_path);

	ani_ctxmeta_rec_t *stats = calloc((size_t)infile_num, sizeof(stats[0]));
	if (!stats)
		err(EXIT_FAILURE, "%s(): calloc ctxmeta", __func__);

	char *line = NULL;
	size_t cap = 0;
	ssize_t len = getline(&line, &cap, fp);
	if (len < 0) {
		free(line);
		fclose(fp);
		free(ctxmeta_path);
		free(stats);
		return NULL;
	}

	int line_no = 1;
	while ((len = getline(&line, &cap, fp)) >= 0) {
		++line_no;
		line[strcspn(line, "\r\n")] = '\0';
		if (line[0] == '\0')
			continue;

		char *fields[13] = {0};
		const int nfields = split_tsv_fields(line, fields, 13);
		if (nfields < 13) {
			warnx("%s(): skipping malformed %s line %d: expected 13 columns, found %d",
				  __func__, ctxmeta_path, line_no, nfields);
			continue;
		}

		uint64_t sample_id = 0, valid = 0, hash_bits = 0, threshold = 0;
		uint64_t sketch_entries = 0, selected_obs = 0, selected_est = 0;
		uint64_t pre_obs = 0, pre_est = 0, post_obs = 0, post_est = 0;
		if (!parse_u64_field(fields[0], &sample_id) ||
			!parse_u64_field(fields[3], &valid) ||
			!parse_u64_field(fields[4], &hash_bits) ||
			!parse_u64_field(fields[5], &threshold) ||
			!parse_u64_field(fields[6], &sketch_entries) ||
			!parse_u64_field(fields[7], &selected_obs) ||
			!parse_u64_field(fields[8], &selected_est) ||
			!parse_u64_field(fields[9], &pre_obs) ||
			!parse_u64_field(fields[10], &pre_est) ||
			!parse_u64_field(fields[11], &post_obs) ||
			!parse_u64_field(fields[12], &post_est)) {
			warnx("%s(): skipping malformed numeric fields in %s line %d",
				  __func__, ctxmeta_path, line_no);
			continue;
		}
		if (sample_id >= (uint64_t)infile_num) {
			warnx("%s(): skipping out-of-range sample_id %" PRIu64 " in %s line %d",
				  __func__, sample_id, ctxmeta_path, line_no);
			continue;
		}

		stats[sample_id] = (ani_ctxmeta_rec_t){
			.valid = (uint8_t)(valid != 0),
			.hash_bits = (uint32_t)hash_bits,
			.threshold = threshold,
			.sketch_entries = sketch_entries,
			.selected_observed_ctx = selected_obs,
			.selected_estimated_unique_ctx = selected_est,
			.preconflict_observed_ctx = pre_obs,
			.preconflict_estimated_unique_ctx = pre_est,
			.postconflict_observed_ctx = post_obs,
			.postconflict_estimated_unique_ctx = post_est,
		};
	}

	free(line);
	if (fclose(fp) != 0)
		err(errno, "%s", ctxmeta_path);
	free(ctxmeta_path);
	return stats;
}

typedef struct ani_density_af
{
	double qry;
	double ref;
	unsigned char available;
} ani_density_af_t;

static inline const ani_ctxmeta_rec_t *ani_ctxmeta_at(const ani_ctxmeta_rec_t *stats, uint32_t idx)
{
	return stats ? &stats[idx] : NULL;
}

static inline bool ani_ctxmeta_usable(const ani_ctxmeta_rec_t *m)
{
	return m && m->valid && m->hash_bits > 0 && m->hash_bits < 64 &&
		   m->selected_estimated_unique_ctx > 0;
}

static inline ani_density_af_t ani_density_af_fallback(double sketch_af_qry, double sketch_af_ref)
{
	return (ani_density_af_t){
		.qry = bounded_align_fraction(sketch_af_qry),
		.ref = bounded_align_fraction(sketch_af_ref),
		.available = 0,
	};
}

static ani_density_af_t ani_estimate_density_af(const ani_ctxmeta_rec_t *qry_meta,
												const ani_ctxmeta_rec_t *ref_meta,
												uint32_t shared_ctx,
												double sketch_af_qry,
												double sketch_af_ref)
{
	ani_density_af_t out = ani_density_af_fallback(sketch_af_qry, sketch_af_ref);
	if (shared_ctx == 0 || !ani_ctxmeta_usable(qry_meta) || !ani_ctxmeta_usable(ref_meta) ||
		qry_meta->hash_bits != ref_meta->hash_bits)
		return out;

	const uint64_t common_threshold =
		qry_meta->threshold < ref_meta->threshold ? qry_meta->threshold : ref_meta->threshold;
	const long double hash_space = ldexpl(1.0L, (int)qry_meta->hash_bits);
	const long double common_fraction = ((long double)common_threshold + 1.0L) / hash_space;
	if (common_fraction <= 0.0L)
		return out;

	const long double intersection = (long double)shared_ctx / common_fraction;
	if (intersection <= 0.0L)
		return out;

	out.qry = bounded_align_fraction((double)(intersection /
		(long double)qry_meta->selected_estimated_unique_ctx));
	out.ref = bounded_align_fraction((double)(intersection /
		(long double)ref_meta->selected_estimated_unique_ctx));
	out.available = 1;
	return out;
}

static uint32_t ani_density_af_needed_ctx(const ani_opt_t *ani_opt,
										  uint32_t qry_ctx,
										  uint32_t ref_ctx,
										  const ani_ctxmeta_rec_t *qry_meta,
										  const ani_ctxmeta_rec_t *ref_meta)
{
	uint32_t need = ani_report_af_needed_ctx(ani_opt, qry_ctx, ref_ctx);
	if (!ani_opt || ani_opt->afcut <= 0.0f ||
		!ani_ctxmeta_usable(qry_meta) || !ani_ctxmeta_usable(ref_meta) ||
		qry_meta->hash_bits != ref_meta->hash_bits)
		return need;

	const uint64_t common_threshold =
		qry_meta->threshold < ref_meta->threshold ? qry_meta->threshold : ref_meta->threshold;
	const long double hash_space = ldexpl(1.0L, (int)qry_meta->hash_bits);
	const long double common_fraction = ((long double)common_threshold + 1.0L) / hash_space;
	if (common_fraction <= 0.0L)
		return need;

	const long double min_total =
		qry_meta->selected_estimated_unique_ctx < ref_meta->selected_estimated_unique_ctx
			? (long double)qry_meta->selected_estimated_unique_ctx
			: (long double)ref_meta->selected_estimated_unique_ctx;
	const long double need_ld =
		(long double)ani_opt->afcut * min_total * common_fraction;
	if (need_ld <= 0.0L)
		return 0;
	if (need_ld >= (long double)UINT32_MAX)
		return UINT32_MAX;
	return (uint32_t)ceill(need_ld);
}

static inline bool ani_best_guard_enabled(const ani_opt_t *ani_opt)
{
	return ani_opt && ani_opt->fmt == 0 && !ani_opt->unassembled && !ani_opt->v;
}

static void load_infile_meta_for_best_guard(unify_sketch_t *sketch, const char *sketch_dir,
											   const ani_opt_t *ani_opt)
{
	if (!ani_best_guard_enabled(ani_opt) || !sketch || sketch->infile_meta)
		return;
	sketch->infile_meta = read_optional_sketch_infile_meta_stats(sketch_dir, sketch->infile_num);
}

static unsigned ani_query_parse_flags(const ani_opt_t *ani_opt)
{
	return ani_best_guard_enabled(ani_opt) ? SKETCH_PARSE_INFILE_META : SKETCH_PARSE_NONE;
}

static unsigned ani_ref_parse_flags(const ani_opt_t *ani_opt)
{
	unsigned flags = SKETCH_PARSE_NONE;
	if (ani_opt && ani_opt->fmt == 0)
		flags |= SKETCH_PARSE_ANNOTATION;
	if (ani_best_guard_enabled(ani_opt))
		flags |= SKETCH_PARSE_INFILE_META;
	return flags;
}

static bool force_ref_index_requested(void)
{
	const char *env = getenv("MINCO_FORCE_REF_INDEX");
	return env && env[0] != '\0' && strcmp(env, "0") != 0;
}

static void *read_mmap_ro_file(const char *file_path, size_t *file_size)
{
	int fd = open(file_path, O_RDONLY);
	if (fd == -1)
		err(EXIT_FAILURE, "%s(): open %s failed", __func__, file_path);

	struct stat st;
	if (fstat(fd, &st) != 0)
	{
		close(fd);
		err(EXIT_FAILURE, "%s(): fstat %s failed", __func__, file_path);
	}
	if (st.st_size <= 0)
	{
		close(fd);
		errx(EXIT_FAILURE, "%s(): file is empty: %s", __func__, file_path);
	}

	void *buffer = mmap(NULL, (size_t)st.st_size, PROT_READ, MAP_SHARED, fd, 0);
	if (buffer == MAP_FAILED)
	{
		close(fd);
		err(EXIT_FAILURE, "%s(): mmap %s failed", __func__, file_path);
	}
	close(fd);
#ifdef MADV_WILLNEED
	(void)madvise(buffer, (size_t)st.st_size, MADV_WILLNEED);
#endif
	if (file_size)
		*file_size = (size_t)st.st_size;
	return buffer;
}

static ctxgidobj_t *read_reference_sorted_index(const char *file_path, size_t *file_size,
												bool *is_mmap)
{
	*is_mmap = force_ref_index_requested();
	if (*is_mmap)
		return read_mmap_ro_file(file_path, file_size);
	return read_from_file(file_path, file_size);
}

static void free_reference_sorted_index(ctxgidobj_t *index, size_t file_size, bool is_mmap)
{
	if (is_mmap)
	{
		if (munmap(index, file_size) == -1)
			err(EXIT_FAILURE, "%s(): munmap failed", __func__);
		return;
	}
	free_read_from_file(index, file_size);
}

static uint32_t count_ctx_runs_sorted_ctxobj64_local(const uint64_t *a, size_t n, bool ignoreconflict)
{
	const uint8_t nobjbits = Bitslen.obj;
	uint32_t runs = 0;
	for (size_t i = 0; i < n; ) {
		const uint64_t ctx = a[i] >> nobjbits;
		const size_t begin = i;
		do { ++i; } while (i < n && (a[i] >> nobjbits) == ctx);
		if (!ignoreconflict || i - begin == 1)
			runs++;
	}
	return runs;
}

static uint32_t *ref_ctx_counts_from_ctxgidobj(const ctxgidobj_t *arr, size_t arrlen, int ref_infile_num,
											   const uint64_t *ref_sketch_index, bool ignoreconflict)
{
	uint32_t *counts = calloc((size_t)ref_infile_num, sizeof(*counts));
	if (!counts)
		err(EXIT_FAILURE, "%s(): calloc counts", __func__);

	for (int gid = 0; gid < ref_infile_num; gid++)
		counts[gid] = (uint32_t)(ref_sketch_index[gid + 1] - ref_sketch_index[gid]);

	for (size_t i = 0; i < arrlen; ) {
		const uint64_t ctxgid = arr[i].ctxgid;
		const uint32_t gid = (uint32_t)(ctxgid & ((1ULL << GID_NBITS) - 1ULL));
		size_t run_len = 0;
		do { ++i; ++run_len; } while (i < arrlen && arr[i].ctxgid == ctxgid);
		if (run_len > 1)
			counts[gid] -= (uint32_t)(ignoreconflict ? run_len : run_len - 1);
	}
	return counts;
}

static void fill_ctx_counts_for_query_block(uint32_t *qry_ctx_count, int offset_gid, int this_block_size,
											const uint64_t *qry_sketch_index, const uint64_t *tmp_ctxobj,
											bool qry_conflict)
{
	const uint64_t block_base = qry_sketch_index[offset_gid];
	for (int i = 0; i < this_block_size; i++) {
		const int gid = offset_gid + i;
		const uint64_t begin = qry_sketch_index[gid] - block_base;
		const uint64_t end = qry_sketch_index[gid + 1] - block_base;
		const size_t len = (size_t)(end - begin);
		if (len == 0) {
			qry_ctx_count[gid] = 0;
			continue;
		}
		qry_ctx_count[gid] = qry_conflict
			? count_ctx_runs_sorted_ctxobj64_local(tmp_ctxobj + begin, len, false)
			: (uint32_t)len;
	}
}

#define BLOCK_SIZE (4096) // #of qry genomes per batch, for mem_eff handling
int mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(ani_opt_t *ani_opt)
{
	
	minco_sketch_stat_t *ref_sketch_stat = read_from_file(test_get_fullpath(ani_opt->refdir, sketch_stat), &file_size);
	int ref_infile_num = ref_sketch_stat->infile_num;
	// read index
	size_t ctxgidobj_arr_fsize;
	uint64_t *ref_sketch_index = read_from_file(test_get_fullpath(ani_opt->refdir, idx_sketch_suffix), &file_size);
	assert(file_size == (ref_infile_num + 1) * sizeof(ref_sketch_index[0]));
	size_t ref_sketch_size = ref_sketch_index[ref_infile_num];
	char *sorted_index_path = test_get_fullpath(ani_opt->refdir, sorted_comb_ctxgid64obj32);
	bool sorted_index_is_mmap = false;
	ctxgidobj_t *sortedcomb_ctxgid64obj32 = read_reference_sorted_index(sorted_index_path, &ctxgidobj_arr_fsize, &sorted_index_is_mmap);
	free(sorted_index_path);
	assert(ctxgidobj_arr_fsize == ref_sketch_size * sizeof(sortedcomb_ctxgid64obj32[0]));
	uint32_t *ref_ctx_count = ref_ctx_counts_from_ctxgidobj(sortedcomb_ctxgid64obj32, ref_sketch_size, ref_infile_num, ref_sketch_index, ani_opt->ignoreconflict);

	minco_sketch_stat_t *qry_sketch_stat = read_from_file(test_get_fullpath(ani_opt->qrydir, sketch_stat), &file_size);
	int qry_infile_num = qry_sketch_stat->infile_num;
	assert(qry_sketch_stat->hash_id == ref_sketch_stat->hash_id);
	uint64_t *qry_sketch_index = read_from_file(test_get_fullpath(ani_opt->qrydir, idx_sketch_suffix), &file_size);
	size_t qry_sketch_size = qry_sketch_index[qry_infile_num];
	uint32_t *qry_ctx_count = calloc((size_t)qry_infile_num, sizeof(*qry_ctx_count));
	if (!qry_ctx_count)
		err(EXIT_FAILURE, "%s(): calloc qry_ctx_count", __func__);

	int block_size = qry_infile_num < BLOCK_SIZE ? qry_infile_num : BLOCK_SIZE;
	int offset_gid = 0;
	const char *qry_comb_sketch_path = test_get_fullpath(ani_opt->qrydir, combined_sketch_suffix);
	FILE *fp = fopen(qry_comb_sketch_path, "rb");
	if (fp == NULL)
		err(errno, "%s", qry_comb_sketch_path);
	uint64_t *tmp_ctxobj = NULL;
	size_t tmp_ctxobj_capacity = 0;
	// uint32_t *ctx = malloc(ref_infile_num * block_size * sizeof(uint32_t));
	// uint32_t *obj = malloc( ref_infile_num * block_size * sizeof(uint32_t))  ;
	ctx_mut2_t *ctx = malloc(ref_infile_num * block_size * sizeof(ctx_mut2_t));
	obj_section_t *obj = malloc(ref_infile_num * block_size * sizeof(obj_section_t));
	// for order id by descending ani
	uint32_t *num_passid_block = malloc(block_size * sizeof(uint32_t));
	idani_t **sort_idani_block = malloc(block_size * sizeof(idani_t *));
	for (int i = 0; i < block_size; i++)
		sort_idani_block[i] = malloc(ref_infile_num * sizeof(idani_t));

	char (*refname)[PATHLEN] = (char (*)[PATHLEN])(ref_sketch_stat + 1);
	char (*qryname)[PATHLEN] = (char (*)[PATHLEN])(qry_sketch_stat + 1);
	char (*refanno)[PATHLEN] = read_optional_sketch_annotations(ani_opt->refdir, ref_infile_num);
	const bool enable_best_guard = ani_best_guard_enabled(ani_opt);
	infile_meta_t *qry_infile_meta =
		enable_best_guard ? read_optional_sketch_infile_meta_stats(ani_opt->qrydir, qry_infile_num) : NULL;
	infile_meta_t *ref_infile_meta =
		enable_best_guard ? read_optional_sketch_infile_meta_stats(ani_opt->refdir, ref_infile_num) : NULL;
	ani_ctxmeta_rec_t *qry_ctxmeta =
		read_optional_ani_ctxmeta_stats(ani_opt->qrydir, qry_infile_num);
	ani_ctxmeta_rec_t *ref_ctxmeta =
		read_optional_ani_ctxmeta_stats(ani_opt->refdir, ref_infile_num);
	ani_opt_t scan_opt = *ani_opt;
	if ((qry_ctxmeta || ref_ctxmeta) && scan_opt.fmt == 0)
		scan_opt.afcut = 0.0f;

	FILE *outfp = ani_opt->outf[0] == '\0' ? stdout : fopen(ani_opt->outf, "w");
	if (outfp == NULL)
		err(errno, "%s", ani_opt->outf);

	/* load model
	if (ani_opt->model[0] == '\0')
		err(EXIT_FAILURE, "%s(): need specify model file using -M ", __func__);
	init_model(ani_opt->model); // f8C9O7_model.xgb
	*/

	// printf header
	if (ani_opt->fmt)
	{ // matrix format
		for (int i = 0; i < ref_infile_num; i++)
			fprintf(outfp, "\t%s", refname[i]);
		fprintf(outfp, "\n");
	}
	else
		print_ani_detail_header(outfp, ani_opt, true);

	while (offset_gid < qry_infile_num)
	{

		int this_block_size = qry_infile_num - offset_gid;
		if (this_block_size > block_size)
			this_block_size = block_size;
		size_t this_sketch_size = (size_t)(qry_sketch_index[offset_gid + this_block_size] - qry_sketch_index[offset_gid]);
		if (this_sketch_size > tmp_ctxobj_capacity)
		{
			uint64_t *new_tmp_ctxobj = realloc(tmp_ctxobj, this_sketch_size * sizeof(uint64_t));
			if (!new_tmp_ctxobj)
				err(EXIT_FAILURE, "%s(): realloc tmp_ctxobj", __func__);
			tmp_ctxobj = new_tmp_ctxobj;
			tmp_ctxobj_capacity = this_sketch_size;
		}
		size_t read_sketch_size = fread(tmp_ctxobj, sizeof(uint64_t), this_sketch_size, fp);
		uint64_t *this_sketch_index = qry_sketch_index + offset_gid;
		assert(this_sketch_size == read_sketch_size);
		fill_ctx_counts_for_query_block(qry_ctx_count, offset_gid, this_block_size, qry_sketch_index, tmp_ctxobj, qry_sketch_stat->conflict);

		memset(ctx, 0, ref_infile_num * block_size * sizeof(ctx_mut2_t));
		memset(obj, 0, ref_infile_num * block_size * sizeof(obj_section_t)); // memset(obj,0,ref_infile_num * block_size * sizeof(uint32_t));
		count_ctx_obj_frm_comb_sketch_section(ctx, obj, sortedcomb_ctxgid64obj32, ref_sketch_size, ref_infile_num, ref_ctx_count, qry_ctx_count + offset_gid, this_block_size, tmp_ctxobj, this_sketch_index, num_passid_block, sort_idani_block, &scan_opt);
		ani_block_print(ref_infile_num, offset_gid, this_block_size, ref_sketch_index, qry_sketch_index, ref_ctx_count, qry_ctx_count, ctx, obj, refname, qryname, refanno, qry_infile_meta, ref_infile_meta, qry_ctxmeta, ref_ctxmeta, num_passid_block, sort_idani_block, outfp, ani_opt, ani_opt->fmt);

		offset_gid += this_block_size;
	}

	for (int i = 0; i < block_size; i++)
		free(sort_idani_block[i]);
	if (refanno)
		free_read_from_file(refanno, (size_t)ref_infile_num * PATHLEN);
	if (qry_infile_meta)
		free_read_from_file(qry_infile_meta, (size_t)qry_infile_num * sizeof(qry_infile_meta[0]));
	if (ref_infile_meta)
		free_read_from_file(ref_infile_meta, (size_t)ref_infile_num * sizeof(ref_infile_meta[0]));
	free(qry_ctxmeta);
	free(ref_ctxmeta);
	free_all(ref_sketch_stat, ref_sketch_index, ref_ctx_count, qry_sketch_stat, qry_sketch_index, qry_ctx_count, tmp_ctxobj, ctx, obj, num_passid_block, sort_idani_block, NULL);
	free_reference_sorted_index(sortedcomb_ctxgid64obj32, ctxgidobj_arr_fsize, sorted_index_is_mmap);
	fclose(fp);
	if (outfp != stdout)
		fclose(outfp);
	// clean xgb model
	// cleanup_model();
	return ctxgidobj_arr_fsize;
}

//
#define DIFF_OBJ_BITS 1
size_t dedup_with_ctxobj_counts(uint32_t *arr, size_t n, co_distance_t **ctxobj_cnt)
{
	*ctxobj_cnt = NULL;
	if (n == 0)
		return 0;

	co_distance_t *tmp_ctxobj_cnt = malloc(n * sizeof(co_distance_t));
	if (!tmp_ctxobj_cnt)
		err(EXIT_FAILURE, "%s(): tmp_ctxobj_cnt malloc failure", __func__);

	size_t j = 0;
	tmp_ctxobj_cnt[0].ctx_ct = 1;
	tmp_ctxobj_cnt[0].diff_obj = arr[0] % 2;
	arr[0] >>= DIFF_OBJ_BITS;

	for (size_t i = 1; i < n; i++)
	{
		if ((arr[i] >> DIFF_OBJ_BITS) == arr[j])
		{
			tmp_ctxobj_cnt[j].ctx_ct++;
			tmp_ctxobj_cnt[j].diff_obj += (arr[i] % 2);
		}
		else
		{
			j++;
			tmp_ctxobj_cnt[j].ctx_ct = 1;
			tmp_ctxobj_cnt[j].diff_obj = arr[i] % 2;
			arr[j] = arr[i] >> DIFF_OBJ_BITS;
		}
	}
	// Trim arrays
	co_distance_t *ctxobj_tmp = realloc(tmp_ctxobj_cnt, (j + 1) * sizeof(co_distance_t));
	*ctxobj_cnt = ctxobj_tmp ? ctxobj_tmp : tmp_ctxobj_cnt;

	return j + 1;
}

ctxgidobj_t *comb_sortedsketch64_2sortedcomb_ctxgid64obj32(unify_sketch_t *ref_result)
{
	// const_comask_init(&ref_result->stats.minco_stat);
	uint64_t sketch_size = ref_result->sketch_index[ref_result->infile_num];
	if (sketch_size > (float)UINT32_MAX * LD_FCTR)
		err(EXIT_FAILURE, "%s():sketch_index maximun %lu exceed UINT32_MAX*LF;%f", __func__, sketch_size, (float)UINT32_MAX * LD_FCTR);
	if (ref_result->infile_num >= (1 << GID_NBITS))
		err(EXIT_FAILURE, "%s(): genome numer %d exceed maximum:%u", __func__, ref_result->infile_num, 1 << GID_NBITS);
	if (GID_NBITS + 4 * hclen > 64)
		err(EXIT_FAILURE, "%s(): context_bits_len(%d)+gid_bits_len(%d) exceed 64", __func__, 4 * hclen, GID_NBITS);
	ctxgidobj_t *ctxgidobj_arr = ctxobj64_2ctxgidobj(ref_result->sketch_index, ref_result->comb_sketch, ref_result->infile_num, sketch_size);
	ctxgidobj_sort_array(ctxgidobj_arr, sketch_size);
	return ctxgidobj_arr;
}

sort_sketch_summary_t *summarize_ctxgidobj_arr(ctxgidobj_t *ctxgidobj_arr, uint64_t *sketch_index, uint32_t arrlen, int infile_num)
{

	uint64_t gidmask = UINT64_MAX >> (64 - GID_NBITS);
	num_ctx_cfltobj_t *num_ctx_cfltobj_arr = malloc(infile_num * sizeof(num_ctx_cfltobj_t));
	for (int i = 0; i < infile_num; i++)
	{
		num_ctx_cfltobj_arr[i].num_ctx = sketch_index[i + 1] - sketch_index[i];
		num_ctx_cfltobj_arr[i].num_conflictobj = 0;
	}

	uint32_t num_ctx = 1;
	for (uint32_t i = 1; i < arrlen; i++)
	{
		if (ctxgidobj_arr[i].ctxgid >> GID_NBITS != ctxgidobj_arr[i - 1].ctxgid >> GID_NBITS)
			num_ctx++;
		if (ctxgidobj_arr[i].ctxgid == ctxgidobj_arr[i - 1].ctxgid)
		{
			uint32_t gid = (uint32_t)ctxgidobj_arr[i].ctxgid & gidmask;
			if (num_ctx_cfltobj_arr[gid].num_ctx == (sketch_index[gid + 1] - sketch_index[gid]))
				num_ctx_cfltobj_arr[gid].num_conflictobj++;

			num_ctx_cfltobj_arr[gid].num_ctx--;
		}
	}
	//	for(int i = 0 ; i < infile_num; i++) num_conflictobj +=  num_ctx_cfltobj_arr[i].num_conflictobj;
	sort_sketch_summary_t *sort_sketch_summary = malloc(sizeof(sort_sketch_summary_t));
	sort_sketch_summary->arrlen = arrlen;
	sort_sketch_summary->numgids_perctx = (double)arrlen / num_ctx;
	sort_sketch_summary->num_ctx = num_ctx;
	sort_sketch_summary->infile_num = infile_num;
	sort_sketch_summary->num_ctx_cfltobj_arr = num_ctx_cfltobj_arr;
	return sort_sketch_summary;
}

void free_sort_sketch_summary(sort_sketch_summary_t *sort_sketch_summary)
{
	free(sort_sketch_summary->num_ctx_cfltobj_arr);
	free(sort_sketch_summary);
}

#define CTX(X, Y) (ctx[(size_t)(((X) * ((X) + 1)) / 2 + (Y))])
#define OBJ(X, Y) (obj[(size_t)(((X) * ((X) + 1)) / 2 + (Y))])

/* use sketch variants to caculate distance */
// 1.global sorted comb_sketch64 (i.e. sorted_ctxgidobj_arr or inverted index):
// slow when dist matrix is large, low memory cache efficient, but may be very fast when matrix is small?
void sorted_ctxgidobj_arr2triangle(ctxgidobj_t *ctxgidobj_arr, sort_sketch_summary_t *sort_sketch_summary)
{
	uint64_t gidmask = UINT64_MAX >> (64 - GID_NBITS);
	uint32_t arrlen = sort_sketch_summary->arrlen;
	int infile_num = sort_sketch_summary->infile_num;
	uint16_t *ctx = calloc((size_t)infile_num * (infile_num + 1) / 2, sizeof(uint16_t));
	uint16_t *obj = calloc((size_t)infile_num * (infile_num + 1) / 2, sizeof(uint16_t));

	for (uint32_t i = 0, j; i < arrlen - 1; i = j)
	{
		j = i + 1;
		// find range i..j;
		for (; j < arrlen && (ctxgidobj_arr[i].ctxgid >> GID_NBITS == ctxgidobj_arr[j].ctxgid >> GID_NBITS); j++)
			;

#pragma omp parallel for num_threads(32) schedule(guided)
		for (uint32_t a = i + 1; a < j; a++)
		{
			if (ctxgidobj_arr[a].ctxgid == ctxgidobj_arr[a - 1].ctxgid || ctxgidobj_arr[a].ctxgid == ctxgidobj_arr[a + 1].ctxgid)
				continue; // with confclit object
			uint32_t x = ctxgidobj_arr[a].ctxgid & gidmask;
#pragma omp parallel for num_threads(32) schedule(guided)
			for (uint32_t b = i; b < a; b++)
			{

				if ((b > 0 && ctxgidobj_arr[b].ctxgid == ctxgidobj_arr[b - 1].ctxgid) || ctxgidobj_arr[b].ctxgid == ctxgidobj_arr[b + 1].ctxgid)
					continue;
				uint32_t y = ctxgidobj_arr[b].ctxgid & gidmask;
				CTX(x, y)
				++;
				if (ctxgidobj_arr[a].obj != ctxgidobj_arr[b].obj)
					OBJ(x, y)
				++;
			}
		}
		printf("\ri=%d", i);
	}
	for (int x = 1; x < infile_num; x++)
	{

		for (int y = 0; y < x; y++)
		{
			if (CTX(x, y) > 0)
			{
				printf("%d\t%d\t%d\t%d\t%f\n", x, y, CTX(x, y), OBJ(x, y), (float)OBJ(x, y) / CTX(x, y));
			}
		}
	}
}

// sparse_mem_eff.. seems slower than mem_eff,
int sparse_mem_eff_sorted_ctxgidobj_arrXcomb_sortedsketch64(ani_opt_t *ani_opt)
{
	// initialize
	minco_sketch_stat_t *ref_sketch_stat = read_from_file(test_get_fullpath(ani_opt->refdir, sketch_stat), &file_size);
	const_comask_init(ref_sketch_stat);

	uint64_t gidmask = UINT64_MAX >> (64 - GID_NBITS);
	uint64_t objmask = (1UL << Bitslen.obj) - 1;

	int ref_infile_num = ref_sketch_stat->infile_num;
	// read index
	size_t ctxgidobj_arr_fsize;
	uint64_t *ref_sketch_index = read_from_file(test_get_fullpath(ani_opt->refdir, idx_sketch_suffix), &file_size);
	assert(file_size == (ref_infile_num + 1) * sizeof(ref_sketch_index[0]));
	size_t ref_sketch_size = ref_sketch_index[ref_infile_num];
	ctxgidobj_t *sortedcomb_ctxgid64obj32 = read_from_file(test_get_fullpath(ani_opt->refdir, sorted_comb_ctxgid64obj32), &ctxgidobj_arr_fsize);
	assert(ctxgidobj_arr_fsize == ref_sketch_size * sizeof(sortedcomb_ctxgid64obj32[0]));

	minco_sketch_stat_t *qry_sketch_stat = read_from_file(test_get_fullpath(ani_opt->qrydir, sketch_stat), &file_size);
	int qry_infile_num = qry_sketch_stat->infile_num;
	assert(qry_sketch_stat->hash_id == ref_sketch_stat->hash_id);
	uint64_t *qry_sketch_index = read_from_file(test_get_fullpath(ani_opt->qrydir, idx_sketch_suffix), &file_size);
	size_t qry_sketch_size = qry_sketch_index[qry_infile_num];

	int block_size = BLOCK_SIZE;
	int offset_gid = 0;
	const char *qry_comb_sketch_path = test_get_fullpath(ani_opt->qrydir, combined_sketch_suffix);
	FILE *fp = fopen(qry_comb_sketch_path, "rb");
	if (fp == NULL)
		err(errno, "%s", qry_comb_sketch_path);
	uint64_t *tmp_ctxobj = malloc(qry_sketch_size * sizeof(uint64_t));
	// arry of dynamic gids array per query, note gid <<= DIFF_OBJ_BITS + (ref_obj==qry_obj?0:1);  >>sparse only
	//    Vector *ref_gids_perqry_arr  = malloc( block_size * sizeof(Vector));
	// ref_gids_perqry_arr[i][0] is capacity, ref_gids_perqry_arr[i][1] is length of the inner array
	uint32_t **ref_gids_perqry_arr = malloc(block_size * sizeof(uint32_t *));
	co_distance_t **ctxobj_cnt_perqry_arr = malloc(block_size * sizeof(co_distance_t *));
	size_t *lens = malloc(block_size * sizeof(size_t));
	for (int i = 0; i < block_size; i++)
	{
		ref_gids_perqry_arr[i] = malloc(6 * sizeof(uint32_t)); // 4 for innitialized capacity
		ref_gids_perqry_arr[i][0] = 4;
		ref_gids_perqry_arr[i][1] = 0;
	}

	for (int b = 0; b <= qry_infile_num / block_size; b++)
	{

		int this_block_size = (b == qry_infile_num / block_size) ? (qry_infile_num % block_size) : block_size;
		int this_sketch_size = qry_sketch_index[offset_gid + this_block_size] - qry_sketch_index[offset_gid];
		int read_sketch_size = fread(tmp_ctxobj, sizeof(uint64_t), this_sketch_size, fp);
		uint64_t *this_sketch_index = qry_sketch_index + offset_gid;
		assert(this_sketch_size == read_sketch_size);

#pragma omp parallel for num_threads(ani_opt->p) schedule(guided)
		for (int i = 0; i < this_block_size; i++)
		{
			uint64_t *a = tmp_ctxobj + (this_sketch_index[i] - this_sketch_index[0]);
			size_t a_size = this_sketch_index[i + 1] - this_sketch_index[i];
			size_t *idx = find_first_occurrences_AT_ctxgidobj_arr(a, a_size, sortedcomb_ctxgid64obj32, ref_sketch_size);

			for (int j = 0; j < a_size; j++)
			{
				if (idx[j] == SIZE_MAX)
					continue;
				// skip conlict object;
				if ((j > 0) && ((a[j] >> Bitslen.obj) == (a[j - 1] >> Bitslen.obj)))
					continue;
				if ((j < a_size - 1) && ((a[j] >> Bitslen.obj) == (a[j + 1] >> Bitslen.obj)))
					continue;

				for (int d = idx[j];; d++)
				{
					if ((sortedcomb_ctxgid64obj32[d].ctxgid >> Bitslen.gid) != (a[j] >> Bitslen.obj))
						break;
					uint32_t gid01 = sortedcomb_ctxgid64obj32[d].ctxgid & gidmask << DIFF_OBJ_BITS;
					if ((a[j] & objmask) != sortedcomb_ctxgid64obj32[d].obj)
						gid01 |= 1;
					if (ref_gids_perqry_arr[i][1] + 2 == ref_gids_perqry_arr[i][0])
					{
						ref_gids_perqry_arr[i][0] += 100;
						ref_gids_perqry_arr[i] = realloc(ref_gids_perqry_arr[i], ref_gids_perqry_arr[i][0] * sizeof(ref_gids_perqry_arr[i][0]));
					}
					ref_gids_perqry_arr[i][2 + ref_gids_perqry_arr[i][1]] = gid01;
					ref_gids_perqry_arr[i][1]++;
					// vector_push(&ref_gids_perqry_arr[i],&gid01);
				}
			}
			free(idx);
			qsort(ref_gids_perqry_arr[i] + 2, ref_gids_perqry_arr[i][1], sizeof(uint32_t), qsort_comparator_uint32);
			lens[i] = dedup_with_ctxobj_counts(ref_gids_perqry_arr[i] + 2, ref_gids_perqry_arr[i][1], &ctxobj_cnt_perqry_arr[i]);
		}

		for (int i = 0; i < this_block_size; i++)
		{
			int qry_gid = b * block_size + i;
			int qry_sketch_size = ref_sketch_index[qry_gid + 1] - ref_sketch_index[qry_gid];
			// sparse only code >>
			for (int l = 0; l < lens[i]; l++)
			{
				double dist = ctxobj_cnt_perqry_arr[i][l].diff_obj / ctxobj_cnt_perqry_arr[i][l].ctx_ct;
				double ani = (1 - dist);
#define VEC_GET_AS(type, vec, idx) (((type *)(vec).data)[idx])
				//			if(ctxobj_cnt_perqry_arr[i][l].ctx_ct > 20)
				//          	printf ("%d|%d|%d|%d|%d|%lf|%lf\t",qry_gid,VEC_GET_AS(uint32_t, ref_gids_perqry_arr[i], l),qry_sketch_size,ctxobj_cnt_perqry_arr[i][l].ctx_ct,ctxobj_cnt_perqry_arr[i][l].diff_obj,dist,ani);
			}
			//        printf("\n");
			//<<
		}
		//	for(int i =0 ; i< this_block_size; i++) vector_free(&ref_gids_perqry_arr[i]);//ref_gids_perqry_arr[i].size = 0;
#pragma omp parallel for num_threads(ani_opt->p) schedule(guided)
		for (int i = 0; i < this_block_size; i++)
		{
			ref_gids_perqry_arr[i][1] = 0;
			//			ref_gids_perqry_arr[i][0] = 4;
			free(ctxobj_cnt_perqry_arr[i]);
		}
		offset_gid += this_block_size;
	} // all blocks loop end

	for (int i = 0; i < block_size; i++)
		free(ref_gids_perqry_arr[i]);
	free_all(ref_gids_perqry_arr, ctxobj_cnt_perqry_arr, NULL); ////sparse only code >>
	free_all(ref_sketch_stat, ref_sketch_index, qry_sketch_stat, qry_sketch_index, tmp_ctxobj, NULL);
	free_read_from_file(sortedcomb_ctxgid64obj32, ctxgidobj_arr_fsize);
	fclose(fp);

	return ctxgidobj_arr_fsize;
}

// 2. inverted index(i.e. global sorted_ctxgidobj_arr) X common index(i.e. genome-wise sorted comb_sortedsketch64 ) :
//  ** the fatest method when dist matrix is sparse
void sorted_ctxgidobj_arrXcomb_sortedsketch64(unify_sketch_t *qry_result, ctxgidobj_t *ctxgidobj_arr, sort_sketch_summary_t *sort_sketch_summary)
{

	uint64_t *qry_sketch_index = qry_result->sketch_index;
	uint64_t *qry_comb_sketch = qry_result->comb_sketch;
	int qry_infile_num = qry_result->infile_num;

	uint64_t gidmask = (1UL << Bitslen.gid) - 1;
	uint64_t objmask = (1UL << Bitslen.obj) - 1;
	uint32_t ref_arrlen = sort_sketch_summary->arrlen;
	int ref_infile_num = sort_sketch_summary->infile_num;
	// for self comparision only
	assert(qry_infile_num == ref_infile_num);

	uint16_t *ctx = calloc((size_t)ref_infile_num * (ref_infile_num + 1) / 2, sizeof(uint16_t));
	uint16_t *obj = calloc((size_t)ref_infile_num * (ref_infile_num + 1) / 2, sizeof(uint16_t));
#pragma omp parallel for num_threads(32) schedule(guided)
	for (uint32_t rn = 0; rn < qry_infile_num; rn++)
	{
		//		if(rn > 0) break;
		uint64_t *a = qry_comb_sketch + qry_sketch_index[rn];
		size_t a_size = qry_sketch_index[rn + 1] - qry_sketch_index[rn];
		size_t *idx = find_first_occurrences_AT_ctxgidobj_arr(a, a_size, ctxgidobj_arr, ref_arrlen);
		for (int i = 0; i < a_size; i++)
		{
			if (idx[i] == SIZE_MAX)
				continue;
			// skip conlict object;
			if ((i > 0) && ((a[i] >> Bitslen.obj) == (a[i - 1] >> Bitslen.obj)))
				continue;
			if ((i < a_size - 1) && ((a[i] >> Bitslen.obj) == (a[i + 1] >> Bitslen.obj)))
				continue;

			for (int d = idx[i];; d++)
			{
				uint32_t gid = ctxgidobj_arr[d].ctxgid & gidmask;
				if (gid > rn || (ctxgidobj_arr[d].ctxgid >> Bitslen.gid) != (a[i] >> Bitslen.obj))
					break;
				CTX(rn, gid)
				++;
				if ((a[i] & objmask) != ctxgidobj_arr[d].obj) // wrong: if(a[i] & objmask != ctxgidobj_arr[d].obj) ...
					OBJ(rn, gid)
				++;
			}
		}
		free(idx);
		/*
		#pragma omp critical
				{
					for(int i = 0 ; i <= rn;i++){
						if( CTX(rn,i) > 0 ) printf("\t%d|%d|%d|%d|%f",i,rn, CTX(rn,i),OBJ(rn,i),(float) OBJ(rn,i)/CTX(rn, i));
					}
					printf("\n");
				}
		*/
	}

	for (uint32_t rn = 0; rn < qry_infile_num; rn++)
	{
		printf("%u", rn);
		for (int i = 0; i < rn; i++)
		{
			if (CTX(rn, i) > 0)
				printf("\t%d|%d|%d|%f", i, CTX(rn, i), OBJ(rn, i), (float)OBJ(rn, i) / CTX(rn, i));
		}
		printf("\n");
	}

} // end

// 3. common index X common index (i.e. genome-wise sorted comb_sortedsketch64 ):
//     code hits: use paris wise small sortted arrays overlapping
/*  Advatages:
		** immediate output (no need precompute dist matrix )
		** most convience and memory efficience, no invert indxeing needed
	** when genome are highly similar (dense dist matrix) speed is even faster than inverted index(i.e. global sorted_ctxgidobj_arr) X common index
		** small sortted arrays overlapping is ~ 2X-3X faster than hashtable lookup based overlapping.
*/

typedef struct {
    uint32_t rn;
    double   ani;          /* raw model ANI used internally */
    double   metric;       /* selected distance */
    double   selected_ani;
    double   moe_dist;
    double   naive_dist;
    double   mash_dist;
    double   aaf_dist;
    double   calibrated_ani;
    double   best_ani;
    unsigned char confidence;
    unsigned char best_guarded;
    unsigned char af_pass;
    uint8_t domain_profile;
    double   af_qry, blastn_af_qry, af_ref, blastn_af_ref;
    double   real_af_qry, real_af_ref;
    unsigned char real_af_available;
    int      XnY_ctx, N_diff_obj, N_diff_obj_section, N_mut2_ctx;
    uint64_t readwise_total_reads;
    uint64_t readwise_reads_with_ctx_match;
    uint64_t readwise_unique_query_ctx;
    uint64_t readwise_unique_query_ctx_hit;
    uint64_t readwise_unique_ref_ctx_hit;
    uint64_t readwise_ref_ctx_total;
    uint64_t readwise_density_block_ctx;
    uint64_t readwise_total_density_blocks;
    uint64_t readwise_blocks_with_ctx_match;
    uint64_t readwise_raw_xny_ctx;
    uint64_t readwise_rejected_ctx;
    uint64_t readwise_rejected_diff_ctx;
    double   readwise_fake_ctx_fraction;
    double   readwise_fake_ctx_prob_mean;
    double   readwise_fake_ctx_prob_weighted;
    double   abundance_ref_breadth;
    double   abundance_ref_mean_depth;
    double   abundance_ref_hit_mean_depth;
    double   abundance_ref_depth_variance;
    double   abundance_ref_depth_cv;
    double   abundance_ref_zero_fraction;
    double   abundance_relative_depth;
    double   abundance_normalized_depth;
    double   abundance_effective_depth;
    double   abundance_normalized_effective_depth;
    double   abundance_ref_zip_af;
    double   abundance_ref_zip_aaf_ani;
    double   reliable_ref_breadth;
    double   reliable_ref_mean_depth;
    uint64_t reliable_ref_hit_ctx;
    double   reliable_ref_hit_mean_depth;
    double   reliable_ref_hit_median_depth;
    double   reliable_ref_hit_depth_variance;
    double   reliable_ref_zip_af;
    uint64_t marker_xny_ctx;
    uint64_t marker_raw_xny_ctx;
    uint64_t marker_n_diff_obj;
    uint64_t marker_n_diff_obj_section;
    uint64_t marker_n_mut2_ctx;
    uint64_t marker_ref_ctx_total;
    double   marker_ref_breadth;
    double   marker_ref_mean_depth;
    double   marker_ref_hit_mean_depth;
    double   marker_ref_depth_variance;
    double   marker_ref_depth_cv;
    double   marker_ref_zero_fraction;
    double   marker_relative_depth;
    double   marker_effective_depth;
    double   marker_ref_zip_af;
    double   marker_ref_zip_aaf_ani;
    double   marker_reliable_ref_breadth;
    double   marker_reliable_ref_mean_depth;
    uint64_t marker_reliable_ref_hit_ctx;
    double   marker_reliable_ref_hit_mean_depth;
    double   marker_reliable_ref_hit_median_depth;
    double   marker_reliable_ref_hit_depth_variance;
    double   marker_reliable_ref_zip_af;
} ani_row_t;

static inline double ani_row_report_af_qry(const ani_row_t *r)
{
	return r && r->real_af_available ? r->real_af_qry : (r ? r->af_qry : 0.0);
}

static inline double ani_row_report_af_ref(const ani_row_t *r)
{
	return r && r->real_af_available ? r->real_af_ref : (r ? r->af_ref : 0.0);
}

static int cmp_ani_desc(const void *pa, const void *pb) {
    const ani_row_t *a = (const ani_row_t*)pa;
    const ani_row_t *b = (const ani_row_t*)pb;
    if (a->selected_ani < b->selected_ani) return  1;
    if (a->selected_ani > b->selected_ani) return -1;
    const double a_af = ani_row_report_af_qry(a);
    const double b_af = ani_row_report_af_qry(b);
    if (a_af < b_af) return  1;
    if (a_af > b_af) return -1;
    return (a->rn > b->rn) - (a->rn < b->rn);
}

/* typedef to avoid anonymous-struct warnings from kvec_t params */
typedef kvec_t(ani_row_t) kv_ani_row_t;

typedef enum ani_readwise_default_call
{
	ANI_READWISE_CALL_WEAK = 0,
	ANI_READWISE_CALL_LOW_ABUNDANCE = 1,
	ANI_READWISE_CALL_SCREENING = 2,
	ANI_READWISE_CALL_MAJOR = 3
} ani_readwise_default_call_t;

static inline uint64_t ani_readwise_default_support_cut_for_ref(uint64_t ref_ctx_total)
{
	uint64_t target = ref_ctx_total ? ref_ctx_total :
		(ani_model_target_sketch_size
			 ? (uint64_t)ani_model_target_sketch_size
			 : (uint64_t)ANI_MODEL_REFERENCE_SKETCH_SIZE);
	if (!target)
		target = (uint64_t)ANI_MODEL_REFERENCE_SKETCH_SIZE;
	uint64_t cut = (target + 99u) / 100u;
	if (cut < 100u)
		cut = 100u;
	if (cut > target)
		cut = target;
	return cut;
}

static inline uint64_t ani_readwise_default_support_cut(const ani_row_t *r)
{
	return ani_readwise_default_support_cut_for_ref(
		r ? r->readwise_ref_ctx_total : 0);
}

static uint64_t ani_readwise_default_min_support_cut(const uint32_t *ref_ctx_total,
													 uint32_t ref_n)
{
	uint64_t min_cut = UINT64_MAX;
	if (ref_ctx_total) {
		for (uint32_t rn = 0; rn < ref_n; ++rn) {
			if (!ref_ctx_total[rn])
				continue;
			const uint64_t cut =
				ani_readwise_default_support_cut_for_ref(ref_ctx_total[rn]);
			if (cut < min_cut)
				min_cut = cut;
		}
	}
	return min_cut == UINT64_MAX
			   ? ani_readwise_default_support_cut_for_ref(0)
			   : min_cut;
}

static inline bool ani_readwise_default_support_pass(const ani_row_t *r,
													 uint64_t support_cut)
{
	return r &&
		   (uint64_t)(r->XnY_ctx > 0 ? r->XnY_ctx : 0) >= support_cut &&
		   r->readwise_unique_ref_ctx_hit >= support_cut;
}

static inline uint64_t ani_readwise_fractional_support_cut(const ani_row_t *r,
														   uint64_t min_support,
														   double support_fraction)
{
	const uint64_t total = r ? r->readwise_ref_ctx_total : 0;
	uint64_t fractional = 0;
	if (total > 0 && support_fraction > 0.0)
		fractional = (uint64_t)ceill((long double)total *
									 (long double)support_fraction);
	return fractional > min_support ? fractional : min_support;
}

static inline bool ani_readwise_small_feature_profile(uint8_t profile)
{
	return profile == MINCO_DOMAIN_PROFILE_AMR ||
		   profile == MINCO_DOMAIN_PROFILE_GENE;
}

static inline ani_readwise_default_call_t ani_readwise_default_call(const ani_row_t *r)
{
	if (!r)
		return ANI_READWISE_CALL_WEAK;
	if (ani_readwise_small_feature_profile(r->domain_profile)) {
		const uint64_t high_support =
			ani_readwise_fractional_support_cut(r, 30, 0.03);
		if (r->abundance_ref_breadth >= 0.60 &&
			ani_readwise_default_support_pass(r, high_support) &&
			r->selected_ani >= 0.96)
			return ANI_READWISE_CALL_MAJOR;

		const uint64_t screen_support =
			ani_readwise_fractional_support_cut(r, 20, 0.02);
		if (r->abundance_ref_breadth >= 0.35 &&
			ani_readwise_default_support_pass(r, screen_support) &&
			r->selected_ani >= 0.95)
			return ANI_READWISE_CALL_SCREENING;
		return ANI_READWISE_CALL_WEAK;
	}
	const uint64_t support_cut = ani_readwise_default_support_cut(r);
	if (r->abundance_ref_breadth >= 0.5 &&
		r->abundance_relative_depth >= 1e-4 &&
		ani_readwise_default_support_pass(r, support_cut) &&
		r->selected_ani >= 0.95)
		return ANI_READWISE_CALL_MAJOR;
	if (r->abundance_ref_breadth >= 0.5 &&
		r->abundance_relative_depth >= 1e-5 &&
		ani_readwise_default_support_pass(r, support_cut) &&
		r->selected_ani >= 0.95)
		return ANI_READWISE_CALL_LOW_ABUNDANCE;
	return ANI_READWISE_CALL_WEAK;
}

static inline const char *ani_readwise_default_call_label(ani_readwise_default_call_t call)
{
	switch (call)
	{
	case ANI_READWISE_CALL_MAJOR:
		return "major";
	case ANI_READWISE_CALL_SCREENING:
		return "screening";
	case ANI_READWISE_CALL_LOW_ABUNDANCE:
		return "low_abundance";
	case ANI_READWISE_CALL_WEAK:
	default:
		return "weak";
	}
}

static inline const char *ani_readwise_default_call_rule(const ani_row_t *r,
														ani_readwise_default_call_t call)
{
	switch (call)
	{
	case ANI_READWISE_CALL_MAJOR:
		if (r && ani_readwise_small_feature_profile(r->domain_profile))
			return "domain=amr/gene;major:breadth>=0.60;support>=max(30,ceil(0.03*S));ANI>=0.96";
		return "breadth>=0.5;rel_depth>=1e-4;support>=min(S,max(100,ceil(S/100)));ANI>=0.95";
	case ANI_READWISE_CALL_SCREENING:
		return "domain=amr/gene;screening:breadth>=0.35;support>=max(20,ceil(0.02*S));ANI>=0.95";
	case ANI_READWISE_CALL_LOW_ABUNDANCE:
		return "breadth>=0.5;rel_depth>=1e-5;support>=min(S,max(100,ceil(S/100)));ANI>=0.95";
	case ANI_READWISE_CALL_WEAK:
	default:
		return "below_default_call_thresholds";
	}
}

static int cmp_readwise_default_call_desc(const void *pa, const void *pb)
{
	const ani_row_t *a = (const ani_row_t *)pa;
	const ani_row_t *b = (const ani_row_t *)pb;
	const int ca = (int)ani_readwise_default_call(a);
	const int cb = (int)ani_readwise_default_call(b);
	if (ca != cb)
		return cb - ca;
	if (a->abundance_relative_depth < b->abundance_relative_depth)
		return 1;
	if (a->abundance_relative_depth > b->abundance_relative_depth)
		return -1;
	if (a->abundance_ref_breadth < b->abundance_ref_breadth)
		return 1;
	if (a->abundance_ref_breadth > b->abundance_ref_breadth)
		return -1;
	if (a->selected_ani < b->selected_ani)
		return 1;
	if (a->selected_ani > b->selected_ani)
		return -1;
	return 0;
}

static void keep_readwise_default_calls(kv_ani_row_t *rows)
{
	if (!rows || kv_size(*rows) == 0)
		return;
	size_t w = 0;
	for (size_t i = 0; i < kv_size(*rows); ++i)
	{
		if (ani_readwise_default_call(&kv_A(*rows, i)) == ANI_READWISE_CALL_WEAK)
			continue;
		if (w != i)
			kv_A(*rows, w) = kv_A(*rows, i);
		++w;
	}
	rows->n = w;
}

static void normalize_readwise_abundance_depth(kv_ani_row_t *rows, size_t out_n)
{
	if (!rows || out_n == 0)
		return;
	if (out_n > kv_size(*rows))
		out_n = kv_size(*rows);
	long double sum = 0.0L;
	long double effective_sum = 0.0L;
	for (size_t i = 0; i < out_n; ++i) {
		const double value = kv_A(*rows, i).abundance_relative_depth;
		if (value > 0.0)
			sum += (long double)value;
		const double effective_value = kv_A(*rows, i).abundance_effective_depth;
		if (effective_value > 0.0)
			effective_sum += (long double)effective_value;
	}
	for (size_t i = 0; i < out_n; ++i) {
		ani_row_t *row = &kv_A(*rows, i);
		row->abundance_normalized_depth =
			sum > 0.0L && row->abundance_relative_depth > 0.0
				? (double)((long double)row->abundance_relative_depth / sum)
				: 0.0;
		row->abundance_normalized_effective_depth =
			effective_sum > 0.0L && row->abundance_effective_depth > 0.0
				? (double)((long double)row->abundance_effective_depth / effective_sum)
				: 0.0;
	}
}

static double ani_readwise_effective_abundance_depth_values(double median_depth,
															double mean_depth,
															double zip_af)
{
	if (median_depth >=
		MINCO_READWISE_EFFECTIVE_MEDIAN_DEPTH_CUTOFF)
		return median_depth;
	if (!isfinite(mean_depth) || mean_depth <= 0.0)
		return 0.0;
	if (!isfinite(zip_af) || zip_af <= 0.0)
		return mean_depth;
	return mean_depth / pow(zip_af, MINCO_READWISE_EFFECTIVE_AF_EXPONENT);
}

static double ani_readwise_effective_abundance_depth(const ani_row_t *row)
{
	if (!row)
		return 0.0;
	return ani_readwise_effective_abundance_depth_values(
		row->reliable_ref_hit_median_depth,
		row->reliable_ref_mean_depth,
		row->reliable_ref_zip_af);
}

typedef struct ani_cami_tax_record
{
	char *key;
	char *taxid;
	char *rank;
	char *taxpath;
	char *taxpathsn;
	char *genome_id;
	char *otu;
} ani_cami_tax_record_t;

typedef kvec_t(ani_cami_tax_record_t) kv_cami_tax_record_t;

typedef struct ani_cami_profile_entry
{
	const ani_cami_tax_record_t *tax;
	long double abundance;
} ani_cami_profile_entry_t;

typedef kvec_t(ani_cami_profile_entry_t) kv_cami_profile_entry_t;

static char *ani_cami_strdup_field(const char *s)
{
	if (!s || s[0] == '\0')
		s = "NA";
	char *out = strdup(s);
	if (!out)
		err(EXIT_FAILURE, "%s(): OOM string copy", __func__);
	return out;
}

static char *ani_cami_trim(char *s)
{
	if (!s)
		return s;
	while (*s && isspace((unsigned char)*s))
		++s;
	char *end = s + strlen(s);
	while (end > s && isspace((unsigned char)end[-1]))
		*--end = '\0';
	return s;
}

static void ani_cami_tax_record_destroy(ani_cami_tax_record_t *rec)
{
	if (!rec)
		return;
	free(rec->key);
	free(rec->taxid);
	free(rec->rank);
	free(rec->taxpath);
	free(rec->taxpathsn);
	free(rec->genome_id);
	free(rec->otu);
	memset(rec, 0, sizeof(*rec));
}

static void ani_cami_tax_records_destroy(kv_cami_tax_record_t *records)
{
	if (!records)
		return;
	for (size_t i = 0; i < kv_size(*records); ++i)
		ani_cami_tax_record_destroy(&kv_A(*records, i));
	kv_destroy(*records);
}

static int ani_cami_tax_record_key_cmp(const void *a, const void *b)
{
	const ani_cami_tax_record_t *ra = (const ani_cami_tax_record_t *)a;
	const ani_cami_tax_record_t *rb = (const ani_cami_tax_record_t *)b;
	return strcmp(ra->key, rb->key);
}

static bool ani_cami_looks_like_header(const char *field0, const char *field1)
{
	if (!field0)
		return false;
	return strcasecmp(field0, "ref") == 0 ||
		   strcasecmp(field0, "key") == 0 ||
		   strcasecmp(field0, "sample") == 0 ||
		   strcasecmp(field0, "ref_key") == 0 ||
		   strcasecmp(field0, "accession") == 0 ||
		   (field1 && strcasecmp(field1, "taxid") == 0);
}

static void ani_cami_load_taxmap(const char *path, kv_cami_tax_record_t *records)
{
	kv_init(*records);
	FILE *fp = fopen(path, "r");
	if (!fp)
		err(errno, "%s(): cannot open --cami-taxmap %s", __func__, path);

	char *line = NULL;
	size_t cap = 0;
	ssize_t len = 0;
	uint64_t line_no = 0;
	while ((len = getline(&line, &cap, fp)) >= 0) {
		++line_no;
		while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r'))
			line[--len] = '\0';
		char *cur = ani_cami_trim(line);
		if (cur[0] == '\0' || cur[0] == '#')
			continue;

		char *fields[7] = {0};
		int nf = 0;
		fields[nf++] = cur;
		for (char *p = cur; *p && nf < 7; ++p) {
			if (*p == '\t') {
				*p = '\0';
				fields[nf++] = p + 1;
			}
		}
		for (int i = 0; i < nf; ++i)
			fields[i] = ani_cami_trim(fields[i]);
		if (ani_cami_looks_like_header(fields[0], nf > 1 ? fields[1] : NULL))
			continue;
		if (nf < 5)
			errx(EINVAL,
				 "%s(): %s line %" PRIu64 " needs at least 5 tab-separated fields: ref_key TAXID RANK TAXPATH TAXPATHSN",
				 __func__, path, line_no);
		if (fields[0][0] == '\0' || fields[1][0] == '\0' ||
			fields[2][0] == '\0' || fields[3][0] == '\0' ||
			fields[4][0] == '\0')
			errx(EINVAL, "%s(): %s line %" PRIu64 " has an empty required CAMI field",
				 __func__, path, line_no);

		ani_cami_tax_record_t rec = {
			.key = ani_cami_strdup_field(fields[0]),
			.taxid = ani_cami_strdup_field(fields[1]),
			.rank = ani_cami_strdup_field(fields[2]),
			.taxpath = ani_cami_strdup_field(fields[3]),
			.taxpathsn = ani_cami_strdup_field(fields[4]),
			.genome_id = ani_cami_strdup_field(nf > 5 ? fields[5] : "NA"),
			.otu = ani_cami_strdup_field(nf > 6 ? fields[6] : "NA"),
		};
		kv_push(ani_cami_tax_record_t, *records, rec);
	}
	free(line);
	if (fclose(fp) != 0)
		err(errno, "%s(): cannot close --cami-taxmap %s", __func__, path);
	if (kv_size(*records) == 0)
		errx(EINVAL, "%s(): --cami-taxmap %s contains no usable records", __func__, path);
	qsort(&kv_A(*records, 0), kv_size(*records), sizeof(kv_A(*records, 0)),
		  ani_cami_tax_record_key_cmp);
}

static const ani_cami_tax_record_t *ani_cami_find_tax_record(
	const kv_cami_tax_record_t *records, const char *key)
{
	if (!records || kv_size(*records) == 0 || !key || key[0] == '\0')
		return NULL;
	ani_cami_tax_record_t needle = {.key = (char *)key};
	return (const ani_cami_tax_record_t *)bsearch(
		&needle, &kv_A(*records, 0), kv_size(*records),
		sizeof(kv_A(*records, 0)), ani_cami_tax_record_key_cmp);
}

static bool ani_cami_find_tax_record_range(
	const kv_cami_tax_record_t *records, const char *key,
	size_t *begin, size_t *end)
{
	if (begin)
		*begin = 0;
	if (end)
		*end = 0;
	const ani_cami_tax_record_t *hit = ani_cami_find_tax_record(records, key);
	if (!hit)
		return false;
	const ani_cami_tax_record_t *base = &kv_A(*records, 0);
	size_t b = (size_t)(hit - base);
	size_t e = b + 1u;
	while (b > 0 && strcmp(kv_A(*records, b - 1u).key, key) == 0)
		--b;
	while (e < kv_size(*records) && strcmp(kv_A(*records, e).key, key) == 0)
		++e;
	if (begin)
		*begin = b;
	if (end)
		*end = e;
	return e > b;
}

static void ani_cami_copy_basename(const char *path, char *out, size_t out_size)
{
	if (!out || out_size == 0)
		return;
	out[0] = '\0';
	if (!path || path[0] == '\0')
		return;
	const char *base = strrchr(path, '/');
	base = base ? base + 1 : path;
	snprintf(out, out_size, "%s", base);
}

static bool ani_cami_extract_accession(const char *s, char *out, size_t out_size)
{
	if (!out || out_size == 0)
		return false;
	out[0] = '\0';
	if (!s)
		return false;

	const char *best = NULL;
	const char *prefixes[] = {"GCF_", "GCA_", "NC_", "NZ_"};
	for (size_t p = 0; p < sizeof(prefixes) / sizeof(prefixes[0]); ++p) {
		const char *hit = strstr(s, prefixes[p]);
		if (hit && (!best || hit < best))
			best = hit;
	}
	if (!best)
		return false;

	size_t n = 0;
	if (strncmp(best, "GCF_", 4) == 0 || strncmp(best, "GCA_", 4) == 0) {
		const char *p = best;
		for (int i = 0; i < 4 && p[i]; ++i)
			out[n++] = p[i];
		p += 4;
		while (*p && isdigit((unsigned char)*p) && n + 1 < out_size)
			out[n++] = *p++;
		if (*p == '.' && n + 1 < out_size) {
			out[n++] = *p++;
			while (*p && isdigit((unsigned char)*p) && n + 1 < out_size)
				out[n++] = *p++;
		}
	} else {
		const char *p = best;
		while (*p && (isalnum((unsigned char)*p) || *p == '_' || *p == '.') &&
			   n + 1 < out_size)
			out[n++] = *p++;
	}
	out[n] = '\0';
	return n > 0;
}

static void ani_cami_profile_add(kv_cami_profile_entry_t *entries,
								 const ani_cami_tax_record_t *tax,
								 long double abundance)
{
	if (!entries || !tax || abundance <= 0.0L)
		return;
	for (size_t i = 0; i < kv_size(*entries); ++i) {
		ani_cami_profile_entry_t *entry = &kv_A(*entries, i);
		if (strcmp(entry->tax->taxid, tax->taxid) == 0 &&
			strcmp(entry->tax->rank, tax->rank) == 0 &&
			strcmp(entry->tax->taxpath, tax->taxpath) == 0) {
			entry->abundance += abundance;
			return;
		}
	}
	ani_cami_profile_entry_t entry = {.tax = tax, .abundance = abundance};
	kv_push(ani_cami_profile_entry_t, *entries, entry);
}

static size_t ani_cami_profile_add_key_matches(
	const kv_cami_tax_record_t *records,
	const char *key,
	kv_cami_profile_entry_t *entries,
	long double abundance)
{
	size_t begin = 0;
	size_t end = 0;
	if (!ani_cami_find_tax_record_range(records, key, &begin, &end))
		return 0;
	for (size_t i = begin; i < end; ++i)
		ani_cami_profile_add(entries, &kv_A(*records, i), abundance);
	return end - begin;
}

static size_t ani_cami_profile_add_ref_tax_matches(
	const kv_cami_tax_record_t *records,
	const char *ref_name,
	const char *annotation,
	kv_cami_profile_entry_t *entries,
	long double abundance)
{
	size_t n = ani_cami_profile_add_key_matches(records, ref_name,
												entries, abundance);
	if (n)
		return n;
	char buf[PATHLEN];
	ani_cami_copy_basename(ref_name, buf, sizeof(buf));
	n = ani_cami_profile_add_key_matches(records, buf, entries, abundance);
	if (n)
		return n;
	if (ani_cami_extract_accession(ref_name, buf, sizeof(buf))) {
		n = ani_cami_profile_add_key_matches(records, buf, entries, abundance);
		if (n)
			return n;
	}
	n = ani_cami_profile_add_key_matches(records, annotation,
										 entries, abundance);
	if (n)
		return n;
	if (ani_cami_extract_accession(annotation, buf, sizeof(buf)))
		return ani_cami_profile_add_key_matches(records, buf, entries, abundance);
	return 0;
}

static int ani_cami_rank_index(const char *rank)
{
	static const char *ranks[] = {
		"superkingdom", "phylum", "class", "order",
		"family", "genus", "species", "strain"};
	for (int i = 0; i < (int)(sizeof(ranks) / sizeof(ranks[0])); ++i)
		if (rank && strcasecmp(rank, ranks[i]) == 0)
			return i;
	return -1;
}

static void ani_cami_write_ranks_header(FILE *fp,
										const kv_cami_profile_entry_t *entries)
{
	static const char *ranks[] = {
		"superkingdom", "phylum", "class", "order",
		"family", "genus", "species", "strain"};
	bool present[8] = {0};
	for (size_t i = 0; entries && i < kv_size(*entries); ++i) {
		int idx = ani_cami_rank_index(kv_A(*entries, i).tax->rank);
		if (idx >= 0)
			present[idx] = true;
	}
	bool any = false;
	fputs("@Ranks:", fp);
	for (int i = 0; i < 8; ++i) {
		if (!present[i])
			continue;
		fprintf(fp, "%s%s", any ? "|" : "", ranks[i]);
		any = true;
	}
	if (!any)
		fputs("species", fp);
	fputc('\n', fp);
}

static void ani_cami_default_sample_id(const char *query_path, char *out,
									   size_t out_size)
{
	if (!out || out_size == 0)
		return;
	if (!query_path || query_path[0] == '\0' || strcmp(query_path, "-") == 0) {
		snprintf(out, out_size, "minco_sample");
		return;
	}
	ani_cami_copy_basename(query_path, out, out_size);
	if (out[0] == '\0')
		snprintf(out, out_size, "minco_sample");
}

static void ani_write_cami_profile(const ani_opt_t *ani_opt,
								   const char *query_path,
								   char (*refname)[PATHLEN],
								   char (*refanno)[PATHLEN],
								   const kv_ani_row_t *rows,
								   size_t out_n)
{
	if (!ani_opt || ani_opt->cami_profile[0] == '\0')
		return;
	if (ani_opt->abundance_model == ANI_ABUNDANCE_NONE)
		errx(EINVAL, "%s(): --cami-profile requires --abundance-est depth", __func__);
	if (ani_opt->cami_taxmap[0] == '\0')
		errx(EINVAL, "%s(): --cami-profile requires --cami-taxmap", __func__);

	kv_cami_tax_record_t taxmap;
	ani_cami_load_taxmap(ani_opt->cami_taxmap, &taxmap);
	kv_cami_profile_entry_t entries;
	kv_init(entries);
	size_t missing = 0;
	for (size_t i = 0; rows && i < out_n && i < kv_size(*rows); ++i) {
		const ani_row_t *row = &kv_A(*rows, i);
		const char *annotation = annotation_at(refanno, row->rn);
		const size_t added = ani_cami_profile_add_ref_tax_matches(
			&taxmap, refname[row->rn], annotation, &entries,
			(long double)row->abundance_normalized_depth);
		if (!added) {
			++missing;
			continue;
		}
	}

	FILE *fp = strcmp(ani_opt->cami_profile, "-") == 0
				   ? stdout
				   : fopen(ani_opt->cami_profile, "w");
	if (!fp)
		err(errno, "%s(): cannot open --cami-profile %s", __func__,
			ani_opt->cami_profile);

	char sample_id[PATHLEN];
	if (ani_opt->cami_sample_id[0] != '\0')
		snprintf(sample_id, sizeof(sample_id), "%s", ani_opt->cami_sample_id);
	else
		ani_cami_default_sample_id(query_path, sample_id, sizeof(sample_id));

	fprintf(fp, "@SampleID:%s\n", sample_id);
	fputs("@Version:0.9.1\n", fp);
	ani_cami_write_ranks_header(fp, &entries);
	fputs("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\t_CAMI_genomeID\t_CAMI_OTU\n", fp);
	for (size_t i = 0; i < kv_size(entries); ++i) {
		const ani_cami_profile_entry_t *entry = &kv_A(entries, i);
		const double pct = (double)(entry->abundance * 100.0L);
		if (pct <= 0.0)
			continue;
		fprintf(fp, "%s\t%s\t%s\t%s\t%.10g\t%s\t%s\n",
				entry->tax->taxid,
				entry->tax->rank,
				entry->tax->taxpath,
				entry->tax->taxpathsn,
				pct,
				entry->tax->genome_id,
				entry->tax->otu);
	}
	if (fp != stdout && fclose(fp) != 0)
		err(errno, "%s(): cannot close --cami-profile %s", __func__,
			ani_opt->cami_profile);
	if (missing)
		warnx("%s(): skipped %zu printed ANI rows missing --cami-taxmap records",
			  __func__, missing);
	if (kv_size(entries) == 0)
		warnx("%s(): CAMI profile contains only headers because no printed rows matched --cami-taxmap",
			  __func__);
	kv_destroy(entries);
	ani_cami_tax_records_destroy(&taxmap);
}

/* --- tiny helpers --- */
static inline const char *ani_best_confidence_label(const ani_row_t *r)
{
	return r->best_guarded ? "guarded_low_confidence" : ani_confidence_label(r->confidence);
}

static inline void fill_row_calibration(ani_row_t *r, bool unassembled,
                                        bool enable_best_guard,
                                        const infile_meta_t *qry_asm,
                                        const infile_meta_t *ref_asm)
{
    if (unassembled) {
        r->calibrated_ani = r->ani;
        r->best_ani = r->ani;
        r->confidence = ANI_CONF_HIGH;
        r->best_guarded = 0;
        return;
    }
    r->calibrated_ani = refaf_hgb_predict_ani(r->ani, r->af_ref, (unsigned int)r->XnY_ctx,
                                              (unsigned int)r->N_diff_obj,
                                              (unsigned int)r->N_diff_obj_section,
                                              (unsigned int)r->N_mut2_ctx);
    r->confidence = ani_confidence_from_values(r->ani, r->calibrated_ani, r->af_ref);
    r->best_ani = r->calibrated_ani;
    r->best_guarded = 0;

    const bool query_assembly = has_assembly_meta_record(qry_asm);
    const double guard_af_qry = ani_row_report_af_qry(r);
    const double guard_af_ref = ani_row_report_af_ref(r);
    const double min_af = guard_af_qry < guard_af_ref ? guard_af_qry : guard_af_ref;
    if (!enable_best_guard ||
        !r->real_af_available ||
        !query_assembly ||
        !infile_meta_complete_like_assembly(qry_asm) ||
        r->ani < 0.958 ||
        r->calibrated_ani < 0.962 ||
        min_af < 0.25 || min_af >= 0.50 ||
        r->XnY_ctx <= 0 || guard_af_qry <= 0.0 || guard_af_ref <= 0.0)
        return;

    const double qry_ctx = (double)r->XnY_ctx / guard_af_qry;
    const double ref_ctx = (double)r->XnY_ctx / guard_af_ref;
    double exact = (double)r->XnY_ctx - (double)r->N_diff_obj;
    if (exact < 0.0)
        exact = 0.0;
    const double ctx_exact_mean =
        (mash_ani_from_counts((double)r->XnY_ctx, qry_ctx, ref_ctx) +
         aaf_ani_from_counts((double)r->XnY_ctx, qry_ctx, ref_ctx) +
         mash_ani_from_counts(exact, qry_ctx, ref_ctx) +
         aaf_ani_from_counts(exact, qry_ctx, ref_ctx)) / 4.0;
    if (r->ani - ctx_exact_mean < 0.015)
        return;

    r->best_ani = aaf_ani_from_counts(exact, qry_ctx, ref_ctx);
    r->best_guarded = 1;

    (void)ref_asm;
}

static inline int selected_metric_request_id(const ani_opt_t *ani_opt)
{
	int s = ani_opt ? abs(ani_opt->s) : 1;
	if (s < 1 || s > ANI_SELECTED_METRIC_COUNT)
		s = 1;
	return s;
}

static inline bool ani_value_available(double value)
{
	return value > 0.0 && value <= 1.0;
}

static inline bool row_best_available(const ani_row_t *row, const ani_opt_t *ani_opt)
{
	return row && ani_opt && !ani_opt->raw_output && !ani_opt->unassembled &&
		   ani_value_available(row->best_ani);
}

static inline bool row_recal_available(const ani_row_t *row, const ani_opt_t *ani_opt)
{
	return row && ani_opt && !ani_opt->raw_output && !ani_opt->unassembled &&
		   ani_value_available(row->calibrated_ani);
}

static inline int selected_metric_effective_id(const ani_row_t *row, const ani_opt_t *ani_opt)
{
	int s = selected_metric_request_id(ani_opt);
	if (ani_opt && ani_opt->unassembled && !ani_opt->unified_metric && s <= 4)
		return 4;
	if (s == 1 && !row_best_available(row, ani_opt))
		s = 2;
	if (s == 2 && !row_recal_available(row, ani_opt))
		s = 3;
	return s;
}

static inline double selected_metric_distance_from_row(const ani_row_t *row, const ani_opt_t *ani_opt)
{
	if (row && ani_opt && ani_opt->readwise_query &&
		ani_opt->readwise_ani_model == ANI_READWISE_ANI_ZIP_AAF &&
		row->abundance_ref_zip_aaf_ani > 0.0)
		return 1.0 - row->abundance_ref_zip_aaf_ani;
	switch (selected_metric_effective_id(row, ani_opt)) {
	case 1:
		return 1.0 - row->best_ani;
	case 2:
		return 1.0 - row->calibrated_ani;
	case 3:
		return row->moe_dist;
	case 4:
		return row->naive_dist;
	case 5:
		return row->mash_dist;
	case 6:
		return row->aaf_dist;
	case 7:
		return row->af_pass ? row->moe_dist : row->mash_dist;
	case 8:
		return row->af_pass ? row->moe_dist : row->aaf_dist;
	default:
		return row->metric;
	}
}

static inline const char *selected_metric_name_for_row(const ani_row_t *row, const ani_opt_t *ani_opt)
{
	if (row && ani_opt && ani_opt->readwise_query &&
		ani_opt->readwise_ani_model == ANI_READWISE_ANI_ZIP_AAF &&
		row->abundance_ref_zip_aaf_ani > 0.0)
		return "ZipAaf";
	return select_metrics_header[selected_metric_effective_id(row, ani_opt) - 1];
}

static inline const char *selected_metric_confidence_label(const ani_row_t *row, const ani_opt_t *ani_opt)
{
	if (!ani_opt || ani_opt->raw_output)
		return "raw";
	if (row && ani_opt->readwise_query &&
		ani_opt->readwise_ani_model == ANI_READWISE_ANI_ZIP_AAF &&
		row->abundance_ref_zip_aaf_ani > 0.0)
		return "readwise_zip";
	if (ani_opt->unassembled && !ani_opt->unified_metric &&
		selected_metric_request_id(ani_opt) <= 4)
		return "unassembled";
	switch (selected_metric_effective_id(row, ani_opt)) {
	case 1:
		return ani_best_confidence_label(row);
	case 2:
		return ani_confidence_label(row->confidence);
	default:
		return "raw";
	}
}

static inline void fill_row_base_distances(ani_row_t *row, const ani_features_t *features,
										   const ani_opt_t *ani_opt,
										   uint32_t qry_ctx, uint32_t ref_ctx,
										   double af_q, double af_r)
{
	row->af_pass = ani_report_af_pass(ani_opt, af_q, af_r);
	if (row->af_pass) {
		ani_features_t tmp = *features;
		row->moe_dist = lm3ways_dist_from_features(&tmp);
		tmp = *features;
		row->naive_dist = get_naive_dist(&tmp);
	} else {
		row->moe_dist = ani_opt->e;
		row->naive_dist = ani_opt->e;
	}
	row->mash_dist = get_mashD(Bitslen.ctx / 2, ref_ctx, qry_ctx, features->XnY_ctx);
	row->aaf_dist = get_aafD(Bitslen.ctx / 2, ref_ctx, qry_ctx, features->XnY_ctx);
}

static inline void finalize_row_selected_metric(ani_row_t *row, const ani_opt_t *ani_opt)
{
	row->metric = selected_metric_distance_from_row(row, ani_opt);
	row->selected_ani = 1.0 - row->metric;
}

static inline ani_row_t make_selected_output_row(uint32_t rn,
												 const ani_features_t *features,
												 const ani_opt_t *ani_opt,
												 uint32_t qry_ctx,
												 uint32_t ref_ctx,
												 double af_qry,
												 double blastn_af_qry,
												 double af_ref,
												 double blastn_af_ref,
												 ani_density_af_t density_af,
												 const infile_meta_t *qry_asm,
												 const infile_meta_t *ref_asm)
{
	ani_row_t row;
	memset(&row, 0, sizeof(row));
	row.rn = rn;
	row.af_qry = af_qry;
	row.blastn_af_qry = blastn_af_qry;
	row.af_ref = af_ref;
	row.blastn_af_ref = blastn_af_ref;
	row.real_af_qry = density_af.qry;
	row.real_af_ref = density_af.ref;
	row.real_af_available = density_af.available;
	row.XnY_ctx = features->XnY_ctx;
	row.N_diff_obj = features->N_diff_obj;
	row.N_diff_obj_section = features->N_diff_obj_section;
	row.N_mut2_ctx = features->N_mut2_ctx;
	fill_row_base_distances(&row, features, ani_opt, qry_ctx, ref_ctx,
							ani_row_report_af_qry(&row),
							ani_row_report_af_ref(&row));
	row.ani = 1.0 - (ani_opt->v ? row.naive_dist : row.moe_dist);
	if (!ani_opt->raw_output && row.af_pass)
		fill_row_calibration(&row, ani_opt->unassembled, ani_best_guard_enabled(ani_opt),
							 qry_asm, ref_asm);
	finalize_row_selected_metric(&row, ani_opt);
	return row;
}

static inline void append_unified_detail_row(kstring_t *ks_out,
                                             const ani_opt_t *ani_opt,
                                             const char *qry_name,
                                             const char *ref_name,
                                             const ani_row_t *r,
                                             const char *annotation)
{
    const double selected_distance = r->metric;
    const double selected_similarity = r->selected_ani;
    const char *metric_name = selected_metric_name_for_row(r, ani_opt);
    const char *confidence = selected_metric_confidence_label(r, ani_opt);
    const char *ref_annotation = annotation_or_na(annotation);
    const double blastn_af_qry = bounded_align_fraction(r->blastn_af_qry);
    const double blastn_af_ref = bounded_align_fraction(r->blastn_af_ref);
    const double real_af_qry = bounded_align_fraction(r->real_af_qry);
    const double real_af_ref = bounded_align_fraction(r->real_af_ref);
    const double real_min_af = real_af_qry < real_af_ref ? real_af_qry : real_af_ref;
    const char *af_source = (ani_opt && ani_opt->readwise_query)
        ? "readwise_coverage"
        : (r->real_af_available ? "density" : "sketch");

    ksprintf(ks_out, "%s\t%s\t%lf\t%lf\t%s\t%s\t%d\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%s\t%f\t%f\t%f\t%s",
             qry_name, ref_name, selected_similarity, selected_distance, confidence, metric_name,
             r->XnY_ctx, r->af_qry, blastn_af_qry, r->af_ref, blastn_af_ref,
             r->N_diff_obj, r->N_diff_obj_section, r->N_mut2_ctx, ref_annotation,
             real_af_qry, real_af_ref, real_min_af, af_source);
    if (ani_opt && ani_opt->readwise_query) {
        const double read_match_fraction = r->readwise_total_reads
            ? (double)r->readwise_reads_with_ctx_match / (double)r->readwise_total_reads
            : 0.0;
        const double block_match_fraction = r->readwise_total_density_blocks
            ? (double)r->readwise_blocks_with_ctx_match / (double)r->readwise_total_density_blocks
            : 0.0;
        ksprintf(ks_out, "\t%" PRIu64 "\t%" PRIu64 "\t%f\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%f\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%f\t%f\t%f",
                 r->readwise_reads_with_ctx_match,
                 r->readwise_total_reads,
                 read_match_fraction,
                 r->readwise_unique_query_ctx,
                 r->readwise_unique_query_ctx_hit,
                 r->readwise_unique_ref_ctx_hit,
                 r->readwise_density_block_ctx,
                 r->readwise_total_density_blocks,
                 r->readwise_blocks_with_ctx_match,
                 block_match_fraction,
                 r->readwise_raw_xny_ctx,
                 r->readwise_rejected_ctx,
                 r->readwise_rejected_diff_ctx,
                 r->readwise_fake_ctx_fraction,
                 r->readwise_fake_ctx_prob_mean,
                 r->readwise_fake_ctx_prob_weighted);
        if (ani_opt->abundance_model != ANI_ABUNDANCE_NONE) {
			const ani_readwise_default_call_t default_call =
				ani_readwise_default_call(r);
            ksprintf(ks_out, "\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%" PRIu64 "\t%f\t%f\t%f\t%f\t%s\t%s",
                     r->abundance_ref_breadth,
                     r->abundance_ref_mean_depth,
                     r->abundance_ref_hit_mean_depth,
                     r->abundance_ref_depth_variance,
                     r->abundance_ref_depth_cv,
                     r->abundance_ref_zero_fraction,
                     r->abundance_relative_depth,
                     r->abundance_normalized_depth,
                     r->abundance_effective_depth,
                     r->abundance_normalized_effective_depth,
					 r->abundance_ref_zip_af,
					 r->abundance_ref_zip_aaf_ani,
                     r->reliable_ref_breadth,
                     r->reliable_ref_mean_depth,
                     r->reliable_ref_hit_ctx,
                     r->reliable_ref_hit_mean_depth,
					 r->reliable_ref_hit_median_depth,
                     r->reliable_ref_hit_depth_variance,
                     r->reliable_ref_zip_af,
					 ani_readwise_default_call_label(default_call),
					 ani_readwise_default_call_rule(r, default_call));
			if (ani_opt->readwise_dual_evidence) {
				ksprintf(ks_out, "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%" PRIu64 "\t%f\t%f\t%f\t%f",
						 r->marker_xny_ctx,
						 r->marker_raw_xny_ctx,
						 r->marker_n_diff_obj,
						 r->marker_n_diff_obj_section,
						 r->marker_n_mut2_ctx,
						 r->marker_ref_ctx_total,
						 r->marker_ref_ctx_total,
						 r->marker_ref_breadth,
						 r->marker_ref_mean_depth,
						 r->marker_ref_hit_mean_depth,
						 r->marker_ref_depth_variance,
						 r->marker_ref_depth_cv,
						 r->marker_ref_zero_fraction,
						 r->marker_relative_depth,
						 r->marker_effective_depth,
						 r->marker_ref_zip_af,
						 r->marker_ref_zip_aaf_ani,
						 r->marker_reliable_ref_breadth,
						 r->marker_reliable_ref_mean_depth,
						 r->marker_reliable_ref_hit_ctx,
						 r->marker_reliable_ref_hit_mean_depth,
						 r->marker_reliable_ref_hit_median_depth,
						 r->marker_reliable_ref_hit_depth_variance,
						 r->marker_reliable_ref_zip_af);
			}
        }
    }
    kputc('\n', ks_out);
}

static inline void print_unified_detail_row(FILE *outfp,
                                            const ani_opt_t *ani_opt,
                                            const char *qry_name,
                                            const char *ref_name,
                                            const ani_row_t *r,
                                            const char *annotation)
{
    kstring_t ks = (kstring_t){0, 0, 0};
    append_unified_detail_row(&ks, ani_opt, qry_name, ref_name, r, annotation);
    if (ks.l)
        fwrite(ks.s, 1, ks.l, outfp);
	free(ks.s);
}

static inline double ani_matrix_exception_value(const ani_opt_t *ani_opt)
{
	return ani_opt->s < 0 ? 1.0 - ani_opt->e : ani_opt->e;
}

static inline double ani_matrix_diagonal_value(const ani_opt_t *ani_opt)
{
	return ani_opt->s < 0 ? 1.0 : 0.0;
}

static double ani_matrix_pair_value(const unify_sketch_t *qry_result, uint32_t qn,
									const unify_sketch_t *ref_result, uint32_t rn,
									const ani_opt_t *ani_opt,
									const ani_ctxmeta_rec_t *qry_ctxmeta,
									const ani_ctxmeta_rec_t *ref_ctxmeta)
{
	uint64_t *arr_qry = qry_result->comb_sketch + qry_result->sketch_index[qn];
	size_t len_qry = qry_result->sketch_index[qn + 1] - qry_result->sketch_index[qn];
	const uint32_t qry_ctx = qry_result->conflict
		? count_ctx_runs_sorted_ctxobj64_local(arr_qry, len_qry, false)
		: (uint32_t)len_qry;
	uint64_t *arr_ref = ref_result->comb_sketch + ref_result->sketch_index[rn];
	size_t len_ref = ref_result->sketch_index[rn + 1] - ref_result->sketch_index[rn];
	const uint32_t ref_ctx = ref_result->conflict
		? count_ctx_runs_sorted_ctxobj64_local(arr_ref, len_ref, ani_opt->ignoreconflict)
		: (uint32_t)len_ref;

	if (qry_ctx == 0 || ref_ctx == 0)
		return ani_matrix_exception_value(ani_opt);

	ani_features_t ani_features;
	get_ani_features_ctx_min_over_conflicts_both_filtered(
		arr_qry, len_qry, arr_ref, len_ref,
		ani_opt->ignoreconflict && ref_result->conflict,
		&ani_features);
	const double af_qry = (double)ani_features.XnY_ctx / (double)qry_ctx;
	const double af_ref = (double)ani_features.XnY_ctx / (double)ref_ctx;
	ani_features_t tmp = ani_features;
	tmp.X_ctx = qry_ctx;
	const double blastn_af_qry = lm3ways_af_ANIb_from_features(&tmp);
	tmp = ani_features;
	tmp.X_ctx = ref_ctx;
	const double blastn_af_ref = lm3ways_af_ANIb_from_features(&tmp);
	const ani_density_af_t density_af = ani_estimate_density_af(
		ani_ctxmeta_at(qry_ctxmeta, qn), ani_ctxmeta_at(ref_ctxmeta, rn),
		(uint32_t)ani_features.XnY_ctx, af_qry, af_ref);
	ani_row_t outrow = make_selected_output_row(rn, &ani_features, ani_opt,
												qry_ctx, ref_ctx,
												af_qry, blastn_af_qry,
												af_ref, blastn_af_ref,
												density_af,
												infile_meta_at(qry_result, qn),
												infile_meta_at(ref_result, rn));
	return ani_opt->s < 0 ? outrow.selected_ani : outrow.metric;
}

static int ani_choose_dense_matrix_block_size(int ref_n, int qry_n)
{
	if (ref_n <= 0 || qry_n <= 0)
		return 0;
	const size_t bytes_per_query =
		(size_t)ref_n * (sizeof(ctx_mut2_t) + sizeof(obj_section_t));
	if (bytes_per_query == 0)
		return 0;
	size_t budget = GetAvailableMemory() / 4;
	if (budget < bytes_per_query * 64)
		budget = bytes_per_query * 64;
	size_t block = budget / bytes_per_query;
	if (block < 64)
		block = 64;
	if (block > 64)
		block = 64;
	if (block > (size_t)qry_n)
		block = (size_t)qry_n;

	const char *env = getenv("MINCO_ANI_MATRIX_BLOCK_SIZE");
	if (!env || env[0] == '\0')
		env = getenv("MINCO_MATRIX_BLOCK_SIZE");
	if (env && env[0] != '\0') {
		char *end = NULL;
		errno = 0;
		long requested = strtol(env, &end, 10);
		if (errno == 0 && end != env && *end == '\0' && requested > 0) {
			if (requested > qry_n)
				requested = qry_n;
			block = (size_t)requested;
		}
	}
	return (int)block;
}

static uint32_t *ani_ctx_counts_for_sketch(const unify_sketch_t *sketch, bool ignoreconflict)
{
	uint32_t *counts = calloc((size_t)sketch->infile_num, sizeof(counts[0]));
	if (!counts)
		err(EXIT_FAILURE, "%s(): calloc ctx counts", __func__);
	for (int i = 0; i < sketch->infile_num; ++i) {
		const uint64_t begin = sketch->sketch_index[i];
		const uint64_t end = sketch->sketch_index[i + 1];
		const uint64_t *arr = sketch->comb_sketch + begin;
		const size_t len = (size_t)(end - begin);
		counts[i] = sketch->conflict
						? count_ctx_runs_sorted_ctxobj64_local(arr, len, ignoreconflict)
						: (uint32_t)len;
	}
	return counts;
}

static ctxgidobj_t *load_self_sorted_index_or_build(const unify_sketch_t *sketch,
													const char *sketch_dir,
													size_t *index_bytes,
													bool *from_file,
													bool *is_mmap)
{
	*index_bytes = 0;
	*from_file = false;
	*is_mmap = false;
	const size_t total_entries = (size_t)sketch->sketch_index[sketch->infile_num];
	if (sketch_dir && sketch_dir[0] != '\0' &&
		file_exists_in_folder((char *)sketch_dir, (char *)sorted_comb_ctxgid64obj32)) {
		char *index_path = test_get_fullpath(sketch_dir, sorted_comb_ctxgid64obj32);
		ctxgidobj_t *index = read_reference_sorted_index(index_path, index_bytes, is_mmap);
		free(index_path);
		if (*index_bytes != total_entries * sizeof(index[0]))
			errx(EINVAL, "%s(): sorted index size mismatch", __func__);
		*from_file = true;
		return index;
	}
	if (force_ref_index_requested())
		errx(EXIT_FAILURE, "MINCO_FORCE_REF_INDEX requires reference index '%s/%s'",
			 sketch_dir, sorted_comb_ctxgid64obj32);
	*index_bytes = total_entries * sizeof(ctxgidobj_t);
	return comb_sortedsketch64_2sortedcomb_ctxgid64obj32((unify_sketch_t *)sketch);
}

static void free_self_sorted_index(ctxgidobj_t *index, size_t index_bytes,
								   bool from_file, bool is_mmap)
{
	if (!index)
		return;
	if (from_file)
		free_reference_sorted_index(index, index_bytes, is_mmap);
	else
		free(index);
}

static inline size_t ani_lower_offset(int row, int col)
{
	return (size_t)row * (size_t)(row - 1) / 2u + (size_t)col;
}

static int32_t ani_scale_fixed6(double value)
{
	if (!isfinite(value))
		value = 1.0;
	if (value > 2147.0)
		value = 2147.0;
	if (value < -2147.0)
		value = -2147.0;
	return (int32_t)(value >= 0.0
					 ? value * 1000000.0 + 0.5
					 : value * 1000000.0 - 0.5);
}

static void ani_kput_scaled6(kstring_t *s, int32_t scaled)
{
	uint32_t mag;
	if (scaled < 0) {
		kputc('-', s);
		mag = (uint32_t)(-scaled);
	} else {
		mag = (uint32_t)scaled;
	}
	const uint32_t whole = mag / 1000000U;
	uint32_t frac = mag % 1000000U;
	kputuw(whole, s);
	kputc('.', s);
	char frac_buf[6];
	for (int i = 5; i >= 0; --i) {
		frac_buf[i] = (char)('0' + (frac % 10));
		frac /= 10;
	}
	kputsn(frac_buf, 6, s);
}

static double ani_indexed_self_pair_metric(const unify_sketch_t *sketch,
										   const ani_opt_t *ani_opt,
										   int qgid, int rgid,
										   const ani_features_t *features,
										   const uint32_t *qry_ctx_count,
										   const uint32_t *ref_ctx_count,
										   const ani_ctxmeta_rec_t *ctxmeta)
{
	const uint32_t qry_ctx = qry_ctx_count[qgid];
	const uint32_t ref_ctx = ref_ctx_count[rgid];
	double metric = ani_matrix_exception_value(ani_opt);
	if (qry_ctx > 0 && ref_ctx > 0) {
		const double af_qry = (double)features->XnY_ctx / (double)qry_ctx;
		const double af_ref = (double)features->XnY_ctx / (double)ref_ctx;
		ani_features_t tmp = *features;
		tmp.X_ctx = qry_ctx;
		const double blastn_af_qry = lm3ways_af_ANIb_from_features(&tmp);
		tmp = *features;
		tmp.X_ctx = ref_ctx;
		const double blastn_af_ref = lm3ways_af_ANIb_from_features(&tmp);
		const ani_density_af_t density_af = ani_estimate_density_af(
			ani_ctxmeta_at(ctxmeta, (uint32_t)qgid), ani_ctxmeta_at(ctxmeta, (uint32_t)rgid),
			(uint32_t)features->XnY_ctx, af_qry, af_ref);
		ani_row_t outrow = make_selected_output_row(
			(uint32_t)rgid, features, ani_opt,
			qry_ctx, ref_ctx, af_qry, blastn_af_qry,
			af_ref, blastn_af_ref,
			density_af,
			infile_meta_at(sketch, (uint32_t)qgid),
			infile_meta_at(sketch, (uint32_t)rgid));
		metric = ani_opt->s < 0 ? outrow.selected_ani : outrow.metric;
	}
	return metric;
}

bool comb_sortedsketch64_indexed_self_full(ani_opt_t *ani_opt)
{
	unify_sketch_t *sketch = generic_sketch_parse(ani_opt->qrydir, SKETCH_PARSE_NONE);
	if (sketch->stat_type != 2) {
		free_unify_sketch(sketch);
		return false;
	}

	const int n = sketch->infile_num;
	size_t index_bytes = 0;
	bool index_from_file = false;
	bool index_is_mmap = false;
	ctxgidobj_t *sorted_index =
		load_self_sorted_index_or_build(sketch, ani_opt->qrydir,
										&index_bytes, &index_from_file,
										&index_is_mmap);
	const size_t ref_sksize = (size_t)sketch->sketch_index[n];
	uint32_t *qry_ctx_count = ani_ctx_counts_for_sketch(sketch, false);
	uint32_t *ref_ctx_count = ani_opt->ignoreconflict
								  ? ani_ctx_counts_for_sketch(sketch, true)
								  : qry_ctx_count;
	ani_ctxmeta_rec_t *ctxmeta =
		read_optional_ani_ctxmeta_stats(ani_opt->qrydir, sketch->infile_num);

	const int block_size = ani_choose_dense_matrix_block_size(n, n);
	if (block_size <= 0)
		errx(EINVAL, "%s(): invalid ANI full-matrix block size", __func__);
	ctx_mut2_t *ctx = malloc((size_t)n * (size_t)block_size * sizeof(ctx[0]));
	obj_section_t *obj = malloc((size_t)n * (size_t)block_size * sizeof(obj[0]));
	if (!ctx || !obj)
		err(EXIT_FAILURE, "%s(): OOM indexed ANI full-matrix block", __func__);

	const size_t lower_n = (size_t)n * (size_t)(n - 1) / 2u;
	int32_t *lower = lower_n > 0 ? malloc(lower_n * sizeof(lower[0])) : NULL;
	if (lower_n > 0 && !lower)
		err(EXIT_FAILURE, "%s(): OOM indexed ANI full lower matrix", __func__);

	ani_opt_t scan_opt = *ani_opt;
	scan_opt.p = ani_opt->p > 0 ? ani_opt->p : 1;
	for (int offset = 0; offset < n; offset += block_size) {
		const int this_block = (offset + block_size > n) ? n - offset : block_size;
		memset(ctx, 0, (size_t)n * (size_t)this_block * sizeof(ctx[0]));
		memset(obj, 0, (size_t)n * (size_t)this_block * sizeof(obj[0]));
		count_ctx_obj_frm_comb_sketch_section_lower(
			ctx, obj, sorted_index, ref_sksize, n, offset, this_block,
			sketch->comb_sketch + sketch->sketch_index[offset],
			sketch->sketch_index + offset, &scan_opt);

		for (int i = 0; i < this_block; ++i) {
			const int qgid = offset + i;
			for (int rgid = 0; rgid < qgid; ++rgid) {
				ani_features_t features = {
					.XnY_ctx = MCTX(n, i, rgid).num_ctx,
					.N_diff_obj = MOBJ(n, i, rgid).diff_obj,
					.N_diff_obj_section = MOBJ(n, i, rgid).diff_obj_section,
					.N_mut2_ctx = MCTX(n, i, rgid).num_mut2_ctx,
				};
				const double metric = ani_indexed_self_pair_metric(
					sketch, ani_opt, qgid, rgid, &features,
					qry_ctx_count, ref_ctx_count, ctxmeta);
				lower[ani_lower_offset(qgid, rgid)] = ani_scale_fixed6(metric);
			}
		}
	}

	FILE *outfp = ani_opt->outf[0] == '\0' ? stdout : fopen(ani_opt->outf, "w");
	if (outfp == NULL)
		err(errno, "%s", ani_opt->outf);
	kstring_t row = {0, 0, NULL};
	for (int rn = 0; rn < n; ++rn) {
		kputc('\t', &row);
		kputs(sketch->gname[rn], &row);
	}
	kputc('\n', &row);
	fwrite(row.s, 1, row.l, outfp);

	const int32_t diag_scaled = ani_scale_fixed6(ani_matrix_diagonal_value(ani_opt));
	for (int qgid = 0; qgid < n; ++qgid) {
		row.l = 0;
		kputs(sketch->gname[qgid], &row);
		for (int rgid = 0; rgid < n; ++rgid) {
			int32_t scaled = diag_scaled;
			if (qgid > rgid)
				scaled = lower[ani_lower_offset(qgid, rgid)];
			else if (qgid < rgid)
				scaled = lower[ani_lower_offset(rgid, qgid)];
			kputc('\t', &row);
			ani_kput_scaled6(&row, scaled);
		}
		kputc('\n', &row);
		fwrite(row.s, 1, row.l, outfp);
	}

	if (outfp != stdout)
		fclose(outfp);
	free(row.s);
	free(lower);
	free(ctx);
	free(obj);
	free(qry_ctx_count);
	if (ref_ctx_count != qry_ctx_count)
		free(ref_ctx_count);
	free(ctxmeta);
	free_self_sorted_index(sorted_index, index_bytes, index_from_file, index_is_mmap);
	free_unify_sketch(sketch);
	return true;
}

bool comb_sortedsketch64_indexed_self_triangle(ani_opt_t *ani_opt)
{
	unify_sketch_t *sketch = generic_sketch_parse(ani_opt->qrydir, SKETCH_PARSE_NONE);
	if (sketch->stat_type != 2) {
		free_unify_sketch(sketch);
		return false;
	}

	const int n = sketch->infile_num;
	size_t index_bytes = 0;
	bool index_from_file = false;
	bool index_is_mmap = false;
	ctxgidobj_t *sorted_index =
		load_self_sorted_index_or_build(sketch, ani_opt->qrydir,
										&index_bytes, &index_from_file,
										&index_is_mmap);
	const size_t ref_sksize = (size_t)sketch->sketch_index[n];
	uint32_t *qry_ctx_count = ani_ctx_counts_for_sketch(sketch, false);
	uint32_t *ref_ctx_count = ani_opt->ignoreconflict
								  ? ani_ctx_counts_for_sketch(sketch, true)
								  : qry_ctx_count;
	ani_ctxmeta_rec_t *ctxmeta =
		read_optional_ani_ctxmeta_stats(ani_opt->qrydir, sketch->infile_num);

	const int block_size = ani_choose_dense_matrix_block_size(n, n);
	if (block_size <= 0)
		errx(EINVAL, "%s(): invalid ANI triangle block size", __func__);
	ctx_mut2_t *ctx = malloc((size_t)n * (size_t)block_size * sizeof(ctx[0]));
	obj_section_t *obj = malloc((size_t)n * (size_t)block_size * sizeof(obj[0]));
	if (!ctx || !obj)
		err(EXIT_FAILURE, "%s(): OOM indexed ANI triangle block", __func__);

	ani_opt_t scan_opt = *ani_opt;
	scan_opt.p = ani_opt->p > 0 ? ani_opt->p : 1;
	kstring_t row = {0, 0, NULL};
	FILE *outfp = ani_opt->outf[0] == '\0' ? stdout : fopen(ani_opt->outf, "w");
	if (outfp == NULL)
		err(errno, "%s", ani_opt->outf);

	for (int offset = 0; offset < n; offset += block_size) {
		const int this_block = (offset + block_size > n) ? n - offset : block_size;
		memset(ctx, 0, (size_t)n * (size_t)this_block * sizeof(ctx[0]));
		memset(obj, 0, (size_t)n * (size_t)this_block * sizeof(obj[0]));
		count_ctx_obj_frm_comb_sketch_section_lower(
			ctx, obj, sorted_index, ref_sksize, n, offset, this_block,
			sketch->comb_sketch + sketch->sketch_index[offset],
			sketch->sketch_index + offset, &scan_opt);

		for (int i = 0; i < this_block; ++i) {
			const int qgid = offset + i;
			row.l = 0;
			kputs(sketch->gname[qgid], &row);
			for (int rgid = 0; rgid < qgid; ++rgid) {
				ani_features_t features = {
					.XnY_ctx = MCTX(n, i, rgid).num_ctx,
					.N_diff_obj = MOBJ(n, i, rgid).diff_obj,
					.N_diff_obj_section = MOBJ(n, i, rgid).diff_obj_section,
					.N_mut2_ctx = MCTX(n, i, rgid).num_mut2_ctx,
				};
				const double metric = ani_indexed_self_pair_metric(
					sketch, ani_opt, qgid, rgid, &features,
					qry_ctx_count, ref_ctx_count, ctxmeta);
				ksprintf(&row, "\t%lf", metric);
			}
			if (ani_opt->d)
				ksprintf(&row, "\t%lf", ani_matrix_diagonal_value(ani_opt));
			kputc('\n', &row);
			fwrite(row.s, 1, row.l, outfp);
		}
	}

	if (outfp != stdout)
		fclose(outfp);
	free(row.s);
	free(ctx);
	free(obj);
	free(qry_ctx_count);
	if (ref_ctx_count != qry_ctx_count)
		free(ref_ctx_count);
	free(ctxmeta);
	free_self_sorted_index(sorted_index, index_bytes, index_from_file, index_is_mmap);
	free_unify_sketch(sketch);
	return true;
}

static inline void kv_append_rows(kv_ani_row_t *dst, const kv_ani_row_t *src)
{
    const size_t n = kv_size(*src);
    if (!n) return;
    const size_t old = kv_size(*dst);
    kv_resize(ani_row_t, *dst, old + n);
    memcpy(&kv_A(*dst, old), &kv_A(*src, 0), n * sizeof(ani_row_t));
    kv_size(*dst) = old + n;
}

static inline void format_rows_to_kstr(const unify_sketch_t *qry, const unify_sketch_t *ref,
                                       uint32_t qn, const kv_ani_row_t *rows,
                                       const ani_opt_t *ani_opt, kstring_t *ks_out)
{
    for (size_t i = 0; i < kv_size(*rows); ++i) {
        const ani_row_t *r = &kv_A(*rows, i);
        append_unified_detail_row(ks_out, ani_opt, qry->gname[qn], ref->gname[r->rn],
                                  r, unify_annotation_at(ref, r->rn));
    }
}


/* Compute (qn,rn); return 1 if it passes filters and fill *row_out */
// called by void comb_sortedsketch64Xcomb_sortedsketch64_filter_and_sort_survivors(ani_opt_t *ani_opt);
static inline int compute_row_if_survivor(
    const unify_sketch_t *qry, const unify_sketch_t *ref,
    uint32_t qn, uint32_t rn, const ani_opt_t *opt,
    ani_row_t *row_out)
{
    uint64_t *arr_qry = qry->comb_sketch + qry->sketch_index[qn];
    size_t    len_qry = qry->sketch_index[qn + 1] - qry->sketch_index[qn];

    uint64_t *arr_ref = ref->comb_sketch + ref->sketch_index[rn];
    size_t    len_ref = ref->sketch_index[rn + 1] - ref->sketch_index[rn];

    ani_features_t f;
    /* Compare one shared context once, choosing the closest object pair when either side keeps conflicts. */
    get_ani_features_ctx_min_over_conflicts_both_filtered(arr_qry, len_qry, arr_ref, len_ref, opt->ignoreconflict && ref->conflict, &f);

    const uint32_t qry_ctx = qry->conflict ? count_ctx_runs_sorted_ctxobj64_local(arr_qry, len_qry, false) : (uint32_t)len_qry;
    const uint32_t ref_ctx = ref->conflict ? count_ctx_runs_sorted_ctxobj64_local(arr_ref, len_ref, opt->ignoreconflict) : (uint32_t)len_ref;
    if (qry_ctx == 0 || ref_ctx == 0) return 0;

    const double af_q = (double)f.XnY_ctx / (double)qry_ctx;
    const double af_r = (double)f.XnY_ctx / (double)ref_ctx;
	const ani_density_af_t density_af = ani_estimate_density_af(NULL, NULL,
		(uint32_t)f.XnY_ctx, af_q, af_r);

	if (!ani_report_af_pass(opt, density_af.qry, density_af.ref)) return 0;

    ani_features_t tmp = f;
    tmp.X_ctx = qry_ctx;
    const double blastn_af_q = lm3ways_af_ANIb_from_features(&tmp);
    tmp = f;
    tmp.X_ctx = ref_ctx;
    const double blastn_af_r = lm3ways_af_ANIb_from_features(&tmp);

    ani_row_t r = make_selected_output_row(rn, &f, opt, qry_ctx, ref_ctx,
                                           af_q, blastn_af_q, af_r, blastn_af_r,
                                           density_af,
                                           infile_meta_at(qry, qn), infile_meta_at(ref, rn));
    if (r.selected_ani <= opt->anicut) return 0;

    *row_out = r;
    return 1;
}

typedef struct {
    uint64_t ctx;
    size_t beg;
    size_t end;
} qry_ctx_run_t;

typedef struct {
    qry_ctx_run_t *runs;
    size_t n_runs;
    size_t *fence;
    int fence_k;
} qry_ctx_lookup_t;

typedef struct {
    uint64_t *keys;  /* ctx + 1; 0 means empty */
    uint64_t *vals;  /* packed beg/end */
    uint32_t mask;
    uint32_t n_runs;
} qry_ctx_hash_lookup_t;

typedef enum {
    QRAW_LOOKUP_SORTED = 0,
    QRAW_LOOKUP_HASH = 1
} qraw_lookup_mode_t;

static inline uint32_t next_pow2_u32_local(uint32_t x)
{
    if (x <= 1)
        return 1;
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return x + 1;
}

static inline uint64_t pack_qry_run_be(size_t beg, size_t end)
{
    return ((uint64_t)(uint32_t)beg << 32) | (uint64_t)(uint32_t)end;
}

static inline size_t unpack_qry_run_beg(uint64_t be) { return (size_t)(uint32_t)(be >> 32); }
static inline size_t unpack_qry_run_end(uint64_t be) { return (size_t)(uint32_t)be; }

static inline uint32_t top_k_bits_ctx_local(uint64_t ctx, int k)
{
    if (k <= 0)
        return 0;
    const int ctx_bits = Bitslen.ctx;
    if (k >= ctx_bits)
        return (uint32_t)ctx;
    return (uint32_t)(ctx >> (ctx_bits - k));
}

static int choose_qry_fence_k(size_t n_runs)
{
    int max_k = Bitslen.ctx < 20 ? Bitslen.ctx : 20;
    if (max_k < 0)
        max_k = 0;
    int k = max_k < 12 ? max_k : 12;
    while (k > 4 && ((size_t)1u << k) > n_runs * 4)
        k--;
    return k;
}

static void build_qry_ctx_lookup(const uint64_t *qry, size_t n, qry_ctx_lookup_t *lookup)
{
    memset(lookup, 0, sizeof(*lookup));
    if (n == 0)
        return;

    lookup->runs = malloc(n * sizeof(*lookup->runs));
    if (!lookup->runs)
        err(EXIT_FAILURE, "%s(): OOM query context runs", __func__);

    const uint8_t nobjbits = Bitslen.obj;
    size_t i = 0;
    while (i < n) {
        const uint64_t ctx = qry[i] >> nobjbits;
        const size_t beg = i;
        do { ++i; } while (i < n && (qry[i] >> nobjbits) == ctx);
        lookup->runs[lookup->n_runs++] = (qry_ctx_run_t){ctx, beg, i};
    }

    lookup->fence_k = choose_qry_fence_k(lookup->n_runs);
    const size_t buckets = (size_t)1u << lookup->fence_k;
    lookup->fence = malloc((buckets + 1) * sizeof(*lookup->fence));
    if (!lookup->fence)
        err(EXIT_FAILURE, "%s(): OOM query context fence", __func__);

    for (size_t t = 0; t <= buckets; ++t)
        lookup->fence[t] = lookup->n_runs;

    size_t next = 0;
    for (size_t r = 0; r < lookup->n_runs; ++r) {
        const uint32_t topk = top_k_bits_ctx_local(lookup->runs[r].ctx, lookup->fence_k);
        while (next <= topk && next <= buckets)
            lookup->fence[next++] = r;
        if (next > buckets)
            break;
    }
    while (next <= buckets)
        lookup->fence[next++] = lookup->n_runs;
}

static void free_qry_ctx_lookup(qry_ctx_lookup_t *lookup)
{
    if (!lookup)
        return;
    free(lookup->runs);
    free(lookup->fence);
    memset(lookup, 0, sizeof(*lookup));
}

static void build_qry_ctx_hash_lookup(const uint64_t *qry, size_t n, qry_ctx_hash_lookup_t *lookup)
{
    memset(lookup, 0, sizeof(*lookup));
    if (n == 0)
        return;
    if (n > UINT32_MAX)
        errx(EXIT_FAILURE, "%s(): query sketch has too many entries for packed hash lookup", __func__);

    const uint32_t cap = next_pow2_u32_local((uint32_t)n * 2u);
    lookup->keys = calloc(cap, sizeof(*lookup->keys));
    lookup->vals = malloc((size_t)cap * sizeof(*lookup->vals));
    if (!lookup->keys || !lookup->vals)
        err(EXIT_FAILURE, "%s(): OOM query context hash", __func__);
    lookup->mask = cap - 1;

    const uint8_t nobjbits = Bitslen.obj;
    size_t i = 0;
    while (i < n) {
        const uint64_t ctx = qry[i] >> nobjbits;
        const size_t beg = i;
        do { ++i; } while (i < n && (qry[i] >> nobjbits) == ctx);
        const uint64_t key = ctx + 1;
        uint32_t pos = (uint32_t)mix64(key) & lookup->mask;
        while (lookup->keys[pos] && lookup->keys[pos] != key)
            pos = (pos + 1) & lookup->mask;
        lookup->keys[pos] = key;
        lookup->vals[pos] = pack_qry_run_be(beg, i);
        lookup->n_runs++;
    }
}

static void free_qry_ctx_hash_lookup(qry_ctx_hash_lookup_t *lookup)
{
    if (!lookup)
        return;
    free(lookup->keys);
    free(lookup->vals);
    memset(lookup, 0, sizeof(*lookup));
}

static inline int lookup_qry_ctx_hash_run(const qry_ctx_hash_lookup_t *lookup, uint64_t ctx,
                                          size_t *beg, size_t *end)
{
    if (!lookup || lookup->n_runs == 0)
        return 0;
    const uint64_t key = ctx + 1;
    uint32_t pos = (uint32_t)mix64(key) & lookup->mask;
    for (;;) {
        const uint64_t cur = lookup->keys[pos];
        if (!cur)
            return 0;
        if (cur == key) {
            const uint64_t be = lookup->vals[pos];
            *beg = unpack_qry_run_beg(be);
            *end = unpack_qry_run_end(be);
            return 1;
        }
        pos = (pos + 1) & lookup->mask;
    }
}

static inline const qry_ctx_run_t *lookup_qry_ctx_run(const qry_ctx_lookup_t *lookup, uint64_t ctx)
{
    if (!lookup || lookup->n_runs == 0)
        return NULL;

    const size_t buckets = (size_t)1u << lookup->fence_k;
    uint32_t topk = top_k_bits_ctx_local(ctx, lookup->fence_k);
    if (topk >= buckets)
        topk = (uint32_t)buckets - 1;

    size_t lo = lookup->fence[topk];
    size_t hi = lookup->fence[topk + 1];
    while (lo < hi) {
        const size_t mid = lo + ((hi - lo) >> 1);
        if (lookup->runs[mid].ctx < ctx)
            lo = mid + 1;
        else
            hi = mid;
    }

    if (lo < lookup->n_runs && lookup->runs[lo].ctx == ctx)
        return &lookup->runs[lo];
    return NULL;
}

static int compute_streamed_ref_row_one_qry(
    const uint64_t *qry, const qry_ctx_lookup_t *qry_lookup,
    const qry_ctx_hash_lookup_t *qry_hash_lookup, qraw_lookup_mode_t lookup_mode,
    const uint64_t *ref, size_t ref_len,
    uint32_t rn, const ani_opt_t *opt, bool ref_conflict,
    const infile_meta_t *qry_asm,
    const infile_meta_t *ref_asm,
    const ani_ctxmeta_rec_t *qry_meta,
    const ani_ctxmeta_rec_t *ref_meta,
    ani_row_t *row_out)
{
    const uint32_t qry_ctx = lookup_mode == QRAW_LOOKUP_HASH
        ? qry_hash_lookup->n_runs
        : (uint32_t)qry_lookup->n_runs;
    if (qry_ctx == 0 || ref_len == 0)
        return 0;
    if (ref_len > UINT32_MAX)
        errx(EXIT_FAILURE, "%s(): reference sketch for genome %u has too many entries", __func__, rn);
    const uint32_t ref_ctx_total = ref_conflict
        ? count_ctx_runs_sorted_ctxobj64_local(ref, ref_len, opt->ignoreconflict)
        : (uint32_t)ref_len;
    if (ref_ctx_total == 0)
        return 0;
    const uint32_t need_X = ani_density_af_needed_ctx(opt, qry_ctx, ref_ctx_total,
                                                      qry_meta, ref_meta);

    const uint8_t nobjbits = Bitslen.obj;
    const uint64_t objmask = (nobjbits == 64) ? UINT64_MAX : ((1ULL << nobjbits) - 1ULL);
    ani_features_t f;
    memset(&f, 0, sizeof(f));

    uint32_t ref_ctx = 0;
    for (size_t j = 0; j < ref_len; ) {
        const uint64_t ctx = ref[j] >> nobjbits;
        const size_t ref_beg = j;
        do { ++j; } while (j < ref_len && (ref[j] >> nobjbits) == ctx);
        const size_t ref_end = j;

        if (opt->ignoreconflict && ref_end - ref_beg > 1)
            continue;
        ref_ctx++;
        if ((uint32_t)f.XnY_ctx + (ref_ctx_total - ref_ctx) < need_X)
            return 0;

        size_t qry_beg = 0, qry_end = 0;
        if (lookup_mode == QRAW_LOOKUP_HASH) {
            if (!lookup_qry_ctx_hash_run(qry_hash_lookup, ctx, &qry_beg, &qry_end))
                continue;
        } else {
            const qry_ctx_run_t *qry_run = lookup_qry_ctx_run(qry_lookup, ctx);
            if (!qry_run)
                continue;
            qry_beg = qry_run->beg;
            qry_end = qry_run->end;
        }
        f.XnY_ctx++;
        const int min_diff = min_diff_sections_ctxobj64_runs(qry, qry_beg, qry_end,
                                                             ref, ref_beg, ref_end, objmask);
        if (min_diff > 0) {
            f.N_diff_obj++;
            f.N_diff_obj_section += min_diff;
            if (min_diff > 1)
                f.N_mut2_ctx++;
        }
    }

    if (ref_ctx == 0 || (uint32_t)f.XnY_ctx < (uint32_t)opt->ctxcut ||
        (uint32_t)f.XnY_ctx < need_X)
        return 0;

    const double af_q = (double)f.XnY_ctx / (double)qry_ctx;
    const double af_r = (double)f.XnY_ctx / (double)ref_ctx;
	const ani_density_af_t density_af = ani_estimate_density_af(qry_meta, ref_meta,
		(uint32_t)f.XnY_ctx, af_q, af_r);
	if (!ani_report_af_pass(opt, density_af.qry, density_af.ref))
		return 0;

    ani_features_t tmp = f;
    tmp.X_ctx = qry_ctx;
    const double blastn_af_q = lm3ways_af_ANIb_from_features(&tmp);
    tmp = f;
    tmp.X_ctx = ref_ctx;
    const double blastn_af_r = lm3ways_af_ANIb_from_features(&tmp);

    ani_row_t row = make_selected_output_row(rn, &f, opt, qry_ctx, ref_ctx,
                                             af_q, blastn_af_q, af_r, blastn_af_r,
                                             density_af,
                                             qry_asm, ref_asm);
    if (row.selected_ani < opt->anicut)
        return 0;
    *row_out = row;
    return 1;
}

int stream_ref_sketches_one_qraw_lookup(ani_opt_t *ani_opt)
{
    unify_sketch_t *qry = generic_sketch_parse(ani_opt->qrydir, ani_query_parse_flags(ani_opt));
    load_infile_meta_for_best_guard(qry, ani_opt->qrydir, ani_opt);
    if (qry->stat_type != 2)
        errx(EXIT_FAILURE, "%s(): small-query streaming requires a 64-bit query sketch", __func__);
    if (qry->infile_num != 1)
        errx(EXIT_FAILURE, "%s(): expected exactly one query sample, found %d",
             __func__, qry->infile_num);
    const_comask_init(&qry->stats.minco_stat);

    size_t ref_stat_size = 0;
    char *ref_stat_path = test_get_fullpath(ani_opt->refdir, sketch_stat);
    minco_sketch_stat_t *ref_stat = read_from_file(ref_stat_path, &ref_stat_size);
    free(ref_stat_path);
    if (ref_stat->hash_id != qry->hash_id)
        errx(EXIT_FAILURE, "%s(): hash_id mismatch between query and reference sketches", __func__);
    const uint32_t ref_n = (uint32_t)ref_stat->infile_num;
    char (*refname)[PATHLEN] = (char (*)[PATHLEN])(ref_stat + 1);

    size_t ref_idx_size = 0;
    char *ref_idx_path = test_get_fullpath(ani_opt->refdir, idx_sketch_suffix);
    uint64_t *ref_idx = read_from_file(ref_idx_path, &ref_idx_size);
    free(ref_idx_path);
    if (ref_idx_size != ((size_t)ref_n + 1) * sizeof(ref_idx[0]))
        errx(EXIT_FAILURE, "%s(): %s/%s has %zu bytes, expected %zu",
             __func__, ani_opt->refdir, idx_sketch_suffix, ref_idx_size,
             ((size_t)ref_n + 1) * sizeof(ref_idx[0]));

    char (*refanno)[PATHLEN] = read_optional_sketch_annotations(ani_opt->refdir, (int)ref_n);
    infile_meta_t *ref_infile_meta =
        ani_best_guard_enabled(ani_opt) ? read_optional_sketch_infile_meta_stats(ani_opt->refdir, (int)ref_n) : NULL;
    ani_ctxmeta_rec_t *qry_ctxmeta =
        read_optional_ani_ctxmeta_stats(ani_opt->qrydir, qry->infile_num);
    ani_ctxmeta_rec_t *ref_ctxmeta =
        read_optional_ani_ctxmeta_stats(ani_opt->refdir, (int)ref_n);

    const uint64_t *qry_arr = qry->comb_sketch + qry->sketch_index[0];
    const size_t qry_len = (size_t)(qry->sketch_index[1] - qry->sketch_index[0]);
    qraw_lookup_mode_t lookup_mode = QRAW_LOOKUP_HASH;
    const char *lookup_env = getenv("MINCO_QRAW_LOOKUP");
    if (lookup_env && lookup_env[0]) {
        if (strcmp(lookup_env, "sorted") == 0 || strcmp(lookup_env, "fence") == 0 ||
            strcmp(lookup_env, "fencepost") == 0) {
            lookup_mode = QRAW_LOOKUP_SORTED;
        } else if (strcmp(lookup_env, "hash") == 0) {
            lookup_mode = QRAW_LOOKUP_HASH;
        } else {
            errx(EXIT_FAILURE, "%s(): invalid MINCO_QRAW_LOOKUP='%s' (use hash or sorted)",
                 __func__, lookup_env);
        }
    }
    qry_ctx_lookup_t qry_lookup;
    qry_ctx_hash_lookup_t qry_hash_lookup;
    memset(&qry_lookup, 0, sizeof(qry_lookup));
    memset(&qry_hash_lookup, 0, sizeof(qry_hash_lookup));
    if (lookup_mode == QRAW_LOOKUP_HASH)
        build_qry_ctx_hash_lookup(qry_arr, qry_len, &qry_hash_lookup);
    else
        build_qry_ctx_lookup(qry_arr, qry_len, &qry_lookup);

    FILE *outfp = ani_opt->outf[0] == '\0' ? stdout : fopen(ani_opt->outf, "w");
    if (!outfp)
        err(errno, "%s", ani_opt->outf);
    print_ani_detail_header(outfp, ani_opt, true);

    char *ref_comb_path = test_get_fullpath(ani_opt->refdir, combined_sketch_suffix);
    FILE *ref_fp = fopen(ref_comb_path, "rb");
    if (!ref_fp)
        err(errno, "%s", ref_comb_path);
    free(ref_comb_path);
    if (ref_idx[0] != 0) {
        const off_t offset = (off_t)(ref_idx[0] * sizeof(uint64_t));
        if (fseeko(ref_fp, offset, SEEK_SET) != 0)
            err(errno, "%s(): failed to seek reference sketch", __func__);
    }

    const int P = ani_opt->p > 0 ? ani_opt->p : 1;
    const size_t max_block_entries = (size_t)1u << 22; /* 32 MiB of uint64_t sketch data */
    uint64_t *ref_block = NULL;
    size_t ref_block_cap = 0;
    kv_ani_row_t survivors;
    kv_init(survivors);

    for (uint32_t block_start = 0; block_start < ref_n; ) {
        const uint64_t block_begin = ref_idx[block_start];
        uint32_t block_end_gid = block_start;
        uint64_t block_end = block_begin;
        while (block_end_gid < ref_n) {
            const uint64_t next_end = ref_idx[block_end_gid + 1];
            if (block_end_gid > block_start && next_end - block_begin > max_block_entries)
                break;
            block_end_gid++;
            block_end = next_end;
        }
        const size_t block_entries = (size_t)(block_end - block_begin);
        if (block_entries > ref_block_cap) {
            uint64_t *new_block = realloc(ref_block, block_entries * sizeof(*ref_block));
            if (!new_block)
                err(EXIT_FAILURE, "%s(): OOM reference stream block", __func__);
            ref_block = new_block;
            ref_block_cap = block_entries;
        }
        if (block_entries > 0) {
            const size_t got = fread(ref_block, sizeof(*ref_block), block_entries, ref_fp);
            if (got != block_entries)
                errx(EXIT_FAILURE, "%s(): short read from reference sketch: got %zu, expected %zu",
                     __func__, got, block_entries);
        }

        kv_ani_row_t *tls = calloc((size_t)P, sizeof(*tls));
        if (!tls)
            err(EXIT_FAILURE, "%s(): OOM survivor buffers", __func__);
        for (int t = 0; t < P; ++t)
            kv_init(tls[t]);

        const uint32_t block_count = block_end_gid - block_start;
#pragma omp parallel num_threads(P)
        {
            const int tid = omp_get_thread_num();
            ani_row_t row;
#pragma omp for schedule(dynamic, 64)
            for (uint32_t bi = 0; bi < block_count; ++bi) {
                const uint32_t rn = block_start + bi;
                const size_t ref_off = (size_t)(ref_idx[rn] - block_begin);
                const size_t ref_len = (size_t)(ref_idx[rn + 1] - ref_idx[rn]);
                const uint64_t *ref_arr = ref_len ? ref_block + ref_off : NULL;
                if (compute_streamed_ref_row_one_qry(qry_arr, &qry_lookup, &qry_hash_lookup,
                                                     lookup_mode, ref_arr, ref_len,
                                                     rn, ani_opt, ref_stat->conflict,
                                                     infile_meta_at(qry, 0),
                                                     ref_infile_meta ? &ref_infile_meta[rn] : NULL,
                                                     ani_ctxmeta_at(qry_ctxmeta, 0),
                                                     ani_ctxmeta_at(ref_ctxmeta, rn),
                                                     &row))
                    kv_push(ani_row_t, tls[tid], row);
            }
        }

        for (int t = 0; t < P; ++t) {
            if (kv_size(tls[t]))
                kv_append_rows(&survivors, &tls[t]);
            kv_destroy(tls[t]);
        }
        free(tls);
        block_start = block_end_gid;
    }

    if (kv_size(survivors))
        qsort(&kv_A(survivors, 0), kv_size(survivors), sizeof(ani_row_t), cmp_ani_desc);

    size_t out_n = kv_size(survivors);
    if (ani_opt->ntop > 0 && (size_t)ani_opt->ntop < out_n)
        out_n = (size_t)ani_opt->ntop;
    for (size_t i = 0; i < out_n; ++i) {
        const ani_row_t *r = &kv_A(survivors, i);
        print_unified_detail_row(outfp, ani_opt, qry->gname[0], refname[r->rn],
                                 r, annotation_at(refanno, r->rn));
    }

    kv_destroy(survivors);
    free(ref_block);
    fclose(ref_fp);
    if (outfp != stdout)
        fclose(outfp);
    free_qry_ctx_lookup(&qry_lookup);
    free_qry_ctx_hash_lookup(&qry_hash_lookup);
    if (refanno)
        free_read_from_file(refanno, (size_t)ref_n * PATHLEN);
    if (ref_infile_meta)
        free_read_from_file(ref_infile_meta, (size_t)ref_n * sizeof(ref_infile_meta[0]));
    free(qry_ctxmeta);
    free(ref_ctxmeta);
    free_read_from_file(ref_idx, ref_idx_size);
    free_read_from_file(ref_stat, ref_stat_size);
    free_unify_sketch(qry);
    return 0;
}

typedef struct {
    uint64_t ctx;
    uint32_t qid;
    uint32_t beg;
    uint32_t end;
} qraw_multi_qctx_run_t;

typedef struct {
    qraw_multi_qctx_run_t *runs;
    size_t n_runs;
    size_t *fence;
    int fence_k;
    uint32_t *qry_ctx_counts;
    uint32_t qry_n;
} qraw_multi_lookup_t;

static int cmp_qraw_multi_qctx_run(const void *pa, const void *pb)
{
    const qraw_multi_qctx_run_t *a = (const qraw_multi_qctx_run_t *)pa;
    const qraw_multi_qctx_run_t *b = (const qraw_multi_qctx_run_t *)pb;
    if (a->ctx < b->ctx) return -1;
    if (a->ctx > b->ctx) return 1;
    return (a->qid > b->qid) - (a->qid < b->qid);
}

static void build_qraw_multi_lookup(const unify_sketch_t *qry, qraw_multi_lookup_t *lookup)
{
    memset(lookup, 0, sizeof(*lookup));
    const uint32_t Q = (uint32_t)qry->infile_num;
    lookup->qry_n = Q;
    lookup->qry_ctx_counts = calloc(Q, sizeof(*lookup->qry_ctx_counts));
    if (!lookup->qry_ctx_counts)
        err(EXIT_FAILURE, "%s(): OOM query context counts", __func__);

    const uint64_t total_entries = qry->sketch_index[Q];
    if (total_entries == 0)
        return;
    if (total_entries > UINT32_MAX)
        errx(EXIT_FAILURE, "%s(): query sketch has too many entries for packed multi-query lookup", __func__);

    lookup->runs = malloc((size_t)total_entries * sizeof(*lookup->runs));
    if (!lookup->runs)
        err(EXIT_FAILURE, "%s(): OOM multi-query context runs", __func__);

    const uint8_t nobjbits = Bitslen.obj;
    for (uint32_t q = 0; q < Q; ++q) {
        size_t i = (size_t)qry->sketch_index[q];
        const size_t end_q = (size_t)qry->sketch_index[q + 1];
        while (i < end_q) {
            const uint64_t ctx = qry->comb_sketch[i] >> nobjbits;
            const size_t beg = i;
            do { ++i; } while (i < end_q && (qry->comb_sketch[i] >> nobjbits) == ctx);
            lookup->runs[lookup->n_runs++] = (qraw_multi_qctx_run_t){
                .ctx = ctx, .qid = q, .beg = (uint32_t)beg, .end = (uint32_t)i
            };
            lookup->qry_ctx_counts[q]++;
        }
    }

    if (lookup->n_runs > 1)
        qsort(lookup->runs, lookup->n_runs, sizeof(*lookup->runs), cmp_qraw_multi_qctx_run);

    lookup->fence_k = choose_qry_fence_k(lookup->n_runs);
    const size_t buckets = (size_t)1u << lookup->fence_k;
    lookup->fence = malloc((buckets + 1) * sizeof(*lookup->fence));
    if (!lookup->fence)
        err(EXIT_FAILURE, "%s(): OOM multi-query context fence", __func__);

    for (size_t t = 0; t <= buckets; ++t)
        lookup->fence[t] = lookup->n_runs;

    size_t next = 0;
    for (size_t r = 0; r < lookup->n_runs; ++r) {
        const uint32_t topk = top_k_bits_ctx_local(lookup->runs[r].ctx, lookup->fence_k);
        while (next <= topk && next <= buckets)
            lookup->fence[next++] = r;
        if (next > buckets)
            break;
    }
    while (next <= buckets)
        lookup->fence[next++] = lookup->n_runs;
}

static void free_qraw_multi_lookup(qraw_multi_lookup_t *lookup)
{
    if (!lookup)
        return;
    free(lookup->runs);
    free(lookup->fence);
    free(lookup->qry_ctx_counts);
    memset(lookup, 0, sizeof(*lookup));
}

static inline void qraw_multi_lookup_ctx_range(const qraw_multi_lookup_t *lookup, uint64_t ctx,
                                               size_t *lo_out, size_t *hi_out)
{
    *lo_out = *hi_out = 0;
    if (!lookup || lookup->n_runs == 0)
        return;

    const size_t buckets = (size_t)1u << lookup->fence_k;
    uint32_t topk = top_k_bits_ctx_local(ctx, lookup->fence_k);
    if (topk >= buckets)
        topk = (uint32_t)buckets - 1;

    size_t lo = lookup->fence[topk];
    size_t hi = lookup->fence[topk + 1];
    while (lo < hi) {
        const size_t mid = lo + ((hi - lo) >> 1);
        if (lookup->runs[mid].ctx < ctx)
            lo = mid + 1;
        else
            hi = mid;
    }
    if (lo >= lookup->n_runs || lookup->runs[lo].ctx != ctx)
        return;

    size_t end = lo + 1;
    while (end < lookup->n_runs && lookup->runs[end].ctx == ctx)
        end++;
    *lo_out = lo;
    *hi_out = end;
}

int stream_ref_sketches_multi_qraw_sortedindex(ani_opt_t *ani_opt)
{
    unify_sketch_t *qry = generic_sketch_parse(ani_opt->qrydir, ani_query_parse_flags(ani_opt));
    load_infile_meta_for_best_guard(qry, ani_opt->qrydir, ani_opt);
    if (qry->stat_type != 2)
        errx(EXIT_FAILURE, "%s(): small-query streaming requires 64-bit query sketches", __func__);
    if (qry->infile_num <= 1)
        errx(EXIT_FAILURE, "%s(): expected multiple query samples, found %d",
             __func__, qry->infile_num);
    const_comask_init(&qry->stats.minco_stat);
    const uint32_t Q = (uint32_t)qry->infile_num;

    size_t ref_stat_size = 0;
    char *ref_stat_path = test_get_fullpath(ani_opt->refdir, sketch_stat);
    minco_sketch_stat_t *ref_stat = read_from_file(ref_stat_path, &ref_stat_size);
    free(ref_stat_path);
    if (ref_stat->hash_id != qry->hash_id)
        errx(EXIT_FAILURE, "%s(): hash_id mismatch between query and reference sketches", __func__);
    const uint32_t ref_n = (uint32_t)ref_stat->infile_num;
    char (*refname)[PATHLEN] = (char (*)[PATHLEN])(ref_stat + 1);

    size_t ref_idx_size = 0;
    char *ref_idx_path = test_get_fullpath(ani_opt->refdir, idx_sketch_suffix);
    uint64_t *ref_idx = read_from_file(ref_idx_path, &ref_idx_size);
    free(ref_idx_path);
    if (ref_idx_size != ((size_t)ref_n + 1) * sizeof(ref_idx[0]))
        errx(EXIT_FAILURE, "%s(): %s/%s has %zu bytes, expected %zu",
             __func__, ani_opt->refdir, idx_sketch_suffix, ref_idx_size,
             ((size_t)ref_n + 1) * sizeof(ref_idx[0]));

    char (*refanno)[PATHLEN] = read_optional_sketch_annotations(ani_opt->refdir, (int)ref_n);
    infile_meta_t *ref_infile_meta =
        ani_best_guard_enabled(ani_opt) ? read_optional_sketch_infile_meta_stats(ani_opt->refdir, (int)ref_n) : NULL;
    ani_ctxmeta_rec_t *qry_ctxmeta =
        read_optional_ani_ctxmeta_stats(ani_opt->qrydir, qry->infile_num);
    ani_ctxmeta_rec_t *ref_ctxmeta =
        read_optional_ani_ctxmeta_stats(ani_opt->refdir, (int)ref_n);
    qraw_multi_lookup_t lookup;
    build_qraw_multi_lookup(qry, &lookup);

    FILE *outfp = ani_opt->outf[0] == '\0' ? stdout : fopen(ani_opt->outf, "w");
    if (!outfp)
        err(errno, "%s", ani_opt->outf);
    print_ani_detail_header(outfp, ani_opt, true);

    char *ref_comb_path = test_get_fullpath(ani_opt->refdir, combined_sketch_suffix);
    FILE *ref_fp = fopen(ref_comb_path, "rb");
    if (!ref_fp)
        err(errno, "%s", ref_comb_path);
    free(ref_comb_path);
    if (ref_idx[0] != 0) {
        const off_t offset = (off_t)(ref_idx[0] * sizeof(uint64_t));
        if (fseeko(ref_fp, offset, SEEK_SET) != 0)
            err(errno, "%s(): failed to seek reference sketch", __func__);
    }

    const int P = ani_opt->p > 0 ? ani_opt->p : 1;
    const size_t max_block_entries = (size_t)1u << 22;
    uint64_t *ref_block = NULL;
    size_t ref_block_cap = 0;
    kv_ani_row_t *survivors = calloc(Q, sizeof(*survivors));
    if (!survivors)
        err(EXIT_FAILURE, "%s(): OOM survivor vectors", __func__);
    for (uint32_t q = 0; q < Q; ++q)
        kv_init(survivors[q]);

    const uint8_t nobjbits = Bitslen.obj;
    const uint64_t objmask = (nobjbits == 64) ? UINT64_MAX : ((1ULL << nobjbits) - 1ULL);

    for (uint32_t block_start = 0; block_start < ref_n; ) {
        const uint64_t block_begin = ref_idx[block_start];
        uint32_t block_end_gid = block_start;
        uint64_t block_end = block_begin;
        while (block_end_gid < ref_n) {
            const uint64_t next_end = ref_idx[block_end_gid + 1];
            if (block_end_gid > block_start && next_end - block_begin > max_block_entries)
                break;
            block_end_gid++;
            block_end = next_end;
        }
        const size_t block_entries = (size_t)(block_end - block_begin);
        if (block_entries > ref_block_cap) {
            uint64_t *new_block = realloc(ref_block, block_entries * sizeof(*ref_block));
            if (!new_block)
                err(EXIT_FAILURE, "%s(): OOM reference stream block", __func__);
            ref_block = new_block;
            ref_block_cap = block_entries;
        }
        if (block_entries > 0) {
            const size_t got = fread(ref_block, sizeof(*ref_block), block_entries, ref_fp);
            if (got != block_entries)
                errx(EXIT_FAILURE, "%s(): short read from reference sketch: got %zu, expected %zu",
                     __func__, got, block_entries);
        }

        kv_ani_row_t *tls = calloc((size_t)P * Q, sizeof(*tls));
        if (!tls)
            err(EXIT_FAILURE, "%s(): OOM thread-local survivor vectors", __func__);
        for (size_t i = 0; i < (size_t)P * Q; ++i)
            kv_init(tls[i]);

        const uint32_t block_count = block_end_gid - block_start;
#pragma omp parallel num_threads(P)
        {
            const int tid = omp_get_thread_num();
            ani_features_t *features = calloc(Q, sizeof(*features));
            if (!features)
                err(EXIT_FAILURE, "%s(): OOM per-reference features", __func__);

#pragma omp for schedule(dynamic, 64)
            for (uint32_t bi = 0; bi < block_count; ++bi) {
                const uint32_t rn = block_start + bi;
                const size_t ref_off = (size_t)(ref_idx[rn] - block_begin);
                const size_t ref_len = (size_t)(ref_idx[rn + 1] - ref_idx[rn]);
                const uint64_t *ref = ref_len ? ref_block + ref_off : NULL;
                memset(features, 0, (size_t)Q * sizeof(*features));
                uint32_t ref_ctx = 0;

                for (size_t j = 0; j < ref_len; ) {
                    const uint64_t ctx = ref[j] >> nobjbits;
                    const size_t ref_beg = j;
                    do { ++j; } while (j < ref_len && (ref[j] >> nobjbits) == ctx);
                    const size_t ref_end = j;

                    if (ani_opt->ignoreconflict && ref_end - ref_beg > 1)
                        continue;
                    ref_ctx++;

                    size_t qlo = 0, qhi = 0;
                    qraw_multi_lookup_ctx_range(&lookup, ctx, &qlo, &qhi);
                    for (size_t qi = qlo; qi < qhi; ++qi) {
                        const qraw_multi_qctx_run_t *qr = &lookup.runs[qi];
                        ani_features_t *f = &features[qr->qid];
                        f->XnY_ctx++;
                        const int min_diff = min_diff_sections_ctxobj64_runs(
                            qry->comb_sketch, qr->beg, qr->end,
                            ref, ref_beg, ref_end, objmask);
                        if (min_diff > 0) {
                            f->N_diff_obj++;
                            f->N_diff_obj_section += min_diff;
                            if (min_diff > 1)
                                f->N_mut2_ctx++;
                        }
                    }
                }

                if (ref_ctx == 0)
                    continue;

                for (uint32_t q = 0; q < Q; ++q) {
                    ani_features_t f = features[q];
                    const uint32_t qry_ctx = lookup.qry_ctx_counts[q];
                    if (qry_ctx == 0 || (uint32_t)f.XnY_ctx < (uint32_t)ani_opt->ctxcut)
                        continue;
                    const ani_ctxmeta_rec_t *qmeta = ani_ctxmeta_at(qry_ctxmeta, q);
                    const ani_ctxmeta_rec_t *rmeta = ani_ctxmeta_at(ref_ctxmeta, rn);
                    const uint32_t need_X = ani_density_af_needed_ctx(ani_opt, qry_ctx, ref_ctx,
                                                                      qmeta, rmeta);
                    if ((uint32_t)f.XnY_ctx < need_X)
                        continue;

                    const double af_q = (double)f.XnY_ctx / (double)qry_ctx;
                    const double af_r = (double)f.XnY_ctx / (double)ref_ctx;
                    const ani_density_af_t density_af = ani_estimate_density_af(
                        qmeta, rmeta, (uint32_t)f.XnY_ctx, af_q, af_r);
					if (!ani_report_af_pass(ani_opt, density_af.qry, density_af.ref))
						continue;

                    ani_features_t tmp = f;
                    tmp.X_ctx = qry_ctx;
                    const double blastn_af_q = lm3ways_af_ANIb_from_features(&tmp);
                    tmp = f;
                    tmp.X_ctx = ref_ctx;
                    const double blastn_af_r = lm3ways_af_ANIb_from_features(&tmp);

                    ani_row_t row = make_selected_output_row(rn, &f, ani_opt, qry_ctx, ref_ctx,
                                                             af_q, blastn_af_q, af_r, blastn_af_r,
                                                             density_af,
                                                             infile_meta_at(qry, q),
                                                             ref_infile_meta ? &ref_infile_meta[rn] : NULL);
                    if (row.selected_ani < ani_opt->anicut)
                        continue;
                    kv_push(ani_row_t, tls[(size_t)tid * Q + q], row);
                }
            }

            free(features);
        }

        for (int t = 0; t < P; ++t) {
            for (uint32_t q = 0; q < Q; ++q) {
                kv_ani_row_t *src = &tls[(size_t)t * Q + q];
                if (kv_size(*src))
                    kv_append_rows(&survivors[q], src);
                kv_destroy(*src);
            }
        }
        free(tls);
        block_start = block_end_gid;
    }

    for (uint32_t q = 0; q < Q; ++q) {
        if (kv_size(survivors[q]))
            qsort(&kv_A(survivors[q], 0), kv_size(survivors[q]), sizeof(ani_row_t), cmp_ani_desc);
        size_t out_n = kv_size(survivors[q]);
        if (ani_opt->ntop > 0 && (size_t)ani_opt->ntop < out_n)
            out_n = (size_t)ani_opt->ntop;
        for (size_t i = 0; i < out_n; ++i) {
            const ani_row_t *r = &kv_A(survivors[q], i);
            print_unified_detail_row(outfp, ani_opt, qry->gname[q], refname[r->rn],
                                     r, annotation_at(refanno, r->rn));
        }
        kv_destroy(survivors[q]);
    }

    free(survivors);
    free(ref_block);
    fclose(ref_fp);
    if (outfp != stdout)
        fclose(outfp);
    free_qraw_multi_lookup(&lookup);
    if (refanno)
        free_read_from_file(refanno, (size_t)ref_n * PATHLEN);
    if (ref_infile_meta)
        free_read_from_file(ref_infile_meta, (size_t)ref_n * sizeof(ref_infile_meta[0]));
    free(qry_ctxmeta);
    free(ref_ctxmeta);
    free_read_from_file(ref_idx, ref_idx_size);
    free_read_from_file(ref_stat, ref_stat_size);
    free_unify_sketch(qry);
    return 0;
}

/* --- main --- */
/** 
void comb_sortedsketch64Xcomb_sortedsketch64_filter_and_sort_survivors(ani_opt_t *ani_opt)
{
    unify_sketch_t *qry = generic_sketch_parse(ani_opt->qrydir, ani_query_parse_flags(ani_opt));
    unify_sketch_t *ref = generic_sketch_parse(ani_opt->refdir, ani_ref_parse_flags(ani_opt));
    load_infile_meta_for_best_guard(qry, ani_opt->qrydir, ani_opt);
    load_infile_meta_for_best_guard(ref, ani_opt->refdir, ani_opt);
    ani_ctxmeta_rec_t *qry_ctxmeta =
        read_optional_ani_ctxmeta_stats(ani_opt->qrydir, qry->infile_num);
    ani_ctxmeta_rec_t *ref_ctxmeta =
        read_optional_ani_ctxmeta_stats(ani_opt->refdir, ref->infile_num);
    if (ref->conflict)
        errx(EXIT_FAILURE, "%s(): ref '%s' contains conflicting objects!", __func__, ani_opt->refdir);

    const uint32_t Q = qry->infile_num;
    const uint32_t R = ref->infile_num;

    FILE *outfp = (ani_opt->outf[0] == '\0') ? stdout : fopen(ani_opt->outf, "w");
    if (!outfp) err(errno, "%s", ani_opt->outf);
    print_ani_detail_header(outfp, ani_opt, false);

    const int P = (ani_opt->p > 0) ? ani_opt->p : 1;

    // ---- P==1: fast serial path (logic unchanged) ---- //
    if (P == 1) {
        for (uint32_t qn = 0; qn < Q; ++qn) {
            kv_ani_row_t surv; kv_init(surv);
            kstring_t ks_out = (kstring_t){0,0,0};

            ani_row_t tmp;
            for (uint32_t rn = 0; rn < R; ++rn) {
                if (compute_row_if_survivor(qry, ref, qn, rn, ani_opt, &tmp))
                    kv_push(ani_row_t, surv, tmp);
            }

            if (kv_size(surv))
                qsort(&kv_A(surv,0), kv_size(surv), sizeof(ani_row_t), cmp_ani_desc);

            format_rows_to_kstr(qry, ref, qn, &surv, ani_opt, &ks_out);
            if (ks_out.l) fwrite(ks_out.s, 1, ks_out.l, outfp);

            kv_destroy(surv);
            free(ks_out.s);
        }
        if (outfp != stdout) fclose(outfp);
        return;
    }

#ifdef PRINT_UNORDERED
    // ---- Parallel, unordered print (as-ready) ---- //
    #pragma omp parallel num_threads(P)
    {
        #pragma omp single nowait
        {
            const int many_queries = (Q >= (uint32_t)P);

            for (uint32_t qn = 0; qn < Q; ++qn) {
                #pragma omp task firstprivate(qn, many_queries) shared(qry,ref,ani_opt,outfp)
                {
                    kv_ani_row_t surv; kv_init(surv);

                    if (many_queries) {
                        ani_row_t tmp;
                        for (uint32_t rn = 0; rn < R; ++rn) {
                            if (compute_row_if_survivor(qry, ref, qn, rn, ani_opt, &tmp))
                                kv_push(ani_row_t, surv, tmp);
                        }
                    } else {
                        #pragma omp taskgroup
                        {
                            const uint32_t CHUNK = 256;
                            for (uint32_t start = 0; start < R; start += CHUNK) {
                                const uint32_t end = (start + CHUNK < R) ? start + CHUNK : R;
                                #pragma omp task firstprivate(start,end,qn) shared(qry,ref,ani_opt,surv)
                                {
                                    kv_ani_row_t local; kv_init(local);
                                    ani_row_t tmp;
                                    for (uint32_t rn = start; rn < end; ++rn) {
                                        if (compute_row_if_survivor(qry, ref, qn, rn, ani_opt, &tmp))
                                            kv_push(ani_row_t, local, tmp);
                                    }
                                    if (kv_size(local)) {
                                        #pragma omp critical (append_survivors)
                                        kv_append_rows(&surv, &local);
                                    }
                                    kv_destroy(local);
                                }
                            }
                        }
                    }

                    if (kv_size(surv))
                        qsort(&kv_A(surv,0), kv_size(surv), sizeof(ani_row_t), cmp_ani_desc);

                    kstring_t ks_out = (kstring_t){0,0,0};
                    format_rows_to_kstr(qry, ref, qn, &surv, ani_opt, &ks_out);

                    if (ks_out.l) {
                        #pragma omp critical (print)
                        fwrite(ks_out.s, 1, ks_out.l, outfp);
                    }

                    free(ks_out.s);
                    kv_destroy(surv);
                } // compute+print task //
            }
        } // single //
    } // parallel //
#else
    // ---- Parallel, ordered print: compute in parallel; only print is serialized ---- //

    // tokens: comp[qn] signals compute done; order[i] forms a print chain i->i+1 //
    char *comp  = (char*)calloc(Q,   1);
    char *order = (char*)calloc(Q+1, 1);
    if (!comp || !order) err(EXIT_FAILURE, "calloc tokens");

    // per-qn output buffers produced by compute, consumed by print //
    kstring_t *ks_arr = (kstring_t*)calloc(Q, sizeof(kstring_t));
    if (!ks_arr) err(EXIT_FAILURE, "calloc ks_arr");

    #pragma omp parallel num_threads(P)
    {
        #pragma omp single nowait
        {
            const int many_queries = (Q >= (uint32_t)P);

            for (uint32_t qn = 0; qn < Q; ++qn) {

                // 1) compute task: fully parallel across qn //
                #pragma omp task firstprivate(qn, many_queries) shared(qry,ref,ani_opt,ks_arr,comp) depend(out: comp[qn])
                {
                    kv_ani_row_t surv; kv_init(surv);

                    if (many_queries) {
                        ani_row_t tmp;
                        for (uint32_t rn = 0; rn < R; ++rn) {
                            if (compute_row_if_survivor(qry, ref, qn, rn, ani_opt, &tmp))
                                kv_push(ani_row_t, surv, tmp);
                        }
                    } else {
                        #pragma omp taskgroup
                        {
                            const uint32_t CHUNK = 256;
                            for (uint32_t start = 0; start < R; start += CHUNK) {
                                const uint32_t end = (start + CHUNK < R) ? start + CHUNK : R;
                                #pragma omp task firstprivate(start,end,qn) shared(qry,ref,ani_opt,surv)
                                {
                                    kv_ani_row_t local; kv_init(local);
                                    ani_row_t tmp;
                                    for (uint32_t rn = start; rn < end; ++rn) {
                                        if (compute_row_if_survivor(qry, ref, qn, rn, ani_opt, &tmp))
                                            kv_push(ani_row_t, local, tmp);
                                    }
                                    if (kv_size(local)) {
                                        #pragma omp critical (append_survivors)
                                        kv_append_rows(&surv, &local);
                                    }
                                    kv_destroy(local);
                                }
                            }
                        }
                    }

                    if (kv_size(surv))
                        qsort(&kv_A(surv,0), kv_size(surv), sizeof(ani_row_t), cmp_ani_desc);

                    kstring_t ks_out = (kstring_t){0,0,0};
                    format_rows_to_kstr(qry, ref, qn, &surv, ani_opt, &ks_out);

                    ks_arr[qn] = ks_out; // move: keep buffer for print task //
                    kv_destroy(surv);
                }

                // 2) print task: serialized by order[] but can run as soon as qn computed //
                #pragma omp task firstprivate(qn) shared(outfp, ks_arr, comp, order) depend(in: comp[qn]) depend(in: order[qn]) depend(out: order[qn+1])
                {
                    kstring_t *ksp = &ks_arr[qn];
                    if (ksp->l)
                        fwrite(ksp->s, 1, ksp->l, outfp);
                    free(ksp->s);
                    ksp->s = NULL; ksp->l = ksp->m = 0;
                }
            } // for qn //
        } // single //
    } // parallel //

    free(ks_arr);
    free(order);
    free(comp);
#endif // PRINT_UNORDERED //

    if (outfp != stdout) fclose(outfp);
}
*/

static void comb_sortedsketch96Xcomb_sortedsketch96_detail(ani_opt_t *ani_opt,
														   unify_sketch_t *qry_result,
														   unify_sketch_t *ref_result)
{
	if (ani_opt->fmt != 0)
		errx(EXIT_FAILURE, "ctxobj96 ANI currently supports detail output (-m0) only");
	pairwise_prepare_minco_model(ref_result);
	load_infile_meta_for_best_guard(qry_result, ani_opt->qrydir, ani_opt);
	load_infile_meta_for_best_guard(ref_result, ani_opt->refdir, ani_opt);
	ani_ctxmeta_rec_t *qry_ctxmeta =
		read_optional_ani_ctxmeta_stats(ani_opt->qrydir, qry_result->infile_num);
	const bool same_sketch = strcmp(ani_opt->qrydir, ani_opt->refdir) == 0;
	ani_ctxmeta_rec_t *ref_ctxmeta = same_sketch
		? qry_ctxmeta
		: read_optional_ani_ctxmeta_stats(ani_opt->refdir, ref_result->infile_num);

	FILE *outfp = ani_opt->outf[0] == '\0' ? stdout : fopen(ani_opt->outf, "w");
	if (outfp == NULL)
		err(errno, "%s", ani_opt->outf);
	print_ani_detail_header(outfp, ani_opt, false);
	for (uint32_t rn = 0; rn < (uint32_t)ref_result->infile_num; rn++)
	{
		ctxobj96_t *arr_ref = ref_result->comb_sketch96 + ref_result->sketch_index[rn];
		size_t len_ref = ref_result->sketch_index[rn + 1] - ref_result->sketch_index[rn];
		for (uint32_t qn = 0; qn < (uint32_t)qry_result->infile_num; qn++)
		{
			ctxobj96_t *arr_qry = qry_result->comb_sketch96 + qry_result->sketch_index[qn];
			size_t len_qry = qry_result->sketch_index[qn + 1] - qry_result->sketch_index[qn];
			if (len_ref == 0 || len_qry == 0)
				continue;
			ani_features_t ani_features;
			get_ani_features_from_two_sorted_ctxobj96(arr_ref, len_ref, arr_qry, len_qry,
													  &ani_features);
			double af_qry = (double)ani_features.XnY_ctx / (double)len_qry;
			double af_ref = (double)ani_features.XnY_ctx / (double)len_ref;
			const ani_density_af_t density_af = ani_estimate_density_af(
				ani_ctxmeta_at(qry_ctxmeta, qn), ani_ctxmeta_at(ref_ctxmeta, rn),
				(uint32_t)ani_features.XnY_ctx, af_qry, af_ref);
			if (!ani_report_af_pass(ani_opt, density_af.qry, density_af.ref))
				continue;
			ani_features.X_ctx = len_qry;
			double blastn_af_qry = lm3ways_af_ANIb_from_features(&ani_features);
			ani_features.X_ctx = len_ref;
			double blastn_af_ref = lm3ways_af_ANIb_from_features(&ani_features);
			ani_row_t outrow = make_selected_output_row(rn, &ani_features, ani_opt,
														(uint32_t)len_qry,
														(uint32_t)len_ref,
														af_qry, blastn_af_qry,
														af_ref, blastn_af_ref,
														density_af,
														infile_meta_at(qry_result, qn),
														infile_meta_at(ref_result, rn));
			if (outrow.selected_ani < ani_opt->anicut)
				continue;
			print_unified_detail_row(outfp, ani_opt, qry_result->gname[qn], ref_result->gname[rn],
									 &outrow, unify_annotation_at(ref_result, rn));
		}
	}
	if (outfp != stdout)
		fclose(outfp);
	if (ref_ctxmeta != qry_ctxmeta)
		free(ref_ctxmeta);
	free(qry_ctxmeta);
}

void comb_sortedsketch64Xcomb_sortedsketch64(ani_opt_t *ani_opt)
{
	unify_sketch_t *qry_result = generic_sketch_parse(ani_opt->qrydir, ani_query_parse_flags(ani_opt));
	unify_sketch_t *ref_result = generic_sketch_parse(ani_opt->refdir, ani_ref_parse_flags(ani_opt));
	pairwise_check_compatible(ref_result, qry_result);
	if (qry_result->payload_layout == MINCO_PAYLOAD_CTXOBJ96 ||
		ref_result->payload_layout == MINCO_PAYLOAD_CTXOBJ96)
	{
		if (qry_result->payload_layout != MINCO_PAYLOAD_CTXOBJ96 ||
			ref_result->payload_layout != MINCO_PAYLOAD_CTXOBJ96)
			errx(EXIT_FAILURE, "%s(): cannot mix ctxobj64 and ctxobj96 payload sketches", __func__);
		comb_sortedsketch96Xcomb_sortedsketch96_detail(ani_opt, qry_result, ref_result);
		free_unify_sketch(qry_result);
		free_unify_sketch(ref_result);
		return;
	}
	pairwise_prepare_minco_model(ref_result);
	load_infile_meta_for_best_guard(qry_result, ani_opt->qrydir, ani_opt);
	load_infile_meta_for_best_guard(ref_result, ani_opt->refdir, ani_opt);
	const bool same_sketch = strcmp(ani_opt->qrydir, ani_opt->refdir) == 0;
	ani_ctxmeta_rec_t *qry_ctxmeta =
		read_optional_ani_ctxmeta_stats(ani_opt->qrydir, qry_result->infile_num);
	ani_ctxmeta_rec_t *ref_ctxmeta = same_sketch
		? qry_ctxmeta
		: read_optional_ani_ctxmeta_stats(ani_opt->refdir, ref_result->infile_num);

	FILE *outfp = ani_opt->outf[0] == '\0' ? stdout : fopen(ani_opt->outf, "w");
	if (outfp == NULL)
		err(errno, "%s", ani_opt->outf);
	if (ani_opt->fmt == 1)
	{
		for (uint32_t rn = 0; rn < ref_result->infile_num; rn++)
			fprintf(outfp, "\t%s", ref_result->gname[rn]);
		fprintf(outfp, "\n");
		for (uint32_t qn = 0; qn < qry_result->infile_num; qn++)
		{
			fprintf(outfp, "%s", qry_result->gname[qn]);
			for (uint32_t rn = 0; rn < ref_result->infile_num; rn++)
			{
				double metric;
				if (same_sketch && qn == rn)
					metric = ani_matrix_diagonal_value(ani_opt);
				else if (same_sketch)
				{
					const uint32_t q = qn > rn ? qn : rn;
					const uint32_t r = qn > rn ? rn : qn;
					metric = ani_matrix_pair_value(qry_result, q, ref_result, r, ani_opt,
												   qry_ctxmeta, ref_ctxmeta);
				}
				else
					metric = ani_matrix_pair_value(qry_result, qn, ref_result, rn, ani_opt,
												   qry_ctxmeta, ref_ctxmeta);
				fprintf(outfp, "\t%lf", metric);
			}
			fprintf(outfp, "\n");
		}
		if (outfp != stdout)
			fclose(outfp);
		if (ref_ctxmeta != qry_ctxmeta)
			free(ref_ctxmeta);
		free(qry_ctxmeta);
		free_unify_sketch(qry_result);
		free_unify_sketch(ref_result);
		return;
	}
	if (ani_opt->fmt == 2)
	{
		if (!same_sketch)
			errx(EXIT_FAILURE, "%s(): triangle output requires one sketch or identical -r/-q sketches", __func__);
		for (uint32_t qn = 0; qn < qry_result->infile_num; qn++)
		{
			fprintf(outfp, "%s", qry_result->gname[qn]);
			for (uint32_t rn = 0; rn < qn; rn++)
				fprintf(outfp, "\t%lf", ani_matrix_pair_value(qry_result, qn, ref_result, rn, ani_opt,
															  qry_ctxmeta, ref_ctxmeta));
			if (ani_opt->d)
				fprintf(outfp, "\t%lf", ani_matrix_diagonal_value(ani_opt));
			fprintf(outfp, "\n");
		}
		if (outfp != stdout)
			fclose(outfp);
		if (ref_ctxmeta != qry_ctxmeta)
			free(ref_ctxmeta);
		free(qry_ctxmeta);
		free_unify_sketch(qry_result);
		free_unify_sketch(ref_result);
		return;
	}

	print_ani_detail_header(outfp, ani_opt, false);
	// #pragma omp parallel for num_threads(32) schedule(guided)
	for (uint32_t rn = 0; rn < ref_result->infile_num; rn++)
	{

		uint64_t *arr_ref = ref_result->comb_sketch + ref_result->sketch_index[rn];
		size_t len_ref = ref_result->sketch_index[rn + 1] - ref_result->sketch_index[rn];
		// printf(outfp,"%s", ref_result->gname[rn]);

		for (uint32_t qn = 0; qn < qry_result->infile_num; qn++)
		{
			ani_features_t ani_features;
			uint64_t *arr_qry = qry_result->comb_sketch + qry_result->sketch_index[qn];
			size_t len_qry = qry_result->sketch_index[qn + 1] - qry_result->sketch_index[qn];
			get_ani_features_from_two_sorted_ctxobj64(arr_ref, len_ref, arr_qry, len_qry, &ani_features);
			double af_qry = (double)ani_features.XnY_ctx / len_qry;
			double af_ref = (double)ani_features.XnY_ctx / len_ref;
			const ani_density_af_t density_af = ani_estimate_density_af(
				ani_ctxmeta_at(qry_ctxmeta, qn), ani_ctxmeta_at(ref_ctxmeta, rn),
				(uint32_t)ani_features.XnY_ctx, af_qry, af_ref);

			if (!ani_report_af_pass(ani_opt, density_af.qry, density_af.ref))
				continue;
			ani_features.X_ctx = len_qry;
			double blastn_af_qry = lm3ways_af_ANIb_from_features(&ani_features);
			ani_features.X_ctx = len_ref;
			double blastn_af_ref = lm3ways_af_ANIb_from_features(&ani_features);
			ani_row_t outrow = make_selected_output_row(rn, &ani_features, ani_opt,
														(uint32_t)len_qry,
														(uint32_t)len_ref,
														af_qry, blastn_af_qry,
														af_ref, blastn_af_ref,
														density_af,
														infile_meta_at(qry_result, qn),
														infile_meta_at(ref_result, rn));
			if (outrow.selected_ani < ani_opt->anicut)
				continue;
			print_unified_detail_row(outfp, ani_opt, qry_result->gname[qn], ref_result->gname[rn],
									 &outrow, unify_annotation_at(ref_result, rn));
		}
	}
	if (outfp != stdout)
		fclose(outfp);
	if (ref_ctxmeta != qry_ctxmeta)
		free(ref_ctxmeta);
	free(qry_ctxmeta);
	free_unify_sketch(qry_result);
	free_unify_sketch(ref_result);
}

void comb_sortedsketch64_self_matrix(ani_opt_t *ani_opt)
{
	unify_sketch_t *sketch = generic_sketch_parse(ani_opt->qrydir, SKETCH_PARSE_NONE);
	load_infile_meta_for_best_guard(sketch, ani_opt->qrydir, ani_opt);
	ani_ctxmeta_rec_t *ctxmeta =
		read_optional_ani_ctxmeta_stats(ani_opt->qrydir, sketch->infile_num);

	FILE *outfp = ani_opt->outf[0] == '\0' ? stdout : fopen(ani_opt->outf, "w");
	if (outfp == NULL)
		err(errno, "%s", ani_opt->outf);

	if (ani_opt->fmt == 1)
	{
		for (uint32_t rn = 0; rn < sketch->infile_num; rn++)
			fprintf(outfp, "\t%s", sketch->gname[rn]);
		fprintf(outfp, "\n");
		for (uint32_t qn = 0; qn < sketch->infile_num; qn++)
		{
			fprintf(outfp, "%s", sketch->gname[qn]);
			for (uint32_t rn = 0; rn < sketch->infile_num; rn++)
			{
				double metric;
				if (qn == rn)
					metric = ani_matrix_diagonal_value(ani_opt);
				else
				{
					const uint32_t q = qn > rn ? qn : rn;
					const uint32_t r = qn > rn ? rn : qn;
					metric = ani_matrix_pair_value(sketch, q, sketch, r, ani_opt,
												   ctxmeta, ctxmeta);
				}
				fprintf(outfp, "\t%lf", metric);
			}
			fprintf(outfp, "\n");
		}
	}
	else if (ani_opt->fmt == 2)
	{
		for (uint32_t qn = 0; qn < sketch->infile_num; qn++)
		{
			fprintf(outfp, "%s", sketch->gname[qn]);
			for (uint32_t rn = 0; rn < qn; rn++)
				fprintf(outfp, "\t%lf", ani_matrix_pair_value(sketch, qn, sketch, rn, ani_opt,
															  ctxmeta, ctxmeta));
			if (ani_opt->d)
				fprintf(outfp, "\t%lf", ani_matrix_diagonal_value(ani_opt));
			fprintf(outfp, "\n");
		}
	}
	else
		errx(EXIT_FAILURE, "%s(): one-sketch ANI requires -m1 full matrix or -m2 triangle", __func__);

	if (outfp != stdout)
		fclose(outfp);
	free(ctxmeta);
	free_unify_sketch(sketch);
}

void check_comb_sortedsketch64(unify_sketch_t *result)
{
	for (uint32_t rn = 0; rn < result->infile_num; rn++)
	{
		uint64_t *arr = result->comb_sketch + result->sketch_index[rn];
		size_t len = result->sketch_index[rn + 1] - result->sketch_index[rn];
		for (uint32_t i = 1; i < len; i++)
		{
			if (arr[i] < arr[i - 1])
				err(EXIT_FAILURE, "%s(): %dth genome %dth kmer < %dth kmer (%lx<%lx)", __func__, rn, i, i - 1, arr[i], arr[i - 1]);
		}
	}
}

//

/**
 * Finds the first occurrence index in large sorted array `b` for each element in a small sorted array `a`.
 *
 * @param a       Sorted array of uint64_t elements (ascending order)
 * @param a_size  Number of elements in array `a`
 * @param b       Sorted array of uint64_t elements (ascending order)
 * @param b_size  Number of elements in array `b`
 * @return        Array of indices (size_t*) where indices[i] = first occurrence of a[i] in b,
 *                SIZE_MAX if not found. Caller must free() the returned array.
 */
size_t *find_first_occurrences_AT_ctxgidobj_arr(const uint64_t *a, size_t a_size,
												const ctxgidobj_t *b, size_t b_size)
{
	size_t *indices = malloc(a_size * sizeof(size_t));
	if (!indices)
		return NULL;
	size_t low = 0; // Track lower bound for binary search
	int nobjbits = Bitslen.obj;

	for (size_t i = 0; i < a_size; ++i)
	{
		const uint64_t target = a[i] >> nobjbits; // (2*(2*holen + iolen ));
		// when conflict objects are kept, skip searching if target[i+1] == target[i].
		if (i > 0 && target == a[i - 1] >> nobjbits)
		{
			indices[i] = indices[i - 1]; // Use the previous index if the current target is the same
			continue;
		}

		size_t high = b_size;
		// Leftmost binary search within [low, high)
		while (low < high)
		{
			size_t mid = low + (high - low) / 2;
			if ((b[mid].ctxgid >> GID_NBITS) < target)
			{
				low = mid + 1;
			}
			else
			{
				high = mid;
			}
		}

		// Check if target was found: modified from: if (low < b_size && b[low] == target) {
		if (low < b_size && (b[low].ctxgid >> GID_NBITS == target))
			indices[i] = low;
		else
			indices[i] = SIZE_MAX; // Not found
	}

	return indices;
}

void ani_block_print(
	int ref_infile_num, int qry_gid_offset, int this_block_size,
	uint64_t *ref_sketch_index, uint64_t *qry_sketch_index,
	const uint32_t *ref_ctx_count, const uint32_t *qry_ctx_count,
	ctx_mut2_t *ctx, obj_section_t *obj,
	char (*refname)[PATHLEN], char (*qryfname)[PATHLEN],
	char (*refanno)[PATHLEN],
	const infile_meta_t *qry_infile_meta,
	const infile_meta_t *ref_infile_meta,
	const ani_ctxmeta_rec_t *qry_ctxmeta,
	const ani_ctxmeta_rec_t *ref_ctxmeta,
	uint32_t *num_passid_block, idani_t **sort_idani_block,
	FILE *outfp, ani_opt_t *ani_opt, int matrix_mode)
{
	ani_features_t ani_features;

	for (int i = 0; i < this_block_size; i++)
	{
		int qry_gid = qry_gid_offset + i;
		int qry_sketch_size = qry_ctx_count[qry_gid];

		if (matrix_mode)
		{
			fprintf(outfp, "%s", qryfname[qry_gid]);
		}

		int loop_j = matrix_mode ? ref_infile_num : num_passid_block[i];
		for (int n = 0; n < loop_j; n++)
		{
			int j = matrix_mode ? n : sort_idani_block[i][n].id;
			ani_features.XnY_ctx = MCTX(ref_infile_num, i, j).num_ctx;
			ani_features.N_diff_obj_section = MOBJ(ref_infile_num, i, j).diff_obj_section;
			ani_features.N_mut2_ctx = MCTX(ref_infile_num, i, j).num_mut2_ctx;
			ani_features.N_diff_obj = MOBJ(ref_infile_num, i, j).diff_obj;

			int ref_sketch_size = ref_ctx_count[j];
				if (qry_sketch_size == 0 || ref_sketch_size == 0)
				{
					if (matrix_mode)
					{
						double metric = ani_opt->s < 0 ? 1.0 - ani_opt->e : ani_opt->e;
						fprintf(outfp, "\t%lf", metric);
					}
					continue;
				}
			float af_qry = (float)ani_features.XnY_ctx / qry_sketch_size;
			float af_ref = (float)ani_features.XnY_ctx / ref_sketch_size;
			const ani_density_af_t density_af = ani_estimate_density_af(
				ani_ctxmeta_at(qry_ctxmeta, (uint32_t)qry_gid),
				ani_ctxmeta_at(ref_ctxmeta, (uint32_t)j),
				(uint32_t)ani_features.XnY_ctx,
				(double)af_qry, (double)af_ref);
			ani_features.X_ctx = qry_sketch_size;
			float blastn_af_qry = lm3ways_af_ANIb_from_features(&ani_features);
			ani_features.X_ctx = ref_sketch_size;
			float blastn_af_ref = lm3ways_af_ANIb_from_features(&ani_features);

				ani_row_t outrow = make_selected_output_row((uint32_t)j, &ani_features, ani_opt,
															(uint32_t)qry_sketch_size,
															(uint32_t)ref_sketch_size,
															af_qry, blastn_af_qry,
															af_ref, blastn_af_ref,
															density_af,
															qry_infile_meta ? &qry_infile_meta[qry_gid] : NULL,
															ref_infile_meta ? &ref_infile_meta[j] : NULL);
				double metric = ani_opt->s < 0 ? outrow.selected_ani : outrow.metric;

				if (matrix_mode)
				{
					fprintf(outfp, "\t%lf", metric);
				}
				else if (outrow.af_pass && outrow.selected_ani >= ani_opt->anicut)
				{
					print_unified_detail_row(outfp, ani_opt, qryfname[qry_gid], refname[j],
											 &outrow, annotation_at(refanno, (uint32_t)j));
				}
		}
		if (matrix_mode)
		{
			fprintf(outfp, "\n");
		}
	}
}

void simple_sortedsketch64Xcomb_sortedsketch64(simple_sketch_t *simple_sketch, infile_tab_t *genomes_infiletab, ani_opt_t *ani_opt)
{
	uint64_t *arr_ref = simple_sketch->comb_sketch;
	if (arr_ref == NULL)
		err(EXIT_FAILURE, "%s(): simple_sketch->comb_sketch is NULL", __func__);

	if (simple_sketch->infile_num != genomes_infiletab->infile_num)
		err(EXIT_FAILURE, "%s(): infile_num mismatch: %d vs %d", __func__, simple_sketch->infile_num, genomes_infiletab->infile_num);
	size_t ref_sketch_size = simple_sketch->sketch_index[1] - simple_sketch->sketch_index[0];
	assert(simple_sketch->infile_num == genomes_infiletab->infile_num);

	FILE *outfp = ani_opt->outf[0] == '\0' ? stdout : fopen(ani_opt->outf, "w");
	if (outfp == NULL)
		err(errno, "%s", ani_opt->outf);
	print_ani_detail_header(outfp, ani_opt, true);

	for (uint32_t qn = 1; qn < simple_sketch->infile_num; qn++)
	{
		ani_features_t ani_features;
		uint64_t *arr_qry = simple_sketch->comb_sketch + simple_sketch->sketch_index[qn];
		size_t qry_sketch_size = simple_sketch->sketch_index[qn + 1] - simple_sketch->sketch_index[qn];
		if (qry_sketch_size == 0 || ref_sketch_size == 0)
			continue;

		get_ani_features_from_two_sorted_ctxobj64(arr_ref, ref_sketch_size, arr_qry, qry_sketch_size, &ani_features);

		float af_qry = (double)ani_features.XnY_ctx / qry_sketch_size;

		float af_ref = (double)ani_features.XnY_ctx / ref_sketch_size;
		const ani_density_af_t density_af =
			ani_density_af_fallback((double)af_qry, (double)af_ref);
		ani_features.X_ctx = qry_sketch_size;
		float blastn_af_qry = lm3ways_af_ANIb_from_features(&ani_features);
		ani_features.X_ctx = ref_sketch_size;
		float blastn_af_ref = lm3ways_af_ANIb_from_features(&ani_features);

		if (!ani_report_af_pass(ani_opt, density_af.qry, density_af.ref))
			continue;
		ani_row_t outrow = make_selected_output_row(0, &ani_features, ani_opt,
													(uint32_t)qry_sketch_size,
													(uint32_t)ref_sketch_size,
													af_qry, blastn_af_qry,
													af_ref, blastn_af_ref,
													density_af,
													NULL, NULL);
		if (outrow.selected_ani >= ani_opt->anicut)
		{
			print_unified_detail_row(outfp, ani_opt,
									 genomes_infiletab->organized_infile_tab[qn].fpath,
									 genomes_infiletab->organized_infile_tab[0].fpath,
									 &outrow, NULL);
		}
	}
	if (outfp != stdout)
		fclose(outfp);
}

// fencepost search method
#include <limits.h>

static inline uint32_t top_k_bits_u64(uint64_t x, int k)
{
	if (k <= 0)
		return 0u;
	if (k >= 64)
		return (uint32_t)x; /* defensive; we clamp k well below 64 */
	return (uint32_t)(x >> (64 - k));
}

int minco_choose_k_fenceposts(size_t b_size, size_t a_size)
{
	if (a_size == 0 || b_size == 0)
		return 14;
	double ratio = (double)b_size / (double)a_size;			   /* ≈ m = b/a */
	int k = (ratio > 1.0) ? (int)floor(log2(ratio) + 0.5) : 0; /* ~log2(m) */
	if (k < 12)
		k = 12; /* 4K buckets  (~32 KB of fenceposts) */
	if (k > 20)
		k = 20; /* 1M buckets  (~8  MB of fenceposts, 64-bit size_t) */
	return k;
}

int minco_build_fenceposts_ctxgid(const ctxgidobj_t *b, size_t b_size, int k, size_t *F)
{
	if (!b || !F || k < 0 || k > 32)
		return -1;
	const size_t buckets = (size_t)1u << k;

	for (size_t t = 0; t <= buckets; ++t)
		F[t] = b_size;

	size_t next = 0;
	for (size_t i = 0; i < b_size; ++i)
	{
		uint64_t key = (uint64_t)b[i].ctxgid >> GID_NBITS; /* effective key */
		uint32_t topk = top_k_bits_u64(key, k);
		while (next <= topk)
			F[next++] = i;
		if (next > buckets)
			break;
	}
	while (next <= buckets)
		F[next++] = b_size;
	return 0;
}

static inline size_t lb_in_bucket_ctxgid(const ctxgidobj_t *b, const size_t *F, int k,
										 uint64_t target_key)
{
	const size_t buckets = (size_t)1u << k;
	uint32_t t = top_k_bits_u64(target_key, k);
	if (t > buckets)
		t = (uint32_t)buckets; /* defensive */

	size_t lo = F[t];
	size_t hi = F[t + 1];

	while (lo < hi)
	{
		size_t mid = lo + ((hi - lo) >> 1);
		uint64_t mid_key = (uint64_t)b[mid].ctxgid >> GID_NBITS;
		if (mid_key < target_key)
			lo = mid + 1;
		else
			hi = mid;
	}
	return lo; /* first position with key >= target_key (may be hi) */
}

static int minco_build_fenceposts_ctxgid128(const ctxgidobj128_t *b, size_t b_size,
											int k, size_t *F)
{
	if (!b || !F || k < 0 || k > 32)
		return -1;
	const size_t buckets = (size_t)1u << k;
	for (size_t t = 0; t <= buckets; ++t)
		F[t] = b_size;

	size_t next = 0;
	for (size_t i = 0; i < b_size; ++i) {
		uint32_t topk = top_k_bits_u64(b[i].ctx, k);
		while (next <= topk)
			F[next++] = i;
		if (next > buckets)
			break;
	}
	while (next <= buckets)
		F[next++] = b_size;
	return 0;
}

static inline size_t lb_in_bucket_ctxgid128(const ctxgidobj128_t *b,
											const size_t *F,
											int k,
											uint64_t target_key)
{
	const size_t buckets = (size_t)1u << k;
	uint32_t t = top_k_bits_u64(target_key, k);
	if (t > buckets)
		t = (uint32_t)buckets;

	size_t lo = F[t];
	size_t hi = F[t + 1];
	while (lo < hi) {
		size_t mid = lo + ((hi - lo) >> 1);
		if (b[mid].ctx < target_key)
			lo = mid + 1;
		else
			hi = mid;
	}
	return lo;
}

size_t *minco_find_first_occurrences_fenceposts(const uint64_t *a, size_t a_size,
											   const ctxgidobj_t *b, size_t b_size,
											   const size_t *F, int k, unsigned nobjbits)
{
	if (!a || !b || !F)
		return NULL;

	size_t *idx = (size_t *)malloc(a_size * sizeof *idx);
	if (!idx)
		return NULL;

	uint64_t prev_key = UINT64_MAX;
	size_t prev_idx = SIZE_MAX;

	for (size_t i = 0; i < a_size; ++i)
	{
		uint64_t target_key = a[i] >> nobjbits;

		if (i && target_key == prev_key)
		{ /* reuse for duplicates in a */
			idx[i] = prev_idx;
			continue;
		}

		size_t pos = lb_in_bucket_ctxgid(b, F, k, target_key);

		if (pos < b_size && (((uint64_t)b[pos].ctxgid >> GID_NBITS) == target_key))
			idx[i] = pos;
		else
			idx[i] = SIZE_MAX;

		prev_key = target_key;
		prev_idx = idx[i];
	}
	return idx;
}

#ifndef MINCO_HASH_BOTTOMK
#define MINCO_HASH_BOTTOMK 0
#endif
#ifndef MINCO_KEEP_SOURCE_FILTER
#define MINCO_KEEP_SOURCE_FILTER 1
#endif
#ifndef MINCO_SEED
#define MINCO_SEED 0x9e3779b97f4a7c15ULL
#endif
#ifndef MINCO_HASH_SPARSE_CTX
#define MINCO_HASH_SPARSE_CTX 0
#endif
#ifndef likely
#define likely(x) __builtin_expect(!!(x), 1)
#endif
#ifndef unlikely
#define unlikely(x) __builtin_expect(!!(x), 0)
#endif

#define ANI_APPLY_SOURCE_FILTER (!MINCO_HASH_BOTTOMK || MINCO_KEEP_SOURCE_FILTER)
#define ANI_U64SET_MAX_LOAD_NUM 7u
#define ANI_U64SET_MAX_LOAD_DEN 10u
#define ANI_READWISE_BATCH_READS 16384u

#ifndef MINCO_REPORT_FILTERED_READWISE_ANI
#define MINCO_REPORT_FILTERED_READWISE_ANI 0
#endif

typedef struct {
	uint64_t *keys;
	uint8_t *used;
	size_t cap;
	size_t n;
} ani_u64_set_t;

typedef struct {
	uint64_t XnY_ctx;
	uint64_t N_diff_obj;
	uint64_t N_diff_obj_section;
	uint64_t N_mut2_ctx;
	uint64_t qry_ctx_hit;
	uint64_t ref_ctx_hit;
	uint64_t reads_with_ctx_match;
	uint64_t blocks_with_ctx_match;
} ani_readwise_acc_t;

typedef struct {
	double ref_breadth;
	double ref_mean_depth;
	double ref_hit_mean_depth;
	double ref_depth_variance;
	double ref_depth_cv;
	double ref_zero_fraction;
	double relative_depth;
} ani_readwise_abundance_t;

typedef struct {
	double ref_breadth;
	double ref_mean_depth;
	uint64_t ref_hit_ctx;
	double ref_hit_mean_depth;
	double ref_hit_median_depth;
	double ref_hit_depth_variance;
	double ref_zip_af;
} ani_readwise_reliable_abundance_t;

typedef struct {
	uint64_t raw_xny_ctx;
	uint64_t rejected_ctx;
	uint64_t rejected_diff_ctx;
	long double fake_prob_sum;
	long double fake_prob_weighted_sum;
	long double fake_prob_weight_sum;
} ani_readwise_filter_stats_t;

typedef struct {
	gzFile gz;
	FILE *pipe_fp;
	char *cmd;
	uint64_t input_size;
	bool input_size_known;
} ani_fastx_stream_t;

typedef struct {
	bool enabled;
	bool input_size_known;
	uint64_t input_size;
	uint64_t next_reads;
	time_t started_at;
	time_t last_at;
} ani_readwise_progress_t;

typedef kvec_t(size_t) kv_size_t;
typedef kvec_t(uint64_t) kv_u64_t;
typedef kvec_t(uint32_t) kv_u32_t;

typedef struct {
	uint64_t packed;
	uint32_t offset;
} ani_read_ctxobj64_t;

typedef struct {
	ctxobj96_t rec;
	uint32_t offset;
} ani_read_ctxobj96_t;

typedef kvec_t(ani_read_ctxobj64_t) kv_read_ctxobj64_t;
typedef kvec_t(ani_read_ctxobj96_t) kv_read_ctxobj96_t;

typedef struct {
	ctxobj96_t *a;
	size_t n;
	size_t m;
} ani_ctxobj96_vec_t;

typedef struct {
	ani_ctxobj96_vec_t vec;
	uint64_t id;
	uint64_t reads_with_density_ctx;
	char *read_name;
} ani_density_unit96_t;

typedef struct {
	double *a;
	uint32_t n;
	uint32_t cap;
} ani_double_vec_t;

typedef struct {
	double product_threshold;
	uint32_t median_diff;
	double median_cov;
} ani_readwise_product_topfrac_median_t;

typedef struct {
	size_t ref_begin;
	uint32_t gid;
	uint32_t diff;
} ani_readwise_candidate_t;

typedef kvec_t(ani_readwise_candidate_t) kv_readwise_candidate_t;

typedef struct {
	char *seq;
	char *name;
	int len;
} ani_readwise_seq_rec_t;

typedef struct {
	u64vec vec;
	uint64_t id;
	uint64_t reads_with_density_ctx;
	char *read_name;
} ani_density_unit_t;

typedef kvec_t(ani_density_unit_t) kv_density_unit_t;
typedef kvec_t(ani_density_unit96_t) kv_density_unit96_t;

static void ani_ctxobj96_vec_free(ani_ctxobj96_vec_t *v);

typedef struct {
	FILE *fp;
	uint32_t gid;
	const char *ref_name;
} ani_readwise_trace_t;

typedef struct {
	FILE *fp;
	uint64_t emitted_edges;
	uint64_t max_edges;
	uint32_t max_candidate_refs;
	bool ambiguous_only;
	bool selected_only;
} ani_readwise_edge_trace_t;

typedef struct {
	bool enabled;
	kv_cami_tax_record_t records;
	const ani_cami_tax_record_t **ref_tax;
	const char *label;
} ani_readwise_tax_namespace_t;

typedef struct {
	bool enabled;
	FILE *fp;
	char summary_path[PATHLEN];
	ani_readwise_taxonomy_mode_t taxonomy_mode;
	ani_readwise_tax_namespace_t gtdb;
	ani_readwise_tax_namespace_t ncbi;
	const char *query_path;
	char (*refname)[PATHLEN];
	char (*refanno)[PATHLEN];
	uint32_t ref_n;
	const ani_ctxmeta_rec_t *ref_ctxmeta;
	long double density_probability;
	bool sketch_corrected_available;
	uint64_t total_reads;
	uint64_t reads_with_density_ctx;
	uint64_t reads_with_ref_hit;
	uint64_t reads_with_multi_ref_hit;
	uint64_t reads_without_ref_hit_but_density_ctx;
	uint64_t total_possible_ctx;
	uint64_t total_density_ctx;
	uint64_t total_matched_ctx;
	uint64_t total_selected_ctx;
	uint64_t total_selected_ref_events;
	long double expected_trackable_reads;
	uint64_t sketch_corrected_observed_ctx;
	uint64_t sketch_corrected_missing_meta_ctx;
	long double sketch_corrected_present_ctx;
	long double sketch_corrected_capture_prob_sum;
} ani_readwise_tracker_t;

typedef struct {
	ani_readwise_acc_t *acc;
	uint64_t *read_marks;
	uint8_t *ref_hit_bits;
	ani_u64_set_t qry_ctx_seen;
	ani_u64_set_t qry_ref_ctx_seen;
	kv_size_t ref_ctx_hit_positions;
	kv_u64_t ref_ctx_cov_hits;
} ani_readwise_thread_state_t;

#define ANI_READWISE_PROGRESS_READ_INTERVAL 1000000ULL
#define ANI_READWISE_PROGRESS_TIME_INTERVAL 30

static inline size_t ani_next_pow2_size(size_t x)
{
	if (x <= 2)
		return 2;
	--x;
	for (size_t s = 1; s < sizeof(size_t) * CHAR_BIT; s <<= 1)
		x |= x >> s;
	return x + 1;
}

static void ani_u64_set_init(ani_u64_set_t *set, size_t initial)
{
	memset(set, 0, sizeof(*set));
	set->cap = ani_next_pow2_size(initial < 1024 ? 1024 : initial);
	set->keys = calloc(set->cap, sizeof(set->keys[0]));
	set->used = calloc(set->cap, sizeof(set->used[0]));
	if (!set->keys || !set->used)
		err(EXIT_FAILURE, "%s(): OOM hash set", __func__);
}

static void ani_u64_set_destroy(ani_u64_set_t *set)
{
	if (!set)
		return;
	free(set->keys);
	free(set->used);
	memset(set, 0, sizeof(*set));
}

static void ani_u64_set_clear(ani_u64_set_t *set)
{
	if (!set || !set->used)
		return;
	memset(set->used, 0, set->cap * sizeof(set->used[0]));
	set->n = 0;
}

static size_t ani_u64_set_capacity_for_expected(size_t expected)
{
	if (!expected)
		return 1024;
	if (expected > (SIZE_MAX - ANI_U64SET_MAX_LOAD_NUM) / ANI_U64SET_MAX_LOAD_DEN)
		errx(EXIT_FAILURE, "%s(): hash set too large", __func__);
	const size_t min_cap =
		(expected * ANI_U64SET_MAX_LOAD_DEN + ANI_U64SET_MAX_LOAD_NUM - 1) /
		ANI_U64SET_MAX_LOAD_NUM;
	return ani_next_pow2_size(min_cap < 1024 ? 1024 : min_cap);
}

static void ani_u64_set_reset_for_expected(ani_u64_set_t *set, size_t expected)
{
	if (!set)
		return;
	free(set->keys);
	free(set->used);
	memset(set, 0, sizeof(*set));
	set->cap = ani_u64_set_capacity_for_expected(expected);
	set->keys = calloc(set->cap, sizeof(set->keys[0]));
	set->used = calloc(set->cap, sizeof(set->used[0]));
	if (!set->keys || !set->used)
		err(EXIT_FAILURE, "%s(): OOM hash set", __func__);
}

static void ani_u64_set_grow(ani_u64_set_t *set)
{
	ani_u64_set_t next;
	ani_u64_set_init(&next, set->cap << 1);
	for (size_t i = 0; i < set->cap; ++i) {
		if (!set->used[i])
			continue;
		uint64_t key = set->keys[i];
		size_t pos = (size_t)mix64(key) & (next.cap - 1);
		while (next.used[pos])
			pos = (pos + 1) & (next.cap - 1);
		next.used[pos] = 1;
		next.keys[pos] = key;
		next.n++;
	}
	free(set->keys);
	free(set->used);
	*set = next;
}

static bool ani_u64_set_insert(ani_u64_set_t *set, uint64_t key)
{
	if ((set->n + 1) * ANI_U64SET_MAX_LOAD_DEN >
		set->cap * ANI_U64SET_MAX_LOAD_NUM)
		ani_u64_set_grow(set);
	size_t pos = (size_t)mix64(key) & (set->cap - 1);
	for (;;) {
		if (!set->used[pos]) {
			set->used[pos] = 1;
			set->keys[pos] = key;
			set->n++;
			return true;
		}
		if (set->keys[pos] == key)
			return false;
		pos = (pos + 1) & (set->cap - 1);
	}
}

static inline uint32_t ani_clamp_u64_to_u32(uint64_t x)
{
	return x > (uint64_t)UINT32_MAX ? UINT32_MAX : (uint32_t)x;
}

static inline int ani_clamp_u64_to_int(uint64_t x)
{
	return x > (uint64_t)INT_MAX ? INT_MAX : (int)x;
}

static inline bool ani_bitset_test_set(uint8_t *bits, size_t idx)
{
	uint8_t *byte = &bits[idx >> 3];
	const uint8_t mask = (uint8_t)(1u << (idx & 7u));
	const bool was_set = (*byte & mask) != 0;
	*byte |= mask;
	return was_set;
}

#define ANI_REFCOV_COV_MASK 0x0fffffffU
#define ANI_REFCOV_DIFF_SHIFT 28u
#define ANI_REFCOV_DIFF_MAX_STORED 14u
#define ANI_REFCOV_HIT_INC_BITS 16u
#define ANI_REFCOV_HIT_INC_MASK ((1ULL << ANI_REFCOV_HIT_INC_BITS) - 1ULL)
#define ANI_REFCOV_HIT_SHIFT (ANI_REFCOV_HIT_INC_BITS + 4u)
#define ANI_REFCOV_SPLIT_SCALE 1024u
#define ANI_READWISE_PROB_SCALE 1024u

static inline uint32_t ani_ref_covdiff_coverage(uint32_t x)
{
	return x & ANI_REFCOV_COV_MASK;
}

static inline uint32_t ani_ref_covdiff_clip_diff(uint32_t diff)
{
	return diff > ANI_REFCOV_DIFF_MAX_STORED
			   ? ANI_REFCOV_DIFF_MAX_STORED
			   : diff;
}

static inline uint32_t ani_ref_covdiff_encode_diff(uint32_t diff)
{
	return ani_ref_covdiff_clip_diff(diff) + 1u;
}

static inline uint32_t ani_ref_covdiff_decode_diff(uint32_t code)
{
	return code ? code - 1u : 0u;
}

static inline void ani_ref_covdiff_add_scaled_hit(uint32_t *x, uint32_t diff, uint32_t inc)
{
	if (!inc)
		return;
	uint32_t cov = ani_ref_covdiff_coverage(*x);
	if (cov > ANI_REFCOV_COV_MASK - inc)
		cov = ANI_REFCOV_COV_MASK;
	else
		cov += inc;
	uint32_t code = *x >> ANI_REFCOV_DIFF_SHIFT;
	const uint32_t new_code = ani_ref_covdiff_encode_diff(diff);
	if (code == 0u || new_code < code)
		code = new_code;
	*x = (code << ANI_REFCOV_DIFF_SHIFT) | cov;
}

static inline void ani_ref_covmindiff_add_scaled_hit(uint32_t *x, uint32_t diff, uint32_t inc)
{
	if (!inc)
		return;
	const uint32_t new_code = ani_ref_covdiff_encode_diff(diff);
	uint32_t code = *x >> ANI_REFCOV_DIFF_SHIFT;
	uint32_t cov = ani_ref_covdiff_coverage(*x);
	if (code == 0u || new_code < code) {
		code = new_code;
		cov = inc;
	} else if (new_code == code) {
		if (cov > ANI_REFCOV_COV_MASK - inc)
			cov = ANI_REFCOV_COV_MASK;
		else
			cov += inc;
	} else {
		return;
	}
	*x = (code << ANI_REFCOV_DIFF_SHIFT) | cov;
}

static inline void ani_ref_covdiff_add_hit(uint32_t *x, uint32_t diff)
{
	ani_ref_covdiff_add_scaled_hit(x, diff, 1u);
}

static inline uint64_t ani_ref_covdiff_pack_hit(size_t idx, uint32_t diff, uint32_t inc)
{
	if (inc > (uint32_t)ANI_REFCOV_HIT_INC_MASK)
		inc = (uint32_t)ANI_REFCOV_HIT_INC_MASK;
	return ((uint64_t)idx << ANI_REFCOV_HIT_SHIFT) |
		   ((uint64_t)inc << 4) |
		   (uint64_t)ani_ref_covdiff_clip_diff(diff);
}

static double ani_poisson_tail_neglog10(uint32_t k, double lambda)
{
	if (k == 0)
		return 0.0;
	if (!isfinite(lambda) || lambda <= 0.0)
		return 16.0;
	if (k > 1024u)
		return 16.0;
	long double term = expl(-(long double)lambda);
	long double cdf = term;
	for (uint32_t i = 1; i < k; ++i) {
		term *= (long double)lambda / (long double)i;
		cdf += term;
		if (cdf >= 1.0L)
			return 0.0;
	}
	long double tail = 1.0L - cdf;
	if (tail <= 1e-16L)
		return 16.0;
	return (double)(-log10l(tail));
}

static double ani_readwise_poisson_diff_score(uint32_t diff,
											  double cov,
											  const ani_readwise_abundance_t *abund)
{
	if (!abund)
		return (double)diff;
	double zero = abund->ref_zero_fraction;
	if (!isfinite(zero))
		zero = 0.0;
	if (zero <= 0.0)
		zero = 1e-12;
	if (zero >= 1.0)
		zero = 1.0 - 1e-12;
	const double lambda_breadth = -log(zero);
	const double lambda = lambda_breadth > 1e-12
							  ? lambda_breadth
							  : (abund->ref_mean_depth > 1e-12 ? abund->ref_mean_depth : 1e-12);
	const uint32_t k = cov > 0.0 ? (uint32_t)ceil(cov) : 0u;
	const double depth_surprise = ani_poisson_tail_neglog10(k, lambda);
	double pressure = 0.0;
	if (lambda_breadth > 1e-12 && abund->ref_mean_depth > 0.0)
		pressure = log1p(abund->ref_mean_depth / lambda_breadth);
	return (double)diff + 0.30 * depth_surprise + 0.25 * pressure;
}

static inline double ani_logistic(double x)
{
	if (x >= 40.0)
		return 1.0;
	if (x <= -40.0)
		return 0.0;
	return 1.0 / (1.0 + exp(-x));
}

static double ani_readwise_fake_ctx_probability(uint32_t diff,
											   double cov,
											   const ani_readwise_abundance_t *abund,
											   double threshold)
{
	if (!abund || threshold <= 0.0)
		return 0.0;
	const double score = ani_readwise_poisson_diff_score(diff, cov, abund);
	return ani_logistic(score - threshold);
}

static bool ani_readwise_reject_poisson_diff_ctx(uint32_t diff,
												 double cov,
												 const ani_readwise_abundance_t *abund,
												 double threshold)
{
	if (diff == 0)
		return false;
	if (!abund || threshold <= 0.0)
		return false;
	const double score = ani_readwise_poisson_diff_score(diff, cov, abund);
	return score >= threshold;
}

static bool ani_readwise_reject_poisson_depth_ctx(uint32_t diff,
												  double cov,
												  const ani_readwise_abundance_t *abund,
												  uint32_t ref_ctx_total,
												  double threshold)
{
	if (diff == 0)
		return false;
	if (!abund || threshold <= 0.0)
		return false;
	double zero = abund->ref_zero_fraction;
	if (!isfinite(zero))
		zero = 0.0;
	if (zero <= 0.0)
		zero = 1e-12;
	if (zero >= 1.0)
		zero = 1.0 - 1e-12;
	const double lambda_breadth = -log(zero);
	const double lambda = lambda_breadth > 1e-12
							  ? lambda_breadth
							  : (abund->ref_mean_depth > 1e-12 ? abund->ref_mean_depth : 1e-12);
	if (!isfinite(cov) || cov <= lambda)
		return false;
	const uint32_t k = cov > 0.0 ? (uint32_t)ceil(cov) : 0u;
	double score = ani_poisson_tail_neglog10(k, lambda);
	if (ref_ctx_total > 1u)
		score -= log10((double)ref_ctx_total);
	return score >= threshold;
}

static bool ani_readwise_reject_poisson_product_ctx(uint32_t diff,
													double cov,
													double lambda,
													uint32_t ref_ctx_total,
													double threshold)
{
	if (diff == 0)
		return false;
	if (threshold <= 0.0)
		return false;
	if (!isfinite(cov) || cov <= 0.0)
		return false;
	if (!isfinite(lambda) || lambda <= 0.0)
		lambda = 1e-12;
	const double product = (double)diff * cov;
	if (product <= lambda)
		return false;
	const uint32_t k = (uint32_t)ceil(product);
	double score = ani_poisson_tail_neglog10(k, lambda);
	if (ref_ctx_total > 1u)
		score -= log10((double)ref_ctx_total);
	return score >= threshold;
}

#define ANI_READWISE_PRODUCT_NB_MIN_NONZERO 10u

static inline bool ani_readwise_ctx_filter_uses_product(ani_readwise_ctx_filter_model_t filter_model)
{
	return filter_model == ANI_READWISE_CTX_FILTER_POISSON_PRODUCT ||
		   filter_model == ANI_READWISE_CTX_FILTER_PRODUCT_NB;
}

static inline double ani_logaddexp(double a, double b)
{
	if (!isfinite(a))
		return b;
	if (!isfinite(b))
		return a;
	if (a < b) {
		const double tmp = a;
		a = b;
		b = tmp;
	}
	return a + log1p(exp(b - a));
}

static double ani_nb_tail_neglog10(uint32_t k, double mean, double variance)
{
	if (k == 0)
		return 0.0;
	if (!isfinite(mean) || mean <= 0.0)
		mean = 1e-12;
	if (k > 1024u)
		return 16.0;
	if (!isfinite(variance) || variance <= mean)
		return ani_poisson_tail_neglog10(k, mean);

	const double r = (mean * mean) / (variance - mean);
	if (!isfinite(r) || r <= 0.0)
		return ani_poisson_tail_neglog10(k, mean);
	const double p = r / (r + mean);
	if (!isfinite(p) || p <= 0.0 || p >= 1.0)
		return ani_poisson_tail_neglog10(k, mean);

	const double log_p = log(p);
	const double log_q = log1p(-p);
	double log_pmf = r * log_p; /* NB count of failures before r successes. */
	double log_cdf = log_pmf;
	for (uint32_t i = 1; i < k; ++i) {
		log_pmf += log(((double)i - 1.0 + r) / (double)i) + log_q;
		log_cdf = ani_logaddexp(log_cdf, log_pmf);
	}
	if (log_cdf >= 0.0)
		return 16.0;
	const double cdf = exp(log_cdf);
	if (cdf >= 1.0)
		return 16.0;
	const double tail = 1.0 - cdf;
	if (tail <= 1e-16)
		return 16.0;

	const double p0 = exp(r * log_p);
	if (p0 > 0.0 && p0 < 1.0) {
		const double pos = 1.0 - p0;
		if (pos > 0.0) {
			const double conditional_tail = tail / pos;
			if (conditional_tail > 0.0 && conditional_tail <= 1.0)
				return -log10(conditional_tail);
		}
	}
	return -log10(tail);
}

static bool ani_readwise_reject_product_nb_ctx(uint32_t diff,
											   double cov,
											   double mean,
											   double variance,
											   uint32_t product_nonzero,
											   uint32_t ref_ctx_total,
											   double threshold)
{
	if (diff == 0)
		return false;
	if (threshold <= 0.0)
		return false;
	if (product_nonzero < ANI_READWISE_PRODUCT_NB_MIN_NONZERO)
		return false;
	if (!isfinite(cov) || cov <= 0.0)
		return false;
	if (!isfinite(mean) || mean <= 0.0)
		return false;
	const double product = (double)diff * cov;
	if (!isfinite(product) || product <= mean)
		return false;
	const uint32_t k = product >= (double)UINT32_MAX
						   ? UINT32_MAX
						   : (uint32_t)ceil(product);
	double score = ani_nb_tail_neglog10(k, mean, variance);
	if (ref_ctx_total > 1u)
		score -= log10((double)ref_ctx_total);
	return score >= threshold;
}

static int ani_cmp_double_asc(const void *pa, const void *pb)
{
	const double a = *(const double *)pa;
	const double b = *(const double *)pb;
	return (a > b) - (a < b);
}

static void ani_double_vec_push(ani_double_vec_t *v, double x)
{
	if (!v)
		return;
	if (v->n == v->cap) {
		uint32_t new_cap = v->cap ? v->cap * 2u : 4u;
		double *next = realloc(v->a, (size_t)new_cap * sizeof(next[0]));
		if (!next)
			err(EXIT_FAILURE, "%s(): OOM product vector", __func__);
		v->a = next;
		v->cap = new_cap;
	}
	v->a[v->n++] = x;
}

static void ani_double_vecs_destroy(ani_double_vec_t *vecs, uint32_t n)
{
	if (!vecs)
		return;
	for (uint32_t i = 0; i < n; ++i)
		free(vecs[i].a);
	free(vecs);
}

static double ani_readwise_product_topfrac_fraction(double threshold)
{
	if (!isfinite(threshold) || threshold <= 0.0)
		return 0.0;
	if (threshold > 1.0)
		return 1.0;
	return threshold;
}

static inline bool ani_ctxgid_group_marker_unique64(const ctxgidobj_t *index,
													size_t index_n,
													size_t begin,
													size_t end)
{
	if (!index || begin >= index_n || end <= begin || end > index_n)
		return false;
	const uint64_t ctx = index[begin].ctxgid >> GID_NBITS;
	if (begin > 0 && (index[begin - 1].ctxgid >> GID_NBITS) == ctx)
		return false;
	if (end < index_n && (index[end].ctxgid >> GID_NBITS) == ctx)
		return false;
	return true;
}

static inline bool ani_ctxgid_group_marker_unique128(const ctxgidobj128_t *index,
													 size_t index_n,
													 size_t begin,
													 size_t end)
{
	if (!index || begin >= index_n || end <= begin || end > index_n)
		return false;
	const uint64_t ctx = index[begin].ctx;
	if (begin > 0 && index[begin - 1].ctx == ctx)
		return false;
	if (end < index_n && index[end].ctx == ctx)
		return false;
	return true;
}

static double *ani_readwise_product_topfrac_thresholds(
	const ctxgidobj_t *index,
	size_t index_n,
	uint32_t ref_n,
	bool ignoreconflict,
	const uint32_t *ref_ctx_cov,
	double coverage_scale,
	double fraction,
	bool marker_only)
{
	if (!index || !ref_ctx_cov || ref_n == 0)
		return NULL;
	fraction = ani_readwise_product_topfrac_fraction(fraction);
	if (fraction <= 0.0)
		return NULL;
	if (!isfinite(coverage_scale) || coverage_scale <= 0.0)
		coverage_scale = 1.0;

	ani_double_vec_t *products = calloc((size_t)ref_n, sizeof(products[0]));
	double *thresholds = malloc((size_t)ref_n * sizeof(thresholds[0]));
	if (!products || !thresholds)
		err(EXIT_FAILURE, "%s(): OOM product top-fraction thresholds", __func__);
	for (uint32_t rn = 0; rn < ref_n; ++rn)
		thresholds[rn] = NAN;

	const uint64_t gidmask_local = (1ULL << GID_NBITS) - 1ULL;
	for (size_t i = 0; i < index_n; ) {
		const uint64_t ctxgid = index[i].ctxgid;
		const uint32_t gid = (uint32_t)(ctxgid & gidmask_local);
		const size_t begin = i;
		do { ++i; } while (i < index_n && index[i].ctxgid == ctxgid);
		if (gid >= ref_n || (ignoreconflict && i - begin > 1) ||
			(marker_only && !ani_ctxgid_group_marker_unique64(index, index_n, begin, i)))
			continue;
		const uint32_t packed = ref_ctx_cov[begin];
		const uint32_t coverage_raw = ani_ref_covdiff_coverage(packed);
		const uint32_t code = packed >> ANI_REFCOV_DIFF_SHIFT;
		if (!code || !coverage_raw)
			continue;
		const uint32_t diff = ani_ref_covdiff_decode_diff(code);
		const double cov = (double)coverage_raw / coverage_scale;
		const double product = (double)diff * cov;
		if (!isfinite(product) || product < 0.0)
			continue;
		ani_double_vec_push(&products[gid], product);
	}

	for (uint32_t rn = 0; rn < ref_n; ++rn) {
		ani_double_vec_t *v = &products[rn];
		if (!v->n)
			continue;
		qsort(v->a, v->n, sizeof(v->a[0]), ani_cmp_double_asc);
		uint32_t drop_n = (uint32_t)ceil((double)v->n * fraction);
		if (drop_n < 1u)
			drop_n = 1u;
		if (drop_n > v->n)
			drop_n = v->n;
		thresholds[rn] = v->a[v->n - drop_n];
	}

	ani_double_vecs_destroy(products, ref_n);
	return thresholds;
}

static double ani_double_vec_lower_median(ani_double_vec_t *v)
{
	if (!v || !v->n)
		return NAN;
	qsort(v->a, v->n, sizeof(v->a[0]), ani_cmp_double_asc);
	return v->a[(v->n - 1u) / 2u];
}

static inline double ani_readwise_topfrac_product_score(uint32_t diff,
														double cov,
														bool plus_one)
{
	if (!isfinite(cov) || cov <= 0.0)
		return NAN;
	const double factor = plus_one ? (double)diff + 1.0 : (double)diff;
	return factor * cov;
}

static ani_readwise_product_topfrac_median_t *ani_readwise_product_topfrac_median_stats(
	const ctxgidobj_t *index,
	size_t index_n,
	uint32_t ref_n,
	bool ignoreconflict,
	const uint32_t *ref_ctx_cov,
	double coverage_scale,
	double fraction,
	bool plus_one,
	bool marker_only)
{
	if (!index || !ref_ctx_cov || ref_n == 0)
		return NULL;
	fraction = ani_readwise_product_topfrac_fraction(fraction);
	if (fraction <= 0.0)
		return NULL;
	if (!isfinite(coverage_scale) || coverage_scale <= 0.0)
		coverage_scale = 1.0;

	ani_double_vec_t *products = calloc((size_t)ref_n, sizeof(products[0]));
	ani_double_vec_t *diffs = calloc((size_t)ref_n, sizeof(diffs[0]));
	ani_double_vec_t *covs = calloc((size_t)ref_n, sizeof(covs[0]));
	ani_readwise_product_topfrac_median_t *stats =
		malloc((size_t)ref_n * sizeof(stats[0]));
	if (!products || !diffs || !covs || !stats)
		err(EXIT_FAILURE, "%s(): OOM product top-fraction median stats", __func__);
	for (uint32_t rn = 0; rn < ref_n; ++rn) {
		stats[rn].product_threshold = NAN;
		stats[rn].median_diff = 0u;
		stats[rn].median_cov = NAN;
	}

	const uint64_t gidmask_local = (1ULL << GID_NBITS) - 1ULL;
	for (size_t i = 0; i < index_n; ) {
		const uint64_t ctxgid = index[i].ctxgid;
		const uint32_t gid = (uint32_t)(ctxgid & gidmask_local);
		const size_t begin = i;
		do { ++i; } while (i < index_n && index[i].ctxgid == ctxgid);
		if (gid >= ref_n || (ignoreconflict && i - begin > 1) ||
			(marker_only && !ani_ctxgid_group_marker_unique64(index, index_n, begin, i)))
			continue;
		const uint32_t packed = ref_ctx_cov[begin];
		const uint32_t coverage_raw = ani_ref_covdiff_coverage(packed);
		const uint32_t code = packed >> ANI_REFCOV_DIFF_SHIFT;
		if (!code || !coverage_raw)
			continue;
		const uint32_t diff = ani_ref_covdiff_decode_diff(code);
		const double cov = (double)coverage_raw / coverage_scale;
		const double product = ani_readwise_topfrac_product_score(diff, cov, plus_one);
		if (!isfinite(product) || product < 0.0 || !isfinite(cov) || cov <= 0.0)
			continue;
		ani_double_vec_push(&products[gid], product);
		ani_double_vec_push(&diffs[gid], (double)diff);
		ani_double_vec_push(&covs[gid], cov);
	}

	for (uint32_t rn = 0; rn < ref_n; ++rn) {
		ani_double_vec_t *pv = &products[rn];
		if (!pv->n)
			continue;
		qsort(pv->a, pv->n, sizeof(pv->a[0]), ani_cmp_double_asc);
		uint32_t drop_n = (uint32_t)ceil((double)pv->n * fraction);
		if (drop_n < 1u)
			drop_n = 1u;
		if (drop_n > pv->n)
			drop_n = pv->n;
		stats[rn].product_threshold = pv->a[pv->n - drop_n];

		const double median_diff = ani_double_vec_lower_median(&diffs[rn]);
		const double median_cov = ani_double_vec_lower_median(&covs[rn]);
		if (isfinite(median_diff) && median_diff > 0.0)
			stats[rn].median_diff = median_diff >= (double)UINT32_MAX
										 ? UINT32_MAX
										 : (uint32_t)llround(median_diff);
		if (isfinite(median_cov) && median_cov > 0.0)
			stats[rn].median_cov = median_cov;
	}

	ani_double_vecs_destroy(products, ref_n);
	ani_double_vecs_destroy(diffs, ref_n);
	ani_double_vecs_destroy(covs, ref_n);
	return stats;
}

static double *ani_readwise_product_topfrac_thresholds128(
	const ctxgidobj128_t *index,
	size_t index_n,
	uint32_t ref_n,
	bool ignoreconflict,
	const uint32_t *ref_ctx_cov,
	double coverage_scale,
	double fraction,
	bool marker_only)
{
	if (!index || !ref_ctx_cov || ref_n == 0)
		return NULL;
	fraction = ani_readwise_product_topfrac_fraction(fraction);
	if (fraction <= 0.0)
		return NULL;
	if (!isfinite(coverage_scale) || coverage_scale <= 0.0)
		coverage_scale = 1.0;

	ani_double_vec_t *products = calloc((size_t)ref_n, sizeof(products[0]));
	double *thresholds = malloc((size_t)ref_n * sizeof(thresholds[0]));
	if (!products || !thresholds)
		err(EXIT_FAILURE, "%s(): OOM product top-fraction thresholds", __func__);
	for (uint32_t rn = 0; rn < ref_n; ++rn)
		thresholds[rn] = NAN;

	for (size_t i = 0; i < index_n; ) {
		const uint64_t ctx = index[i].ctx;
		const uint32_t gid = index[i].gid;
		const size_t begin = i;
		do { ++i; } while (i < index_n &&
							index[i].ctx == ctx &&
							index[i].gid == gid);
		if (gid >= ref_n || (ignoreconflict && i - begin > 1) ||
			(marker_only && !ani_ctxgid_group_marker_unique128(index, index_n, begin, i)))
			continue;
		const uint32_t packed = ref_ctx_cov[begin];
		const uint32_t coverage_raw = ani_ref_covdiff_coverage(packed);
		const uint32_t code = packed >> ANI_REFCOV_DIFF_SHIFT;
		if (!code || !coverage_raw)
			continue;
		const uint32_t diff = ani_ref_covdiff_decode_diff(code);
		const double cov = (double)coverage_raw / coverage_scale;
		const double product = (double)diff * cov;
		if (!isfinite(product) || product < 0.0)
			continue;
		ani_double_vec_push(&products[gid], product);
	}

	for (uint32_t rn = 0; rn < ref_n; ++rn) {
		ani_double_vec_t *v = &products[rn];
		if (!v->n)
			continue;
		qsort(v->a, v->n, sizeof(v->a[0]), ani_cmp_double_asc);
		uint32_t drop_n = (uint32_t)ceil((double)v->n * fraction);
		if (drop_n < 1u)
			drop_n = 1u;
		if (drop_n > v->n)
			drop_n = v->n;
		thresholds[rn] = v->a[v->n - drop_n];
	}

	ani_double_vecs_destroy(products, ref_n);
	return thresholds;
}

static ani_readwise_product_topfrac_median_t *ani_readwise_product_topfrac_median_stats128(
	const ctxgidobj128_t *index,
	size_t index_n,
	uint32_t ref_n,
	bool ignoreconflict,
	const uint32_t *ref_ctx_cov,
	double coverage_scale,
	double fraction,
	bool plus_one,
	bool marker_only)
{
	if (!index || !ref_ctx_cov || ref_n == 0)
		return NULL;
	fraction = ani_readwise_product_topfrac_fraction(fraction);
	if (fraction <= 0.0)
		return NULL;
	if (!isfinite(coverage_scale) || coverage_scale <= 0.0)
		coverage_scale = 1.0;

	ani_double_vec_t *products = calloc((size_t)ref_n, sizeof(products[0]));
	ani_double_vec_t *diffs = calloc((size_t)ref_n, sizeof(diffs[0]));
	ani_double_vec_t *covs = calloc((size_t)ref_n, sizeof(covs[0]));
	ani_readwise_product_topfrac_median_t *stats =
		malloc((size_t)ref_n * sizeof(stats[0]));
	if (!products || !diffs || !covs || !stats)
		err(EXIT_FAILURE, "%s(): OOM product top-fraction median stats", __func__);
	for (uint32_t rn = 0; rn < ref_n; ++rn) {
		stats[rn].product_threshold = NAN;
		stats[rn].median_diff = 0u;
		stats[rn].median_cov = NAN;
	}

	for (size_t i = 0; i < index_n; ) {
		const uint64_t ctx = index[i].ctx;
		const uint32_t gid = index[i].gid;
		const size_t begin = i;
		do { ++i; } while (i < index_n &&
							index[i].ctx == ctx &&
							index[i].gid == gid);
		if (gid >= ref_n || (ignoreconflict && i - begin > 1) ||
			(marker_only && !ani_ctxgid_group_marker_unique128(index, index_n, begin, i)))
			continue;
		const uint32_t packed = ref_ctx_cov[begin];
		const uint32_t coverage_raw = ani_ref_covdiff_coverage(packed);
		const uint32_t code = packed >> ANI_REFCOV_DIFF_SHIFT;
		if (!code || !coverage_raw)
			continue;
		const uint32_t diff = ani_ref_covdiff_decode_diff(code);
		const double cov = (double)coverage_raw / coverage_scale;
		const double product = ani_readwise_topfrac_product_score(diff, cov, plus_one);
		if (!isfinite(product) || product < 0.0 || !isfinite(cov) || cov <= 0.0)
			continue;
		ani_double_vec_push(&products[gid], product);
		ani_double_vec_push(&diffs[gid], (double)diff);
		ani_double_vec_push(&covs[gid], cov);
	}

	for (uint32_t rn = 0; rn < ref_n; ++rn) {
		ani_double_vec_t *pv = &products[rn];
		if (!pv->n)
			continue;
		qsort(pv->a, pv->n, sizeof(pv->a[0]), ani_cmp_double_asc);
		uint32_t drop_n = (uint32_t)ceil((double)pv->n * fraction);
		if (drop_n < 1u)
			drop_n = 1u;
		if (drop_n > pv->n)
			drop_n = pv->n;
		stats[rn].product_threshold = pv->a[pv->n - drop_n];

		const double median_diff = ani_double_vec_lower_median(&diffs[rn]);
		const double median_cov = ani_double_vec_lower_median(&covs[rn]);
		if (isfinite(median_diff) && median_diff > 0.0)
			stats[rn].median_diff = median_diff >= (double)UINT32_MAX
										 ? UINT32_MAX
										 : (uint32_t)llround(median_diff);
		if (isfinite(median_cov) && median_cov > 0.0)
			stats[rn].median_cov = median_cov;
	}

	ani_double_vecs_destroy(products, ref_n);
	ani_double_vecs_destroy(diffs, ref_n);
	ani_double_vecs_destroy(covs, ref_n);
	return stats;
}

static bool ani_readwise_reject_product_topfrac_ctx(uint32_t diff,
													double cov,
													double product_threshold)
{
	if (diff == 0)
		return false;
	if (!isfinite(cov) || cov <= 0.0)
		return false;
	if (!isfinite(product_threshold) || product_threshold < 0.0)
		return false;
	const double product = (double)diff * cov;
	return product > 0.0 && product >= product_threshold;
}

static bool ani_readwise_apply_product_topfrac_median_ctx(uint32_t *diff,
														  double *cov,
														  const ani_readwise_product_topfrac_median_t *stats,
														  bool plus_one)
{
	if (!diff || !cov || !stats)
		return false;
	if (*diff == 0 && !plus_one)
		return false;
	if (!isfinite(*cov) || *cov <= 0.0)
		return false;
	if (!isfinite(stats->product_threshold) || stats->product_threshold < 0.0)
		return false;
	if (!isfinite(stats->median_cov) || stats->median_cov <= 0.0)
		return false;
	const double product = ani_readwise_topfrac_product_score(*diff, *cov, plus_one);
	if (!(product > 0.0 && product >= stats->product_threshold))
		return false;
	*diff = stats->median_diff;
	*cov = stats->median_cov;
	return true;
}


static void ani_readwise_trace_emit(const ani_readwise_trace_t *trace,
									const char *read_name,
									uint64_t unit_id,
									uint64_t qctx,
									size_t ref_begin,
									uint32_t diff,
									uint32_t best_diff,
									size_t candidate_n,
									size_t selected_n,
									uint32_t cov_inc)
{
	if (!trace || !trace->fp)
		return;
#ifdef _OPENMP
#pragma omp critical(minco_readwise_trace)
#endif
	{
		fprintf(trace->fp,
				"%s\t%" PRIu64 "\t%" PRIu64 "\t%zu\t%u\t%u\t%zu\t%zu\t%u\t%s\n",
				read_name ? read_name : "",
				unit_id,
				qctx,
				ref_begin,
				diff,
				best_diff,
				candidate_n,
				selected_n,
				cov_inc,
				trace->ref_name ? trace->ref_name : "");
	}
}

static bool ani_readwise_candidate_selected_by_mode(ani_readwise_assign_mode_t assign_mode,
													uint32_t diff,
													uint32_t best_diff,
													size_t selected_n)
{
	if (assign_mode == ANI_READWISE_ASSIGN_ALL)
		return true;
	if (diff != best_diff)
		return false;
	if (assign_mode == ANI_READWISE_ASSIGN_BEST_DIFF_UNIQUE)
		return selected_n == 1;
	return selected_n > 0;
}

static void ani_readwise_edge_trace_emit_group(ani_readwise_edge_trace_t *edge,
											   const char *read_name,
											   uint64_t unit_id,
											   uint64_t qctx,
											   const kv_readwise_candidate_t *candidates,
											   uint32_t best_diff,
											   size_t selected_n,
											   uint32_t cov_inc,
											   ani_readwise_assign_mode_t assign_mode)
{
	if (!edge || !edge->fp || !candidates)
		return;
	const size_t candidate_n = kv_size(*candidates);
	if (!candidate_n)
		return;
	if (edge->ambiguous_only && candidate_n <= 1)
		return;
	if (edge->max_candidate_refs && candidate_n > edge->max_candidate_refs)
		return;
#ifdef _OPENMP
#pragma omp critical(minco_readwise_edge_trace)
#endif
	{
		for (size_t ci = 0; ci < candidate_n; ++ci) {
			if (edge->max_edges && edge->emitted_edges >= edge->max_edges)
				break;
			const ani_readwise_candidate_t *cand = &kv_A(*candidates, ci);
			const bool selected = ani_readwise_candidate_selected_by_mode(
				assign_mode, cand->diff, best_diff, selected_n);
			if (edge->selected_only && !selected)
				continue;
			++edge->emitted_edges;
			fprintf(edge->fp,
					"%" PRIu64 "\t%s\t%" PRIu64 "\t%" PRIu64 "\t%zu\t%zu\t%u\t%u\t%u\t%zu\t%zu\t%u\t%u\n",
					edge->emitted_edges,
					read_name ? read_name : "",
					unit_id,
					qctx,
					ci,
					cand->ref_begin,
					cand->gid,
					cand->diff,
					best_diff,
					candidate_n,
					selected_n,
					selected ? 1u : 0u,
					cov_inc);
		}
	}
}

static inline uint64_t ani_make_hashed_ctxobj(uint64_t unituple,
											  uint32_t n_obj_bits,
											  uint64_t density_threshold);
static inline __uint128_t ani_mask128_bits(unsigned bits);
static inline ctxobj96_t ani_coden128_to_ctxobj96(__uint128_t tuple,
												  int codens);

#define ANI_READWISE_TRACK_MAX_LIST_ITEMS 64u
#define ANI_READWISE_TRACK_MAX_OFFSETS 256u

static int ani_read_ctxobj64_cmp(const void *pa, const void *pb)
{
	const ani_read_ctxobj64_t *a = (const ani_read_ctxobj64_t *)pa;
	const ani_read_ctxobj64_t *b = (const ani_read_ctxobj64_t *)pb;
	if (a->packed != b->packed)
		return (a->packed > b->packed) - (a->packed < b->packed);
	return (a->offset > b->offset) - (a->offset < b->offset);
}

static int ani_read_ctxobj96_cmp(const void *pa, const void *pb)
{
	const ani_read_ctxobj96_t *a = (const ani_read_ctxobj96_t *)pa;
	const ani_read_ctxobj96_t *b = (const ani_read_ctxobj96_t *)pb;
	if (a->rec.ctx != b->rec.ctx)
		return (a->rec.ctx > b->rec.ctx) - (a->rec.ctx < b->rec.ctx);
	if (a->rec.obj != b->rec.obj)
		return (a->rec.obj > b->rec.obj) - (a->rec.obj < b->rec.obj);
	return (a->offset > b->offset) - (a->offset < b->offset);
}

static void ani_read_ctxobj64_push(kv_read_ctxobj64_t *vec,
								   uint64_t packed,
								   uint32_t offset)
{
	ani_read_ctxobj64_t rec = {.packed = packed, .offset = offset};
	kv_push(ani_read_ctxobj64_t, *vec, rec);
}

static void ani_read_ctxobj96_push(kv_read_ctxobj96_t *vec,
								   ctxobj96_t packed,
								   uint32_t offset)
{
	ani_read_ctxobj96_t rec = {.rec = packed, .offset = offset};
	kv_push(ani_read_ctxobj96_t, *vec, rec);
}

static void ani_extract_read_density_ctxobjs_with_offsets(
	const char *s,
	int len,
	kv_read_ctxobj64_t *vec,
	uint32_t n_obj_bits,
	uint64_t density_threshold)
{
	if (len < (int)klen)
		return;
	const uint32_t len_mv = (uint32_t)(2 * klen - 2);
	uint64_t tuple = 0, crv = 0;
	int base = 0;

	for (int pos = 0; pos < len; ++pos) {
		const int bmap = Basemap[(unsigned char)s[pos]];
		if (unlikely(bmap == DEFAULT)) {
			base = 0;
			tuple = 0;
			crv = 0;
			continue;
		}
		const uint64_t b2 = (uint64_t)bmap;
		tuple = (tuple << 2) | b2;
		crv = (crv >> 2) | ((b2 ^ 3ull) << len_mv);
		if (unlikely(++base < (int)klen))
			continue;

		const uint64_t t_ctx = tuple & ctxmask;
		const uint64_t r_ctx = crv & ctxmask;
		const uint64_t unictx = t_ctx < r_ctx ? t_ctx : r_ctx;
#if ANI_APPLY_SOURCE_FILTER
		if (unlikely((uint32_t)mix64(unictx) > FILTER))
			continue;
#endif
		const uint64_t unituple = (t_ctx < r_ctx ? tuple : crv) & tupmask;
		const uint64_t packed = ani_make_hashed_ctxobj(unituple, n_obj_bits,
													   density_threshold);
		if (packed != UINT64_MAX)
			ani_read_ctxobj64_push(vec, packed, (uint32_t)(pos + 1 - (int)klen));
	}
}

static void ani_extract_read_density_ctxobjs96_with_offsets(
	const char *s,
	int len,
	kv_read_ctxobj96_t *vec,
	uint64_t density_threshold)
{
	if (len < (int)klen)
		return;
	const unsigned tuple_bits = 2u * klen;
	const __uint128_t tuple_mask = ani_mask128_bits(tuple_bits);
	const unsigned rev_shift = 2u * (klen - 1u);
	__uint128_t tuple = 0, crv = 0;
	int base = 0;

	for (int pos = 0; pos < len; ++pos) {
		const int bmap = Basemap[(unsigned char)s[pos]];
		if (unlikely(bmap == DEFAULT)) {
			base = 0;
			tuple = 0;
			crv = 0;
			continue;
		}
		const __uint128_t b2 = (uint64_t)bmap;
		tuple = ((tuple << 2) | b2) & tuple_mask;
		crv = (crv >> 2) | ((b2 ^ 3u) << rev_shift);
		if (unlikely(++base < (int)klen))
			continue;

		const ctxobj96_t fwd = ani_coden128_to_ctxobj96(tuple, NUM_CODENS);
		const ctxobj96_t rev = ani_coden128_to_ctxobj96(crv, NUM_CODENS);
		const ctxobj96_t rec =
			(fwd.ctx < rev.ctx || (fwd.ctx == rev.ctx && fwd.obj <= rev.obj))
				? fwd
				: rev;
#if ANI_APPLY_SOURCE_FILTER
		if (unlikely((uint32_t)mix64(rec.ctx) > FILTER))
			continue;
#endif
		if (mix64(rec.ctx ^ (uint64_t)MINCO_SEED) > density_threshold)
			continue;
		ani_read_ctxobj96_push(vec, rec, (uint32_t)(pos + 1 - (int)klen));
	}
}

static int ani_min_diff_sections_tracked64_vs_ref_index(
	const ani_read_ctxobj64_t *qry,
	size_t qry_begin,
	size_t qry_end,
	const ctxgidobj_t *ref,
	size_t ref_begin,
	size_t ref_end,
	uint64_t objmask)
{
	int min_diff_sections = NUM_CODENS + 1;
	for (size_t qi = qry_begin; qi < qry_end; ++qi) {
		const uint32_t obj_q = (uint32_t)(qry[qi].packed & objmask);
		for (size_t ri = ref_begin; ri < ref_end; ++ri) {
			const uint32_t diff = obj_q ^ ref[ri].obj;
			if (diff == 0)
				return 0;
			const int d = dna_popcount(diff);
			if (d < min_diff_sections)
				min_diff_sections = d;
		}
	}
	return min_diff_sections;
}

static int ani_min_diff_sections_tracked96_vs_ref_index(
	const ani_read_ctxobj96_t *qry,
	size_t qry_begin,
	size_t qry_end,
	const ctxgidobj128_t *ref,
	size_t ref_begin,
	size_t ref_end)
{
	int min_diff_sections = NUM_CODENS + 1;
	for (size_t qi = qry_begin; qi < qry_end; ++qi) {
		const uint32_t obj_q = qry[qi].rec.obj;
		for (size_t ri = ref_begin; ri < ref_end; ++ri) {
			const uint32_t diff = obj_q ^ ref[ri].obj;
			if (diff == 0)
				return 0;
			const int d = dna_popcount(diff);
			if (d < min_diff_sections)
				min_diff_sections = d;
		}
	}
	return min_diff_sections;
}

static void ani_track_u32_push_unique(kv_u32_t *vec, uint32_t value)
{
	for (size_t i = 0; i < kv_size(*vec); ++i)
		if (kv_A(*vec, i) == value)
			return;
	kv_push(uint32_t, *vec, value);
}

static void ani_track_offset_push(kv_u32_t *vec, uint32_t value, bool *truncated)
{
	if (kv_size(*vec) >= ANI_READWISE_TRACK_MAX_OFFSETS) {
		if (truncated)
			*truncated = true;
		return;
	}
	kv_push(uint32_t, *vec, value);
}

static int ani_track_tax_record_depth(const ani_cami_tax_record_t *tax)
{
	if (!tax)
		return -1;
	const int rank_idx = ani_cami_rank_index(tax->rank);
	if (rank_idx >= 0)
		return rank_idx + 1;
	int depth = 0;
	const char *p = tax->taxpathsn;
	if (!p || p[0] == '\0' || strcmp(p, "NA") == 0)
		return 0;
	for (;;) {
		++depth;
		const char *bar = strchr(p, '|');
		if (!bar)
			break;
		p = bar + 1;
	}
	return depth;
}

static const ani_cami_tax_record_t *ani_track_best_tax_for_key(
	const kv_cami_tax_record_t *records,
	const char *key)
{
	size_t begin = 0;
	size_t end = 0;
	if (!ani_cami_find_tax_record_range(records, key, &begin, &end))
		return NULL;
	const ani_cami_tax_record_t *best = NULL;
	int best_depth = -1;
	for (size_t i = begin; i < end; ++i) {
		const ani_cami_tax_record_t *tax = &kv_A(*records, i);
		const int depth = ani_track_tax_record_depth(tax);
		if (!best || depth > best_depth) {
			best = tax;
			best_depth = depth;
		}
	}
	return best;
}

static const ani_cami_tax_record_t *ani_track_best_ref_tax_match(
	const kv_cami_tax_record_t *records,
	const char *ref_name,
	const char *annotation)
{
	const ani_cami_tax_record_t *hit = ani_track_best_tax_for_key(records, ref_name);
	if (hit)
		return hit;
	char buf[PATHLEN];
	ani_cami_copy_basename(ref_name, buf, sizeof(buf));
	hit = ani_track_best_tax_for_key(records, buf);
	if (hit)
		return hit;
	if (ani_cami_extract_accession(ref_name, buf, sizeof(buf))) {
		hit = ani_track_best_tax_for_key(records, buf);
		if (hit)
			return hit;
	}
	hit = ani_track_best_tax_for_key(records, annotation);
	if (hit)
		return hit;
	if (ani_cami_extract_accession(annotation, buf, sizeof(buf)))
		return ani_track_best_tax_for_key(records, buf);
	return NULL;
}

static void ani_readwise_tax_namespace_load(
	ani_readwise_tax_namespace_t *ns,
	const char *label,
	const char *path,
	char (*refname)[PATHLEN],
	char (*refanno)[PATHLEN],
	uint32_t ref_n)
{
	memset(ns, 0, sizeof(*ns));
	ns->enabled = true;
	ns->label = label;
	ani_cami_load_taxmap(path, &ns->records);
	ns->ref_tax = calloc((size_t)ref_n, sizeof(ns->ref_tax[0]));
	if (!ns->ref_tax)
		err(EXIT_FAILURE, "%s(): OOM %s taxonomy refs", __func__, label ? label : "readwise");
	size_t missing = 0;
	for (uint32_t rn = 0; rn < ref_n; ++rn) {
		ns->ref_tax[rn] = ani_track_best_ref_tax_match(
			&ns->records, refname[rn], annotation_at(refanno, rn));
		if (!ns->ref_tax[rn])
			++missing;
	}
	if (missing)
		warnx("%s(): %s taxonomy missing for %zu/%u reference entries",
			  __func__, label ? label : "readwise", missing, ref_n);
}

static void ani_readwise_tax_namespace_destroy(ani_readwise_tax_namespace_t *ns)
{
	if (!ns)
		return;
	free(ns->ref_tax);
	ns->ref_tax = NULL;
	if (ns->enabled)
		ani_cami_tax_records_destroy(&ns->records);
	memset(ns, 0, sizeof(*ns));
}

static size_t ani_track_common_taxpath_prefix_len(const char *a,
												  size_t a_len,
												  const char *b)
{
	if (!a || !b || !a_len)
		return 0;
	size_t pa = 0;
	size_t pb = 0;
	size_t last_end = 0;
	while (pa < a_len && b[pb]) {
		size_t ea = pa;
		while (ea < a_len && a[ea] != '|')
			++ea;
		size_t eb = pb;
		while (b[eb] && b[eb] != '|')
			++eb;
		const size_t la = ea - pa;
		const size_t lb = eb - pb;
		if (la != lb || strncmp(a + pa, b + pb, la) != 0)
			break;
		last_end = ea;
		if (ea >= a_len || b[eb] == '\0')
			break;
		pa = ea + 1u;
		pb = eb + 1u;
	}
	return last_end;
}

static int ani_track_taxpath_token_count(const char *path, size_t prefix_len)
{
	if (!path || prefix_len == 0)
		return 0;
	int count = 1;
	for (size_t i = 0; i < prefix_len; ++i)
		if (path[i] == '|')
			++count;
	return count;
}

static void ani_track_taxpath_last_token(const char *path,
										 size_t prefix_len,
										 char *out,
										 size_t out_size)
{
	if (!out || out_size == 0)
		return;
	out[0] = '\0';
	if (!path || prefix_len == 0) {
		snprintf(out, out_size, "root");
		return;
	}
	size_t begin = 0;
	for (size_t i = 0; i < prefix_len; ++i)
		if (path[i] == '|')
			begin = i + 1u;
	const size_t len = prefix_len > begin ? prefix_len - begin : 0;
	const size_t copy = len + 1u < out_size ? len : out_size - 1u;
	memcpy(out, path + begin, copy);
	out[copy] = '\0';
	if (out[0] == '\0')
		snprintf(out, out_size, "root");
}

static const char *ani_track_rank_for_depth(int depth)
{
	static const char *ranks[] = {
		"root", "superkingdom", "phylum", "class", "order",
		"family", "genus", "species", "strain"};
	if (depth < 0)
		depth = 0;
	if (depth >= (int)(sizeof(ranks) / sizeof(ranks[0])))
		depth = (int)(sizeof(ranks) / sizeof(ranks[0])) - 1;
	return ranks[depth];
}

static void ani_track_lca_namespace(
	const ani_readwise_tax_namespace_t *ns,
	const kv_u32_t *gids,
	char *rank_out,
	size_t rank_size,
	char *name_out,
	size_t name_size)
{
	if (rank_out && rank_size)
		snprintf(rank_out, rank_size, "NA");
	if (name_out && name_size)
		snprintf(name_out, name_size, "NA");
	if (!ns || !ns->enabled || !ns->ref_tax || !gids || kv_size(*gids) == 0)
		return;

	const ani_cami_tax_record_t *first = NULL;
	for (size_t i = 0; i < kv_size(*gids); ++i) {
		const uint32_t gid = kv_A(*gids, i);
		const ani_cami_tax_record_t *tax = ns->ref_tax[gid];
		if (tax && tax->taxpathsn && tax->taxpathsn[0] &&
			strcmp(tax->taxpathsn, "NA") != 0) {
			first = tax;
			break;
		}
	}
	if (!first)
		return;

	size_t common_len = strlen(first->taxpathsn);
	for (size_t i = 0; i < kv_size(*gids); ++i) {
		const uint32_t gid = kv_A(*gids, i);
		const ani_cami_tax_record_t *tax = ns->ref_tax[gid];
		if (!tax || !tax->taxpathsn || !tax->taxpathsn[0] ||
			strcmp(tax->taxpathsn, "NA") == 0)
			continue;
		common_len = ani_track_common_taxpath_prefix_len(
			first->taxpathsn, common_len, tax->taxpathsn);
		if (common_len == 0)
			break;
	}

	const int depth = ani_track_taxpath_token_count(first->taxpathsn, common_len);
	if (rank_out && rank_size)
		snprintf(rank_out, rank_size, "%s", ani_track_rank_for_depth(depth));
	if (name_out && name_size)
		ani_track_taxpath_last_token(first->taxpathsn, common_len, name_out, name_size);
}

static long double ani_readwise_density_probability(uint64_t density_threshold,
													bool ref_uses_ctxobj96,
													uint8_t nobjbits)
{
	if (density_threshold == UINT64_MAX)
		return 1.0L;
	if (ref_uses_ctxobj96) {
		return ((long double)density_threshold + 1.0L) / ldexpl(1.0L, 64);
	}
	const int hash_bits = 64 - (int)nobjbits;
	if (hash_bits <= 0)
		return 1.0L;
	return ((long double)density_threshold + 1.0L) / ldexpl(1.0L, hash_bits);
}

static long double ani_readwise_threshold_density(uint64_t threshold,
												  uint32_t hash_bits)
{
	if (hash_bits == 0 || hash_bits > 64)
		return 0.0L;
	if (hash_bits < 64) {
		const uint64_t full_threshold = (1ULL << hash_bits) - 1ULL;
		if (threshold >= full_threshold)
			return 1.0L;
	}
	if (hash_bits == 64 && threshold == UINT64_MAX)
		return 1.0L;
	long double density = ((long double)threshold + 1.0L) /
						  ldexpl(1.0L, (int)hash_bits);
	if (density < 0.0L)
		density = 0.0L;
	if (density > 1.0L)
		density = 1.0L;
	return density;
}

static bool ani_readwise_tracker_ref_density(const ani_readwise_tracker_t *tracker,
											 uint32_t gid,
											 long double *density_out)
{
	if (density_out)
		*density_out = 0.0L;
	if (!tracker || !tracker->ref_ctxmeta || gid >= tracker->ref_n)
		return false;
	const ani_ctxmeta_rec_t *meta = &tracker->ref_ctxmeta[gid];
	if (!meta->valid)
		return false;
	const long double density =
		ani_readwise_threshold_density(meta->threshold, meta->hash_bits);
	if (density <= 0.0L)
		return false;
	if (density_out)
		*density_out = density;
	return true;
}

static void ani_readwise_tracker_add_sketch_corrected_hit(
	ani_readwise_tracker_t *tracker,
	uint64_t occurrences,
	long double max_ref_density,
	bool have_ref_density)
{
	if (!tracker || !tracker->enabled || !tracker->sketch_corrected_available ||
		occurrences == 0)
		return;
	if (!have_ref_density || max_ref_density <= 0.0L ||
		tracker->density_probability <= 0.0L) {
		tracker->sketch_corrected_missing_meta_ctx += occurrences;
		return;
	}
	long double capture_probability =
		max_ref_density / tracker->density_probability;
	if (capture_probability > 1.0L)
		capture_probability = 1.0L;
	if (capture_probability <= 0.0L) {
		tracker->sketch_corrected_missing_meta_ctx += occurrences;
		return;
	}
	tracker->sketch_corrected_observed_ctx += occurrences;
	tracker->sketch_corrected_capture_prob_sum +=
		(long double)occurrences * capture_probability;
	tracker->sketch_corrected_present_ctx +=
		(long double)occurrences / capture_probability;
}

static long double ani_readwise_trackable_probability(size_t possible_ctx,
													  long double density_p)
{
	if (possible_ctx == 0 || density_p <= 0.0L)
		return 0.0L;
	if (density_p >= 1.0L)
		return 1.0L;
	return 1.0L - powl(1.0L - density_p, (long double)possible_ctx);
}

static void ani_track_fprint_text(FILE *fp, const char *s)
{
	if (!fp)
		return;
	if (!s || s[0] == '\0') {
		fputs("NA", fp);
		return;
	}
	for (const char *p = s; *p; ++p) {
		const unsigned char c = (unsigned char)*p;
		fputc((c == '\t' || c == '\n' || c == '\r') ? ' ' : (int)c, fp);
	}
}

static void ani_track_fprint_u32_list(FILE *fp, const kv_u32_t *values)
{
	if (!fp || !values || kv_size(*values) == 0) {
		fputs("NA", fp);
		return;
	}
	const size_t n = kv_size(*values);
	const size_t cap = n < ANI_READWISE_TRACK_MAX_LIST_ITEMS
						   ? n
						   : ANI_READWISE_TRACK_MAX_LIST_ITEMS;
	for (size_t i = 0; i < cap; ++i) {
		if (i)
			fputc(',', fp);
		fprintf(fp, "%u", kv_A(*values, i));
	}
	if (cap < n)
		fputs(",...", fp);
}

static void ani_track_fprint_ref_list(FILE *fp,
									  const kv_u32_t *gids,
									  char (*refname)[PATHLEN])
{
	if (!fp || !gids || kv_size(*gids) == 0) {
		fputs("NA", fp);
		return;
	}
	const size_t n = kv_size(*gids);
	const size_t cap = n < ANI_READWISE_TRACK_MAX_LIST_ITEMS
						   ? n
						   : ANI_READWISE_TRACK_MAX_LIST_ITEMS;
	for (size_t i = 0; i < cap; ++i) {
		if (i)
			fputc(',', fp);
		ani_track_fprint_text(fp, refname ? refname[kv_A(*gids, i)] : "NA");
	}
	if (cap < n)
		fputs(",...", fp);
}

static void ani_readwise_tracker_init(
	ani_readwise_tracker_t *tracker,
	const ani_opt_t *ani_opt,
	const char *query_path,
	char (*refname)[PATHLEN],
	char (*refanno)[PATHLEN],
	uint32_t ref_n,
	const ani_ctxmeta_rec_t *ref_ctxmeta,
	long double density_probability)
{
	memset(tracker, 0, sizeof(*tracker));
	if (!ani_opt || ani_opt->readwise_track[0] == '\0')
		return;
	tracker->enabled = true;
	tracker->query_path = query_path;
	tracker->refname = refname;
	tracker->refanno = refanno;
	tracker->ref_n = ref_n;
	tracker->ref_ctxmeta = ref_ctxmeta;
	tracker->taxonomy_mode = ani_opt->readwise_taxonomy_mode;
	tracker->density_probability = density_probability;
	tracker->sketch_corrected_available =
		ref_ctxmeta != NULL && density_probability > 0.0L;
	if (ani_opt->readwise_track_summary[0] != '\0') {
		snprintf(tracker->summary_path, sizeof(tracker->summary_path), "%s",
				 ani_opt->readwise_track_summary);
	} else if (strcmp(ani_opt->readwise_track, "-") == 0) {
		snprintf(tracker->summary_path, sizeof(tracker->summary_path),
				 "minco.readwise_track.summary.tsv");
	} else {
		snprintf(tracker->summary_path, sizeof(tracker->summary_path),
				 "%s.summary.tsv", ani_opt->readwise_track);
	}

	if (tracker->taxonomy_mode == ANI_READWISE_TAXONOMY_GTDB ||
		tracker->taxonomy_mode == ANI_READWISE_TAXONOMY_BOTH)
		ani_readwise_tax_namespace_load(&tracker->gtdb, "GTDB",
										ani_opt->gtdb_taxmap,
										refname, refanno, ref_n);
	if (tracker->taxonomy_mode == ANI_READWISE_TAXONOMY_NCBI ||
		tracker->taxonomy_mode == ANI_READWISE_TAXONOMY_BOTH)
		ani_readwise_tax_namespace_load(&tracker->ncbi, "NCBI",
										ani_opt->ncbi_taxmap,
										refname, refanno, ref_n);

	tracker->fp = strcmp(ani_opt->readwise_track, "-") == 0
					  ? stdout
					  : fopen(ani_opt->readwise_track, "w");
	if (!tracker->fp)
		err(errno, "%s(): cannot open --readwise-track %s", __func__,
			ani_opt->readwise_track);
	fputs("read_id\tread_ord\tread_len\tpossible_ctx\tdensity_ctx\tmatched_ctx\tselected_ctx\ttarget_ref_count\ttarget_refs\ttarget_ref_ids\tgtdb_rank\tgtdb_name\tncbi_rank\tncbi_name\tctx_offsets\tctx_offsets_truncated\tassignment_status\n",
		  tracker->fp);
	fprintf(stderr,
			"minco readwise: read tracking active; out=%s summary=%s taxonomy=%d density_p=%.6Lg\n",
			ani_opt->readwise_track, tracker->summary_path,
			(int)tracker->taxonomy_mode, tracker->density_probability);
}

static void ani_readwise_tracker_write_summary(const ani_readwise_tracker_t *tracker)
{
	if (!tracker || !tracker->enabled || tracker->summary_path[0] == '\0')
		return;
	FILE *fp = fopen(tracker->summary_path, "w");
	if (!fp)
		err(errno, "%s(): cannot open readwise tracking summary %s", __func__,
			tracker->summary_path);
	const long double total = tracker->total_reads ? (long double)tracker->total_reads : 1.0L;
	const long double tracked_pct = 100.0L * (long double)tracker->reads_with_ref_hit / total;
	const long double density_pct = 100.0L * (long double)tracker->reads_with_density_ctx / total;
	const long double density_nohit_pct =
		100.0L * (long double)tracker->reads_without_ref_hit_but_density_ctx / total;
	long double density_trackable_nohit_pct = 0.0L;
	if (tracker->expected_trackable_reads > 0.0L) {
		density_trackable_nohit_pct =
			100.0L * (long double)tracker->reads_without_ref_hit_but_density_ctx /
			tracker->expected_trackable_reads;
		if (density_trackable_nohit_pct < 0.0L)
			density_trackable_nohit_pct = 0.0L;
		if (density_trackable_nohit_pct > 100.0L)
			density_trackable_nohit_pct = 100.0L;
	}
	const long double density_ctx_total =
		tracker->total_density_ctx ? (long double)tracker->total_density_ctx : 1.0L;
	long double ref_present_ctx_pct =
		100.0L * (long double)tracker->total_matched_ctx / density_ctx_total;
	if (ref_present_ctx_pct < 0.0L)
		ref_present_ctx_pct = 0.0L;
	if (ref_present_ctx_pct > 100.0L)
		ref_present_ctx_pct = 100.0L;
	fputs("metric\tvalue\n", fp);
	fprintf(fp, "query\t%s\n", tracker->query_path ? tracker->query_path : "NA");
	fprintf(fp, "taxonomy_mode\t%d\n", (int)tracker->taxonomy_mode);
	fprintf(fp, "density_probability\t%.12Lg\n", tracker->density_probability);
	fprintf(fp, "total_reads\t%" PRIu64 "\n", tracker->total_reads);
	fprintf(fp, "reads_with_density_ctx\t%" PRIu64 "\n", tracker->reads_with_density_ctx);
	fprintf(fp, "reads_with_ref_hit\t%" PRIu64 "\n", tracker->reads_with_ref_hit);
	fprintf(fp, "reads_with_multi_ref_hit\t%" PRIu64 "\n", tracker->reads_with_multi_ref_hit);
	fprintf(fp, "reads_without_ref_hit_but_density_ctx\t%" PRIu64 "\n",
			tracker->reads_without_ref_hit_but_density_ctx);
	fprintf(fp, "tracked_read_pct\t%.10Lg\n", tracked_pct);
	fprintf(fp, "density_positive_pct\t%.10Lg\n", density_pct);
	fprintf(fp, "density_positive_no_ref_hit_read_pct\t%.10Lg\n", density_nohit_pct);
	fprintf(fp, "expected_trackable_reads\t%.10Lg\n", tracker->expected_trackable_reads);
	fprintf(fp, "density_trackable_no_ref_hit_read_pct\t%.10Lg\n",
			density_trackable_nohit_pct);
	fprintf(fp, "total_possible_ctx\t%" PRIu64 "\n", tracker->total_possible_ctx);
	fprintf(fp, "total_density_ctx\t%" PRIu64 "\n", tracker->total_density_ctx);
	fprintf(fp, "total_matched_ctx\t%" PRIu64 "\n", tracker->total_matched_ctx);
	fprintf(fp, "total_selected_ctx\t%" PRIu64 "\n", tracker->total_selected_ctx);
	fprintf(fp, "total_selected_ref_events\t%" PRIu64 "\n", tracker->total_selected_ref_events);
	fprintf(fp, "sampled_ctx_ref_hit_pct\t%.10Lg\n", ref_present_ctx_pct);
	fprintf(fp, "sketch_corrected_available\t%u\n",
			tracker->sketch_corrected_available ? 1u : 0u);
	fprintf(fp, "sketch_corrected_observed_ctx\t%" PRIu64 "\n",
			tracker->sketch_corrected_observed_ctx);
	fprintf(fp, "sketch_corrected_missing_meta_ctx\t%" PRIu64 "\n",
			tracker->sketch_corrected_missing_meta_ctx);
	if (tracker->sketch_corrected_available) {
		long double corrected_present = tracker->sketch_corrected_present_ctx;
		if (corrected_present < 0.0L)
			corrected_present = 0.0L;
		long double corrected_pct =
			100.0L * corrected_present / density_ctx_total;
		if (corrected_pct < 0.0L)
			corrected_pct = 0.0L;
		if (corrected_pct > 100.0L)
			corrected_pct = 100.0L;
		const long double corrected_absent_pct = 100.0L - corrected_pct;
		if (tracker->sketch_corrected_observed_ctx > 0) {
			const long double mean_capture_probability =
				tracker->sketch_corrected_capture_prob_sum /
				(long double)tracker->sketch_corrected_observed_ctx;
			fprintf(fp, "sketch_corrected_mean_capture_probability\t%.12Lg\n",
					mean_capture_probability);
		} else {
			fprintf(fp, "sketch_corrected_mean_capture_probability\tNA\n");
		}
		fprintf(fp, "sketch_corrected_estimated_ref_present_ctx\t%.10Lg\n",
				tracker->sketch_corrected_present_ctx);
		fprintf(fp, "sketch_corrected_ref_present_ctx_pct\t%.10Lg\n",
				corrected_pct);
		fprintf(fp, "estimated_unknown_reads_pct\t%.10Lg\n",
				corrected_absent_pct);
		fprintf(fp, "estimated_unknown_reads_basis\t%s\n",
				"sketch_corrected_ref_absent_contexts;ref_ctxmeta_density_horvitz_thompson;whole_genome_estimate_only_for_full_sketch_refdb;shared_ctx_uses_max_ref_density");
	} else {
		fprintf(fp, "sketch_corrected_mean_capture_probability\tNA\n");
		fprintf(fp, "sketch_corrected_estimated_ref_present_ctx\tNA\n");
		fprintf(fp, "sketch_corrected_ref_present_ctx_pct\tNA\n");
		fprintf(fp, "estimated_unknown_reads_pct\tNA\n");
		fprintf(fp, "estimated_unknown_reads_basis\t%s\n",
				"unavailable:no_ref_ctxmeta_or_query_density");
	}
	if (fclose(fp) != 0)
		err(errno, "%s(): cannot close readwise tracking summary %s", __func__,
			tracker->summary_path);
}

static void ani_readwise_tracker_destroy(ani_readwise_tracker_t *tracker)
{
	if (!tracker || !tracker->enabled)
		return;
	ani_readwise_tracker_write_summary(tracker);
	if (tracker->fp && tracker->fp != stdout)
		fclose(tracker->fp);
	tracker->fp = NULL;
	ani_readwise_tax_namespace_destroy(&tracker->gtdb);
	ani_readwise_tax_namespace_destroy(&tracker->ncbi);
	memset(tracker, 0, sizeof(*tracker));
}

static void ani_readwise_track_emit_row(
	ani_readwise_tracker_t *tracker,
	const char *read_name,
	uint64_t read_ord,
	int read_len,
	uint64_t possible_ctx,
	uint64_t density_ctx,
	uint64_t matched_ctx,
	uint64_t selected_ctx,
	const kv_u32_t *target_gids,
	const kv_u32_t *offsets,
	bool offsets_truncated)
{
	if (!tracker || !tracker->enabled || !tracker->fp || !target_gids ||
		kv_size(*target_gids) == 0)
		return;
	char gtdb_rank[64], gtdb_name[PATHLEN];
	char ncbi_rank[64], ncbi_name[PATHLEN];
	snprintf(gtdb_rank, sizeof(gtdb_rank), "NA");
	snprintf(gtdb_name, sizeof(gtdb_name), "NA");
	snprintf(ncbi_rank, sizeof(ncbi_rank), "NA");
	snprintf(ncbi_name, sizeof(ncbi_name), "NA");
	if (tracker->taxonomy_mode == ANI_READWISE_TAXONOMY_GTDB ||
		tracker->taxonomy_mode == ANI_READWISE_TAXONOMY_BOTH)
		ani_track_lca_namespace(&tracker->gtdb, target_gids,
								gtdb_rank, sizeof(gtdb_rank),
								gtdb_name, sizeof(gtdb_name));
	if (tracker->taxonomy_mode == ANI_READWISE_TAXONOMY_NCBI ||
		tracker->taxonomy_mode == ANI_READWISE_TAXONOMY_BOTH)
		ani_track_lca_namespace(&tracker->ncbi, target_gids,
								ncbi_rank, sizeof(ncbi_rank),
								ncbi_name, sizeof(ncbi_name));

	ani_track_fprint_text(tracker->fp, read_name);
	fprintf(tracker->fp,
			"\t%" PRIu64 "\t%d\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%zu\t",
			read_ord, read_len, possible_ctx, density_ctx, matched_ctx,
			selected_ctx, kv_size(*target_gids));
	ani_track_fprint_ref_list(tracker->fp, target_gids, tracker->refname);
	fputc('\t', tracker->fp);
	ani_track_fprint_u32_list(tracker->fp, target_gids);
	fputc('\t', tracker->fp);
	ani_track_fprint_text(tracker->fp, gtdb_rank);
	fputc('\t', tracker->fp);
	ani_track_fprint_text(tracker->fp, gtdb_name);
	fputc('\t', tracker->fp);
	ani_track_fprint_text(tracker->fp, ncbi_rank);
	fputc('\t', tracker->fp);
	ani_track_fprint_text(tracker->fp, ncbi_name);
	fputc('\t', tracker->fp);
	ani_track_fprint_u32_list(tracker->fp, offsets);
	fprintf(tracker->fp, "\t%u\t%s\n",
			offsets_truncated ? 1u : 0u,
			kv_size(*target_gids) > 1 ? "multi_ref_lca" : "single_ref");
}

static void ani_readwise_track_read64(
	ani_readwise_tracker_t *tracker,
	const char *read_name,
	uint64_t read_ord,
	const char *seq,
	int len,
	const ctxgidobj_t *index,
	size_t index_n,
	const size_t *fence,
	int fence_k,
	uint32_t ref_n,
	bool ignoreconflict,
	uint8_t nobjbits,
	uint64_t gidmask_local,
	uint64_t objmask,
	uint64_t density_threshold,
	ani_readwise_assign_mode_t assign_mode)
{
	if (!tracker || !tracker->enabled)
		return;
	tracker->total_reads++;
	const uint64_t possible_ctx = len >= (int)klen ? (uint64_t)(len - (int)klen + 1) : 0;
	tracker->total_possible_ctx += possible_ctx;
	tracker->expected_trackable_reads +=
		ani_readwise_trackable_probability((size_t)possible_ctx,
										   tracker->density_probability);

	kv_read_ctxobj64_t vec;
	kv_init(vec);
	ani_extract_read_density_ctxobjs_with_offsets(seq, len, &vec, nobjbits,
												  density_threshold);
	const uint64_t density_ctx = (uint64_t)kv_size(vec);
	if (density_ctx == 0) {
		kv_destroy(vec);
		return;
	}
	tracker->reads_with_density_ctx++;
	tracker->total_density_ctx += density_ctx;
	qsort(&kv_A(vec, 0), kv_size(vec), sizeof(kv_A(vec, 0)),
		  ani_read_ctxobj64_cmp);

	kv_readwise_candidate_t candidates;
	kv_u32_t gids;
	kv_u32_t offsets;
	kv_init(candidates);
	kv_init(gids);
	kv_init(offsets);
	bool offsets_truncated = false;
	uint64_t matched_ctx = 0;
	uint64_t selected_ctx = 0;
	uint64_t selected_ref_events = 0;

	for (size_t q = 0; q < kv_size(vec); ) {
		const uint64_t qctx = kv_A(vec, q).packed >> nobjbits;
		const size_t qbeg = q;
		do { ++q; } while (q < kv_size(vec) &&
							(kv_A(vec, q).packed >> nobjbits) == qctx);
		const size_t qend = q;

		kv_size(candidates) = 0;
		uint32_t best_diff = UINT32_MAX;
		long double max_ref_density = 0.0L;
		bool have_ref_density = false;
		size_t pos = lb_in_bucket_ctxgid(index, fence, fence_k, qctx);
		while (pos < index_n && (index[pos].ctxgid >> GID_NBITS) == qctx) {
			const uint64_t ctxgid = index[pos].ctxgid;
			const uint32_t gid = (uint32_t)(ctxgid & gidmask_local);
			const size_t ref_begin = pos;
			do { ++pos; } while (pos < index_n && index[pos].ctxgid == ctxgid);
			const size_t ref_end = pos;
			if (gid >= ref_n)
				continue;
			if (ignoreconflict && ref_end - ref_begin > 1)
				continue;
			long double ref_density = 0.0L;
			if (ani_readwise_tracker_ref_density(tracker, gid, &ref_density) &&
				(!have_ref_density || ref_density > max_ref_density)) {
				max_ref_density = ref_density;
				have_ref_density = true;
			}
			const int min_diff = ani_min_diff_sections_tracked64_vs_ref_index(
				&kv_A(vec, 0), qbeg, qend, index, ref_begin, ref_end, objmask);
			ani_readwise_candidate_t cand = {
				.ref_begin = ref_begin,
				.gid = gid,
				.diff = (uint32_t)min_diff,
			};
			kv_push(ani_readwise_candidate_t, candidates, cand);
			if ((uint32_t)min_diff < best_diff)
				best_diff = (uint32_t)min_diff;
		}
		if (kv_size(candidates) == 0)
			continue;
		matched_ctx += (uint64_t)(qend - qbeg);
		ani_readwise_tracker_add_sketch_corrected_hit(
			tracker, (uint64_t)(qend - qbeg), max_ref_density, have_ref_density);

		size_t selected_n = kv_size(candidates);
		if (assign_mode != ANI_READWISE_ASSIGN_ALL) {
			selected_n = 0;
			for (size_t ci = 0; ci < kv_size(candidates); ++ci)
				if (kv_A(candidates, ci).diff == best_diff)
					++selected_n;
			if (!selected_n)
				continue;
			if (assign_mode == ANI_READWISE_ASSIGN_BEST_DIFF_UNIQUE &&
				selected_n != 1)
				continue;
		}
		selected_ctx += (uint64_t)(qend - qbeg);
		for (size_t oi = qbeg; oi < qend; ++oi)
			ani_track_offset_push(&offsets, kv_A(vec, oi).offset,
								  &offsets_truncated);
		for (size_t ci = 0; ci < kv_size(candidates); ++ci) {
			const ani_readwise_candidate_t *cand = &kv_A(candidates, ci);
			if (assign_mode != ANI_READWISE_ASSIGN_ALL && cand->diff != best_diff)
				continue;
			ani_track_u32_push_unique(&gids, cand->gid);
			selected_ref_events += (uint64_t)(qend - qbeg);
		}
	}

	tracker->total_matched_ctx += matched_ctx;
	tracker->total_selected_ctx += selected_ctx;
	tracker->total_selected_ref_events += selected_ref_events;
	if (kv_size(gids) > 0) {
		tracker->reads_with_ref_hit++;
		if (kv_size(gids) > 1)
			tracker->reads_with_multi_ref_hit++;
		ani_readwise_track_emit_row(tracker, read_name, read_ord, len,
									possible_ctx, density_ctx, matched_ctx,
									selected_ctx, &gids, &offsets,
									offsets_truncated);
	} else {
		tracker->reads_without_ref_hit_but_density_ctx++;
	}
	kv_destroy(offsets);
	kv_destroy(gids);
	kv_destroy(candidates);
	kv_destroy(vec);
}

static void ani_readwise_track_read96(
	ani_readwise_tracker_t *tracker,
	const char *read_name,
	uint64_t read_ord,
	const char *seq,
	int len,
	const ctxgidobj128_t *index,
	size_t index_n,
	const size_t *fence,
	int fence_k,
	uint32_t ref_n,
	bool ignoreconflict,
	uint64_t density_threshold,
	ani_readwise_assign_mode_t assign_mode)
{
	if (!tracker || !tracker->enabled)
		return;
	tracker->total_reads++;
	const uint64_t possible_ctx = len >= (int)klen ? (uint64_t)(len - (int)klen + 1) : 0;
	tracker->total_possible_ctx += possible_ctx;
	tracker->expected_trackable_reads +=
		ani_readwise_trackable_probability((size_t)possible_ctx,
										   tracker->density_probability);

	kv_read_ctxobj96_t vec;
	kv_init(vec);
	ani_extract_read_density_ctxobjs96_with_offsets(seq, len, &vec,
													density_threshold);
	const uint64_t density_ctx = (uint64_t)kv_size(vec);
	if (density_ctx == 0) {
		kv_destroy(vec);
		return;
	}
	tracker->reads_with_density_ctx++;
	tracker->total_density_ctx += density_ctx;
	qsort(&kv_A(vec, 0), kv_size(vec), sizeof(kv_A(vec, 0)),
		  ani_read_ctxobj96_cmp);

	kv_readwise_candidate_t candidates;
	kv_u32_t gids;
	kv_u32_t offsets;
	kv_init(candidates);
	kv_init(gids);
	kv_init(offsets);
	bool offsets_truncated = false;
	uint64_t matched_ctx = 0;
	uint64_t selected_ctx = 0;
	uint64_t selected_ref_events = 0;

	for (size_t q = 0; q < kv_size(vec); ) {
		const uint64_t qctx = kv_A(vec, q).rec.ctx;
		const size_t qbeg = q;
		do { ++q; } while (q < kv_size(vec) && kv_A(vec, q).rec.ctx == qctx);
		const size_t qend = q;

		kv_size(candidates) = 0;
		uint32_t best_diff = UINT32_MAX;
		long double max_ref_density = 0.0L;
		bool have_ref_density = false;
		size_t pos = lb_in_bucket_ctxgid128(index, fence, fence_k, qctx);
		while (pos < index_n && index[pos].ctx == qctx) {
			const uint64_t ctx = index[pos].ctx;
			const uint32_t gid = index[pos].gid;
			const size_t ref_begin = pos;
			do { ++pos; } while (pos < index_n &&
								  index[pos].ctx == ctx &&
								  index[pos].gid == gid);
			const size_t ref_end = pos;
			if (gid >= ref_n)
				continue;
			if (ignoreconflict && ref_end - ref_begin > 1)
				continue;
			long double ref_density = 0.0L;
			if (ani_readwise_tracker_ref_density(tracker, gid, &ref_density) &&
				(!have_ref_density || ref_density > max_ref_density)) {
				max_ref_density = ref_density;
				have_ref_density = true;
			}
			const int min_diff = ani_min_diff_sections_tracked96_vs_ref_index(
				&kv_A(vec, 0), qbeg, qend, index, ref_begin, ref_end);
			ani_readwise_candidate_t cand = {
				.ref_begin = ref_begin,
				.gid = gid,
				.diff = (uint32_t)min_diff,
			};
			kv_push(ani_readwise_candidate_t, candidates, cand);
			if ((uint32_t)min_diff < best_diff)
				best_diff = (uint32_t)min_diff;
		}
		if (kv_size(candidates) == 0)
			continue;
		matched_ctx += (uint64_t)(qend - qbeg);
		ani_readwise_tracker_add_sketch_corrected_hit(
			tracker, (uint64_t)(qend - qbeg), max_ref_density, have_ref_density);

		size_t selected_n = kv_size(candidates);
		if (assign_mode != ANI_READWISE_ASSIGN_ALL) {
			selected_n = 0;
			for (size_t ci = 0; ci < kv_size(candidates); ++ci)
				if (kv_A(candidates, ci).diff == best_diff)
					++selected_n;
			if (!selected_n)
				continue;
			if (assign_mode == ANI_READWISE_ASSIGN_BEST_DIFF_UNIQUE &&
				selected_n != 1)
				continue;
		}
		selected_ctx += (uint64_t)(qend - qbeg);
		for (size_t oi = qbeg; oi < qend; ++oi)
			ani_track_offset_push(&offsets, kv_A(vec, oi).offset,
								  &offsets_truncated);
		for (size_t ci = 0; ci < kv_size(candidates); ++ci) {
			const ani_readwise_candidate_t *cand = &kv_A(candidates, ci);
			if (assign_mode != ANI_READWISE_ASSIGN_ALL && cand->diff != best_diff)
				continue;
			ani_track_u32_push_unique(&gids, cand->gid);
			selected_ref_events += (uint64_t)(qend - qbeg);
		}
	}

	tracker->total_matched_ctx += matched_ctx;
	tracker->total_selected_ctx += selected_ctx;
	tracker->total_selected_ref_events += selected_ref_events;
	if (kv_size(gids) > 0) {
		tracker->reads_with_ref_hit++;
		if (kv_size(gids) > 1)
			tracker->reads_with_multi_ref_hit++;
		ani_readwise_track_emit_row(tracker, read_name, read_ord, len,
									possible_ctx, density_ctx, matched_ctx,
									selected_ctx, &gids, &offsets,
									offsets_truncated);
	} else {
		tracker->reads_without_ref_hit_but_density_ctx++;
	}
	kv_destroy(offsets);
	kv_destroy(gids);
	kv_destroy(candidates);
	kv_destroy(vec);
}

static void ani_density_units_destroy(kv_density_unit_t *units)
{
	if (!units)
		return;
	for (size_t i = 0; i < kv_size(*units); ++i) {
		v_free(&kv_A(*units, i).vec);
		free(kv_A(*units, i).read_name);
	}
	kv_destroy(*units);
}

static void ani_density_units96_destroy(kv_density_unit96_t *units)
{
	if (!units)
		return;
	for (size_t i = 0; i < kv_size(*units); ++i) {
		ani_ctxobj96_vec_free(&kv_A(*units, i).vec);
		free(kv_A(*units, i).read_name);
	}
	kv_destroy(*units);
}

static void ani_readwise_seq_batch_destroy(ani_readwise_seq_rec_t *batch,
										   size_t n)
{
	if (!batch)
		return;
	for (size_t i = 0; i < n; ++i)
		free(batch[i].seq);
	for (size_t i = 0; i < n; ++i)
		free(batch[i].name);
}

static void ani_readwise_thread_state_init(ani_readwise_thread_state_t *st,
										   uint32_t ref_n,
										   size_t index_n,
										   bool need_ref_cov_hits,
										   bool track_query_sets)
{
	memset(st, 0, sizeof(*st));
	st->acc = calloc((size_t)ref_n, sizeof(st->acc[0]));
	st->read_marks = calloc((size_t)ref_n, sizeof(st->read_marks[0]));
	st->ref_hit_bits = calloc((index_n + 7u) / 8u, 1);
	if (!st->acc || !st->read_marks || !st->ref_hit_bits)
		err(EXIT_FAILURE, "%s(): OOM readwise thread state", __func__);
	if (track_query_sets) {
		ani_u64_set_init(&st->qry_ctx_seen, 1u << 14);
		ani_u64_set_init(&st->qry_ref_ctx_seen, 1u << 15);
	}
	kv_init(st->ref_ctx_hit_positions);
	if (need_ref_cov_hits)
		kv_init(st->ref_ctx_cov_hits);
}

static void ani_readwise_thread_state_destroy(ani_readwise_thread_state_t *st)
{
	if (!st)
		return;
	free(st->acc);
	free(st->read_marks);
	free(st->ref_hit_bits);
	ani_u64_set_destroy(&st->qry_ctx_seen);
	ani_u64_set_destroy(&st->qry_ref_ctx_seen);
	kv_destroy(st->ref_ctx_hit_positions);
	kv_destroy(st->ref_ctx_cov_hits);
	memset(st, 0, sizeof(*st));
}

static void ani_readwise_thread_state_clear_batch(ani_readwise_thread_state_t *st,
												  uint32_t ref_n)
{
	if (!st)
		return;
	if (st->acc)
		memset(st->acc, 0, (size_t)ref_n * sizeof(st->acc[0]));
	if (st->ref_hit_bits) {
		for (size_t i = 0; i < kv_size(st->ref_ctx_hit_positions); ++i) {
			const size_t idx = kv_A(st->ref_ctx_hit_positions, i);
			st->ref_hit_bits[idx >> 3] &= (uint8_t)~(1u << (idx & 7u));
		}
	}
	kv_size(st->ref_ctx_hit_positions) = 0;
	kv_size(st->ref_ctx_cov_hits) = 0;
	ani_u64_set_clear(&st->qry_ctx_seen);
	ani_u64_set_clear(&st->qry_ref_ctx_seen);
}

static void ani_u64_set_merge(ani_u64_set_t *dst,
							  const ani_u64_set_t *src)
{
	if (!dst || !src)
		return;
	for (size_t i = 0; i < src->cap; ++i)
		if (src->used[i])
			(void)ani_u64_set_insert(dst, src->keys[i]);
}

static void ani_qry_ref_set_merge_into_acc(ani_u64_set_t *dst,
										   const ani_u64_set_t *src,
										   ani_readwise_acc_t *acc,
										   uint64_t gidmask_local,
										   uint32_t ref_n)
{
	if (!dst || !src || !acc)
		return;
	for (size_t i = 0; i < src->cap; ++i) {
		if (!src->used[i])
			continue;
		const uint64_t key = src->keys[i];
		if (!ani_u64_set_insert(dst, key))
			continue;
		const uint32_t gid = (uint32_t)(key & gidmask_local);
		if (gid < ref_n)
			acc[gid].qry_ctx_hit++;
	}
}

static void ani_merge_readwise_thread_batch(
	ani_readwise_thread_state_t *states,
	int n_states,
	ani_readwise_acc_t *acc,
	uint8_t *ref_hit_bits,
	uint32_t *ref_ctx_cov,
	uint32_t *ref_ctx_mindiff_cov,
	uint16_t *ref_ctx_hit_weight,
	ani_u64_set_t *qry_ctx_seen,
	ani_u64_set_t *qry_ref_ctx_seen,
	const ctxgidobj_t *index,
	size_t index_n,
	uint32_t ref_n,
	uint64_t gidmask_local)
{
	if (!states || n_states < 1)
		return;
	for (int t = 0; t < n_states; ++t) {
		ani_readwise_thread_state_t *st = &states[t];
		for (uint32_t rn = 0; rn < ref_n; ++rn) {
			acc[rn].XnY_ctx += st->acc[rn].XnY_ctx;
			acc[rn].N_diff_obj += st->acc[rn].N_diff_obj;
			acc[rn].N_diff_obj_section += st->acc[rn].N_diff_obj_section;
			acc[rn].N_mut2_ctx += st->acc[rn].N_mut2_ctx;
			acc[rn].reads_with_ctx_match += st->acc[rn].reads_with_ctx_match;
			acc[rn].blocks_with_ctx_match += st->acc[rn].blocks_with_ctx_match;
		}
		if (qry_ctx_seen)
			ani_u64_set_merge(qry_ctx_seen, &st->qry_ctx_seen);
		if (qry_ref_ctx_seen)
			ani_qry_ref_set_merge_into_acc(qry_ref_ctx_seen, &st->qry_ref_ctx_seen,
											acc, gidmask_local, ref_n);
		for (size_t i = 0; i < kv_size(st->ref_ctx_hit_positions); ++i) {
			const size_t idx = kv_A(st->ref_ctx_hit_positions, i);
			if (idx >= index_n)
				continue;
			if (!ani_bitset_test_set(ref_hit_bits, idx)) {
				const uint32_t gid = (uint32_t)(index[idx].ctxgid & gidmask_local);
				if (gid < ref_n)
					acc[gid].ref_ctx_hit++;
			}
		}
		if (ref_ctx_cov) {
			for (size_t i = 0; i < kv_size(st->ref_ctx_cov_hits); ++i) {
				const uint64_t hit = kv_A(st->ref_ctx_cov_hits, i);
				const size_t idx = (size_t)(hit >> ANI_REFCOV_HIT_SHIFT);
				const uint32_t inc = (uint32_t)((hit >> 4) & ANI_REFCOV_HIT_INC_MASK);
				const uint32_t diff = (uint32_t)(hit & 0xfu);
				if (idx < index_n) {
					ani_ref_covdiff_add_scaled_hit(&ref_ctx_cov[idx], diff, inc);
					if (ref_ctx_mindiff_cov)
						ani_ref_covmindiff_add_scaled_hit(&ref_ctx_mindiff_cov[idx], diff, inc);
					if (ref_ctx_hit_weight && inc > ref_ctx_hit_weight[idx])
						ref_ctx_hit_weight[idx] = (uint16_t)inc;
				}
			}
		}
		ani_readwise_thread_state_clear_batch(st, ref_n);
	}
}

static void ani_merge_readwise_thread_batch128(
	ani_readwise_thread_state_t *states,
	int n_states,
	ani_readwise_acc_t *acc,
	uint8_t *ref_hit_bits,
	uint32_t *ref_ctx_cov,
	uint32_t *ref_ctx_mindiff_cov,
	uint16_t *ref_ctx_hit_weight,
	const ctxgidobj128_t *index,
	size_t index_n,
	uint32_t ref_n)
{
	if (!states || n_states < 1)
		return;
	for (int t = 0; t < n_states; ++t) {
		ani_readwise_thread_state_t *st = &states[t];
		for (uint32_t rn = 0; rn < ref_n; ++rn) {
			acc[rn].XnY_ctx += st->acc[rn].XnY_ctx;
			acc[rn].N_diff_obj += st->acc[rn].N_diff_obj;
			acc[rn].N_diff_obj_section += st->acc[rn].N_diff_obj_section;
			acc[rn].N_mut2_ctx += st->acc[rn].N_mut2_ctx;
			acc[rn].reads_with_ctx_match += st->acc[rn].reads_with_ctx_match;
			acc[rn].blocks_with_ctx_match += st->acc[rn].blocks_with_ctx_match;
		}
		for (size_t i = 0; i < kv_size(st->ref_ctx_hit_positions); ++i) {
			const size_t idx = kv_A(st->ref_ctx_hit_positions, i);
			if (idx >= index_n)
				continue;
			if (!ani_bitset_test_set(ref_hit_bits, idx)) {
				const uint32_t gid = index[idx].gid;
				if (gid < ref_n)
					acc[gid].ref_ctx_hit++;
			}
		}
		if (ref_ctx_cov) {
			for (size_t i = 0; i < kv_size(st->ref_ctx_cov_hits); ++i) {
				const uint64_t hit = kv_A(st->ref_ctx_cov_hits, i);
				const size_t idx = (size_t)(hit >> ANI_REFCOV_HIT_SHIFT);
				const uint32_t inc = (uint32_t)((hit >> 4) & ANI_REFCOV_HIT_INC_MASK);
				const uint32_t diff = (uint32_t)(hit & 0xfu);
				if (idx < index_n) {
					ani_ref_covdiff_add_scaled_hit(&ref_ctx_cov[idx], diff, inc);
					if (ref_ctx_mindiff_cov)
						ani_ref_covmindiff_add_scaled_hit(&ref_ctx_mindiff_cov[idx], diff, inc);
					if (ref_ctx_hit_weight && inc > ref_ctx_hit_weight[idx])
						ref_ctx_hit_weight[idx] = (uint16_t)inc;
				}
			}
		}
		ani_readwise_thread_state_clear_batch(st, ref_n);
	}
}

static char *ani_shell_quote_arg(const char *arg)
{
	size_t len = 2;
	for (const char *p = arg; *p; ++p)
		len += (*p == '\'') ? 4 : 1;
	char *out = malloc(len + 1);
	if (!out)
		err(errno, "%s(): OOM shell quote", __func__);
	char *w = out;
	*w++ = '\'';
	for (const char *p = arg; *p; ++p) {
		if (*p == '\'') {
			memcpy(w, "'\\''", 4);
			w += 4;
		} else {
			*w++ = *p;
		}
	}
	*w++ = '\'';
	*w = '\0';
	return out;
}

static char *ani_pipe_command_for_path(const char *pipecmd, const char *path)
{
	char *quoted = ani_shell_quote_arg(path);
	const char *placeholder = strstr(pipecmd, "{}");
	if (!placeholder) {
		char *cmd = format_string("%s %s", pipecmd, quoted);
		free(quoted);
		return cmd;
	}

	size_t placeholders = 0;
	for (const char *p = pipecmd; (p = strstr(p, "{}")) != NULL; p += 2)
		placeholders++;
	const size_t pipecmd_len = strlen(pipecmd);
	const size_t quoted_len = strlen(quoted);
	const size_t out_len = pipecmd_len - placeholders * 2 + placeholders * quoted_len;
	char *cmd = malloc(out_len + 1);
	if (!cmd)
		err(errno, "%s(): OOM pipe command", __func__);

	const char *src = pipecmd;
	char *dst = cmd;
	while ((placeholder = strstr(src, "{}")) != NULL) {
		const size_t chunk = (size_t)(placeholder - src);
		memcpy(dst, src, chunk);
		dst += chunk;
		memcpy(dst, quoted, quoted_len);
		dst += quoted_len;
		src = placeholder + 2;
	}
	strcpy(dst, src);
	free(quoted);
	return cmd;
}

static ani_fastx_stream_t ani_open_fastx_stream(const char *path, const char *pipecmd)
{
	ani_fastx_stream_t stream = {0};
	if (pipecmd && pipecmd[0]) {
		stream.cmd = ani_pipe_command_for_path(pipecmd, path);
		stream.pipe_fp = popen(stream.cmd, "r");
		if (!stream.pipe_fp)
			err(errno, "%s(): popen %s", __func__, stream.cmd);
		int fd = dup(fileno(stream.pipe_fp));
		if (fd < 0)
			err(errno, "%s(): dup pipe fd", __func__);
		stream.gz = gzdopen(fd, "rb");
		if (!stream.gz)
			err(errno, "%s(): gzdopen pipe %s", __func__, stream.cmd);
		return stream;
	}
	if (strcmp(path, "-") == 0) {
		int fd = dup(STDIN_FILENO);
		if (fd < 0)
			err(errno, "%s(): dup stdin", __func__);
		stream.gz = gzdopen(fd, "rb");
		if (!stream.gz)
			err(errno, "%s(): gzdopen stdin", __func__);
		return stream;
	}
	stream.gz = gzopen(path, "r");
	if (!stream.gz)
		err(errno, "%s(): Cannot open %s", __func__, path);
	struct stat st;
	if (stat(path, &st) == 0 && S_ISREG(st.st_mode) && st.st_size > 0) {
		stream.input_size = (uint64_t)st.st_size;
		stream.input_size_known = true;
	}
	return stream;
}

static bool ani_close_fastx_stream(ani_fastx_stream_t *stream, bool warn_on_gzip_close)
{
	if (!stream)
		return false;
	bool had_warning = false;
	if (stream->gz) {
		const int rc = gzclose(stream->gz);
		stream->gz = NULL;
		if (rc != Z_OK) {
			if (!warn_on_gzip_close)
				errx(EXIT_FAILURE, "%s(): gzclose failed", __func__);
			warnx("%s(): gzclose failed after sequence parsing completed; "
				  "input gzip may be truncated/corrupt, keeping completed readwise result",
				  __func__);
			had_warning = true;
		}
	}
	if (stream->pipe_fp) {
		const int rc = pclose(stream->pipe_fp);
		stream->pipe_fp = NULL;
		if (rc == -1)
			err(errno, "%s(): pclose %s", __func__, stream->cmd ? stream->cmd : "pipe");
		if (rc != 0)
			errx(EXIT_FAILURE, "%s(): pipe command failed with status %d: %s",
				 __func__, rc, stream->cmd ? stream->cmd : "pipe");
	}
	free(stream->cmd);
	stream->cmd = NULL;
	return had_warning;
}

static void ani_progress_format_duration(uint64_t seconds, char *buf, size_t buf_size)
{
	const uint64_t hours = seconds / 3600;
	const uint64_t minutes = (seconds % 3600) / 60;
	const uint64_t secs = seconds % 60;
	if (hours > 9999) {
		snprintf(buf, buf_size, ">9999h");
	} else if (hours > 0) {
		snprintf(buf, buf_size, "%02" PRIu64 ":%02" PRIu64 ":%02" PRIu64,
				 hours, minutes, secs);
	} else {
		snprintf(buf, buf_size, "%02" PRIu64 ":%02" PRIu64, minutes, secs);
	}
}

static bool ani_fastx_stream_offset(const ani_fastx_stream_t *stream, uint64_t *offset)
{
	if (!stream || !stream->gz || !offset)
		return false;
#if defined(ZLIB_VERNUM) && ZLIB_VERNUM >= 0x1240
	const z_off_t pos = gzoffset(stream->gz);
	if (pos < 0)
		return false;
	*offset = (uint64_t)pos;
	return true;
#else
	(void)stream;
	(void)offset;
	return false;
#endif
}

static ani_readwise_progress_t ani_readwise_progress_start(const ani_fastx_stream_t *stream)
{
	ani_readwise_progress_t progress = {
		.enabled = true,
		.input_size_known = stream && stream->input_size_known,
		.input_size = stream ? stream->input_size : 0,
		.next_reads = ANI_READWISE_PROGRESS_READ_INTERVAL,
		.started_at = time(NULL),
		.last_at = 0,
	};
	if (progress.input_size_known) {
		fprintf(stderr,
				"minco readwise: started; input=%" PRIu64
				" bytes; progress every %" PRIu64 " reads or %d seconds\n",
				progress.input_size, (uint64_t)ANI_READWISE_PROGRESS_READ_INTERVAL,
				ANI_READWISE_PROGRESS_TIME_INTERVAL);
	} else {
		fprintf(stderr,
				"minco readwise: started; input size unknown; progress every %" PRIu64
				" reads or %d seconds\n",
				(uint64_t)ANI_READWISE_PROGRESS_READ_INTERVAL,
				ANI_READWISE_PROGRESS_TIME_INTERVAL);
	}
	return progress;
}

static void ani_readwise_progress_update(const ani_fastx_stream_t *stream,
										 ani_readwise_progress_t *progress,
										 uint64_t reads_done, bool force)
{
	if (!progress || !progress->enabled)
		return;
	if (!force && reads_done < progress->next_reads &&
		(reads_done & 0xffffULL) != 0)
		return;

	const time_t now = time(NULL);
	if (!force && reads_done < progress->next_reads &&
		now - progress->last_at < ANI_READWISE_PROGRESS_TIME_INTERVAL)
		return;
	progress->last_at = now;
	while (progress->next_reads <= reads_done)
		progress->next_reads += ANI_READWISE_PROGRESS_READ_INTERVAL;

	const uint64_t elapsed = now >= progress->started_at
								 ? (uint64_t)(now - progress->started_at)
								 : 0;
	const double rate = (double)reads_done / (double)(elapsed > 0 ? elapsed : 1);
	char elapsed_buf[32];
	ani_progress_format_duration(elapsed, elapsed_buf, sizeof(elapsed_buf));

	uint64_t offset = 0;
	if (progress->input_size_known && ani_fastx_stream_offset(stream, &offset)) {
		if (offset > progress->input_size)
			offset = progress->input_size;
		const double pct = progress->input_size > 0
							   ? 100.0 * (double)offset / (double)progress->input_size
							   : 100.0;
		fprintf(stderr,
				"minco readwise: %" PRIu64 " reads processed; %.2f%% input; "
				"elapsed %s; %.0f reads/s\n",
				reads_done, pct, elapsed_buf, rate);
	} else {
		fprintf(stderr,
				"minco readwise: %" PRIu64 " reads processed; elapsed %s; "
				"%.0f reads/s\n",
				reads_done, elapsed_buf, rate);
	}
}

static void ani_readwise_progress_done(const ani_fastx_stream_t *stream,
									   ani_readwise_progress_t *progress,
									   uint64_t reads_done)
{
	ani_readwise_progress_update(stream, progress, reads_done, true);
	if (progress && progress->enabled)
		fprintf(stderr, "minco readwise: complete\n");
}

static inline uint64_t ani_pext_portable(uint64_t value, uint64_t mask)
{
	uint64_t out = 0;
	uint64_t bit = 1;
	while (mask) {
		const uint64_t low = mask & (~mask + 1);
		if (value & low)
			out |= bit;
		mask ^= low;
		bit <<= 1;
	}
	return out;
}

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
__attribute__((target("bmi2")))
static inline uint64_t ani_pext_bmi2(uint64_t value, uint64_t mask)
{
	return _pext_u64(value, mask);
}
#endif

static inline uint64_t ani_pext_u64(uint64_t value, uint64_t mask)
{
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
	static int inited = 0;
	static int has_bmi2 = 0;
	if (!inited) {
		has_bmi2 = __builtin_cpu_supports("bmi2");
		__atomic_store_n(&inited, 1, __ATOMIC_RELAXED);
	}
	if (has_bmi2)
		return ani_pext_bmi2(value, mask);
#endif
	return ani_pext_portable(value, mask);
}

static inline uint64_t ani_hash_ctx(uint64_t ctx, uint32_t n_obj_bits)
{
	return mix64(ctx ^ (uint64_t)MINCO_SEED) >> n_obj_bits;
}

static inline uint64_t ani_hash_sparse_ctx(uint64_t sparse_ctx, uint32_t n_obj_bits)
{
	return mix64(sparse_ctx ^ (uint64_t)MINCO_SEED) >> n_obj_bits;
}

static inline uint64_t ani_make_hashed_ctxobj(uint64_t unituple, uint32_t n_obj_bits,
											  uint64_t density_threshold)
{
#if MINCO_HASH_SPARSE_CTX
	const uint64_t hctx = ani_hash_sparse_ctx(unituple & ctxmask, n_obj_bits);
#else
	const uint64_t ctx = ani_pext_u64(unituple, ctxmask);
	const uint64_t hctx = ani_hash_ctx(ctx, n_obj_bits);
#endif
	if (hctx > density_threshold)
		return UINT64_MAX;
	const uint64_t obj = ani_pext_u64(unituple, tupmask & ~ctxmask);
	return (hctx << n_obj_bits) | obj;
}

static void ani_extract_read_density_ctxobjs(const char *s, int len, u64vec *vec,
											 uint32_t n_obj_bits,
											 uint64_t density_threshold)
{
	if (len < (int)klen)
		return;
	const uint32_t len_mv = (uint32_t)(2 * klen - 2);
	uint64_t tuple = 0, crv = 0;
	int base = 0;

	for (int pos = 0; pos < len; ++pos) {
		const int bmap = Basemap[(unsigned char)s[pos]];
		if (unlikely(bmap == DEFAULT)) {
			base = 0;
			tuple = 0;
			crv = 0;
			continue;
		}
		const uint64_t b2 = (uint64_t)bmap;
		tuple = (tuple << 2) | b2;
		crv = (crv >> 2) | ((b2 ^ 3ull) << len_mv);
		if (unlikely(++base < (int)klen))
			continue;

		const uint64_t t_ctx = tuple & ctxmask;
		const uint64_t r_ctx = crv & ctxmask;
		const uint64_t unictx = t_ctx < r_ctx ? t_ctx : r_ctx;
#if ANI_APPLY_SOURCE_FILTER
		if (unlikely((uint32_t)mix64(unictx) > FILTER))
			continue;
#endif
		const uint64_t unituple = (t_ctx < r_ctx ? tuple : crv) & tupmask;
		const uint64_t packed = ani_make_hashed_ctxobj(unituple, n_obj_bits,
													   density_threshold);
		if (packed != UINT64_MAX)
			v_push(vec, packed);
	}
}

static void ani_ctxobj96_vec_init(ani_ctxobj96_vec_t *v, size_t cap)
{
	v->n = 0;
	v->m = cap;
	v->a = cap ? (ctxobj96_t *)malloc(cap * sizeof(v->a[0])) : NULL;
	if (cap && !v->a)
		err(EXIT_FAILURE, "%s(): OOM ctxobj96 read vector", __func__);
}

static void ani_ctxobj96_vec_free(ani_ctxobj96_vec_t *v)
{
	if (!v)
		return;
	free(v->a);
	v->a = NULL;
	v->n = v->m = 0;
}

static void ani_ctxobj96_vec_reserve(ani_ctxobj96_vec_t *v, size_t need)
{
	if (need <= v->m)
		return;
	size_t cap = v->m ? v->m : 128u;
	while (cap < need)
		cap <<= 1;
	ctxobj96_t *next = (ctxobj96_t *)realloc(v->a, cap * sizeof(v->a[0]));
	if (!next)
		err(EXIT_FAILURE, "%s(): OOM ctxobj96 read vector", __func__);
	v->a = next;
	v->m = cap;
}

static inline void ani_ctxobj96_vec_push(ani_ctxobj96_vec_t *v, ctxobj96_t rec)
{
	if (v->n == v->m)
		ani_ctxobj96_vec_reserve(v, v->m ? (v->m << 1) : 128u);
	v->a[v->n++] = rec;
}

#define ANI_DENSITY_CACHE_MAGIC "MNCDENS1"
#define ANI_DENSITY_CACHE_VERSION 1u
#define ANI_DENSITY_CACHE_STORAGE64 1u
#define ANI_DENSITY_CACHE_STORAGE96 2u

typedef struct {
	char magic[8];
	uint32_t version;
	uint32_t storage;
	uint32_t item_size;
	uint32_t nobjbits;
	uint32_t ctx_bits;
	uint32_t obj_bits;
	uint32_t reserved;
	uint64_t density_threshold;
	uint64_t total_reads;
	uint64_t density_reads;
} ani_density_cache_header_t;

typedef struct {
	FILE *fp;
	const char *path;
	ani_density_cache_header_t hdr;
} ani_density_cache_writer_t;

typedef struct {
	FILE *fp;
	const char *path;
	ani_density_cache_header_t hdr;
	uint64_t records_read;
} ani_density_cache_reader_t;

static void ani_density_cache_fwrite(FILE *fp, const void *ptr, size_t size,
									 size_t n, const char *path)
{
	if (n && fwrite(ptr, size, n, fp) != n)
		err(errno, "%s(): write failed: %s", __func__, path ? path : "density cache");
}

static bool ani_density_cache_fread_record(FILE *fp, void *ptr, size_t size,
										   size_t n, const char *path,
										   bool allow_eof)
{
	const size_t got = fread(ptr, size, n, fp);
	if (got == n)
		return true;
	if (allow_eof && got == 0 && feof(fp))
		return false;
	if (ferror(fp))
		err(errno, "%s(): read failed: %s", __func__, path ? path : "density cache");
	errx(EXIT_FAILURE, "%s(): truncated density cache: %s",
		 __func__, path ? path : "density cache");
}

static void ani_density_cache_writer_open(ani_density_cache_writer_t *w,
										  const char *path,
										  bool use_ctxobj96,
										  uint8_t nobjbits,
										  uint64_t density_threshold)
{
	if (!w || !path || path[0] == '\0')
		return;
	memset(w, 0, sizeof(*w));
	w->path = path;
	w->fp = fopen(path, "wb+");
	if (!w->fp)
		err(errno, "%s(): cannot create %s", __func__, path);
	memcpy(w->hdr.magic, ANI_DENSITY_CACHE_MAGIC, sizeof(w->hdr.magic));
	w->hdr.version = ANI_DENSITY_CACHE_VERSION;
	w->hdr.storage = use_ctxobj96 ? ANI_DENSITY_CACHE_STORAGE96 : ANI_DENSITY_CACHE_STORAGE64;
	w->hdr.item_size = use_ctxobj96 ? (uint32_t)sizeof(ctxobj96_t) : (uint32_t)sizeof(uint64_t);
	w->hdr.nobjbits = nobjbits;
	w->hdr.ctx_bits = Bitslen.ctx;
	w->hdr.obj_bits = Bitslen.obj;
	w->hdr.density_threshold = density_threshold;
	ani_density_cache_fwrite(w->fp, &w->hdr, sizeof(w->hdr), 1, path);
	fprintf(stderr, "minco readwise: density cache write active; out=%s\n", path);
}

static void ani_density_cache_writer_write64(ani_density_cache_writer_t *w,
											 const u64vec *vec)
{
	if (!w || !w->fp || !vec || vec->n == 0)
		return;
	if (vec->n > UINT32_MAX)
		errx(EXIT_FAILURE, "%s(): density cache record too large", __func__);
	const uint32_t n = (uint32_t)vec->n;
	ani_density_cache_fwrite(w->fp, &n, sizeof(n), 1, w->path);
	ani_density_cache_fwrite(w->fp, vec->a, sizeof(vec->a[0]), vec->n, w->path);
	w->hdr.density_reads++;
}

static void ani_density_cache_writer_write96(ani_density_cache_writer_t *w,
											 const ani_ctxobj96_vec_t *vec)
{
	if (!w || !w->fp || !vec || vec->n == 0)
		return;
	if (vec->n > UINT32_MAX)
		errx(EXIT_FAILURE, "%s(): density cache record too large", __func__);
	const uint32_t n = (uint32_t)vec->n;
	ani_density_cache_fwrite(w->fp, &n, sizeof(n), 1, w->path);
	ani_density_cache_fwrite(w->fp, vec->a, sizeof(vec->a[0]), vec->n, w->path);
	w->hdr.density_reads++;
}

static void ani_density_cache_writer_close(ani_density_cache_writer_t *w,
										   uint64_t total_reads)
{
	if (!w || !w->fp)
		return;
	w->hdr.total_reads = total_reads;
	if (fseeko(w->fp, 0, SEEK_SET) != 0)
		err(errno, "%s(): seek failed: %s", __func__, w->path);
	ani_density_cache_fwrite(w->fp, &w->hdr, sizeof(w->hdr), 1, w->path);
	if (fclose(w->fp) != 0)
		err(errno, "%s(): close failed: %s", __func__, w->path);
	fprintf(stderr,
			"minco readwise: wrote density cache %s; total_reads=%" PRIu64
			"; density_reads=%" PRIu64 "\n",
			w->path, w->hdr.total_reads, w->hdr.density_reads);
	w->fp = NULL;
}

static void ani_density_cache_reader_open(ani_density_cache_reader_t *r,
										  const char *path,
										  bool use_ctxobj96,
										  uint8_t nobjbits,
										  uint64_t density_threshold)
{
	if (!r || !path || path[0] == '\0')
		return;
	memset(r, 0, sizeof(*r));
	r->path = path;
	r->fp = fopen(path, "rb");
	if (!r->fp)
		err(errno, "%s(): cannot open %s", __func__, path);
	ani_density_cache_fread_record(r->fp, &r->hdr, sizeof(r->hdr), 1, path, false);
	if (memcmp(r->hdr.magic, ANI_DENSITY_CACHE_MAGIC, sizeof(r->hdr.magic)) != 0 ||
		r->hdr.version != ANI_DENSITY_CACHE_VERSION)
		errx(EXIT_FAILURE, "%s(): unsupported density cache format: %s", __func__, path);
	const uint32_t expected_storage =
		use_ctxobj96 ? ANI_DENSITY_CACHE_STORAGE96 : ANI_DENSITY_CACHE_STORAGE64;
	const uint32_t expected_item_size =
		use_ctxobj96 ? (uint32_t)sizeof(ctxobj96_t) : (uint32_t)sizeof(uint64_t);
	if (r->hdr.storage != expected_storage || r->hdr.item_size != expected_item_size)
		errx(EXIT_FAILURE, "%s(): density cache storage does not match reference sketch: %s",
			 __func__, path);
	if (r->hdr.nobjbits != nobjbits ||
		r->hdr.ctx_bits != Bitslen.ctx ||
		r->hdr.obj_bits != Bitslen.obj ||
		r->hdr.density_threshold != density_threshold)
		errx(EXIT_FAILURE, "%s(): density cache parameters do not match reference sketch: %s",
			 __func__, path);
	fprintf(stderr,
			"minco readwise: density cache replay active; in=%s total_reads=%" PRIu64
			" density_reads=%" PRIu64 "\n",
			path, r->hdr.total_reads, r->hdr.density_reads);
}

static size_t ani_density_cache_reader_read64(ani_density_cache_reader_t *r,
											  u64vec *read_vecs,
											  size_t max_records)
{
	if (!r || !r->fp || !read_vecs || max_records == 0)
		return 0;
	size_t nread = 0;
	while (nread < max_records && r->records_read < r->hdr.density_reads) {
		uint32_t n = 0;
		if (!ani_density_cache_fread_record(r->fp, &n, sizeof(n), 1, r->path, true))
			break;
		v_init(&read_vecs[nread], n);
		if (n) {
			v_reserve(&read_vecs[nread], n);
			ani_density_cache_fread_record(
				r->fp, read_vecs[nread].a, sizeof(read_vecs[nread].a[0]), n,
				r->path, false);
			read_vecs[nread].n = n;
		}
		++nread;
		++r->records_read;
	}
	return nread;
}

static size_t ani_density_cache_reader_read96(ani_density_cache_reader_t *r,
											  ani_ctxobj96_vec_t *read_vecs,
											  size_t max_records)
{
	if (!r || !r->fp || !read_vecs || max_records == 0)
		return 0;
	size_t nread = 0;
	while (nread < max_records && r->records_read < r->hdr.density_reads) {
		uint32_t n = 0;
		if (!ani_density_cache_fread_record(r->fp, &n, sizeof(n), 1, r->path, true))
			break;
		ani_ctxobj96_vec_init(&read_vecs[nread], n);
		if (n) {
			ani_ctxobj96_vec_reserve(&read_vecs[nread], n);
			ani_density_cache_fread_record(
				r->fp, read_vecs[nread].a, sizeof(read_vecs[nread].a[0]), n,
				r->path, false);
			read_vecs[nread].n = n;
		}
		++nread;
		++r->records_read;
	}
	return nread;
}

static void ani_density_cache_reader_close(ani_density_cache_reader_t *r)
{
	if (!r || !r->fp)
		return;
	if (r->records_read != r->hdr.density_reads)
		errx(EXIT_FAILURE,
			 "%s(): density cache ended after %" PRIu64 " records; expected %" PRIu64
			 ": %s",
			 __func__, r->records_read, r->hdr.density_reads, r->path);
	if (fclose(r->fp) != 0)
		err(errno, "%s(): close failed: %s", __func__, r->path);
	r->fp = NULL;
}

static int ani_ctxobj96_cmp(const void *pa, const void *pb)
{
	const ctxobj96_t *a = (const ctxobj96_t *)pa;
	const ctxobj96_t *b = (const ctxobj96_t *)pb;
	if (a->ctx != b->ctx)
		return (a->ctx > b->ctx) - (a->ctx < b->ctx);
	return (a->obj > b->obj) - (a->obj < b->obj);
}

static size_t ani_ctxobj96_dedup_sorted(ctxobj96_t *a, size_t n)
{
	if (n <= 1)
		return n;
	size_t w = 0;
	for (size_t r = 1; r < n; ++r)
		if (a[r].ctx != a[w].ctx || a[r].obj != a[w].obj)
			a[++w] = a[r];
	return w + 1;
}

static inline __uint128_t ani_mask128_bits(unsigned bits)
{
	if (bits >= 128)
		return ~((__uint128_t)0);
	return (((__uint128_t)1) << bits) - 1u;
}

static inline ctxobj96_t ani_coden128_to_ctxobj96(__uint128_t tuple, int codens)
{
	uint64_t ctx = 0;
	uint32_t obj = (uint32_t)(tuple & 0x3u);
	for (int i = 0; i < codens; ++i) {
		tuple >>= 2;
		ctx |= (uint64_t)(tuple & 0xFu) << (4 * i);
		tuple >>= 4;
		obj |= (uint32_t)(tuple & 0x3u) << (2 * (i + 1));
	}
	return ctxobj96_make(ctx, obj);
}

static void ani_extract_read_density_ctxobjs96(const char *s, int len,
											   ani_ctxobj96_vec_t *vec,
											   uint64_t density_threshold)
{
	if (len < (int)klen)
		return;
	const unsigned tuple_bits = 2u * klen;
	const __uint128_t tuple_mask = ani_mask128_bits(tuple_bits);
	const unsigned rev_shift = 2u * (klen - 1u);
	__uint128_t tuple = 0, crv = 0;
	int base = 0;

	for (int pos = 0; pos < len; ++pos) {
		const int bmap = Basemap[(unsigned char)s[pos]];
		if (unlikely(bmap == DEFAULT)) {
			base = 0;
			tuple = 0;
			crv = 0;
			continue;
		}
		const __uint128_t b2 = (uint64_t)bmap;
		tuple = ((tuple << 2) | b2) & tuple_mask;
		crv = (crv >> 2) | ((b2 ^ 3u) << rev_shift);
		if (unlikely(++base < (int)klen))
			continue;

		const ctxobj96_t fwd = ani_coden128_to_ctxobj96(tuple, NUM_CODENS);
		const ctxobj96_t rev = ani_coden128_to_ctxobj96(crv, NUM_CODENS);
		const ctxobj96_t rec =
			(fwd.ctx < rev.ctx || (fwd.ctx == rev.ctx && fwd.obj <= rev.obj))
				? fwd
				: rev;
#if ANI_APPLY_SOURCE_FILTER
		if (unlikely((uint32_t)mix64(rec.ctx) > FILTER))
			continue;
#endif
		if (mix64(rec.ctx ^ (uint64_t)MINCO_SEED) > density_threshold)
			continue;
		ani_ctxobj96_vec_push(vec, rec);
	}
}

static int ani_min_diff_sections_read_run_vs_ref_index(const uint64_t *qry,
													   size_t qry_begin,
													   size_t qry_end,
													   const ctxgidobj_t *ref,
													   size_t ref_begin,
													   size_t ref_end,
													   uint64_t objmask)
{
	int min_diff_sections = NUM_CODENS + 1;
	for (size_t qi = qry_begin; qi < qry_end; ++qi) {
		const uint32_t obj_q = (uint32_t)(qry[qi] & objmask);
		for (size_t ri = ref_begin; ri < ref_end; ++ri) {
			const uint32_t diff = obj_q ^ ref[ri].obj;
			if (diff == 0)
				return 0;
			const int d = dna_popcount(diff);
			if (d < min_diff_sections)
				min_diff_sections = d;
		}
	}
	return min_diff_sections;
}

static void ani_process_density_ctxobj_unit(
	u64vec *unit_vec,
	uint64_t unit_id,
	uint64_t unit_reads_with_density_ctx,
	const ctxgidobj_t *index,
	size_t index_n,
	const size_t *fence,
	int fence_k,
	uint32_t ref_n,
	bool ignoreconflict,
	uint8_t nobjbits,
	uint64_t gidmask_local,
	uint64_t objmask,
	uint64_t *read_marks,
	uint8_t *ref_hit_bits,
	uint32_t *ref_ctx_cov,
	kv_size_t *ref_ctx_hit_positions,
	kv_u64_t *ref_ctx_cov_hits,
	ani_u64_set_t *qry_ctx_seen,
	ani_u64_set_t *qry_ref_ctx_seen,
	ani_readwise_assign_mode_t assign_mode,
	ani_readwise_acc_t *acc,
	uint64_t *unique_read_marks,
	uint8_t *unique_ref_hit_bits,
	kv_size_t *unique_ref_ctx_hit_positions,
	kv_u64_t *unique_ref_ctx_cov_hits,
	ani_readwise_acc_t *unique_acc,
	const char *unit_read_name,
	const ani_readwise_trace_t *trace,
	ani_readwise_edge_trace_t *edge_trace)
{
	if (!unit_vec || unit_vec->n == 0)
		return;
	radix_sort_u64(unit_vec->a, unit_vec->n);
	unit_vec->n = dedup_sorted_uint64(unit_vec->a, unit_vec->n);
	if (unit_vec->n == 0)
		return;
	if (unit_reads_with_density_ctx == 0)
		unit_reads_with_density_ctx = 1;

	kv_readwise_candidate_t candidates;
	kv_init(candidates);
	for (size_t q = 0; q < unit_vec->n; ) {
		const uint64_t qctx = unit_vec->a[q] >> nobjbits;
		const size_t qbeg = q;
		do { ++q; } while (q < unit_vec->n && (unit_vec->a[q] >> nobjbits) == qctx);
		const size_t qend = q;
		if (qry_ctx_seen)
			(void)ani_u64_set_insert(qry_ctx_seen, qctx);

		kv_size(candidates) = 0;
		uint32_t best_diff = UINT32_MAX;
		size_t pos = lb_in_bucket_ctxgid(index, fence, fence_k, qctx);
		while (pos < index_n && (index[pos].ctxgid >> GID_NBITS) == qctx) {
			const uint64_t ctxgid = index[pos].ctxgid;
			const uint32_t gid = (uint32_t)(ctxgid & gidmask_local);
			const size_t ref_begin = pos;
			do { ++pos; } while (pos < index_n && index[pos].ctxgid == ctxgid);
			const size_t ref_end = pos;
			if (gid >= ref_n)
				continue;
			if (ignoreconflict && ref_end - ref_begin > 1)
				continue;

			const int min_diff = ani_min_diff_sections_read_run_vs_ref_index(
				unit_vec->a, qbeg, qend, index, ref_begin, ref_end, objmask);
			ani_readwise_candidate_t cand = {
				.ref_begin = ref_begin,
				.gid = gid,
				.diff = (uint32_t)min_diff,
			};
			kv_push(ani_readwise_candidate_t, candidates, cand);
			if ((uint32_t)min_diff < best_diff)
				best_diff = (uint32_t)min_diff;
		}
		if (kv_size(candidates) == 0)
			continue;

		size_t best_selected_n = 0;
		for (size_t ci = 0; ci < kv_size(candidates); ++ci)
			if (kv_A(candidates, ci).diff == best_diff)
				++best_selected_n;

		size_t selected_n = kv_size(candidates);
		if (assign_mode != ANI_READWISE_ASSIGN_ALL) {
			selected_n = best_selected_n;
			if (!selected_n)
				continue;
			if (assign_mode == ANI_READWISE_ASSIGN_BEST_DIFF_UNIQUE &&
				selected_n != 1)
				continue;
		}
		uint32_t cov_inc = 1u;
		if (assign_mode == ANI_READWISE_ASSIGN_BEST_DIFF_SPLIT) {
			cov_inc = selected_n
						  ? (uint32_t)((ANI_REFCOV_SPLIT_SCALE + selected_n / 2u) / selected_n)
						  : ANI_REFCOV_SPLIT_SCALE;
			if (!cov_inc)
				cov_inc = 1u;
		}
		ani_readwise_edge_trace_emit_group(edge_trace, unit_read_name, unit_id,
										   qctx, &candidates, best_diff,
										   selected_n, cov_inc, assign_mode);

		for (size_t ci = 0; ci < kv_size(candidates); ++ci) {
			const ani_readwise_candidate_t *cand = &kv_A(candidates, ci);
			if (assign_mode != ANI_READWISE_ASSIGN_ALL && cand->diff != best_diff)
				continue;
			const uint32_t gid = cand->gid;
			const size_t ref_begin = cand->ref_begin;
			if (trace && trace->fp && gid == trace->gid)
				ani_readwise_trace_emit(trace, unit_read_name, unit_id, qctx,
										ref_begin, cand->diff, best_diff,
										kv_size(candidates), selected_n, cov_inc);
			if (read_marks[gid] != unit_id) {
				read_marks[gid] = unit_id;
				acc[gid].reads_with_ctx_match += unit_reads_with_density_ctx;
				acc[gid].blocks_with_ctx_match++;
			}
			acc[gid].XnY_ctx++;
			if (!ani_bitset_test_set(ref_hit_bits, ref_begin)) {
				acc[gid].ref_ctx_hit++;
				if (ref_ctx_hit_positions)
					kv_push(size_t, *ref_ctx_hit_positions, ref_begin);
			}
			if (qry_ref_ctx_seen) {
				const uint64_t pair_key = (qctx << GID_NBITS) | gid;
				if (ani_u64_set_insert(qry_ref_ctx_seen, pair_key))
					acc[gid].qry_ctx_hit++;
			}

			if (ref_ctx_cov)
				ani_ref_covdiff_add_scaled_hit(&ref_ctx_cov[ref_begin], cand->diff, cov_inc);
			else if (ref_ctx_cov_hits)
				kv_push(uint64_t, *ref_ctx_cov_hits,
						ani_ref_covdiff_pack_hit(ref_begin, cand->diff, cov_inc));
			if (cand->diff > 0) {
				acc[gid].N_diff_obj++;
				acc[gid].N_diff_obj_section += (uint64_t)cand->diff;
				if (cand->diff > 1)
					acc[gid].N_mut2_ctx++;
			}
		}
		if (unique_acc && best_selected_n == 1) {
			for (size_t ci = 0; ci < kv_size(candidates); ++ci) {
				const ani_readwise_candidate_t *cand = &kv_A(candidates, ci);
				if (cand->diff != best_diff)
					continue;
				const uint32_t gid = cand->gid;
				const size_t ref_begin = cand->ref_begin;
				if (unique_read_marks && unique_read_marks[gid] != unit_id) {
					unique_read_marks[gid] = unit_id;
					unique_acc[gid].reads_with_ctx_match += unit_reads_with_density_ctx;
					unique_acc[gid].blocks_with_ctx_match++;
				}
				unique_acc[gid].XnY_ctx++;
				if (unique_ref_hit_bits &&
					!ani_bitset_test_set(unique_ref_hit_bits, ref_begin)) {
					unique_acc[gid].ref_ctx_hit++;
					if (unique_ref_ctx_hit_positions)
						kv_push(size_t, *unique_ref_ctx_hit_positions, ref_begin);
				}
				if (unique_ref_ctx_cov_hits)
					kv_push(uint64_t, *unique_ref_ctx_cov_hits,
							ani_ref_covdiff_pack_hit(ref_begin, cand->diff, 1u));
				if (cand->diff > 0) {
					unique_acc[gid].N_diff_obj++;
					unique_acc[gid].N_diff_obj_section += (uint64_t)cand->diff;
					if (cand->diff > 1)
						unique_acc[gid].N_mut2_ctx++;
				}
				break;
			}
		}
	}
	kv_destroy(candidates);
}

static void ani_process_density_units_parallel(
	kv_density_unit_t *units,
	ani_readwise_thread_state_t *states,
	int worker_n,
	const ctxgidobj_t *index,
	size_t index_n,
	const size_t *fence,
	int fence_k,
	uint32_t ref_n,
	bool ignoreconflict,
	uint8_t nobjbits,
	uint64_t gidmask_local,
	uint64_t objmask,
	bool need_ref_ctx_cov,
	bool track_query_sets,
	ani_readwise_assign_mode_t assign_mode,
	const ani_readwise_trace_t *trace,
	ani_readwise_edge_trace_t *edge_trace,
	ani_readwise_thread_state_t *unique_states)
{
	const size_t unit_n = units ? kv_size(*units) : 0;
	if (!unit_n)
		return;
#pragma omp parallel for num_threads(worker_n) schedule(dynamic, 64)
	for (size_t i = 0; i < unit_n; ++i) {
		const int tid = omp_get_thread_num();
		ani_readwise_thread_state_t *st = &states[tid];
		ani_readwise_thread_state_t *ust = unique_states ? &unique_states[tid] : NULL;
		ani_density_unit_t *unit = &kv_A(*units, i);
		ani_process_density_ctxobj_unit(
			&unit->vec, unit->id, unit->reads_with_density_ctx,
			index, index_n, fence, fence_k, ref_n, ignoreconflict,
			nobjbits, gidmask_local, objmask, st->read_marks,
			st->ref_hit_bits, NULL,
			&st->ref_ctx_hit_positions,
			need_ref_ctx_cov ? &st->ref_ctx_cov_hits : NULL,
			track_query_sets ? &st->qry_ctx_seen : NULL,
			track_query_sets ? &st->qry_ref_ctx_seen : NULL,
			assign_mode,
			st->acc,
			ust ? ust->read_marks : NULL,
			ust ? ust->ref_hit_bits : NULL,
			ust ? &ust->ref_ctx_hit_positions : NULL,
			ust ? &ust->ref_ctx_cov_hits : NULL,
			ust ? ust->acc : NULL,
			unit->read_name,
			trace,
			edge_trace);
	}
}

static void ani_density_unit_push_move(kv_density_unit_t *units,
									   u64vec *vec,
									   uint64_t id,
									   uint64_t reads_with_density_ctx,
									   const char *read_name)
{
	if (!vec || vec->n == 0)
		return;
	ani_density_unit_t unit = {
		.vec = *vec,
		.id = id,
		.reads_with_density_ctx = reads_with_density_ctx ? reads_with_density_ctx : 1,
		.read_name = read_name ? strdup(read_name) : NULL,
	};
	if (read_name && !unit.read_name)
		err(EXIT_FAILURE, "%s(): OOM readwise trace read name", __func__);
	kv_push(ani_density_unit_t, *units, unit);
	v_init(vec, 0);
}

static void ani_density_unit_push_copy(kv_density_unit_t *units,
									   const u64vec *vec,
									   uint64_t id,
									   uint64_t reads_with_density_ctx,
									   const char *read_name)
{
	if (!vec || vec->n == 0)
		return;
	u64vec copy = {0};
	v_init(&copy, vec->n);
	v_reserve(&copy, vec->n);
	memcpy(copy.a, vec->a, vec->n * sizeof(copy.a[0]));
	copy.n = vec->n;
	ani_density_unit_t unit = {
		.vec = copy,
		.id = id,
		.reads_with_density_ctx = reads_with_density_ctx ? reads_with_density_ctx : 1,
		.read_name = read_name ? strdup(read_name) : NULL,
	};
	if (read_name && !unit.read_name)
		err(EXIT_FAILURE, "%s(): OOM readwise exact sidecar read name", __func__);
	kv_push(ani_density_unit_t, *units, unit);
}

static int ani_min_diff_sections_read_run_vs_ref_index128(const ctxobj96_t *qry,
														  size_t qry_begin,
														  size_t qry_end,
														  const ctxgidobj128_t *ref,
														  size_t ref_begin,
														  size_t ref_end)
{
	int min_diff_sections = NUM_CODENS + 1;
	for (size_t qi = qry_begin; qi < qry_end; ++qi) {
		const uint32_t obj_q = qry[qi].obj;
		for (size_t ri = ref_begin; ri < ref_end; ++ri) {
			const uint32_t diff = obj_q ^ ref[ri].obj;
			if (diff == 0)
				return 0;
			const int d = dna_popcount(diff);
			if (d < min_diff_sections)
				min_diff_sections = d;
		}
	}
	return min_diff_sections;
}

static void ani_process_density_ctxobj96_unit(
	ani_ctxobj96_vec_t *unit_vec,
	uint64_t unit_id,
	uint64_t unit_reads_with_density_ctx,
	const ctxgidobj128_t *index,
	size_t index_n,
	const size_t *fence,
	int fence_k,
	uint32_t ref_n,
	bool ignoreconflict,
	uint64_t *read_marks,
	uint8_t *ref_hit_bits,
	uint32_t *ref_ctx_cov,
	kv_size_t *ref_ctx_hit_positions,
	kv_u64_t *ref_ctx_cov_hits,
	ani_readwise_assign_mode_t assign_mode,
	ani_readwise_acc_t *acc,
	uint64_t *unique_read_marks,
	uint8_t *unique_ref_hit_bits,
	kv_size_t *unique_ref_ctx_hit_positions,
	kv_u64_t *unique_ref_ctx_cov_hits,
	ani_readwise_acc_t *unique_acc,
	const char *unit_read_name,
	const ani_readwise_trace_t *trace,
	ani_readwise_edge_trace_t *edge_trace)
{
	if (!unit_vec || unit_vec->n == 0)
		return;
	qsort(unit_vec->a, unit_vec->n, sizeof(unit_vec->a[0]), ani_ctxobj96_cmp);
	unit_vec->n = ani_ctxobj96_dedup_sorted(unit_vec->a, unit_vec->n);
	if (unit_vec->n == 0)
		return;
	if (unit_reads_with_density_ctx == 0)
		unit_reads_with_density_ctx = 1;

	kv_readwise_candidate_t candidates;
	kv_init(candidates);
	for (size_t q = 0; q < unit_vec->n; ) {
		const uint64_t qctx = unit_vec->a[q].ctx;
		const size_t qbeg = q;
		do { ++q; } while (q < unit_vec->n && unit_vec->a[q].ctx == qctx);
		const size_t qend = q;

		kv_size(candidates) = 0;
		uint32_t best_diff = UINT32_MAX;
		size_t pos = lb_in_bucket_ctxgid128(index, fence, fence_k, qctx);
		while (pos < index_n && index[pos].ctx == qctx) {
			const uint64_t ctx = index[pos].ctx;
			const uint32_t gid = index[pos].gid;
			const size_t ref_begin = pos;
			do { ++pos; } while (pos < index_n &&
								  index[pos].ctx == ctx &&
								  index[pos].gid == gid);
			const size_t ref_end = pos;
			if (gid >= ref_n)
				continue;
			if (ignoreconflict && ref_end - ref_begin > 1)
				continue;

			const int min_diff = ani_min_diff_sections_read_run_vs_ref_index128(
				unit_vec->a, qbeg, qend, index, ref_begin, ref_end);
			ani_readwise_candidate_t cand = {
				.ref_begin = ref_begin,
				.gid = gid,
				.diff = (uint32_t)min_diff,
			};
			kv_push(ani_readwise_candidate_t, candidates, cand);
			if ((uint32_t)min_diff < best_diff)
				best_diff = (uint32_t)min_diff;
		}
		if (kv_size(candidates) == 0)
			continue;

		size_t best_selected_n = 0;
		for (size_t ci = 0; ci < kv_size(candidates); ++ci)
			if (kv_A(candidates, ci).diff == best_diff)
				++best_selected_n;

		size_t selected_n = kv_size(candidates);
		if (assign_mode != ANI_READWISE_ASSIGN_ALL) {
			selected_n = best_selected_n;
			if (!selected_n)
				continue;
			if (assign_mode == ANI_READWISE_ASSIGN_BEST_DIFF_UNIQUE &&
				selected_n != 1)
				continue;
		}
		uint32_t cov_inc = 1u;
		if (assign_mode == ANI_READWISE_ASSIGN_BEST_DIFF_SPLIT) {
			cov_inc = selected_n
						  ? (uint32_t)((ANI_REFCOV_SPLIT_SCALE + selected_n / 2u) / selected_n)
						  : ANI_REFCOV_SPLIT_SCALE;
			if (!cov_inc)
				cov_inc = 1u;
		}
		ani_readwise_edge_trace_emit_group(edge_trace, unit_read_name, unit_id,
										   qctx, &candidates, best_diff,
										   selected_n, cov_inc, assign_mode);

		for (size_t ci = 0; ci < kv_size(candidates); ++ci) {
			const ani_readwise_candidate_t *cand = &kv_A(candidates, ci);
			if (assign_mode != ANI_READWISE_ASSIGN_ALL && cand->diff != best_diff)
				continue;
			const uint32_t gid = cand->gid;
			const size_t ref_begin = cand->ref_begin;
			if (trace && trace->fp && gid == trace->gid)
				ani_readwise_trace_emit(trace, unit_read_name, unit_id, qctx,
										ref_begin, cand->diff, best_diff,
										kv_size(candidates), selected_n, cov_inc);
			if (read_marks[gid] != unit_id) {
				read_marks[gid] = unit_id;
				acc[gid].reads_with_ctx_match += unit_reads_with_density_ctx;
				acc[gid].blocks_with_ctx_match++;
			}
			acc[gid].XnY_ctx++;
			if (!ani_bitset_test_set(ref_hit_bits, ref_begin)) {
				acc[gid].ref_ctx_hit++;
				if (ref_ctx_hit_positions)
					kv_push(size_t, *ref_ctx_hit_positions, ref_begin);
			}
			if (ref_ctx_cov)
				ani_ref_covdiff_add_scaled_hit(&ref_ctx_cov[ref_begin], cand->diff, cov_inc);
			else if (ref_ctx_cov_hits)
				kv_push(uint64_t, *ref_ctx_cov_hits,
						ani_ref_covdiff_pack_hit(ref_begin, cand->diff, cov_inc));
			if (cand->diff > 0) {
				acc[gid].N_diff_obj++;
				acc[gid].N_diff_obj_section += (uint64_t)cand->diff;
				if (cand->diff > 1)
					acc[gid].N_mut2_ctx++;
			}
		}
		if (unique_acc && best_selected_n == 1) {
			for (size_t ci = 0; ci < kv_size(candidates); ++ci) {
				const ani_readwise_candidate_t *cand = &kv_A(candidates, ci);
				if (cand->diff != best_diff)
					continue;
				const uint32_t gid = cand->gid;
				const size_t ref_begin = cand->ref_begin;
				if (unique_read_marks && unique_read_marks[gid] != unit_id) {
					unique_read_marks[gid] = unit_id;
					unique_acc[gid].reads_with_ctx_match += unit_reads_with_density_ctx;
					unique_acc[gid].blocks_with_ctx_match++;
				}
				unique_acc[gid].XnY_ctx++;
				if (unique_ref_hit_bits &&
					!ani_bitset_test_set(unique_ref_hit_bits, ref_begin)) {
					unique_acc[gid].ref_ctx_hit++;
					if (unique_ref_ctx_hit_positions)
						kv_push(size_t, *unique_ref_ctx_hit_positions, ref_begin);
				}
				if (unique_ref_ctx_cov_hits)
					kv_push(uint64_t, *unique_ref_ctx_cov_hits,
							ani_ref_covdiff_pack_hit(ref_begin, cand->diff, 1u));
				if (cand->diff > 0) {
					unique_acc[gid].N_diff_obj++;
					unique_acc[gid].N_diff_obj_section += (uint64_t)cand->diff;
					if (cand->diff > 1)
						unique_acc[gid].N_mut2_ctx++;
				}
				break;
			}
		}
	}
	kv_destroy(candidates);
}

static void ani_process_density_units96_parallel(
	kv_density_unit96_t *units,
	ani_readwise_thread_state_t *states,
	int worker_n,
	const ctxgidobj128_t *index,
	size_t index_n,
	const size_t *fence,
	int fence_k,
	uint32_t ref_n,
	bool ignoreconflict,
	bool need_ref_ctx_cov,
	ani_readwise_assign_mode_t assign_mode,
	const ani_readwise_trace_t *trace,
	ani_readwise_edge_trace_t *edge_trace,
	ani_readwise_thread_state_t *unique_states)
{
	const size_t unit_n = units ? kv_size(*units) : 0;
	if (!unit_n)
		return;
#pragma omp parallel for num_threads(worker_n) schedule(dynamic, 64)
	for (size_t i = 0; i < unit_n; ++i) {
		const int tid = omp_get_thread_num();
		ani_readwise_thread_state_t *st = &states[tid];
		ani_readwise_thread_state_t *ust = unique_states ? &unique_states[tid] : NULL;
		ani_density_unit96_t *unit = &kv_A(*units, i);
		ani_process_density_ctxobj96_unit(
			&unit->vec, unit->id, unit->reads_with_density_ctx,
			index, index_n, fence, fence_k, ref_n, ignoreconflict,
			st->read_marks, st->ref_hit_bits, NULL,
			&st->ref_ctx_hit_positions,
			need_ref_ctx_cov ? &st->ref_ctx_cov_hits : NULL,
			assign_mode,
			st->acc,
			ust ? ust->read_marks : NULL,
			ust ? ust->ref_hit_bits : NULL,
			ust ? &ust->ref_ctx_hit_positions : NULL,
			ust ? &ust->ref_ctx_cov_hits : NULL,
			ust ? ust->acc : NULL,
			unit->read_name,
			trace,
			edge_trace);
	}
}

static void ani_density_unit96_push_move(kv_density_unit96_t *units,
										ani_ctxobj96_vec_t *vec,
										uint64_t id,
										uint64_t reads_with_density_ctx,
										const char *read_name)
{
	if (!vec || vec->n == 0)
		return;
	ani_density_unit96_t unit = {
		.vec = *vec,
		.id = id,
		.reads_with_density_ctx = reads_with_density_ctx ? reads_with_density_ctx : 1,
		.read_name = read_name ? strdup(read_name) : NULL,
	};
	if (read_name && !unit.read_name)
		err(EXIT_FAILURE, "%s(): OOM readwise trace read name", __func__);
	kv_push(ani_density_unit96_t, *units, unit);
	ani_ctxobj96_vec_init(vec, 0);
}

static void ani_density_unit96_push_copy(kv_density_unit96_t *units,
										const ani_ctxobj96_vec_t *vec,
										uint64_t id,
										uint64_t reads_with_density_ctx,
										const char *read_name)
{
	if (!vec || vec->n == 0)
		return;
	ani_ctxobj96_vec_t copy;
	ani_ctxobj96_vec_init(&copy, vec->n);
	ani_ctxobj96_vec_reserve(&copy, vec->n);
	memcpy(copy.a, vec->a, vec->n * sizeof(copy.a[0]));
	copy.n = vec->n;
	ani_density_unit96_t unit = {
		.vec = copy,
		.id = id,
		.reads_with_density_ctx = reads_with_density_ctx ? reads_with_density_ctx : 1,
		.read_name = read_name ? strdup(read_name) : NULL,
	};
	if (read_name && !unit.read_name)
		err(EXIT_FAILURE, "%s(): OOM readwise exact sidecar read name", __func__);
	kv_push(ani_density_unit96_t, *units, unit);
}

static uint32_t *ani_ref_ctx_counts_from_sorted_index(const ctxgidobj_t *index,
													  size_t index_n,
													  uint32_t ref_n,
													  bool ignoreconflict,
													  bool marker_only)
{
	uint32_t *counts = calloc((size_t)ref_n, sizeof(counts[0]));
	if (!counts)
		err(EXIT_FAILURE, "%s(): OOM reference context counts", __func__);
	const uint64_t gidmask_local = (1ULL << GID_NBITS) - 1ULL;
	for (size_t i = 0; i < index_n; ) {
		const uint64_t ctxgid = index[i].ctxgid;
		const uint32_t gid = (uint32_t)(ctxgid & gidmask_local);
		const size_t begin = i;
		do { ++i; } while (i < index_n && index[i].ctxgid == ctxgid);
		if (gid < ref_n &&
			(!ignoreconflict || i - begin == 1) &&
			(!marker_only || ani_ctxgid_group_marker_unique64(index, index_n, begin, i)))
			counts[gid]++;
	}
	return counts;
}

static uint32_t *ani_ref_ctx_counts_from_sorted_index128(const ctxgidobj128_t *index,
														 size_t index_n,
														 uint32_t ref_n,
														 bool ignoreconflict,
														 bool marker_only)
{
	uint32_t *counts = calloc((size_t)ref_n, sizeof(counts[0]));
	if (!counts)
		err(EXIT_FAILURE, "%s(): OOM reference context counts", __func__);
	for (size_t i = 0; i < index_n; ) {
		const uint64_t ctx = index[i].ctx;
		const uint32_t gid = index[i].gid;
		const size_t begin = i;
		do { ++i; } while (i < index_n &&
							index[i].ctx == ctx &&
							index[i].gid == gid);
		if (gid < ref_n &&
			(!ignoreconflict || i - begin == 1) &&
			(!marker_only || ani_ctxgid_group_marker_unique128(index, index_n, begin, i)))
			counts[gid]++;
	}
	return counts;
}

static ani_readwise_abundance_t *ani_depth_abundance_from_ref_coverage(
	const ctxgidobj_t *index,
	size_t index_n,
	uint32_t ref_n,
	bool ignoreconflict,
	const uint32_t *ref_ctx_cov,
	const uint16_t *ref_ctx_hit_weight,
	const uint32_t *ref_ctx_total,
	const ani_readwise_acc_t *acc,
	double coverage_scale,
	bool marker_only)
{
	ani_readwise_abundance_t *stats = calloc((size_t)ref_n, sizeof(stats[0]));
	long double *sum = calloc((size_t)ref_n, sizeof(sum[0]));
	long double *sumsq = calloc((size_t)ref_n, sizeof(sumsq[0]));
	long double *hit_weight = calloc((size_t)ref_n, sizeof(hit_weight[0]));
	if (!stats || !sum || !sumsq || !hit_weight)
		err(EXIT_FAILURE, "%s(): OOM abundance stats", __func__);
	if (!isfinite(coverage_scale) || coverage_scale <= 0.0)
		coverage_scale = 1.0;
	const long double inv_coverage_scale = 1.0L / (long double)coverage_scale;

	const uint64_t gidmask_local = (1ULL << GID_NBITS) - 1ULL;
	for (size_t i = 0; i < index_n; ) {
		const uint64_t ctxgid = index[i].ctxgid;
		const uint32_t gid = (uint32_t)(ctxgid & gidmask_local);
		const size_t begin = i;
		do { ++i; } while (i < index_n && index[i].ctxgid == ctxgid);
		if (gid >= ref_n || (ignoreconflict && i - begin > 1) ||
			(marker_only && !ani_ctxgid_group_marker_unique64(index, index_n, begin, i)))
			continue;
		const long double cov = ref_ctx_cov
									? (long double)ani_ref_covdiff_coverage(ref_ctx_cov[begin]) * inv_coverage_scale
									: 0.0L;
		sum[gid] += cov;
		sumsq[gid] += cov * cov;
		if (ref_ctx_hit_weight && ref_ctx_hit_weight[begin] > 0) {
			long double hit = (long double)ref_ctx_hit_weight[begin] * inv_coverage_scale;
			if (hit > 1.0L)
				hit = 1.0L;
			hit_weight[gid] += hit;
		}
	}

	long double total_mean_depth = 0.0L;
	for (uint32_t rn = 0; rn < ref_n; ++rn) {
		const uint32_t total = ref_ctx_total ? ref_ctx_total[rn] : 0;
		if (!total)
			continue;
		const long double denom = (long double)total;
		const long double mean = sum[rn] / denom;
		long double var = sumsq[rn] / denom - mean * mean;
		if (var < 0.0L && var > -1e-12L)
			var = 0.0L;
		if (var < 0.0L)
			var = 0.0L;
		const long double hit = ref_ctx_hit_weight
									? hit_weight[rn]
									: (long double)(acc ? acc[rn].ref_ctx_hit : 0);
		const long double breadth = hit > 0.0L ? hit / denom : 0.0L;
		stats[rn].ref_breadth = (double)breadth;
		stats[rn].ref_mean_depth = (double)mean;
		stats[rn].ref_hit_mean_depth = hit > 0.0L ? (double)(sum[rn] / hit) : 0.0;
		stats[rn].ref_depth_variance = (double)var;
		stats[rn].ref_depth_cv = mean > 0.0L ? (double)(sqrtl(var) / mean) : 0.0;
		stats[rn].ref_zero_fraction = (double)(1.0L - breadth);
		total_mean_depth += mean;
	}
	if (total_mean_depth > 0.0L) {
		for (uint32_t rn = 0; rn < ref_n; ++rn)
			stats[rn].relative_depth =
				(double)((long double)stats[rn].ref_mean_depth / total_mean_depth);
	}

	free(sum);
	free(sumsq);
	free(hit_weight);
	return stats;
}

static ani_readwise_abundance_t *ani_depth_abundance_from_ref_coverage128(
	const ctxgidobj128_t *index,
	size_t index_n,
	uint32_t ref_n,
	bool ignoreconflict,
	const uint32_t *ref_ctx_cov,
	const uint16_t *ref_ctx_hit_weight,
	const uint32_t *ref_ctx_total,
	const ani_readwise_acc_t *acc,
	double coverage_scale,
	bool marker_only)
{
	ani_readwise_abundance_t *stats = calloc((size_t)ref_n, sizeof(stats[0]));
	long double *sum = calloc((size_t)ref_n, sizeof(sum[0]));
	long double *sumsq = calloc((size_t)ref_n, sizeof(sumsq[0]));
	long double *hit_weight = calloc((size_t)ref_n, sizeof(hit_weight[0]));
	if (!stats || !sum || !sumsq || !hit_weight)
		err(EXIT_FAILURE, "%s(): OOM abundance stats", __func__);
	if (!isfinite(coverage_scale) || coverage_scale <= 0.0)
		coverage_scale = 1.0;
	const long double inv_coverage_scale = 1.0L / (long double)coverage_scale;

	for (size_t i = 0; i < index_n; ) {
		const uint64_t ctx = index[i].ctx;
		const uint32_t gid = index[i].gid;
		const size_t begin = i;
		do { ++i; } while (i < index_n &&
							index[i].ctx == ctx &&
							index[i].gid == gid);
		if (gid >= ref_n || (ignoreconflict && i - begin > 1) ||
			(marker_only && !ani_ctxgid_group_marker_unique128(index, index_n, begin, i)))
			continue;
		const long double cov = ref_ctx_cov
									? (long double)ani_ref_covdiff_coverage(ref_ctx_cov[begin]) * inv_coverage_scale
									: 0.0L;
		sum[gid] += cov;
		sumsq[gid] += cov * cov;
		if (ref_ctx_hit_weight && ref_ctx_hit_weight[begin] > 0) {
			long double hit = (long double)ref_ctx_hit_weight[begin] * inv_coverage_scale;
			if (hit > 1.0L)
				hit = 1.0L;
			hit_weight[gid] += hit;
		}
	}

	long double total_mean_depth = 0.0L;
	for (uint32_t rn = 0; rn < ref_n; ++rn) {
		const uint32_t total = ref_ctx_total ? ref_ctx_total[rn] : 0;
		if (!total)
			continue;
		const long double denom = (long double)total;
		const long double mean = sum[rn] / denom;
		long double var = sumsq[rn] / denom - mean * mean;
		if (var < 0.0L && var > -1e-12L)
			var = 0.0L;
		if (var < 0.0L)
			var = 0.0L;
		const long double hit = ref_ctx_hit_weight
									? hit_weight[rn]
									: (long double)(acc ? acc[rn].ref_ctx_hit : 0);
		const long double breadth = hit > 0.0L ? hit / denom : 0.0L;
		stats[rn].ref_breadth = (double)breadth;
		stats[rn].ref_mean_depth = (double)mean;
		stats[rn].ref_hit_mean_depth = hit > 0.0L ? (double)(sum[rn] / hit) : 0.0;
		stats[rn].ref_depth_variance = (double)var;
		stats[rn].ref_depth_cv = mean > 0.0L ? (double)(sqrtl(var) / mean) : 0.0;
		stats[rn].ref_zero_fraction = (double)(1.0L - breadth);
		total_mean_depth += mean;
	}
	if (total_mean_depth > 0.0L) {
		for (uint32_t rn = 0; rn < ref_n; ++rn)
			stats[rn].relative_depth =
				(double)((long double)stats[rn].ref_mean_depth / total_mean_depth);
	}

	free(sum);
	free(sumsq);
	free(hit_weight);
	return stats;
}

static ani_readwise_acc_t *ani_unique_best_features_from_ref_covdiff(
	const ctxgidobj_t *index,
	size_t index_n,
	uint32_t ref_n,
	bool ignoreconflict,
	const uint32_t *ref_ctx_cov,
	bool marker_only)
{
	if (!ref_ctx_cov)
		return NULL;
	ani_readwise_acc_t *features = calloc((size_t)ref_n, sizeof(features[0]));
	if (!features)
		err(EXIT_FAILURE, "%s(): OOM unique readwise ANI features", __func__);

	const uint64_t gidmask_local = (1ULL << GID_NBITS) - 1ULL;
	for (size_t i = 0; i < index_n; ) {
		const uint64_t ctxgid = index[i].ctxgid;
		const uint32_t gid = (uint32_t)(ctxgid & gidmask_local);
		const size_t begin = i;
		do { ++i; } while (i < index_n && index[i].ctxgid == ctxgid);
		if (gid >= ref_n || (ignoreconflict && i - begin > 1) ||
			(marker_only && !ani_ctxgid_group_marker_unique64(index, index_n, begin, i)))
			continue;
		const uint32_t packed = ref_ctx_cov[begin];
		const uint32_t code = packed >> ANI_REFCOV_DIFF_SHIFT;
		if (!code || !ani_ref_covdiff_coverage(packed))
			continue;
		const uint32_t diff = ani_ref_covdiff_decode_diff(code);
		features[gid].XnY_ctx++;
		if (diff > 0) {
			features[gid].N_diff_obj++;
			features[gid].N_diff_obj_section += (uint64_t)diff;
			if (diff > 1)
				features[gid].N_mut2_ctx++;
		}
	}
	return features;
}

static ani_readwise_acc_t *ani_unique_best_features_from_ref_covdiff128(
	const ctxgidobj128_t *index,
	size_t index_n,
	uint32_t ref_n,
	bool ignoreconflict,
	const uint32_t *ref_ctx_cov,
	bool marker_only)
{
	if (!ref_ctx_cov)
		return NULL;
	ani_readwise_acc_t *features = calloc((size_t)ref_n, sizeof(features[0]));
	if (!features)
		err(EXIT_FAILURE, "%s(): OOM unique readwise ANI features", __func__);

	for (size_t i = 0; i < index_n; ) {
		const uint64_t ctx = index[i].ctx;
		const uint32_t gid = index[i].gid;
		const size_t begin = i;
		do { ++i; } while (i < index_n &&
							index[i].ctx == ctx &&
							index[i].gid == gid);
		if (gid >= ref_n || (ignoreconflict && i - begin > 1) ||
			(marker_only && !ani_ctxgid_group_marker_unique128(index, index_n, begin, i)))
			continue;
		const uint32_t packed = ref_ctx_cov[begin];
		const uint32_t code = packed >> ANI_REFCOV_DIFF_SHIFT;
		if (!code || !ani_ref_covdiff_coverage(packed))
			continue;
		const uint32_t diff = ani_ref_covdiff_decode_diff(code);
		features[gid].XnY_ctx++;
		if (diff > 0) {
			features[gid].N_diff_obj++;
			features[gid].N_diff_obj_section += (uint64_t)diff;
			if (diff > 1)
				features[gid].N_mut2_ctx++;
		}
	}
	return features;
}

static inline uint64_t ani_readwise_scaled_round_u64(uint64_t x, uint32_t scale)
{
	if (scale <= 1u)
		return x;
	return (x + (uint64_t)scale / 2u) / (uint64_t)scale;
}

static inline uint32_t ani_readwise_scaled_round_u32(uint64_t x, uint32_t scale)
{
	return ani_clamp_u64_to_u32(ani_readwise_scaled_round_u64(x, scale));
}

static inline int ani_readwise_scaled_round_int(uint64_t x, uint32_t scale)
{
	return ani_clamp_u64_to_int(ani_readwise_scaled_round_u64(x, scale));
}

static inline uint64_t ani_readwise_probability_weight(double keep_probability)
{
	if (!isfinite(keep_probability) || keep_probability <= 0.0)
		return 0;
	if (keep_probability >= 1.0)
		return ANI_READWISE_PROB_SCALE;
	return (uint64_t)llround(keep_probability * (double)ANI_READWISE_PROB_SCALE);
}

static ani_readwise_acc_t *ani_reliable_features_from_ref_covdiff(
	const ctxgidobj_t *index,
	size_t index_n,
	uint32_t ref_n,
	bool ignoreconflict,
	const uint32_t *ref_ctx_cov,
	const ani_readwise_abundance_t *abundance_stats,
	const uint32_t *ref_ctx_total,
	ani_readwise_ctx_filter_model_t filter_model,
	double coverage_scale,
	double fake_threshold,
	ani_readwise_filter_stats_t *filter_stats,
	bool marker_only)
{
	if (!ref_ctx_cov || !abundance_stats)
		return NULL;
	if (!isfinite(coverage_scale) || coverage_scale <= 0.0)
		coverage_scale = 1.0;
	if (!isfinite(fake_threshold) || fake_threshold <= 0.0)
		fake_threshold = MINCO_DEFAULT_READWISE_FAKE_CTX_THRESHOLD;

	ani_readwise_acc_t *features = calloc((size_t)ref_n, sizeof(features[0]));
	if (!features)
		err(EXIT_FAILURE, "%s(): OOM reliable readwise ANI features", __func__);
	if (filter_stats)
		memset(filter_stats, 0, (size_t)ref_n * sizeof(filter_stats[0]));

	long double *product_sum = NULL;
	long double *product_sumsq = NULL;
	double *product_max = NULL;
	uint32_t *product_nonzero = NULL;
	if (ani_readwise_ctx_filter_uses_product(filter_model)) {
		product_sum = calloc((size_t)ref_n, sizeof(product_sum[0]));
		product_sumsq = calloc((size_t)ref_n, sizeof(product_sumsq[0]));
		product_max = calloc((size_t)ref_n, sizeof(product_max[0]));
		product_nonzero = calloc((size_t)ref_n, sizeof(product_nonzero[0]));
		if (!product_sum || !product_sumsq || !product_max || !product_nonzero)
			err(EXIT_FAILURE, "%s(): OOM readwise product filter stats", __func__);
	}

	const uint64_t gidmask_local = (1ULL << GID_NBITS) - 1ULL;
	if (ani_readwise_ctx_filter_uses_product(filter_model)) {
		for (size_t i = 0; i < index_n; ) {
			const uint64_t ctxgid = index[i].ctxgid;
			const uint32_t gid = (uint32_t)(ctxgid & gidmask_local);
			const size_t begin = i;
			do { ++i; } while (i < index_n && index[i].ctxgid == ctxgid);
			if (gid >= ref_n || (ignoreconflict && i - begin > 1) ||
				(marker_only && !ani_ctxgid_group_marker_unique64(index, index_n, begin, i)))
				continue;
			const uint32_t packed = ref_ctx_cov[begin];
			const uint32_t coverage_raw = ani_ref_covdiff_coverage(packed);
			const uint32_t code = packed >> ANI_REFCOV_DIFF_SHIFT;
			if (!code || !coverage_raw)
				continue;
			const uint32_t diff = ani_ref_covdiff_decode_diff(code);
			if (!diff)
				continue;
			const double cov = (double)coverage_raw / coverage_scale;
			const double product = (double)diff * cov;
			if (!isfinite(product) || product <= 0.0)
				continue;
			product_sum[gid] += (long double)product;
			product_sumsq[gid] += (long double)product * (long double)product;
			if (product > product_max[gid])
				product_max[gid] = product;
			product_nonzero[gid]++;
		}
	}
	double *product_topfrac_thresholds = NULL;
	if (filter_model == ANI_READWISE_CTX_FILTER_PRODUCT_TOPFRAC) {
		product_topfrac_thresholds = ani_readwise_product_topfrac_thresholds(
			index, index_n, ref_n, ignoreconflict, ref_ctx_cov,
			coverage_scale, fake_threshold, marker_only);
	}
	ani_readwise_product_topfrac_median_t *product_topfrac_medians = NULL;
	if (filter_model == ANI_READWISE_CTX_FILTER_PRODUCT_TOPFRAC_MEDIAN ||
		filter_model == ANI_READWISE_CTX_FILTER_PRODUCT1_TOPFRAC_MEDIAN) {
		product_topfrac_medians = ani_readwise_product_topfrac_median_stats(
			index, index_n, ref_n, ignoreconflict, ref_ctx_cov,
			coverage_scale, fake_threshold,
			filter_model == ANI_READWISE_CTX_FILTER_PRODUCT1_TOPFRAC_MEDIAN,
			marker_only);
	}
	for (size_t i = 0; i < index_n; ) {
		const uint64_t ctxgid = index[i].ctxgid;
		const uint32_t gid = (uint32_t)(ctxgid & gidmask_local);
		const size_t begin = i;
		do { ++i; } while (i < index_n && index[i].ctxgid == ctxgid);
		if (gid >= ref_n || (ignoreconflict && i - begin > 1) ||
			(marker_only && !ani_ctxgid_group_marker_unique64(index, index_n, begin, i)))
			continue;
		const uint32_t packed = ref_ctx_cov[begin];
		const uint32_t coverage_raw = ani_ref_covdiff_coverage(packed);
		const uint32_t code = packed >> ANI_REFCOV_DIFF_SHIFT;
		if (!code || !coverage_raw)
			continue;
		const uint32_t raw_diff = ani_ref_covdiff_decode_diff(code);
		if (filter_stats)
			filter_stats[gid].raw_xny_ctx++;
		uint32_t diff = raw_diff;
		double cov = (double)coverage_raw / coverage_scale;
		const double fake_prob = ani_readwise_fake_ctx_probability(
			raw_diff, cov, &abundance_stats[gid], fake_threshold);
		if (filter_stats) {
			filter_stats[gid].fake_prob_sum += (long double)fake_prob;
			filter_stats[gid].fake_prob_weighted_sum += (long double)fake_prob * (long double)cov;
			filter_stats[gid].fake_prob_weight_sum += (long double)cov;
		}
		double product_lambda = 0.0;
		double product_nb_mean = 0.0;
		double product_nb_variance = 0.0;
		uint32_t product_nb_nonzero = 0u;
		if (ani_readwise_ctx_filter_uses_product(filter_model) &&
			ref_ctx_total && ref_ctx_total[gid]) {
			long double robust_sum = product_sum[gid];
			long double robust_sumsq = product_sumsq[gid];
			uint32_t robust_n = product_nonzero[gid];
			if (product_nonzero[gid] > 1u)
				robust_sum -= (long double)product_max[gid];
			if (product_nonzero[gid] > 1u) {
				robust_sumsq -= (long double)product_max[gid] * (long double)product_max[gid];
				robust_n--;
			}
			if (robust_sum < 0.0L)
				robust_sum = 0.0L;
			if (robust_sumsq < 0.0L)
				robust_sumsq = 0.0L;
			product_lambda = (double)(robust_sum / (long double)ref_ctx_total[gid]);
			if (robust_n > 0u) {
				product_nb_nonzero = robust_n;
				product_nb_mean = (double)(robust_sum / (long double)robust_n);
				if (robust_n > 1u) {
					const long double mean_ld = robust_sum / (long double)robust_n;
					long double ss = robust_sumsq - robust_sum * mean_ld;
					if (ss < 0.0L)
						ss = 0.0L;
					product_nb_variance = (double)(ss / (long double)(robust_n - 1u));
				} else {
					product_nb_variance = product_nb_mean;
				}
			}
		}
		const bool adjusted =
			((filter_model == ANI_READWISE_CTX_FILTER_PRODUCT_TOPFRAC_MEDIAN ||
			  filter_model == ANI_READWISE_CTX_FILTER_PRODUCT1_TOPFRAC_MEDIAN) &&
			 ani_readwise_apply_product_topfrac_median_ctx(
				 &diff, &cov,
				 product_topfrac_medians ? &product_topfrac_medians[gid] : NULL,
				 filter_model == ANI_READWISE_CTX_FILTER_PRODUCT1_TOPFRAC_MEDIAN));
		const bool reject =
			(filter_model == ANI_READWISE_CTX_FILTER_POISSON_DIFF &&
			 ani_readwise_reject_poisson_diff_ctx(
				 raw_diff, cov, &abundance_stats[gid], fake_threshold)) ||
			(filter_model == ANI_READWISE_CTX_FILTER_POISSON_DEPTH &&
			 ani_readwise_reject_poisson_depth_ctx(
				 raw_diff, cov, &abundance_stats[gid],
				 ref_ctx_total ? ref_ctx_total[gid] : 0u,
				 fake_threshold)) ||
			(filter_model == ANI_READWISE_CTX_FILTER_POISSON_PRODUCT &&
			 ani_readwise_reject_poisson_product_ctx(
				 raw_diff, cov, product_lambda,
				 ref_ctx_total ? ref_ctx_total[gid] : 0u,
				 fake_threshold)) ||
			(filter_model == ANI_READWISE_CTX_FILTER_PRODUCT_NB &&
			 ani_readwise_reject_product_nb_ctx(
				 raw_diff, cov, product_nb_mean, product_nb_variance,
				 product_nb_nonzero,
				 ref_ctx_total ? ref_ctx_total[gid] : 0u,
				 fake_threshold)) ||
			(filter_model == ANI_READWISE_CTX_FILTER_PRODUCT_TOPFRAC &&
			 ani_readwise_reject_product_topfrac_ctx(
				 raw_diff, cov,
				 product_topfrac_thresholds ? product_topfrac_thresholds[gid] : NAN));
		if (reject) {
			if (filter_stats) {
				filter_stats[gid].rejected_ctx++;
				if (raw_diff > 0)
					filter_stats[gid].rejected_diff_ctx++;
			}
			continue;
		}
		if (adjusted && filter_stats) {
			filter_stats[gid].rejected_ctx++;
			if (raw_diff > 0)
				filter_stats[gid].rejected_diff_ctx++;
		}
		uint64_t weight = 1u;
		if (filter_model == ANI_READWISE_CTX_FILTER_FAKE_PROB)
			weight = raw_diff == 0
						 ? (uint64_t)ANI_READWISE_PROB_SCALE
						 : ani_readwise_probability_weight(1.0 - fake_prob);
		if (!weight)
			continue;
		features[gid].XnY_ctx += weight;
		if (diff > 0) {
			features[gid].N_diff_obj += weight;
			features[gid].N_diff_obj_section += weight * (uint64_t)diff;
			if (diff > 1)
				features[gid].N_mut2_ctx += weight;
		}
	}
	free(product_sum);
	free(product_sumsq);
	free(product_max);
	free(product_nonzero);
	free(product_topfrac_thresholds);
	free(product_topfrac_medians);
	return features;
}

static ani_readwise_acc_t *ani_reliable_features_from_ref_covdiff128(
	const ctxgidobj128_t *index,
	size_t index_n,
	uint32_t ref_n,
	bool ignoreconflict,
	const uint32_t *ref_ctx_cov,
	const ani_readwise_abundance_t *abundance_stats,
	const uint32_t *ref_ctx_total,
	ani_readwise_ctx_filter_model_t filter_model,
	double coverage_scale,
	double fake_threshold,
	ani_readwise_filter_stats_t *filter_stats,
	bool marker_only)
{
	if (!ref_ctx_cov || !abundance_stats)
		return NULL;
	if (!isfinite(coverage_scale) || coverage_scale <= 0.0)
		coverage_scale = 1.0;
	if (!isfinite(fake_threshold) || fake_threshold <= 0.0)
		fake_threshold = MINCO_DEFAULT_READWISE_FAKE_CTX_THRESHOLD;

	ani_readwise_acc_t *features = calloc((size_t)ref_n, sizeof(features[0]));
	if (!features)
		err(EXIT_FAILURE, "%s(): OOM reliable readwise ANI features", __func__);
	if (filter_stats)
		memset(filter_stats, 0, (size_t)ref_n * sizeof(filter_stats[0]));

	long double *product_sum = NULL;
	long double *product_sumsq = NULL;
	double *product_max = NULL;
	uint32_t *product_nonzero = NULL;
	if (ani_readwise_ctx_filter_uses_product(filter_model)) {
		product_sum = calloc((size_t)ref_n, sizeof(product_sum[0]));
		product_sumsq = calloc((size_t)ref_n, sizeof(product_sumsq[0]));
		product_max = calloc((size_t)ref_n, sizeof(product_max[0]));
		product_nonzero = calloc((size_t)ref_n, sizeof(product_nonzero[0]));
		if (!product_sum || !product_sumsq || !product_max || !product_nonzero)
			err(EXIT_FAILURE, "%s(): OOM readwise product filter stats", __func__);
	}

	if (ani_readwise_ctx_filter_uses_product(filter_model)) {
		for (size_t i = 0; i < index_n; ) {
			const uint64_t ctx = index[i].ctx;
			const uint32_t gid = index[i].gid;
			const size_t begin = i;
			do { ++i; } while (i < index_n &&
								index[i].ctx == ctx &&
								index[i].gid == gid);
			if (gid >= ref_n || (ignoreconflict && i - begin > 1) ||
				(marker_only && !ani_ctxgid_group_marker_unique128(index, index_n, begin, i)))
				continue;
			const uint32_t packed = ref_ctx_cov[begin];
			const uint32_t coverage_raw = ani_ref_covdiff_coverage(packed);
			const uint32_t code = packed >> ANI_REFCOV_DIFF_SHIFT;
			if (!code || !coverage_raw)
				continue;
			const uint32_t diff = ani_ref_covdiff_decode_diff(code);
			if (!diff)
				continue;
			const double cov = (double)coverage_raw / coverage_scale;
			const double product = (double)diff * cov;
			if (!isfinite(product) || product <= 0.0)
				continue;
			product_sum[gid] += (long double)product;
			product_sumsq[gid] += (long double)product * (long double)product;
			if (product > product_max[gid])
				product_max[gid] = product;
			product_nonzero[gid]++;
		}
	}
	double *product_topfrac_thresholds = NULL;
	if (filter_model == ANI_READWISE_CTX_FILTER_PRODUCT_TOPFRAC) {
		product_topfrac_thresholds = ani_readwise_product_topfrac_thresholds128(
			index, index_n, ref_n, ignoreconflict, ref_ctx_cov,
			coverage_scale, fake_threshold, marker_only);
	}
	ani_readwise_product_topfrac_median_t *product_topfrac_medians = NULL;
	if (filter_model == ANI_READWISE_CTX_FILTER_PRODUCT_TOPFRAC_MEDIAN ||
		filter_model == ANI_READWISE_CTX_FILTER_PRODUCT1_TOPFRAC_MEDIAN) {
		product_topfrac_medians = ani_readwise_product_topfrac_median_stats128(
			index, index_n, ref_n, ignoreconflict, ref_ctx_cov,
			coverage_scale, fake_threshold,
			filter_model == ANI_READWISE_CTX_FILTER_PRODUCT1_TOPFRAC_MEDIAN,
			marker_only);
	}
	for (size_t i = 0; i < index_n; ) {
		const uint64_t ctx = index[i].ctx;
		const uint32_t gid = index[i].gid;
		const size_t begin = i;
		do { ++i; } while (i < index_n &&
							index[i].ctx == ctx &&
							index[i].gid == gid);
		if (gid >= ref_n || (ignoreconflict && i - begin > 1) ||
			(marker_only && !ani_ctxgid_group_marker_unique128(index, index_n, begin, i)))
			continue;
		const uint32_t packed = ref_ctx_cov[begin];
		const uint32_t coverage_raw = ani_ref_covdiff_coverage(packed);
		const uint32_t code = packed >> ANI_REFCOV_DIFF_SHIFT;
		if (!code || !coverage_raw)
			continue;
		const uint32_t raw_diff = ani_ref_covdiff_decode_diff(code);
		if (filter_stats)
			filter_stats[gid].raw_xny_ctx++;
		uint32_t diff = raw_diff;
		double cov = (double)coverage_raw / coverage_scale;
		const double fake_prob = ani_readwise_fake_ctx_probability(
			raw_diff, cov, &abundance_stats[gid], fake_threshold);
		if (filter_stats) {
			filter_stats[gid].fake_prob_sum += (long double)fake_prob;
			filter_stats[gid].fake_prob_weighted_sum += (long double)fake_prob * (long double)cov;
			filter_stats[gid].fake_prob_weight_sum += (long double)cov;
		}
		double product_lambda = 0.0;
		double product_nb_mean = 0.0;
		double product_nb_variance = 0.0;
		uint32_t product_nb_nonzero = 0u;
		if (ani_readwise_ctx_filter_uses_product(filter_model) &&
			ref_ctx_total && ref_ctx_total[gid]) {
			long double robust_sum = product_sum[gid];
			long double robust_sumsq = product_sumsq[gid];
			uint32_t robust_n = product_nonzero[gid];
			if (product_nonzero[gid] > 1u)
				robust_sum -= (long double)product_max[gid];
			if (product_nonzero[gid] > 1u) {
				robust_sumsq -= (long double)product_max[gid] * (long double)product_max[gid];
				robust_n--;
			}
			if (robust_sum < 0.0L)
				robust_sum = 0.0L;
			if (robust_sumsq < 0.0L)
				robust_sumsq = 0.0L;
			product_lambda = (double)(robust_sum / (long double)ref_ctx_total[gid]);
			if (robust_n > 0u) {
				product_nb_nonzero = robust_n;
				product_nb_mean = (double)(robust_sum / (long double)robust_n);
				if (robust_n > 1u) {
					const long double mean_ld = robust_sum / (long double)robust_n;
					long double ss = robust_sumsq - robust_sum * mean_ld;
					if (ss < 0.0L)
						ss = 0.0L;
					product_nb_variance = (double)(ss / (long double)(robust_n - 1u));
				} else {
					product_nb_variance = product_nb_mean;
				}
			}
		}
		const bool adjusted =
			((filter_model == ANI_READWISE_CTX_FILTER_PRODUCT_TOPFRAC_MEDIAN ||
			  filter_model == ANI_READWISE_CTX_FILTER_PRODUCT1_TOPFRAC_MEDIAN) &&
			 ani_readwise_apply_product_topfrac_median_ctx(
				 &diff, &cov,
				 product_topfrac_medians ? &product_topfrac_medians[gid] : NULL,
				 filter_model == ANI_READWISE_CTX_FILTER_PRODUCT1_TOPFRAC_MEDIAN));
		const bool reject =
			(filter_model == ANI_READWISE_CTX_FILTER_POISSON_DIFF &&
			 ani_readwise_reject_poisson_diff_ctx(
				 raw_diff, cov, &abundance_stats[gid], fake_threshold)) ||
			(filter_model == ANI_READWISE_CTX_FILTER_POISSON_DEPTH &&
			 ani_readwise_reject_poisson_depth_ctx(
				 raw_diff, cov, &abundance_stats[gid],
				 ref_ctx_total ? ref_ctx_total[gid] : 0u,
				 fake_threshold)) ||
			(filter_model == ANI_READWISE_CTX_FILTER_POISSON_PRODUCT &&
			 ani_readwise_reject_poisson_product_ctx(
				 raw_diff, cov, product_lambda,
				 ref_ctx_total ? ref_ctx_total[gid] : 0u,
				 fake_threshold)) ||
			(filter_model == ANI_READWISE_CTX_FILTER_PRODUCT_NB &&
			 ani_readwise_reject_product_nb_ctx(
				 raw_diff, cov, product_nb_mean, product_nb_variance,
				 product_nb_nonzero,
				 ref_ctx_total ? ref_ctx_total[gid] : 0u,
				 fake_threshold)) ||
			(filter_model == ANI_READWISE_CTX_FILTER_PRODUCT_TOPFRAC &&
			 ani_readwise_reject_product_topfrac_ctx(
				 raw_diff, cov,
				 product_topfrac_thresholds ? product_topfrac_thresholds[gid] : NAN));
		if (reject) {
			if (filter_stats) {
				filter_stats[gid].rejected_ctx++;
				if (raw_diff > 0)
					filter_stats[gid].rejected_diff_ctx++;
			}
			continue;
		}
		if (adjusted && filter_stats) {
			filter_stats[gid].rejected_ctx++;
			if (raw_diff > 0)
				filter_stats[gid].rejected_diff_ctx++;
		}
		uint64_t weight = 1u;
		if (filter_model == ANI_READWISE_CTX_FILTER_FAKE_PROB)
			weight = raw_diff == 0
						 ? (uint64_t)ANI_READWISE_PROB_SCALE
						 : ani_readwise_probability_weight(1.0 - fake_prob);
		if (!weight)
			continue;
		features[gid].XnY_ctx += weight;
		if (diff > 0) {
			features[gid].N_diff_obj += weight;
			features[gid].N_diff_obj_section += weight * (uint64_t)diff;
			if (diff > 1)
				features[gid].N_mut2_ctx += weight;
		}
	}
	free(product_sum);
	free(product_sumsq);
	free(product_max);
	free(product_nonzero);
	free(product_topfrac_thresholds);
	free(product_topfrac_medians);
	return features;
}

static ani_readwise_reliable_abundance_t *ani_reliable_abundance_from_ref_covdiff(
	const ctxgidobj_t *index,
	size_t index_n,
	uint32_t ref_n,
	bool ignoreconflict,
	const uint32_t *ref_ctx_cov,
	const uint16_t *ref_ctx_hit_weight,
	const ani_readwise_abundance_t *abundance_stats,
	const uint32_t *ref_ctx_total,
	ani_readwise_ctx_filter_model_t filter_model,
	double coverage_scale,
	double fake_threshold,
	bool marker_only)
{
	if (!ref_ctx_cov || !abundance_stats || !ref_ctx_total)
		return NULL;
	if (!isfinite(coverage_scale) || coverage_scale <= 0.0)
		coverage_scale = 1.0;
	if (!isfinite(fake_threshold) || fake_threshold <= 0.0)
		fake_threshold = MINCO_DEFAULT_READWISE_FAKE_CTX_THRESHOLD;

	ani_readwise_reliable_abundance_t *stats =
		calloc((size_t)ref_n, sizeof(stats[0]));
	long double *sum = calloc((size_t)ref_n, sizeof(sum[0]));
	long double *sumsq = calloc((size_t)ref_n, sizeof(sumsq[0]));
	long double *hit_weight = calloc((size_t)ref_n, sizeof(hit_weight[0]));
	long double *hit_sum = calloc((size_t)ref_n, sizeof(hit_sum[0]));
	long double *hit_sumsq = calloc((size_t)ref_n, sizeof(hit_sumsq[0]));
	uint64_t *hit_ctx = calloc((size_t)ref_n, sizeof(hit_ctx[0]));
	if (!stats || !sum || !sumsq || !hit_weight ||
		!hit_sum || !hit_sumsq || !hit_ctx)
		err(EXIT_FAILURE, "%s(): OOM reliable abundance stats", __func__);

	long double *product_sum = NULL;
	long double *product_sumsq = NULL;
	double *product_max = NULL;
	uint32_t *product_nonzero = NULL;
	if (ani_readwise_ctx_filter_uses_product(filter_model)) {
		product_sum = calloc((size_t)ref_n, sizeof(product_sum[0]));
		product_sumsq = calloc((size_t)ref_n, sizeof(product_sumsq[0]));
		product_max = calloc((size_t)ref_n, sizeof(product_max[0]));
		product_nonzero = calloc((size_t)ref_n, sizeof(product_nonzero[0]));
		if (!product_sum || !product_sumsq || !product_max || !product_nonzero)
			err(EXIT_FAILURE, "%s(): OOM reliable abundance product stats", __func__);
	}

	const long double inv_coverage_scale = 1.0L / (long double)coverage_scale;
	const uint64_t gidmask_local = (1ULL << GID_NBITS) - 1ULL;
	if (ani_readwise_ctx_filter_uses_product(filter_model)) {
		for (size_t i = 0; i < index_n; ) {
			const uint64_t ctxgid = index[i].ctxgid;
			const uint32_t gid = (uint32_t)(ctxgid & gidmask_local);
			const size_t begin = i;
			do { ++i; } while (i < index_n && index[i].ctxgid == ctxgid);
			if (gid >= ref_n || (ignoreconflict && i - begin > 1) ||
				(marker_only && !ani_ctxgid_group_marker_unique64(index, index_n, begin, i)))
				continue;
			const uint32_t packed = ref_ctx_cov[begin];
			const uint32_t coverage_raw = ani_ref_covdiff_coverage(packed);
			const uint32_t code = packed >> ANI_REFCOV_DIFF_SHIFT;
			if (!code || !coverage_raw)
				continue;
			const uint32_t diff = ani_ref_covdiff_decode_diff(code);
			if (!diff)
				continue;
			const double cov = (double)coverage_raw / coverage_scale;
			const double product = (double)diff * cov;
			if (!isfinite(product) || product <= 0.0)
				continue;
			product_sum[gid] += (long double)product;
			product_sumsq[gid] += (long double)product * (long double)product;
			if (product > product_max[gid])
				product_max[gid] = product;
			product_nonzero[gid]++;
		}
	}
	double *product_topfrac_thresholds = NULL;
	if (filter_model == ANI_READWISE_CTX_FILTER_PRODUCT_TOPFRAC) {
		product_topfrac_thresholds = ani_readwise_product_topfrac_thresholds(
			index, index_n, ref_n, ignoreconflict, ref_ctx_cov,
			coverage_scale, fake_threshold, marker_only);
	}
	ani_readwise_product_topfrac_median_t *product_topfrac_medians = NULL;
	if (filter_model == ANI_READWISE_CTX_FILTER_PRODUCT_TOPFRAC_MEDIAN ||
		filter_model == ANI_READWISE_CTX_FILTER_PRODUCT1_TOPFRAC_MEDIAN) {
		product_topfrac_medians = ani_readwise_product_topfrac_median_stats(
			index, index_n, ref_n, ignoreconflict, ref_ctx_cov,
			coverage_scale, fake_threshold,
			filter_model == ANI_READWISE_CTX_FILTER_PRODUCT1_TOPFRAC_MEDIAN,
			marker_only);
	}

	for (size_t i = 0; i < index_n; ) {
		const uint64_t ctxgid = index[i].ctxgid;
		const uint32_t gid = (uint32_t)(ctxgid & gidmask_local);
		const size_t begin = i;
		do { ++i; } while (i < index_n && index[i].ctxgid == ctxgid);
		if (gid >= ref_n || (ignoreconflict && i - begin > 1) ||
			(marker_only && !ani_ctxgid_group_marker_unique64(index, index_n, begin, i)))
			continue;
		const uint32_t packed = ref_ctx_cov[begin];
		const uint32_t coverage_raw = ani_ref_covdiff_coverage(packed);
		const uint32_t code = packed >> ANI_REFCOV_DIFF_SHIFT;
		if (!code || !coverage_raw)
			continue;
		const uint32_t raw_diff = ani_ref_covdiff_decode_diff(code);
		uint32_t diff = raw_diff;
		double cov = (double)coverage_raw / coverage_scale;

		double product_lambda = 0.0;
		double product_nb_mean = 0.0;
		double product_nb_variance = 0.0;
		uint32_t product_nb_nonzero = 0u;
		if (ani_readwise_ctx_filter_uses_product(filter_model) &&
			ref_ctx_total[gid]) {
			long double robust_sum = product_sum[gid];
			long double robust_sumsq = product_sumsq[gid];
			uint32_t robust_n = product_nonzero[gid];
			if (product_nonzero[gid] > 1u)
				robust_sum -= (long double)product_max[gid];
			if (product_nonzero[gid] > 1u) {
				robust_sumsq -= (long double)product_max[gid] * (long double)product_max[gid];
				robust_n--;
			}
			if (robust_sum < 0.0L)
				robust_sum = 0.0L;
			if (robust_sumsq < 0.0L)
				robust_sumsq = 0.0L;
			product_lambda = (double)(robust_sum / (long double)ref_ctx_total[gid]);
			if (robust_n > 0u) {
				product_nb_nonzero = robust_n;
				product_nb_mean = (double)(robust_sum / (long double)robust_n);
				if (robust_n > 1u) {
					const long double mean_ld = robust_sum / (long double)robust_n;
					long double ss = robust_sumsq - robust_sum * mean_ld;
					if (ss < 0.0L)
						ss = 0.0L;
					product_nb_variance = (double)(ss / (long double)(robust_n - 1u));
				} else {
					product_nb_variance = product_nb_mean;
				}
			}
		}
		if (filter_model == ANI_READWISE_CTX_FILTER_PRODUCT_TOPFRAC_MEDIAN ||
			filter_model == ANI_READWISE_CTX_FILTER_PRODUCT1_TOPFRAC_MEDIAN)
			(void)ani_readwise_apply_product_topfrac_median_ctx(
				&diff, &cov,
				product_topfrac_medians ? &product_topfrac_medians[gid] : NULL,
				filter_model == ANI_READWISE_CTX_FILTER_PRODUCT1_TOPFRAC_MEDIAN);
		const bool reject =
			(filter_model == ANI_READWISE_CTX_FILTER_POISSON_DIFF &&
			 ani_readwise_reject_poisson_diff_ctx(
				 raw_diff, cov, &abundance_stats[gid], fake_threshold)) ||
			(filter_model == ANI_READWISE_CTX_FILTER_POISSON_DEPTH &&
			 ani_readwise_reject_poisson_depth_ctx(
				 raw_diff, cov, &abundance_stats[gid],
				 ref_ctx_total ? ref_ctx_total[gid] : 0u,
				 fake_threshold)) ||
			(filter_model == ANI_READWISE_CTX_FILTER_POISSON_PRODUCT &&
			 ani_readwise_reject_poisson_product_ctx(
				 raw_diff, cov, product_lambda,
				 ref_ctx_total ? ref_ctx_total[gid] : 0u,
				 fake_threshold)) ||
			(filter_model == ANI_READWISE_CTX_FILTER_PRODUCT_NB &&
			 ani_readwise_reject_product_nb_ctx(
				 raw_diff, cov, product_nb_mean, product_nb_variance,
				 product_nb_nonzero,
				 ref_ctx_total ? ref_ctx_total[gid] : 0u,
				 fake_threshold)) ||
			(filter_model == ANI_READWISE_CTX_FILTER_PRODUCT_TOPFRAC &&
			 ani_readwise_reject_product_topfrac_ctx(
				 raw_diff, cov,
				 product_topfrac_thresholds ? product_topfrac_thresholds[gid] : NAN));
		if (reject)
			continue;

		const long double cov_ld = (long double)cov;
		sum[gid] += cov_ld;
		sumsq[gid] += cov_ld * cov_ld;
		hit_sum[gid] += cov_ld;
		hit_sumsq[gid] += cov_ld * cov_ld;
		hit_ctx[gid]++;
		if (ref_ctx_hit_weight && ref_ctx_hit_weight[begin] > 0) {
			long double hit = (long double)ref_ctx_hit_weight[begin] * inv_coverage_scale;
			if (hit > 1.0L)
				hit = 1.0L;
			hit_weight[gid] += hit;
		} else {
			hit_weight[gid] += 1.0L;
		}
	}

	for (uint32_t rn = 0; rn < ref_n; ++rn) {
		const uint32_t total = ref_ctx_total[rn];
		if (!total)
			continue;
		const long double denom = (long double)total;
		const long double breadth = hit_weight[rn] > 0.0L
										? hit_weight[rn] / denom
										: 0.0L;
		const long double mean = sum[rn] / denom;
		long double hit_mean = 0.0L;
		long double hit_var = 0.0L;
		if (hit_ctx[rn]) {
			hit_mean = hit_sum[rn] / (long double)hit_ctx[rn];
			hit_var = hit_sumsq[rn] / (long double)hit_ctx[rn] - hit_mean * hit_mean;
			if (hit_var < 0.0L && hit_var > -1e-12L)
				hit_var = 0.0L;
			if (hit_var < 0.0L)
				hit_var = 0.0L;
		}
		stats[rn].ref_breadth = (double)breadth;
		stats[rn].ref_mean_depth = (double)mean;
		stats[rn].ref_hit_ctx = hit_ctx[rn];
		stats[rn].ref_hit_mean_depth = (double)hit_mean;
		stats[rn].ref_hit_median_depth =
			product_topfrac_medians &&
					isfinite(product_topfrac_medians[rn].median_cov) &&
					product_topfrac_medians[rn].median_cov > 0.0
				? product_topfrac_medians[rn].median_cov
				: 0.0;
		stats[rn].ref_hit_depth_variance = (double)hit_var;
		stats[rn].ref_zip_af = ani_zip_corrected_af(
			stats[rn].ref_breadth, stats[rn].ref_mean_depth);
	}

	free(sum);
	free(sumsq);
	free(hit_weight);
	free(hit_sum);
	free(hit_sumsq);
	free(hit_ctx);
	free(product_sum);
	free(product_sumsq);
	free(product_max);
	free(product_nonzero);
	free(product_topfrac_thresholds);
	free(product_topfrac_medians);
	return stats;
}

static ani_readwise_reliable_abundance_t *ani_reliable_abundance_from_ref_covdiff128(
	const ctxgidobj128_t *index,
	size_t index_n,
	uint32_t ref_n,
	bool ignoreconflict,
	const uint32_t *ref_ctx_cov,
	const uint16_t *ref_ctx_hit_weight,
	const ani_readwise_abundance_t *abundance_stats,
	const uint32_t *ref_ctx_total,
	ani_readwise_ctx_filter_model_t filter_model,
	double coverage_scale,
	double fake_threshold,
	bool marker_only)
{
	if (!ref_ctx_cov || !abundance_stats || !ref_ctx_total)
		return NULL;
	if (!isfinite(coverage_scale) || coverage_scale <= 0.0)
		coverage_scale = 1.0;
	if (!isfinite(fake_threshold) || fake_threshold <= 0.0)
		fake_threshold = MINCO_DEFAULT_READWISE_FAKE_CTX_THRESHOLD;

	ani_readwise_reliable_abundance_t *stats =
		calloc((size_t)ref_n, sizeof(stats[0]));
	long double *sum = calloc((size_t)ref_n, sizeof(sum[0]));
	long double *sumsq = calloc((size_t)ref_n, sizeof(sumsq[0]));
	long double *hit_weight = calloc((size_t)ref_n, sizeof(hit_weight[0]));
	long double *hit_sum = calloc((size_t)ref_n, sizeof(hit_sum[0]));
	long double *hit_sumsq = calloc((size_t)ref_n, sizeof(hit_sumsq[0]));
	uint64_t *hit_ctx = calloc((size_t)ref_n, sizeof(hit_ctx[0]));
	if (!stats || !sum || !sumsq || !hit_weight ||
		!hit_sum || !hit_sumsq || !hit_ctx)
		err(EXIT_FAILURE, "%s(): OOM reliable abundance stats", __func__);

	long double *product_sum = NULL;
	long double *product_sumsq = NULL;
	double *product_max = NULL;
	uint32_t *product_nonzero = NULL;
	if (ani_readwise_ctx_filter_uses_product(filter_model)) {
		product_sum = calloc((size_t)ref_n, sizeof(product_sum[0]));
		product_sumsq = calloc((size_t)ref_n, sizeof(product_sumsq[0]));
		product_max = calloc((size_t)ref_n, sizeof(product_max[0]));
		product_nonzero = calloc((size_t)ref_n, sizeof(product_nonzero[0]));
		if (!product_sum || !product_sumsq || !product_max || !product_nonzero)
			err(EXIT_FAILURE, "%s(): OOM reliable abundance product stats", __func__);
	}

	const long double inv_coverage_scale = 1.0L / (long double)coverage_scale;
	if (ani_readwise_ctx_filter_uses_product(filter_model)) {
		for (size_t i = 0; i < index_n; ) {
			const uint64_t ctx = index[i].ctx;
			const uint32_t gid = index[i].gid;
			const size_t begin = i;
			do { ++i; } while (i < index_n &&
								index[i].ctx == ctx &&
								index[i].gid == gid);
			if (gid >= ref_n || (ignoreconflict && i - begin > 1) ||
				(marker_only && !ani_ctxgid_group_marker_unique128(index, index_n, begin, i)))
				continue;
			const uint32_t packed = ref_ctx_cov[begin];
			const uint32_t coverage_raw = ani_ref_covdiff_coverage(packed);
			const uint32_t code = packed >> ANI_REFCOV_DIFF_SHIFT;
			if (!code || !coverage_raw)
				continue;
			const uint32_t diff = ani_ref_covdiff_decode_diff(code);
			if (!diff)
				continue;
			const double cov = (double)coverage_raw / coverage_scale;
			const double product = (double)diff * cov;
			if (!isfinite(product) || product <= 0.0)
				continue;
			product_sum[gid] += (long double)product;
			product_sumsq[gid] += (long double)product * (long double)product;
			if (product > product_max[gid])
				product_max[gid] = product;
			product_nonzero[gid]++;
		}
	}
	double *product_topfrac_thresholds = NULL;
	if (filter_model == ANI_READWISE_CTX_FILTER_PRODUCT_TOPFRAC) {
		product_topfrac_thresholds = ani_readwise_product_topfrac_thresholds128(
			index, index_n, ref_n, ignoreconflict, ref_ctx_cov,
			coverage_scale, fake_threshold, marker_only);
	}
	ani_readwise_product_topfrac_median_t *product_topfrac_medians = NULL;
	if (filter_model == ANI_READWISE_CTX_FILTER_PRODUCT_TOPFRAC_MEDIAN ||
		filter_model == ANI_READWISE_CTX_FILTER_PRODUCT1_TOPFRAC_MEDIAN) {
		product_topfrac_medians = ani_readwise_product_topfrac_median_stats128(
			index, index_n, ref_n, ignoreconflict, ref_ctx_cov,
			coverage_scale, fake_threshold,
			filter_model == ANI_READWISE_CTX_FILTER_PRODUCT1_TOPFRAC_MEDIAN,
			marker_only);
	}

	for (size_t i = 0; i < index_n; ) {
		const uint64_t ctx = index[i].ctx;
		const uint32_t gid = index[i].gid;
		const size_t begin = i;
		do { ++i; } while (i < index_n &&
							index[i].ctx == ctx &&
							index[i].gid == gid);
		if (gid >= ref_n || (ignoreconflict && i - begin > 1) ||
			(marker_only && !ani_ctxgid_group_marker_unique128(index, index_n, begin, i)))
			continue;
		const uint32_t packed = ref_ctx_cov[begin];
		const uint32_t coverage_raw = ani_ref_covdiff_coverage(packed);
		const uint32_t code = packed >> ANI_REFCOV_DIFF_SHIFT;
		if (!code || !coverage_raw)
			continue;
		const uint32_t raw_diff = ani_ref_covdiff_decode_diff(code);
		uint32_t diff = raw_diff;
		double cov = (double)coverage_raw / coverage_scale;

		double product_lambda = 0.0;
		double product_nb_mean = 0.0;
		double product_nb_variance = 0.0;
		uint32_t product_nb_nonzero = 0u;
		if (ani_readwise_ctx_filter_uses_product(filter_model) &&
			ref_ctx_total[gid]) {
			long double robust_sum = product_sum[gid];
			long double robust_sumsq = product_sumsq[gid];
			uint32_t robust_n = product_nonzero[gid];
			if (product_nonzero[gid] > 1u)
				robust_sum -= (long double)product_max[gid];
			if (product_nonzero[gid] > 1u) {
				robust_sumsq -= (long double)product_max[gid] * (long double)product_max[gid];
				robust_n--;
			}
			if (robust_sum < 0.0L)
				robust_sum = 0.0L;
			if (robust_sumsq < 0.0L)
				robust_sumsq = 0.0L;
			product_lambda = (double)(robust_sum / (long double)ref_ctx_total[gid]);
			if (robust_n > 0u) {
				product_nb_nonzero = robust_n;
				product_nb_mean = (double)(robust_sum / (long double)robust_n);
				if (robust_n > 1u) {
					const long double mean_ld = robust_sum / (long double)robust_n;
					long double ss = robust_sumsq - robust_sum * mean_ld;
					if (ss < 0.0L)
						ss = 0.0L;
					product_nb_variance = (double)(ss / (long double)(robust_n - 1u));
				} else {
					product_nb_variance = product_nb_mean;
				}
			}
		}
		if (filter_model == ANI_READWISE_CTX_FILTER_PRODUCT_TOPFRAC_MEDIAN ||
			filter_model == ANI_READWISE_CTX_FILTER_PRODUCT1_TOPFRAC_MEDIAN)
			(void)ani_readwise_apply_product_topfrac_median_ctx(
				&diff, &cov,
				product_topfrac_medians ? &product_topfrac_medians[gid] : NULL,
				filter_model == ANI_READWISE_CTX_FILTER_PRODUCT1_TOPFRAC_MEDIAN);
		const bool reject =
			(filter_model == ANI_READWISE_CTX_FILTER_POISSON_DIFF &&
			 ani_readwise_reject_poisson_diff_ctx(
				 raw_diff, cov, &abundance_stats[gid], fake_threshold)) ||
			(filter_model == ANI_READWISE_CTX_FILTER_POISSON_DEPTH &&
			 ani_readwise_reject_poisson_depth_ctx(
				 raw_diff, cov, &abundance_stats[gid],
				 ref_ctx_total ? ref_ctx_total[gid] : 0u,
				 fake_threshold)) ||
			(filter_model == ANI_READWISE_CTX_FILTER_POISSON_PRODUCT &&
			 ani_readwise_reject_poisson_product_ctx(
				 raw_diff, cov, product_lambda,
				 ref_ctx_total ? ref_ctx_total[gid] : 0u,
				 fake_threshold)) ||
			(filter_model == ANI_READWISE_CTX_FILTER_PRODUCT_NB &&
			 ani_readwise_reject_product_nb_ctx(
				 raw_diff, cov, product_nb_mean, product_nb_variance,
				 product_nb_nonzero,
				 ref_ctx_total ? ref_ctx_total[gid] : 0u,
				 fake_threshold)) ||
			(filter_model == ANI_READWISE_CTX_FILTER_PRODUCT_TOPFRAC &&
			 ani_readwise_reject_product_topfrac_ctx(
				 raw_diff, cov,
				 product_topfrac_thresholds ? product_topfrac_thresholds[gid] : NAN));
		if (reject)
			continue;

		const long double cov_ld = (long double)cov;
		sum[gid] += cov_ld;
		sumsq[gid] += cov_ld * cov_ld;
		hit_sum[gid] += cov_ld;
		hit_sumsq[gid] += cov_ld * cov_ld;
		hit_ctx[gid]++;
		if (ref_ctx_hit_weight && ref_ctx_hit_weight[begin] > 0) {
			long double hit = (long double)ref_ctx_hit_weight[begin] * inv_coverage_scale;
			if (hit > 1.0L)
				hit = 1.0L;
			hit_weight[gid] += hit;
		} else {
			hit_weight[gid] += 1.0L;
		}
	}

	for (uint32_t rn = 0; rn < ref_n; ++rn) {
		const uint32_t total = ref_ctx_total[rn];
		if (!total)
			continue;
		const long double denom = (long double)total;
		const long double breadth = hit_weight[rn] > 0.0L
										? hit_weight[rn] / denom
										: 0.0L;
		const long double mean = sum[rn] / denom;
		long double hit_mean = 0.0L;
		long double hit_var = 0.0L;
		if (hit_ctx[rn]) {
			hit_mean = hit_sum[rn] / (long double)hit_ctx[rn];
			hit_var = hit_sumsq[rn] / (long double)hit_ctx[rn] - hit_mean * hit_mean;
			if (hit_var < 0.0L && hit_var > -1e-12L)
				hit_var = 0.0L;
			if (hit_var < 0.0L)
				hit_var = 0.0L;
		}
		stats[rn].ref_breadth = (double)breadth;
		stats[rn].ref_mean_depth = (double)mean;
		stats[rn].ref_hit_ctx = hit_ctx[rn];
		stats[rn].ref_hit_mean_depth = (double)hit_mean;
		stats[rn].ref_hit_median_depth =
			product_topfrac_medians &&
					isfinite(product_topfrac_medians[rn].median_cov) &&
					product_topfrac_medians[rn].median_cov > 0.0
				? product_topfrac_medians[rn].median_cov
				: 0.0;
		stats[rn].ref_hit_depth_variance = (double)hit_var;
		stats[rn].ref_zip_af = ani_zip_corrected_af(
			stats[rn].ref_breadth, stats[rn].ref_mean_depth);
	}

	free(sum);
	free(sumsq);
	free(hit_weight);
	free(hit_sum);
	free(hit_sumsq);
	free(hit_ctx);
	free(product_sum);
	free(product_sumsq);
	free(product_max);
	free(product_nonzero);
	free(product_topfrac_thresholds);
	free(product_topfrac_medians);
	return stats;
}

static void ani_write_readwise_unique_sidecar(
	ani_opt_t *ani_opt,
	const char *query_path,
	char (*refname)[PATHLEN],
	char (*refanno)[PATHLEN],
	uint32_t ref_n,
	bool ref_uses_ctxobj96,
	const ctxgidobj_t *index,
	const ctxgidobj128_t *index96,
	size_t index_n,
	const uint32_t *ref_ctx_total,
	const ani_readwise_acc_t *unique_acc,
	const uint32_t *unique_ref_ctx_cov,
	uint64_t total_reads,
	uint64_t total_density_blocks,
	bool density_block_mode,
	uint32_t density_block_ctx,
	const infile_meta_t *ref_infile_meta,
	const uint8_t *ref_domain_profiles,
	bool auto_readwise_abundance_report)
{
	if (!ani_opt || ani_opt->readwise_unique_out[0] == '\0')
		return;
	if (!unique_acc || !unique_ref_ctx_cov || !ref_ctx_total)
		return;

	FILE *fp = strcmp(ani_opt->readwise_unique_out, "-") == 0
				   ? stdout
				   : fopen(ani_opt->readwise_unique_out, "w");
	if (!fp)
		err(errno, "%s", ani_opt->readwise_unique_out);

	ani_opt_t side_opt = *ani_opt;
	side_opt.readwise_dual_evidence = false;
	side_opt.readwise_assign_mode = ANI_READWISE_ASSIGN_BEST_DIFF_UNIQUE;
	print_ani_detail_header(fp, &side_opt, true);

	ani_readwise_abundance_t *abundance_stats = ref_uses_ctxobj96
		? ani_depth_abundance_from_ref_coverage128(
			  index96, index_n, ref_n, side_opt.ignoreconflict,
			  unique_ref_ctx_cov, NULL, ref_ctx_total, unique_acc, 1.0, false)
		: ani_depth_abundance_from_ref_coverage(
			  index, index_n, ref_n, side_opt.ignoreconflict,
			  unique_ref_ctx_cov, NULL, ref_ctx_total, unique_acc, 1.0, false);
	ani_readwise_acc_t *raw_feature_acc = ref_uses_ctxobj96
		? ani_unique_best_features_from_ref_covdiff128(
			  index96, index_n, ref_n, side_opt.ignoreconflict,
			  unique_ref_ctx_cov, false)
		: ani_unique_best_features_from_ref_covdiff(
			  index, index_n, ref_n, side_opt.ignoreconflict,
			  unique_ref_ctx_cov, false);
	if (!abundance_stats || !raw_feature_acc)
		err(EXIT_FAILURE, "%s(): OOM unique sidecar features", __func__);

	kv_ani_row_t survivors;
	kv_init(survivors);
	for (uint32_t rn = 0; rn < ref_n; ++rn) {
		const ani_readwise_acc_t *a = &unique_acc[rn];
		const ani_readwise_acc_t *fa = &raw_feature_acc[rn];
		const uint64_t unique_overlap = a->ref_ctx_hit;
		if (unique_overlap < (uint64_t)side_opt.ctxcut)
			continue;
		if (!fa->XnY_ctx || !ref_ctx_total[rn])
			continue;
		if (fa->XnY_ctx < (uint64_t)side_opt.ctxcut)
			continue;
		const double ref_af = (double)a->ref_ctx_hit / (double)ref_ctx_total[rn];
		ani_density_af_t density_af = {
			.qry = ref_af,
			.ref = ref_af,
			.available = 1,
		};
		if (!ani_report_af_pass(&side_opt, density_af.qry, density_af.ref))
			continue;

		ani_features_t f = {
			.XnY_ctx = ani_clamp_u64_to_u32(fa->XnY_ctx),
			.X_ctx = ref_ctx_total[rn],
			.N_diff_obj_section = ani_clamp_u64_to_u32(fa->N_diff_obj_section),
			.N_mut2_ctx = ani_clamp_u64_to_u32(fa->N_mut2_ctx),
			.N_diff_obj = ani_clamp_u64_to_u32(fa->N_diff_obj),
		};
		ani_row_t row = make_selected_output_row(
			rn, &f, &side_opt,
			ref_ctx_total[rn],
			ref_ctx_total[rn],
			ref_af, ref_af,
			ref_af, ref_af,
			density_af,
			NULL,
			ref_infile_meta ? &ref_infile_meta[rn] : NULL);
		row.domain_profile = ref_domain_profiles
								  ? ref_domain_profiles[rn]
								  : MINCO_DOMAIN_PROFILE_DEFAULT;
		row.XnY_ctx = ani_readwise_scaled_round_int(fa->XnY_ctx, 1u);
		row.N_diff_obj = ani_readwise_scaled_round_int(fa->N_diff_obj, 1u);
		row.N_diff_obj_section =
			ani_readwise_scaled_round_int(fa->N_diff_obj_section, 1u);
		row.N_mut2_ctx = ani_readwise_scaled_round_int(fa->N_mut2_ctx, 1u);
		row.readwise_total_reads = total_reads;
		row.readwise_reads_with_ctx_match = a->reads_with_ctx_match;
		row.readwise_unique_query_ctx = 0;
		row.readwise_unique_query_ctx_hit = a->ref_ctx_hit;
		row.readwise_unique_ref_ctx_hit = a->ref_ctx_hit;
		row.readwise_ref_ctx_total = ref_ctx_total[rn];
		row.readwise_density_block_ctx = density_block_mode ? density_block_ctx : 0u;
		row.readwise_total_density_blocks = total_density_blocks;
		row.readwise_blocks_with_ctx_match = a->blocks_with_ctx_match;
		row.readwise_raw_xny_ctx = fa->XnY_ctx;

		row.abundance_ref_breadth = abundance_stats[rn].ref_breadth;
		row.abundance_ref_mean_depth = abundance_stats[rn].ref_mean_depth;
		row.abundance_ref_hit_mean_depth = abundance_stats[rn].ref_hit_mean_depth;
		row.abundance_ref_depth_variance = abundance_stats[rn].ref_depth_variance;
		row.abundance_ref_depth_cv = abundance_stats[rn].ref_depth_cv;
		row.abundance_ref_zero_fraction = abundance_stats[rn].ref_zero_fraction;
		row.abundance_relative_depth = abundance_stats[rn].relative_depth;
		row.abundance_ref_zip_af = ani_zip_corrected_af(
			row.abundance_ref_breadth, row.abundance_ref_mean_depth);
		row.abundance_ref_zip_aaf_ani =
			aaf_ani_from_containment(row.abundance_ref_zip_af);
		row.reliable_ref_breadth = row.abundance_ref_breadth;
		row.reliable_ref_mean_depth = row.abundance_ref_mean_depth;
		row.reliable_ref_hit_ctx = a->ref_ctx_hit;
		row.reliable_ref_hit_mean_depth = row.abundance_ref_hit_mean_depth;
		row.reliable_ref_hit_median_depth = 0.0;
		row.reliable_ref_hit_depth_variance = row.abundance_ref_depth_variance;
		row.reliable_ref_zip_af = row.abundance_ref_zip_af;
		row.abundance_effective_depth =
			ani_readwise_effective_abundance_depth(&row);

		const uint32_t unique_overlap_u32 = ani_clamp_u64_to_u32(unique_overlap);
		row.mash_dist = get_mashD(Bitslen.ctx / 2, ref_ctx_total[rn],
								  ref_ctx_total[rn], unique_overlap_u32);
		row.aaf_dist = get_aafD(Bitslen.ctx / 2, ref_ctx_total[rn],
								ref_ctx_total[rn], unique_overlap_u32);
		finalize_row_selected_metric(&row, &side_opt);
		if (row.selected_ani < side_opt.anicut)
			continue;
		kv_push(ani_row_t, survivors, row);
	}

	if (auto_readwise_abundance_report) {
		keep_readwise_default_calls(&survivors);
		if (kv_size(survivors))
			qsort(&kv_A(survivors, 0), kv_size(survivors),
				  sizeof(ani_row_t), cmp_readwise_default_call_desc);
	} else if (kv_size(survivors)) {
		qsort(&kv_A(survivors, 0), kv_size(survivors), sizeof(ani_row_t),
			  cmp_ani_desc);
	}
	size_t out_n = kv_size(survivors);
	if (side_opt.ntop > 0 && (size_t)side_opt.ntop < out_n)
		out_n = (size_t)side_opt.ntop;
	normalize_readwise_abundance_depth(&survivors, out_n);
	for (size_t i = 0; i < out_n; ++i) {
		const ani_row_t *r = &kv_A(survivors, i);
		print_unified_detail_row(fp, &side_opt, query_path, refname[r->rn],
								 r, annotation_at(refanno, r->rn));
	}
	if (fp != stdout)
		fclose(fp);
	kv_destroy(survivors);
	free(abundance_stats);
	free(raw_feature_acc);
	fprintf(stderr, "minco readwise: wrote best-diff-unique sidecar %s\n",
			ani_opt->readwise_unique_out);
}

static void ani_write_readwise_exact_split_sidecar(
	ani_opt_t *ani_opt,
	const char *query_path,
	char (*refname)[PATHLEN],
	char (*refanno)[PATHLEN],
	uint32_t ref_n,
	bool ref_uses_ctxobj96,
	const ctxgidobj_t *index,
	const ctxgidobj128_t *index96,
	size_t index_n,
	const uint32_t *ref_ctx_total,
	const ani_readwise_acc_t *exact_acc,
	const uint32_t *exact_ref_ctx_cov,
	const uint16_t *exact_ref_ctx_hit_weight,
	uint64_t total_reads,
	uint64_t exact_total_density_blocks,
	const infile_meta_t *ref_infile_meta,
	const uint8_t *ref_domain_profiles,
	bool auto_readwise_abundance_report)
{
	if (!ani_opt || ani_opt->readwise_exact_split_out[0] == '\0')
		return;
	if (!exact_acc || !exact_ref_ctx_cov || !exact_ref_ctx_hit_weight || !ref_ctx_total)
		return;

	FILE *fp = strcmp(ani_opt->readwise_exact_split_out, "-") == 0
				   ? stdout
				   : fopen(ani_opt->readwise_exact_split_out, "w");
	if (!fp)
		err(errno, "%s", ani_opt->readwise_exact_split_out);

	ani_opt_t side_opt = *ani_opt;
	side_opt.readwise_dual_evidence = false;
	side_opt.readwise_assign_mode = ANI_READWISE_ASSIGN_BEST_DIFF_SPLIT;
	print_ani_detail_header(fp, &side_opt, true);

	const double coverage_scale = (double)ANI_REFCOV_SPLIT_SCALE;
	ani_readwise_abundance_t *abundance_stats = ref_uses_ctxobj96
		? ani_depth_abundance_from_ref_coverage128(
			  index96, index_n, ref_n, side_opt.ignoreconflict,
			  exact_ref_ctx_cov, exact_ref_ctx_hit_weight, ref_ctx_total,
			  exact_acc, coverage_scale, false)
		: ani_depth_abundance_from_ref_coverage(
			  index, index_n, ref_n, side_opt.ignoreconflict,
			  exact_ref_ctx_cov, exact_ref_ctx_hit_weight, ref_ctx_total,
			  exact_acc, coverage_scale, false);
	ani_readwise_acc_t *raw_feature_acc = ref_uses_ctxobj96
		? ani_unique_best_features_from_ref_covdiff128(
			  index96, index_n, ref_n, side_opt.ignoreconflict,
			  exact_ref_ctx_cov, false)
		: ani_unique_best_features_from_ref_covdiff(
			  index, index_n, ref_n, side_opt.ignoreconflict,
			  exact_ref_ctx_cov, false);
	if (!abundance_stats || !raw_feature_acc)
		err(EXIT_FAILURE, "%s(): OOM exact split sidecar features", __func__);

	kv_ani_row_t survivors;
	kv_init(survivors);
	for (uint32_t rn = 0; rn < ref_n; ++rn) {
		const ani_readwise_acc_t *a = &exact_acc[rn];
		const ani_readwise_acc_t *fa = &raw_feature_acc[rn];
		const uint64_t unique_overlap = a->ref_ctx_hit;
		if (unique_overlap < (uint64_t)side_opt.ctxcut)
			continue;
		if (!fa->XnY_ctx || !ref_ctx_total[rn])
			continue;
		if (fa->XnY_ctx < (uint64_t)side_opt.ctxcut)
			continue;
		const double ref_af = (double)a->ref_ctx_hit / (double)ref_ctx_total[rn];
		ani_density_af_t density_af = {
			.qry = ref_af,
			.ref = ref_af,
			.available = 1,
		};
		if (!ani_report_af_pass(&side_opt, density_af.qry, density_af.ref))
			continue;

		ani_features_t f = {
			.XnY_ctx = ani_clamp_u64_to_u32(fa->XnY_ctx),
			.X_ctx = ref_ctx_total[rn],
			.N_diff_obj_section = ani_clamp_u64_to_u32(fa->N_diff_obj_section),
			.N_mut2_ctx = ani_clamp_u64_to_u32(fa->N_mut2_ctx),
			.N_diff_obj = ani_clamp_u64_to_u32(fa->N_diff_obj),
		};
		ani_row_t row = make_selected_output_row(
			rn, &f, &side_opt,
			ref_ctx_total[rn],
			ref_ctx_total[rn],
			ref_af, ref_af,
			ref_af, ref_af,
			density_af,
			NULL,
			ref_infile_meta ? &ref_infile_meta[rn] : NULL);
		row.domain_profile = ref_domain_profiles
								  ? ref_domain_profiles[rn]
								  : MINCO_DOMAIN_PROFILE_DEFAULT;
		row.XnY_ctx = ani_readwise_scaled_round_int(fa->XnY_ctx, 1u);
		row.N_diff_obj = ani_readwise_scaled_round_int(fa->N_diff_obj, 1u);
		row.N_diff_obj_section =
			ani_readwise_scaled_round_int(fa->N_diff_obj_section, 1u);
		row.N_mut2_ctx = ani_readwise_scaled_round_int(fa->N_mut2_ctx, 1u);
		row.readwise_total_reads = total_reads;
		row.readwise_reads_with_ctx_match = a->reads_with_ctx_match;
		row.readwise_unique_query_ctx = 0;
		row.readwise_unique_query_ctx_hit = a->ref_ctx_hit;
		row.readwise_unique_ref_ctx_hit = a->ref_ctx_hit;
		row.readwise_ref_ctx_total = ref_ctx_total[rn];
		row.readwise_density_block_ctx = 0u;
		row.readwise_total_density_blocks = exact_total_density_blocks;
		row.readwise_blocks_with_ctx_match = a->blocks_with_ctx_match;
		row.readwise_raw_xny_ctx = fa->XnY_ctx;

		row.abundance_ref_breadth = abundance_stats[rn].ref_breadth;
		row.abundance_ref_mean_depth = abundance_stats[rn].ref_mean_depth;
		row.abundance_ref_hit_mean_depth = abundance_stats[rn].ref_hit_mean_depth;
		row.abundance_ref_depth_variance = abundance_stats[rn].ref_depth_variance;
		row.abundance_ref_depth_cv = abundance_stats[rn].ref_depth_cv;
		row.abundance_ref_zero_fraction = abundance_stats[rn].ref_zero_fraction;
		row.abundance_relative_depth = abundance_stats[rn].relative_depth;
		row.abundance_ref_zip_af = ani_zip_corrected_af(
			row.abundance_ref_breadth, row.abundance_ref_mean_depth);
		row.abundance_ref_zip_aaf_ani =
			aaf_ani_from_containment(row.abundance_ref_zip_af);
		row.reliable_ref_breadth = row.abundance_ref_breadth;
		row.reliable_ref_mean_depth = row.abundance_ref_mean_depth;
		row.reliable_ref_hit_ctx = a->ref_ctx_hit;
		row.reliable_ref_hit_mean_depth = row.abundance_ref_hit_mean_depth;
		row.reliable_ref_hit_median_depth = 0.0;
		row.reliable_ref_hit_depth_variance = row.abundance_ref_depth_variance;
		row.reliable_ref_zip_af = row.abundance_ref_zip_af;
		row.abundance_effective_depth =
			ani_readwise_effective_abundance_depth(&row);

		const uint32_t unique_overlap_u32 = ani_clamp_u64_to_u32(unique_overlap);
		row.mash_dist = get_mashD(Bitslen.ctx / 2, ref_ctx_total[rn],
								  ref_ctx_total[rn], unique_overlap_u32);
		row.aaf_dist = get_aafD(Bitslen.ctx / 2, ref_ctx_total[rn],
								ref_ctx_total[rn], unique_overlap_u32);
		finalize_row_selected_metric(&row, &side_opt);
		if (row.selected_ani < side_opt.anicut)
			continue;
		kv_push(ani_row_t, survivors, row);
	}

	if (auto_readwise_abundance_report) {
		keep_readwise_default_calls(&survivors);
		if (kv_size(survivors))
			qsort(&kv_A(survivors, 0), kv_size(survivors),
				  sizeof(ani_row_t), cmp_readwise_default_call_desc);
	} else if (kv_size(survivors)) {
		qsort(&kv_A(survivors, 0), kv_size(survivors), sizeof(ani_row_t),
			  cmp_ani_desc);
	}
	size_t out_n = kv_size(survivors);
	if (side_opt.ntop > 0 && (size_t)side_opt.ntop < out_n)
		out_n = (size_t)side_opt.ntop;
	normalize_readwise_abundance_depth(&survivors, out_n);
	for (size_t i = 0; i < out_n; ++i) {
		const ani_row_t *r = &kv_A(survivors, i);
		print_unified_detail_row(fp, &side_opt, query_path, refname[r->rn],
								 r, annotation_at(refanno, r->rn));
	}
	if (fp != stdout)
		fclose(fp);
	kv_destroy(survivors);
	free(abundance_stats);
	free(raw_feature_acc);
	fprintf(stderr, "minco readwise: wrote exact best-diff-split sidecar %s\n",
			ani_opt->readwise_exact_split_out);
}

int stream_fastq_query_readwise_density_ani(ani_opt_t *ani_opt, const char *query_path,
											uint64_t density_threshold)
{
	if (!ani_opt || !query_path)
		errx(EXIT_FAILURE, "%s(): missing ANI options or query path", __func__);
	if (ani_opt->fmt != 0)
		errx(EXIT_FAILURE, "readwise density ANI currently supports detail output (-m0) only");
	if (ani_opt->sketch_reads_qc || ani_opt->sketch_abundance)
		errx(EXIT_FAILURE,
			 "readwise density ANI does not yet support --readsQC or --abundance; "
			 "use --save-query-sketch to force the materialized query-sketch path");
	if (ani_opt->readwise_unique_out[0] != '\0' && !ani_opt->readwise_profile_only)
		errx(EXIT_FAILURE, "--readwise-unique-out requires --readwise-profile-only");
	if (ani_opt->readwise_exact_split_out[0] != '\0' && !ani_opt->readwise_profile_only)
		errx(EXIT_FAILURE, "--readwise-exact-split-out requires --readwise-profile-only");
	const bool use_density_cache_in = ani_opt->readwise_density_cache_in[0] != '\0';
	const bool use_density_cache_out = ani_opt->readwise_density_cache_out[0] != '\0';
	if (use_density_cache_in && use_density_cache_out)
		errx(EXIT_FAILURE,
			 "--readwise-density-cache-in cannot be combined with --readwise-density-cache-out");
	const char *trace_out_env = getenv("MINCO_READWISE_TRACE_OUT");
	const char *edge_out_env = getenv("MINCO_READWISE_EDGE_OUT");
	if (use_density_cache_in &&
		(ani_opt->readwise_track[0] != '\0' ||
		 ani_opt->readwise_edge_out[0] != '\0' ||
		 (trace_out_env && trace_out_env[0] != '\0') ||
		 (edge_out_env && edge_out_env[0] != '\0')))
		errx(EXIT_FAILURE,
			 "--readwise-density-cache-in cannot be combined with read/edge tracing");

	char *ref_stat_path = test_get_fullpath(ani_opt->refdir, sketch_stat);
	size_t ref_stat_size = 0;
	minco_sketch_stat_t *ref_stat = read_from_file(ref_stat_path, &ref_stat_size);
	free(ref_stat_path);
	minco_sketch_info_t ref_info = {0};
	if (!minco_stat_decode_mem(ref_stat, ref_stat_size, ref_stat, &ref_info))
		errx(EXIT_FAILURE, "%s(): malformed %s/%s", __func__, ani_opt->refdir, sketch_stat);
	if (ref_info.has_minco_ext)
		ref_stat->hash_id = ref_info.sketch_id;
	ani_model_target_sketch_size = ref_info.target_sketch_size
										? ref_info.target_sketch_size
										: ANI_MODEL_REFERENCE_SKETCH_SIZE;
	const uint32_t ref_n = (uint32_t)ref_stat->infile_num;
	char (*refname)[PATHLEN] = (char (*)[PATHLEN])(ref_stat + 1);
	ani_readwise_trace_t readwise_trace = {0};
	const char *trace_ref_match = getenv("MINCO_READWISE_TRACE_REF");
	const char *trace_out_path = getenv("MINCO_READWISE_TRACE_OUT");
	if (trace_ref_match && trace_ref_match[0] && trace_out_path && trace_out_path[0]) {
		bool found_trace_ref = false;
		for (uint32_t rn = 0; rn < ref_n; ++rn) {
			if (strstr(refname[rn], trace_ref_match) == NULL)
				continue;
			readwise_trace.gid = rn;
			readwise_trace.ref_name = refname[rn];
			found_trace_ref = true;
			break;
		}
		if (!found_trace_ref)
			errx(EXIT_FAILURE,
				 "MINCO_READWISE_TRACE_REF target not found in reference sketch: %s",
				 trace_ref_match);
		readwise_trace.fp = fopen(trace_out_path, "w");
		if (!readwise_trace.fp)
			err(errno, "Cannot write MINCO_READWISE_TRACE_OUT %s", trace_out_path);
		fprintf(readwise_trace.fp,
				"read_id\tunit_id\tqctx\tref_begin\tdiff\tbest_diff\tcandidate_refs\tselected_refs\tcov_inc\tref\n");
		fprintf(stderr,
				"minco readwise: tracing assignments for ref gid=%u match=%s to %s\n",
				readwise_trace.gid, trace_ref_match, trace_out_path);
	}
	ani_readwise_edge_trace_t edge_trace = {0};
	const bool edge_from_cli = ani_opt->readwise_edge_out[0] != '\0';
	const char *edge_out_path = edge_from_cli
		? ani_opt->readwise_edge_out
		: getenv("MINCO_READWISE_EDGE_OUT");
	if (edge_out_path && edge_out_path[0]) {
		edge_trace.max_edges = edge_from_cli ? ani_opt->readwise_edge_max : 10000000ULL;
		edge_trace.max_candidate_refs = edge_from_cli ? ani_opt->readwise_edge_max_candidates : 64u;
		edge_trace.ambiguous_only = edge_from_cli ? ani_opt->readwise_edge_ambiguous_only : true;
		edge_trace.selected_only = edge_from_cli ? ani_opt->readwise_edge_selected_only : false;
		if (!edge_from_cli) {
			const char *edge_max = getenv("MINCO_READWISE_EDGE_MAX");
			if (edge_max && edge_max[0])
				edge_trace.max_edges = strtoull(edge_max, NULL, 10);
			const char *edge_max_cand = getenv("MINCO_READWISE_EDGE_MAX_CANDIDATES");
			if (edge_max_cand && edge_max_cand[0]) {
				const unsigned long v = strtoul(edge_max_cand, NULL, 10);
				edge_trace.max_candidate_refs = v > UINT32_MAX ? UINT32_MAX : (uint32_t)v;
			}
			const char *edge_ambiguous_only = getenv("MINCO_READWISE_EDGE_AMBIGUOUS_ONLY");
			if (edge_ambiguous_only && edge_ambiguous_only[0] &&
				strcmp(edge_ambiguous_only, "0") == 0)
				edge_trace.ambiguous_only = false;
			const char *edge_selected_only = getenv("MINCO_READWISE_EDGE_SELECTED_ONLY");
			if (edge_selected_only && edge_selected_only[0] &&
				strcmp(edge_selected_only, "0") != 0)
				edge_trace.selected_only = true;
		}
		edge_trace.fp = fopen(edge_out_path, "w");
		if (!edge_trace.fp)
			err(errno, "Cannot write readwise edge output %s", edge_out_path);
		fprintf(edge_trace.fp,
				"edge_id\tread_id\tunit_id\tqctx\tedge_rank\tref_begin\tgid\tdiff\tbest_diff\tcandidate_refs\tselected_refs\tselected_by_mode\tcov_inc\n");
		fprintf(stderr,
				"minco readwise: ambiguity edge dump active; out=%s max_edges=%" PRIu64
				" max_candidate_refs=%u ambiguous_only=%u selected_only=%u\n",
				edge_out_path, edge_trace.max_edges, edge_trace.max_candidate_refs,
				edge_trace.ambiguous_only ? 1u : 0u,
				edge_trace.selected_only ? 1u : 0u);
	}
	FILTER = UINT32_MAX >> ref_stat->compat_filter_shift;
	const_comask_init(ref_stat);
	const bool ref_uses_ctxobj96 =
		minco_stat_needs_ctxobj96(ref_stat) ||
		minco_stat_needs_ctxgidobj128(ref_stat);
	fprintf(stderr,
			"minco readwise: reference sketch coden_len=%d; ctx_bits=%u; obj_bits=%u; storage=%s\n",
			ref_stat->coden_len > 0 ? ref_stat->coden_len : NUM_CODENS,
			minco_stat_ctx_bits(ref_stat),
			minco_stat_obj_bits(ref_stat),
			ref_uses_ctxobj96 ? "ctxobj96" : "ctxobj64");
	if (ref_uses_ctxobj96 && !ani_opt->readwise_profile_only)
		errx(EXIT_FAILURE,
			 "ctxobj96 readwise density ANI currently requires --readwise-profile-only");
	if (!ref_uses_ctxobj96 && Bitslen.ctx + GID_NBITS > 64)
		errx(EXIT_FAILURE, "%s(): context bits (%u) + gid bits (%u) exceed 64",
			 __func__, Bitslen.ctx, GID_NBITS);

	const char *index_suffix = ref_uses_ctxobj96
		? sorted_comb_ctx64gid32obj32
		: sorted_comb_ctxgid64obj32;
	if (!file_exists_in_folder(ani_opt->refdir, index_suffix))
		gen_inverted_index_for_minco(ani_opt->refdir);
	char *index_path = test_get_fullpath(ani_opt->refdir, index_suffix);
	size_t index_bytes = 0;
	bool index_is_mmap = false;
	void *index_mem = read_reference_sorted_index(index_path, &index_bytes, &index_is_mmap);
	free(index_path);
	ctxgidobj_t *index = ref_uses_ctxobj96 ? NULL : (ctxgidobj_t *)index_mem;
	ctxgidobj128_t *index96 = ref_uses_ctxobj96 ? (ctxgidobj128_t *)index_mem : NULL;
	const size_t index_item_size = ref_uses_ctxobj96 ? sizeof(index96[0]) : sizeof(index[0]);
	if (index_bytes % index_item_size != 0)
		errx(EXIT_FAILURE, "%s(): malformed sorted reference index", __func__);
	const size_t index_n = index_bytes / index_item_size;
	uint32_t *ref_ctx_total = ref_uses_ctxobj96
		? ani_ref_ctx_counts_from_sorted_index128(index96, index_n, ref_n,
												  ani_opt->ignoreconflict, false)
		: ani_ref_ctx_counts_from_sorted_index(index, index_n, ref_n,
											   ani_opt->ignoreconflict, false);
	char (*refanno)[PATHLEN] = read_optional_sketch_annotations(ani_opt->refdir, (int)ref_n);
	ani_ctxmeta_rec_t *ref_ctxmeta =
		read_optional_ani_ctxmeta_stats(ani_opt->refdir, (int)ref_n);
	infile_meta_t *ref_infile_meta =
		ani_best_guard_enabled(ani_opt) ? read_optional_sketch_infile_meta_stats(ani_opt->refdir, (int)ref_n) : NULL;
	uint8_t *ref_domain_profiles =
		read_optional_domain_profiles(ani_opt->refdir, (int)ref_n);

	const int fence_k = minco_choose_k_fenceposts(index_n, 1024);
	const size_t fence_buckets = (size_t)1u << fence_k;
	size_t *fence = malloc((fence_buckets + 1u) * sizeof(*fence));
	if (!fence)
		err(EXIT_FAILURE, "%s(): OOM fenceposts", __func__);
	if (ref_uses_ctxobj96
			? minco_build_fenceposts_ctxgid128(index96, index_n, fence_k, fence) != 0
			: minco_build_fenceposts_ctxgid(index, index_n, fence_k, fence) != 0)
		errx(EXIT_FAILURE, "%s(): failed to build reference fenceposts", __func__);

	ani_readwise_acc_t *acc = calloc((size_t)ref_n, sizeof(acc[0]));
	uint8_t *ref_hit_bits = calloc((index_n + 7u) / 8u, 1);
	const bool need_unique_sidecar = ani_opt->readwise_unique_out[0] != '\0';
	const bool need_exact_split_sidecar = ani_opt->readwise_exact_split_out[0] != '\0';
	const bool need_ref_ctx_cov =
		ani_opt->abundance_model != ANI_ABUNDANCE_NONE ||
		ani_opt->readwise_ctx_filter_model != ANI_READWISE_CTX_FILTER_NONE;
	const bool test_mindiff_depth =
		need_ref_ctx_cov && getenv("MINCO_TEST_MINDIFF_DEPTH") &&
		getenv("MINCO_TEST_MINDIFF_DEPTH")[0] != '\0' &&
		strcmp(getenv("MINCO_TEST_MINDIFF_DEPTH"), "0") != 0;
	uint32_t *ref_ctx_cov = need_ref_ctx_cov
		? calloc(index_n, sizeof(ref_ctx_cov[0]))
		: NULL;
	uint32_t *ref_ctx_mindiff_cov = test_mindiff_depth
		? calloc(index_n, sizeof(ref_ctx_mindiff_cov[0]))
		: NULL;
	uint16_t *ref_ctx_hit_weight =
		need_ref_ctx_cov &&
				ani_opt->readwise_assign_mode == ANI_READWISE_ASSIGN_BEST_DIFF_SPLIT
			? calloc(index_n, sizeof(ref_ctx_hit_weight[0]))
			: NULL;
	ani_readwise_acc_t *unique_acc = need_unique_sidecar
		? calloc((size_t)ref_n, sizeof(unique_acc[0]))
		: NULL;
	uint8_t *unique_ref_hit_bits = need_unique_sidecar
		? calloc((index_n + 7u) / 8u, 1)
		: NULL;
	uint32_t *unique_ref_ctx_cov = need_unique_sidecar
		? calloc(index_n, sizeof(unique_ref_ctx_cov[0]))
		: NULL;
	ani_readwise_acc_t *exact_split_acc = need_exact_split_sidecar
		? calloc((size_t)ref_n, sizeof(exact_split_acc[0]))
		: NULL;
	uint8_t *exact_split_ref_hit_bits = need_exact_split_sidecar
		? calloc((index_n + 7u) / 8u, 1)
		: NULL;
	uint32_t *exact_split_ref_ctx_cov = need_exact_split_sidecar
		? calloc(index_n, sizeof(exact_split_ref_ctx_cov[0]))
		: NULL;
	uint16_t *exact_split_ref_ctx_hit_weight = need_exact_split_sidecar
		? calloc(index_n, sizeof(exact_split_ref_ctx_hit_weight[0]))
		: NULL;
	if (!acc || !ref_hit_bits)
		err(EXIT_FAILURE, "%s(): OOM readwise accumulators", __func__);
	if (need_unique_sidecar && (!unique_acc || !unique_ref_hit_bits || !unique_ref_ctx_cov))
		err(EXIT_FAILURE, "%s(): OOM readwise unique sidecar accumulators", __func__);
	if (need_exact_split_sidecar &&
		(!exact_split_acc || !exact_split_ref_hit_bits ||
		 !exact_split_ref_ctx_cov || !exact_split_ref_ctx_hit_weight))
		err(EXIT_FAILURE, "%s(): OOM readwise exact split sidecar accumulators", __func__);
	if (need_ref_ctx_cov && !ref_ctx_cov)
		err(EXIT_FAILURE, "%s(): OOM readwise abundance coverage", __func__);
	if (test_mindiff_depth && !ref_ctx_mindiff_cov)
		err(EXIT_FAILURE, "%s(): OOM readwise min-diff coverage", __func__);
	if (need_ref_ctx_cov &&
		ani_opt->readwise_assign_mode == ANI_READWISE_ASSIGN_BEST_DIFF_SPLIT &&
		!ref_ctx_hit_weight)
		err(EXIT_FAILURE, "%s(): OOM readwise split hit weights", __func__);
	if (test_mindiff_depth)
		fprintf(stderr,
				"minco readwise: experimental min-diff-depth mode active for reliable ctx metrics\n");
	const bool profile_only = ani_opt->readwise_profile_only;
	if (profile_only)
		fprintf(stderr,
				"minco readwise: profile-only mode active; exact query-context AF sets disabled\n");
	ani_u64_set_t qry_ctx_seen = {0};
	ani_u64_set_t qry_ref_ctx_seen = {0};
	if (!profile_only) {
		ani_u64_set_init(&qry_ctx_seen, 1u << 16);
		ani_u64_set_init(&qry_ref_ctx_seen, 1u << 18);
	}

	int worker_n = ani_opt->p > 0 ? ani_opt->p : 1;
	const int omp_max = omp_get_max_threads();
	if (worker_n > omp_max)
		worker_n = omp_max;
	if (worker_n < 1)
		worker_n = 1;
	ani_readwise_thread_state_t *thread_states =
		calloc((size_t)worker_n, sizeof(thread_states[0]));
	if (!thread_states)
		err(EXIT_FAILURE, "%s(): OOM readwise thread states", __func__);
	for (int t = 0; t < worker_n; ++t)
		ani_readwise_thread_state_init(&thread_states[t], ref_n, index_n,
										ref_ctx_cov != NULL, !profile_only);
	ani_readwise_thread_state_t *unique_thread_states = need_unique_sidecar
		? calloc((size_t)worker_n, sizeof(unique_thread_states[0]))
		: NULL;
	if (need_unique_sidecar && !unique_thread_states)
		err(EXIT_FAILURE, "%s(): OOM readwise unique sidecar thread states", __func__);
	for (int t = 0; t < worker_n && unique_thread_states; ++t)
		ani_readwise_thread_state_init(&unique_thread_states[t], ref_n, index_n,
										true, false);
	ani_readwise_thread_state_t *exact_split_thread_states = need_exact_split_sidecar
		? calloc((size_t)worker_n, sizeof(exact_split_thread_states[0]))
		: NULL;
	if (need_exact_split_sidecar && !exact_split_thread_states)
		err(EXIT_FAILURE, "%s(): OOM readwise exact split sidecar thread states", __func__);
	for (int t = 0; t < worker_n && exact_split_thread_states; ++t)
		ani_readwise_thread_state_init(&exact_split_thread_states[t], ref_n, index_n,
										true, false);

	const uint8_t nobjbits = Bitslen.obj;
	ani_density_cache_reader_t density_cache_reader = {0};
	ani_density_cache_writer_t density_cache_writer = {0};
	ani_fastx_stream_t stream = {0};
	kseq_t *seq = NULL;
	if (use_density_cache_in) {
		ani_density_cache_reader_open(
			&density_cache_reader, ani_opt->readwise_density_cache_in,
			ref_uses_ctxobj96, nobjbits, density_threshold);
	} else {
		stream = ani_open_fastx_stream(query_path, ani_opt->sketch_pipecmd);
		(void)gzbuffer(stream.gz, 4u << 20);
		seq = kseq_init(stream.gz);
		if (!seq)
			err(errno, "%s(): kseq_init %s", __func__, query_path);
		if (use_density_cache_out)
			ani_density_cache_writer_open(
				&density_cache_writer, ani_opt->readwise_density_cache_out,
				ref_uses_ctxobj96, nobjbits, density_threshold);
	}
	ani_readwise_progress_t progress = ani_readwise_progress_start(&stream);

	const uint64_t objmask = (nobjbits == 64) ? UINT64_MAX : ((1ULL << nobjbits) - 1ULL);
	const uint64_t gidmask_local = (1ULL << GID_NBITS) - 1ULL;
	ani_readwise_tracker_t read_tracker;
	ani_readwise_tracker_init(
		&read_tracker, ani_opt, query_path, refname, refanno, ref_n,
		ref_ctxmeta,
		ani_readwise_density_probability(density_threshold, ref_uses_ctxobj96,
										 nobjbits));
	uint64_t total_reads = use_density_cache_in ? density_cache_reader.hdr.total_reads : 0;
	uint64_t total_density_blocks = 0;
	uint64_t exact_split_total_density_blocks = 0;
	uint64_t block_reads_with_density_ctx = 0;
	u64vec block_vec = {0};
	ani_ctxobj96_vec_t block_vec96 = {0};
	if (ref_uses_ctxobj96)
		ani_ctxobj96_vec_init(&block_vec96,
							   ani_opt->density_block_ctx > 2048u
								   ? ani_opt->density_block_ctx
								   : 2048u);
	else
		v_init(&block_vec, ani_opt->density_block_ctx > 2048u ? ani_opt->density_block_ctx : 2048u);
	const uint32_t density_block_ctx = ani_opt->density_block_ctx;
	const bool density_block_mode = density_block_ctx > 1u;
	if (density_block_mode)
		fprintf(stderr,
				"minco readwise: density block mode active; density_block_ctx=%u\n",
				density_block_ctx);
	if (ani_opt->readwise_assign_mode == ANI_READWISE_ASSIGN_BEST_DIFF)
		fprintf(stderr,
				"minco readwise: shared-context assignment active; mode=best-diff\n");
	else if (ani_opt->readwise_assign_mode == ANI_READWISE_ASSIGN_BEST_DIFF_SPLIT)
		fprintf(stderr,
				"minco readwise: shared-context assignment active; mode=best-diff-split; coverage_scale=%u\n",
				ANI_REFCOV_SPLIT_SCALE);
	else if (ani_opt->readwise_assign_mode == ANI_READWISE_ASSIGN_BEST_DIFF_UNIQUE)
		fprintf(stderr,
				"minco readwise: shared-context assignment active; mode=best-diff-unique\n");
	if (worker_n > 1)
		fprintf(stderr,
				"minco readwise: parallel density processing active; threads=%d; batch_reads=%u\n",
				worker_n, (unsigned)ANI_READWISE_BATCH_READS);

	int read_status = 0;
	for (;;) {
		ani_readwise_seq_rec_t batch[ANI_READWISE_BATCH_READS];
		memset(batch, 0, sizeof(batch));
		size_t batch_n = 0;
		u64vec *read_vecs = NULL;
		ani_ctxobj96_vec_t *read_vecs96 = ref_uses_ctxobj96
			? calloc(ANI_READWISE_BATCH_READS, sizeof(read_vecs96[0]))
			: NULL;
		if (!ref_uses_ctxobj96)
			read_vecs = calloc(ANI_READWISE_BATCH_READS, sizeof(read_vecs[0]));
		if ((!ref_uses_ctxobj96 && !read_vecs) ||
			(ref_uses_ctxobj96 && !read_vecs96))
			err(EXIT_FAILURE, "%s(): OOM read density vectors", __func__);

		if (use_density_cache_in) {
			batch_n = ref_uses_ctxobj96
				? ani_density_cache_reader_read96(
					  &density_cache_reader, read_vecs96, ANI_READWISE_BATCH_READS)
				: ani_density_cache_reader_read64(
					  &density_cache_reader, read_vecs, ANI_READWISE_BATCH_READS);
		} else {
			int kseq_status = 0;
			while (batch_n < ANI_READWISE_BATCH_READS) {
				kseq_status = kseq_read(seq);
				if (kseq_status < 0) {
					read_status = kseq_status;
					break;
				}
				const size_t len = seq->seq.l;
				char *copy = malloc(len + 1u);
				if (!copy)
					err(EXIT_FAILURE, "%s(): OOM read batch sequence", __func__);
				memcpy(copy, seq->seq.s, len);
				copy[len] = '\0';
				char *name_copy = NULL;
				if (seq->name.s) {
					name_copy = strdup(seq->name.s);
					if (!name_copy)
						err(EXIT_FAILURE, "%s(): OOM read batch name", __func__);
				}
				batch[batch_n].seq = copy;
				batch[batch_n].name = name_copy;
				batch[batch_n].len = (int)len;
				++batch_n;
				++total_reads;
				ani_readwise_progress_update(&stream, &progress, total_reads, false);
			}
			if (!batch_n) {
				free(read_vecs);
				free(read_vecs96);
				break;
			}

#pragma omp parallel for num_threads(worker_n) schedule(dynamic, 256)
			for (size_t i = 0; i < batch_n; ++i) {
				if (ref_uses_ctxobj96) {
					ani_ctxobj96_vec_init(&read_vecs96[i], 128);
					ani_extract_read_density_ctxobjs96(batch[i].seq, batch[i].len,
													   &read_vecs96[i],
													   density_threshold);
				} else {
					v_init(&read_vecs[i], 128);
					ani_extract_read_density_ctxobjs(batch[i].seq, batch[i].len,
													 &read_vecs[i], nobjbits,
													 density_threshold);
				}
			}
			if (use_density_cache_out) {
				for (size_t i = 0; i < batch_n; ++i) {
					if (ref_uses_ctxobj96)
						ani_density_cache_writer_write96(&density_cache_writer,
														 &read_vecs96[i]);
					else
						ani_density_cache_writer_write64(&density_cache_writer,
														 &read_vecs[i]);
				}
			}
		}
		if (!batch_n) {
			free(read_vecs);
			free(read_vecs96);
			break;
		}

		if (ref_uses_ctxobj96) {
			if (read_tracker.enabled) {
				const uint64_t first_read_id = total_reads - (uint64_t)batch_n + 1u;
				for (size_t i = 0; i < batch_n; ++i) {
					ani_readwise_track_read96(
						&read_tracker, batch[i].name,
						first_read_id + (uint64_t)i,
						batch[i].seq, batch[i].len,
						index96, index_n, fence, fence_k, ref_n,
						ani_opt->ignoreconflict, density_threshold,
						ani_opt->readwise_assign_mode);
				}
			}
			kv_density_unit96_t units96;
			kv_init(units96);
			kv_density_unit96_t exact_units96;
			kv_init(exact_units96);
			for (size_t i = 0; i < batch_n; ++i) {
				if (read_vecs96[i].n == 0) {
					ani_ctxobj96_vec_free(&read_vecs96[i]);
					continue;
				}
				if (need_exact_split_sidecar) {
					++exact_split_total_density_blocks;
					ani_density_unit96_push_copy(
						&exact_units96, &read_vecs96[i],
						exact_split_total_density_blocks, 1,
						batch[i].name);
				}
				if (density_block_mode) {
					ani_ctxobj96_vec_reserve(&block_vec96,
											 block_vec96.n + read_vecs96[i].n);
					memcpy(block_vec96.a + block_vec96.n, read_vecs96[i].a,
						   read_vecs96[i].n * sizeof(read_vecs96[i].a[0]));
					block_vec96.n += read_vecs96[i].n;
					block_reads_with_density_ctx++;
					ani_ctxobj96_vec_free(&read_vecs96[i]);
					if (block_vec96.n < (size_t)density_block_ctx)
						continue;
					++total_density_blocks;
					ani_density_unit96_push_move(&units96, &block_vec96,
												 total_density_blocks,
												 block_reads_with_density_ctx,
												 NULL);
					ani_ctxobj96_vec_init(&block_vec96,
										   ani_opt->density_block_ctx > 2048u
											   ? ani_opt->density_block_ctx
											   : 2048u);
					block_reads_with_density_ctx = 0;
				} else {
					++total_density_blocks;
					ani_density_unit96_push_move(&units96, &read_vecs96[i],
												 total_density_blocks, 1,
												 batch[i].name);
				}
			}
			free(read_vecs96);
			ani_readwise_seq_batch_destroy(batch, batch_n);

			if (need_exact_split_sidecar) {
				ani_process_density_units96_parallel(
					&exact_units96, exact_split_thread_states, worker_n,
					index96, index_n, fence, fence_k, ref_n,
					ani_opt->ignoreconflict,
					true, ANI_READWISE_ASSIGN_BEST_DIFF_SPLIT,
					NULL, NULL, NULL);
				ani_merge_readwise_thread_batch128(
					exact_split_thread_states, worker_n, exact_split_acc,
					exact_split_ref_hit_bits, exact_split_ref_ctx_cov,
					NULL,
					exact_split_ref_ctx_hit_weight,
					index96, index_n,
					ref_n);
				ani_density_units96_destroy(&exact_units96);
			}

			ani_process_density_units96_parallel(
				&units96, thread_states, worker_n,
				index96, index_n, fence, fence_k, ref_n,
				ani_opt->ignoreconflict,
				ref_ctx_cov != NULL, ani_opt->readwise_assign_mode,
				readwise_trace.fp ? &readwise_trace : NULL,
				edge_trace.fp ? &edge_trace : NULL,
				need_unique_sidecar ? unique_thread_states : NULL);
			ani_merge_readwise_thread_batch128(
				thread_states, worker_n, acc, ref_hit_bits, ref_ctx_cov,
				ref_ctx_mindiff_cov,
				ref_ctx_hit_weight,
				index96, index_n,
				ref_n);
			if (need_unique_sidecar) {
				ani_merge_readwise_thread_batch128(
					unique_thread_states, worker_n, unique_acc,
					unique_ref_hit_bits, unique_ref_ctx_cov,
					NULL, NULL,
					index96, index_n,
					ref_n);
			}
			ani_density_units96_destroy(&units96);
			continue;
		}

		if (read_tracker.enabled) {
			const uint64_t first_read_id = total_reads - (uint64_t)batch_n + 1u;
			for (size_t i = 0; i < batch_n; ++i) {
				ani_readwise_track_read64(
					&read_tracker, batch[i].name,
					first_read_id + (uint64_t)i,
					batch[i].seq, batch[i].len,
					index, index_n, fence, fence_k, ref_n,
					ani_opt->ignoreconflict, nobjbits, gidmask_local,
					objmask, density_threshold,
					ani_opt->readwise_assign_mode);
			}
		}
		kv_density_unit_t units;
		kv_init(units);
		kv_density_unit_t exact_units;
		kv_init(exact_units);
		for (size_t i = 0; i < batch_n; ++i) {
			if (read_vecs[i].n == 0) {
				v_free(&read_vecs[i]);
				continue;
			}
			if (need_exact_split_sidecar) {
				++exact_split_total_density_blocks;
				ani_density_unit_push_copy(
					&exact_units, &read_vecs[i],
					exact_split_total_density_blocks, 1,
					batch[i].name);
			}
			if (density_block_mode) {
				v_reserve(&block_vec, block_vec.n + read_vecs[i].n);
				memcpy(block_vec.a + block_vec.n, read_vecs[i].a,
					   read_vecs[i].n * sizeof(read_vecs[i].a[0]));
				block_vec.n += read_vecs[i].n;
				block_reads_with_density_ctx++;
				v_free(&read_vecs[i]);
				if (block_vec.n < (size_t)density_block_ctx)
					continue;
				++total_density_blocks;
				ani_density_unit_push_move(&units, &block_vec,
										   total_density_blocks,
										   block_reads_with_density_ctx,
										   NULL);
				v_init(&block_vec, ani_opt->density_block_ctx > 2048u
										? ani_opt->density_block_ctx
										: 2048u);
				block_reads_with_density_ctx = 0;
			} else {
				++total_density_blocks;
				ani_density_unit_push_move(&units, &read_vecs[i],
										   total_density_blocks, 1,
										   batch[i].name);
			}
		}
		free(read_vecs);
		ani_readwise_seq_batch_destroy(batch, batch_n);

		if (need_exact_split_sidecar) {
			ani_process_density_units_parallel(
				&exact_units, exact_split_thread_states, worker_n,
				index, index_n, fence, fence_k, ref_n,
				ani_opt->ignoreconflict, nobjbits, gidmask_local, objmask,
				true, false, ANI_READWISE_ASSIGN_BEST_DIFF_SPLIT,
				NULL, NULL, NULL);
			ani_merge_readwise_thread_batch(
				exact_split_thread_states, worker_n, exact_split_acc,
				exact_split_ref_hit_bits, exact_split_ref_ctx_cov,
				NULL,
				exact_split_ref_ctx_hit_weight,
				NULL, NULL,
				index, index_n,
				ref_n, gidmask_local);
			ani_density_units_destroy(&exact_units);
		}

		ani_process_density_units_parallel(
			&units, thread_states, worker_n,
			index, index_n, fence, fence_k, ref_n,
			ani_opt->ignoreconflict, nobjbits, gidmask_local, objmask,
			ref_ctx_cov != NULL, !profile_only, ani_opt->readwise_assign_mode,
			readwise_trace.fp ? &readwise_trace : NULL,
			edge_trace.fp ? &edge_trace : NULL,
			need_unique_sidecar ? unique_thread_states : NULL);
		ani_merge_readwise_thread_batch(
			thread_states, worker_n, acc, ref_hit_bits, ref_ctx_cov,
			ref_ctx_mindiff_cov,
			ref_ctx_hit_weight,
			profile_only ? NULL : &qry_ctx_seen,
			profile_only ? NULL : &qry_ref_ctx_seen,
			index, index_n,
			ref_n, gidmask_local);
		if (need_unique_sidecar) {
			ani_merge_readwise_thread_batch(
				unique_thread_states, worker_n, unique_acc,
				unique_ref_hit_bits, unique_ref_ctx_cov,
				NULL, NULL, NULL, NULL,
				index, index_n,
				ref_n, gidmask_local);
		}
		ani_density_units_destroy(&units);
		if (read_status < 0)
			break;
	}
	if (read_status < -1)
		errx(EXIT_FAILURE,
			 "%s(): malformed FASTA/FASTQ or stream read error from %s (kseq_read=%d)",
			 __func__, query_path, read_status);
	if (density_block_mode && ref_uses_ctxobj96 && block_vec96.n > 0) {
		kv_density_unit96_t units96;
		kv_init(units96);
		++total_density_blocks;
		ani_density_unit96_push_move(&units96, &block_vec96, total_density_blocks,
									 block_reads_with_density_ctx,
									 NULL);
		ani_process_density_units96_parallel(
			&units96, thread_states, worker_n,
			index96, index_n, fence, fence_k, ref_n,
			ani_opt->ignoreconflict,
			ref_ctx_cov != NULL, ani_opt->readwise_assign_mode,
			readwise_trace.fp ? &readwise_trace : NULL,
			edge_trace.fp ? &edge_trace : NULL,
			need_unique_sidecar ? unique_thread_states : NULL);
		ani_merge_readwise_thread_batch128(
			thread_states, worker_n, acc, ref_hit_bits, ref_ctx_cov,
			ref_ctx_mindiff_cov,
			ref_ctx_hit_weight,
			index96, index_n,
			ref_n);
		if (need_unique_sidecar) {
			ani_merge_readwise_thread_batch128(
				unique_thread_states, worker_n, unique_acc,
				unique_ref_hit_bits, unique_ref_ctx_cov,
				NULL, NULL,
				index96, index_n,
				ref_n);
		}
		ani_density_units96_destroy(&units96);
	} else if (density_block_mode && block_vec.n > 0) {
		kv_density_unit_t units;
		kv_init(units);
		++total_density_blocks;
		ani_density_unit_push_move(&units, &block_vec, total_density_blocks,
								   block_reads_with_density_ctx,
								   NULL);
		ani_process_density_units_parallel(
			&units, thread_states, worker_n,
			index, index_n, fence, fence_k, ref_n,
			ani_opt->ignoreconflict, nobjbits, gidmask_local, objmask,
			ref_ctx_cov != NULL, !profile_only, ani_opt->readwise_assign_mode,
			readwise_trace.fp ? &readwise_trace : NULL,
			edge_trace.fp ? &edge_trace : NULL,
			need_unique_sidecar ? unique_thread_states : NULL);
		ani_merge_readwise_thread_batch(
			thread_states, worker_n, acc, ref_hit_bits, ref_ctx_cov,
			ref_ctx_mindiff_cov,
			ref_ctx_hit_weight,
			profile_only ? NULL : &qry_ctx_seen,
			profile_only ? NULL : &qry_ref_ctx_seen,
			index, index_n,
			ref_n, gidmask_local);
		if (need_unique_sidecar) {
			ani_merge_readwise_thread_batch(
				unique_thread_states, worker_n, unique_acc,
				unique_ref_hit_bits, unique_ref_ctx_cov,
				NULL, NULL, NULL, NULL,
				index, index_n,
				ref_n, gidmask_local);
		}
		ani_density_units_destroy(&units);
	}
	ani_readwise_progress_done(&stream, &progress, total_reads);
	if (ref_uses_ctxobj96)
		ani_ctxobj96_vec_free(&block_vec96);
	else
		v_free(&block_vec);
	if (use_density_cache_out)
		ani_density_cache_writer_close(&density_cache_writer, total_reads);
	if (use_density_cache_in)
		ani_density_cache_reader_close(&density_cache_reader);
	if (seq)
		kseq_destroy(seq);
	const bool close_warning = use_density_cache_in
		? false
		: ani_close_fastx_stream(&stream, true);
	if (close_warning)
		fprintf(stderr,
				"minco readwise: warning: gzip integrity/close error after %" PRIu64
				" complete reads; continuing with completed result\n",
				total_reads);
	if (readwise_trace.fp) {
		fclose(readwise_trace.fp);
		readwise_trace.fp = NULL;
	}
	if (edge_trace.fp) {
		fprintf(stderr,
				"minco readwise: ambiguity edge dump wrote %" PRIu64 " edges\n",
				edge_trace.emitted_edges);
		fclose(edge_trace.fp);
		edge_trace.fp = NULL;
	}
	ani_readwise_tracker_destroy(&read_tracker);
	for (int t = 0; t < worker_n; ++t)
		ani_readwise_thread_state_destroy(&thread_states[t]);
	free(thread_states);
	if (unique_thread_states) {
		for (int t = 0; t < worker_n; ++t)
			ani_readwise_thread_state_destroy(&unique_thread_states[t]);
		free(unique_thread_states);
	}
	if (exact_split_thread_states) {
		for (int t = 0; t < worker_n; ++t)
			ani_readwise_thread_state_destroy(&exact_split_thread_states[t]);
		free(exact_split_thread_states);
	}

	FILE *outfp = ani_opt->outf[0] == '\0' ? stdout : fopen(ani_opt->outf, "w");
	if (!outfp)
		err(errno, "%s", ani_opt->outf);
	const bool old_readwise = ani_opt->readwise_query;
	const bool old_unassembled = ani_opt->unassembled;
	const bool old_v = ani_opt->v;
	const int old_ctxcut = ani_opt->ctxcut;
	const float old_afcut = ani_opt->afcut;
	const float old_anicut = ani_opt->anicut;
	get_generic_dist_from_features_fn old_dist_fn = get_generic_dist_from_features;
	const bool auto_readwise_abundance_report =
		ani_opt->abundance_model != ANI_ABUNDANCE_NONE &&
		!ani_opt->ctxcut_set &&
		!ani_opt->afcut_set &&
		!ani_opt->anicut_set &&
		ani_opt->ntop < 0;
	if (auto_readwise_abundance_report)
	{
		uint64_t auto_ctxcut =
			ani_readwise_default_min_support_cut(ref_ctx_total, ref_n);
		if (auto_ctxcut > (uint64_t)INT_MAX)
			auto_ctxcut = (uint64_t)INT_MAX;
		ani_opt->ctxcut = (int)auto_ctxcut;
		ani_opt->afcut = 0.0f;
		ani_opt->anicut = 0.95f;
		fprintf(stderr,
				"minco readwise: default abundance report active; support_cut=%d; "
				"breadth>=0.5; anicut=0.95; calls=major|low_abundance; "
				"final support is checked per reference marker count\n",
				ani_opt->ctxcut);
	}
	ani_opt->readwise_query = true;
	ani_opt->unassembled = true;
	ani_opt->v = true;
	get_generic_dist_from_features = get_naive_dist;
	if (!ani_opt->afcut_set)
		ani_opt->afcut = 0.2f;
	print_ani_detail_header(outfp, ani_opt, true);

	kv_ani_row_t survivors;
	kv_init(survivors);
	const uint64_t total_unique_qry_ctx = profile_only ? 0 : qry_ctx_seen.n;
	const uint32_t global_qry_ctx_for_dist = ani_clamp_u64_to_u32(total_unique_qry_ctx);
	ani_readwise_abundance_t *abundance_stats = NULL;
	const bool need_abundance_stats =
		ani_opt->abundance_model != ANI_ABUNDANCE_NONE ||
		ani_opt->readwise_ctx_filter_model != ANI_READWISE_CTX_FILTER_NONE;
	double coverage_scale =
		ani_opt->readwise_assign_mode == ANI_READWISE_ASSIGN_BEST_DIFF_SPLIT
			? (double)ANI_REFCOV_SPLIT_SCALE
			: 1.0;
	if (need_abundance_stats)
	{
		abundance_stats = ref_uses_ctxobj96
			? ani_depth_abundance_from_ref_coverage128(
				  index96, index_n, ref_n, ani_opt->ignoreconflict,
				  ref_ctx_cov, ref_ctx_hit_weight, ref_ctx_total, acc, coverage_scale, false)
			: ani_depth_abundance_from_ref_coverage(
				  index, index_n, ref_n, ani_opt->ignoreconflict,
				  ref_ctx_cov, ref_ctx_hit_weight, ref_ctx_total, acc, coverage_scale, false);
	}
	ani_readwise_acc_t *raw_feature_acc =
		ref_uses_ctxobj96
			? ani_unique_best_features_from_ref_covdiff128(
				  index96, index_n, ref_n, ani_opt->ignoreconflict, ref_ctx_cov, false)
			: ani_unique_best_features_from_ref_covdiff(
				  index, index_n, ref_n, ani_opt->ignoreconflict, ref_ctx_cov, false);
	ani_readwise_filter_stats_t *ctx_filter_stats = NULL;
	ani_readwise_acc_t *filtered_feature_acc = NULL;
	ani_readwise_reliable_abundance_t *reliable_abundance_stats = NULL;
	const uint32_t *reliable_ref_ctx_cov =
		ref_ctx_mindiff_cov ? ref_ctx_mindiff_cov : ref_ctx_cov;
	if (ani_opt->readwise_ctx_filter_model != ANI_READWISE_CTX_FILTER_NONE &&
		raw_feature_acc && abundance_stats) {
		ctx_filter_stats = calloc((size_t)ref_n, sizeof(ctx_filter_stats[0]));
		if (!ctx_filter_stats)
			err(EXIT_FAILURE, "%s(): OOM readwise ctx filter stats", __func__);
		filtered_feature_acc = ref_uses_ctxobj96
			? ani_reliable_features_from_ref_covdiff128(
				  index96, index_n, ref_n, ani_opt->ignoreconflict,
				  reliable_ref_ctx_cov, abundance_stats, ref_ctx_total,
				  ani_opt->readwise_ctx_filter_model,
				  coverage_scale,
				  ani_opt->readwise_fake_ctx_threshold, ctx_filter_stats, false)
			: ani_reliable_features_from_ref_covdiff(
				  index, index_n, ref_n, ani_opt->ignoreconflict,
				  reliable_ref_ctx_cov, abundance_stats, ref_ctx_total,
				  ani_opt->readwise_ctx_filter_model,
				  coverage_scale,
				  ani_opt->readwise_fake_ctx_threshold, ctx_filter_stats, false);
		reliable_abundance_stats = ref_uses_ctxobj96
			? ani_reliable_abundance_from_ref_covdiff128(
				  index96, index_n, ref_n, ani_opt->ignoreconflict,
				  reliable_ref_ctx_cov, ref_ctx_hit_weight, abundance_stats, ref_ctx_total,
				  ani_opt->readwise_ctx_filter_model,
				  coverage_scale,
				  ani_opt->readwise_fake_ctx_threshold, false)
			: ani_reliable_abundance_from_ref_covdiff(
				  index, index_n, ref_n, ani_opt->ignoreconflict,
				  reliable_ref_ctx_cov, ref_ctx_hit_weight, abundance_stats, ref_ctx_total,
				  ani_opt->readwise_ctx_filter_model,
				  coverage_scale,
				  ani_opt->readwise_fake_ctx_threshold, false);
	}
	uint32_t *marker_ref_ctx_total = NULL;
	ani_readwise_abundance_t *marker_abundance_stats = NULL;
	ani_readwise_acc_t *marker_raw_feature_acc = NULL;
	ani_readwise_acc_t *marker_filtered_feature_acc = NULL;
	ani_readwise_reliable_abundance_t *marker_reliable_abundance_stats = NULL;
	if (ani_opt->readwise_dual_evidence) {
		marker_ref_ctx_total = ref_uses_ctxobj96
			? ani_ref_ctx_counts_from_sorted_index128(
				  index96, index_n, ref_n, ani_opt->ignoreconflict, true)
			: ani_ref_ctx_counts_from_sorted_index(
				  index, index_n, ref_n, ani_opt->ignoreconflict, true);
		marker_abundance_stats = ref_uses_ctxobj96
			? ani_depth_abundance_from_ref_coverage128(
				  index96, index_n, ref_n, ani_opt->ignoreconflict,
				  ref_ctx_cov, ref_ctx_hit_weight, marker_ref_ctx_total, acc,
				  coverage_scale, true)
			: ani_depth_abundance_from_ref_coverage(
				  index, index_n, ref_n, ani_opt->ignoreconflict,
				  ref_ctx_cov, ref_ctx_hit_weight, marker_ref_ctx_total, acc,
				  coverage_scale, true);
		marker_raw_feature_acc = ref_uses_ctxobj96
			? ani_unique_best_features_from_ref_covdiff128(
				  index96, index_n, ref_n, ani_opt->ignoreconflict,
				  ref_ctx_cov, true)
			: ani_unique_best_features_from_ref_covdiff(
				  index, index_n, ref_n, ani_opt->ignoreconflict,
				  ref_ctx_cov, true);
		if (ani_opt->readwise_ctx_filter_model != ANI_READWISE_CTX_FILTER_NONE &&
			marker_raw_feature_acc && marker_abundance_stats) {
			marker_filtered_feature_acc = ref_uses_ctxobj96
				? ani_reliable_features_from_ref_covdiff128(
					  index96, index_n, ref_n, ani_opt->ignoreconflict,
					  reliable_ref_ctx_cov, marker_abundance_stats, marker_ref_ctx_total,
					  ani_opt->readwise_ctx_filter_model,
					  coverage_scale,
					  ani_opt->readwise_fake_ctx_threshold, NULL, true)
				: ani_reliable_features_from_ref_covdiff(
					  index, index_n, ref_n, ani_opt->ignoreconflict,
					  reliable_ref_ctx_cov, marker_abundance_stats, marker_ref_ctx_total,
					  ani_opt->readwise_ctx_filter_model,
					  coverage_scale,
					  ani_opt->readwise_fake_ctx_threshold, NULL, true);
			marker_reliable_abundance_stats = ref_uses_ctxobj96
				? ani_reliable_abundance_from_ref_covdiff128(
					  index96, index_n, ref_n, ani_opt->ignoreconflict,
					  reliable_ref_ctx_cov, ref_ctx_hit_weight, marker_abundance_stats,
					  marker_ref_ctx_total, ani_opt->readwise_ctx_filter_model,
					  coverage_scale, ani_opt->readwise_fake_ctx_threshold, true)
				: ani_reliable_abundance_from_ref_covdiff(
					  index, index_n, ref_n, ani_opt->ignoreconflict,
					  reliable_ref_ctx_cov, ref_ctx_hit_weight, marker_abundance_stats,
					  marker_ref_ctx_total, ani_opt->readwise_ctx_filter_model,
					  coverage_scale, ani_opt->readwise_fake_ctx_threshold, true);
		}
	}
	const ani_readwise_acc_t *support_feature_acc =
		filtered_feature_acc ? filtered_feature_acc :
		(raw_feature_acc ? raw_feature_acc : acc);
	const ani_readwise_acc_t *ani_feature_acc =
#if MINCO_REPORT_FILTERED_READWISE_ANI
		support_feature_acc;
#else
		raw_feature_acc ? raw_feature_acc : support_feature_acc;
#endif
	const ani_readwise_acc_t *marker_feature_acc =
		marker_filtered_feature_acc ? marker_filtered_feature_acc :
		marker_raw_feature_acc;
	const uint32_t feature_scale =
		filtered_feature_acc &&
				ani_opt->readwise_ctx_filter_model == ANI_READWISE_CTX_FILTER_FAKE_PROB
			? ANI_READWISE_PROB_SCALE
			: 1u;
	for (uint32_t rn = 0; rn < ref_n; ++rn) {
		const ani_readwise_acc_t *a = &acc[rn];
		const ani_readwise_acc_t *support_fa = &support_feature_acc[rn];
		const ani_readwise_acc_t *ani_fa = &ani_feature_acc[rn];
		const uint64_t qry_ctx_hit = profile_only ? a->ref_ctx_hit : a->qry_ctx_hit;
		const uint64_t unique_overlap = qry_ctx_hit < a->ref_ctx_hit ? qry_ctx_hit : a->ref_ctx_hit;
		if (unique_overlap < (uint64_t)ani_opt->ctxcut)
			continue;
		const uint64_t effective_xny_ctx =
			ani_readwise_scaled_round_u64(support_fa->XnY_ctx, feature_scale);
		if (!effective_xny_ctx || !ref_ctx_total[rn])
			continue;
		if (effective_xny_ctx < (uint64_t)ani_opt->ctxcut)
			continue;
		const double ref_af = (double)a->ref_ctx_hit / (double)ref_ctx_total[rn];
		if (!profile_only && !total_unique_qry_ctx)
			continue;
		const double qry_af = profile_only
								  ? ref_af
								  : (double)qry_ctx_hit / (double)total_unique_qry_ctx;
		ani_density_af_t density_af = {
			.qry = qry_af,
			.ref = ref_af,
			.available = 1,
		};
		if (!ani_report_af_pass(ani_opt, density_af.qry, density_af.ref))
			continue;

		ani_features_t f = {
			.XnY_ctx = ani_clamp_u64_to_u32(ani_fa->XnY_ctx),
			.X_ctx = profile_only ? ref_ctx_total[rn] : global_qry_ctx_for_dist,
			.N_diff_obj_section = ani_clamp_u64_to_u32(ani_fa->N_diff_obj_section),
			.N_mut2_ctx = ani_clamp_u64_to_u32(ani_fa->N_mut2_ctx),
			.N_diff_obj = ani_clamp_u64_to_u32(ani_fa->N_diff_obj),
		};
		const uint32_t qry_ctx_for_dist =
			profile_only ? ref_ctx_total[rn] : global_qry_ctx_for_dist;
		ani_row_t row = make_selected_output_row(
			rn, &f, ani_opt,
			qry_ctx_for_dist,
			ref_ctx_total[rn],
			qry_af, qry_af,
			ref_af, ref_af,
			density_af,
			NULL,
			ref_infile_meta ? &ref_infile_meta[rn] : NULL);
		row.domain_profile = ref_domain_profiles
								  ? ref_domain_profiles[rn]
								  : MINCO_DOMAIN_PROFILE_DEFAULT;
		row.XnY_ctx = ani_readwise_scaled_round_int(ani_fa->XnY_ctx, 1u);
		row.N_diff_obj = ani_readwise_scaled_round_int(ani_fa->N_diff_obj, 1u);
		row.N_diff_obj_section = ani_readwise_scaled_round_int(ani_fa->N_diff_obj_section, 1u);
		row.N_mut2_ctx = ani_readwise_scaled_round_int(ani_fa->N_mut2_ctx, 1u);
		row.readwise_total_reads = total_reads;
		row.readwise_reads_with_ctx_match = a->reads_with_ctx_match;
		row.readwise_unique_query_ctx = total_unique_qry_ctx;
		row.readwise_unique_query_ctx_hit = qry_ctx_hit;
		row.readwise_unique_ref_ctx_hit = a->ref_ctx_hit;
		row.readwise_ref_ctx_total = ref_ctx_total[rn];
		row.readwise_density_block_ctx = density_block_mode ? density_block_ctx : 0u;
		row.readwise_total_density_blocks = total_density_blocks;
		row.readwise_blocks_with_ctx_match = a->blocks_with_ctx_match;
		row.readwise_raw_xny_ctx = raw_feature_acc
									   ? raw_feature_acc[rn].XnY_ctx
									   : ani_fa->XnY_ctx;
		if (ctx_filter_stats) {
			row.readwise_rejected_ctx = ctx_filter_stats[rn].rejected_ctx;
			row.readwise_rejected_diff_ctx = ctx_filter_stats[rn].rejected_diff_ctx;
			const double fake_prob_mean = ctx_filter_stats[rn].raw_xny_ctx
											  ? (double)(ctx_filter_stats[rn].fake_prob_sum /
														 (long double)ctx_filter_stats[rn].raw_xny_ctx)
											  : 0.0;
			const double fake_prob_weighted = ctx_filter_stats[rn].fake_prob_weight_sum > 0.0L
												  ? (double)(ctx_filter_stats[rn].fake_prob_weighted_sum /
															 ctx_filter_stats[rn].fake_prob_weight_sum)
												  : 0.0;
			row.readwise_fake_ctx_prob_mean = fake_prob_mean;
			row.readwise_fake_ctx_prob_weighted = fake_prob_weighted;
			row.readwise_fake_ctx_fraction =
				ani_opt->readwise_ctx_filter_model == ANI_READWISE_CTX_FILTER_FAKE_PROB
					? fake_prob_mean
					: (ctx_filter_stats[rn].raw_xny_ctx
						   ? (double)ctx_filter_stats[rn].rejected_ctx /
								 (double)ctx_filter_stats[rn].raw_xny_ctx
						   : 0.0);
		}
		if (abundance_stats) {
			row.abundance_ref_breadth = abundance_stats[rn].ref_breadth;
			row.abundance_ref_mean_depth = abundance_stats[rn].ref_mean_depth;
			row.abundance_ref_hit_mean_depth = abundance_stats[rn].ref_hit_mean_depth;
            row.abundance_ref_depth_variance = abundance_stats[rn].ref_depth_variance;
            row.abundance_ref_depth_cv = abundance_stats[rn].ref_depth_cv;
            row.abundance_ref_zero_fraction = abundance_stats[rn].ref_zero_fraction;
            row.abundance_relative_depth = abundance_stats[rn].relative_depth;
			row.abundance_ref_zip_af = ani_zip_corrected_af(
				row.abundance_ref_breadth, row.abundance_ref_mean_depth);
			row.abundance_ref_zip_aaf_ani =
				aaf_ani_from_containment(row.abundance_ref_zip_af);
			if (reliable_abundance_stats) {
				row.reliable_ref_breadth = reliable_abundance_stats[rn].ref_breadth;
				row.reliable_ref_mean_depth = reliable_abundance_stats[rn].ref_mean_depth;
				row.reliable_ref_hit_ctx = reliable_abundance_stats[rn].ref_hit_ctx;
				row.reliable_ref_hit_mean_depth =
					reliable_abundance_stats[rn].ref_hit_mean_depth;
				row.reliable_ref_hit_median_depth =
					reliable_abundance_stats[rn].ref_hit_median_depth;
				row.reliable_ref_hit_depth_variance =
					reliable_abundance_stats[rn].ref_hit_depth_variance;
				row.reliable_ref_zip_af = reliable_abundance_stats[rn].ref_zip_af;
			} else {
				row.reliable_ref_breadth = row.abundance_ref_breadth;
				row.reliable_ref_mean_depth = row.abundance_ref_mean_depth;
				row.reliable_ref_hit_ctx = a->ref_ctx_hit;
				row.reliable_ref_hit_mean_depth = row.abundance_ref_hit_mean_depth;
				row.reliable_ref_hit_median_depth = 0.0;
				row.reliable_ref_hit_depth_variance = row.abundance_ref_depth_variance;
				row.reliable_ref_zip_af = row.abundance_ref_zip_af;
			}
			row.abundance_effective_depth =
				ani_readwise_effective_abundance_depth(&row);
			if (ani_opt->readwise_dual_evidence) {
				const ani_readwise_acc_t *marker_fa =
					marker_feature_acc ? &marker_feature_acc[rn] : NULL;
				row.marker_xny_ctx = marker_fa
					? ani_readwise_scaled_round_u64(marker_fa->XnY_ctx, feature_scale)
					: 0u;
				row.marker_raw_xny_ctx = marker_raw_feature_acc
					? marker_raw_feature_acc[rn].XnY_ctx
					: row.marker_xny_ctx;
				row.marker_n_diff_obj = marker_fa
					? ani_readwise_scaled_round_u64(marker_fa->N_diff_obj, feature_scale)
					: 0u;
				row.marker_n_diff_obj_section = marker_fa
					? ani_readwise_scaled_round_u64(marker_fa->N_diff_obj_section, feature_scale)
					: 0u;
				row.marker_n_mut2_ctx = marker_fa
					? ani_readwise_scaled_round_u64(marker_fa->N_mut2_ctx, feature_scale)
					: 0u;
				row.marker_ref_ctx_total = marker_ref_ctx_total
					? marker_ref_ctx_total[rn]
					: 0u;
				if (marker_abundance_stats) {
					row.marker_ref_breadth = marker_abundance_stats[rn].ref_breadth;
					row.marker_ref_mean_depth = marker_abundance_stats[rn].ref_mean_depth;
					row.marker_ref_hit_mean_depth =
						marker_abundance_stats[rn].ref_hit_mean_depth;
					row.marker_ref_depth_variance =
						marker_abundance_stats[rn].ref_depth_variance;
					row.marker_ref_depth_cv = marker_abundance_stats[rn].ref_depth_cv;
					row.marker_ref_zero_fraction =
						marker_abundance_stats[rn].ref_zero_fraction;
					row.marker_relative_depth = marker_abundance_stats[rn].relative_depth;
					row.marker_ref_zip_af = ani_zip_corrected_af(
						row.marker_ref_breadth, row.marker_ref_mean_depth);
					row.marker_ref_zip_aaf_ani =
						aaf_ani_from_containment(row.marker_ref_zip_af);
				}
				if (marker_reliable_abundance_stats) {
					row.marker_reliable_ref_breadth =
						marker_reliable_abundance_stats[rn].ref_breadth;
					row.marker_reliable_ref_mean_depth =
						marker_reliable_abundance_stats[rn].ref_mean_depth;
					row.marker_reliable_ref_hit_ctx =
						marker_reliable_abundance_stats[rn].ref_hit_ctx;
					row.marker_reliable_ref_hit_mean_depth =
						marker_reliable_abundance_stats[rn].ref_hit_mean_depth;
					row.marker_reliable_ref_hit_median_depth =
						marker_reliable_abundance_stats[rn].ref_hit_median_depth;
					row.marker_reliable_ref_hit_depth_variance =
						marker_reliable_abundance_stats[rn].ref_hit_depth_variance;
					row.marker_reliable_ref_zip_af =
						marker_reliable_abundance_stats[rn].ref_zip_af;
				} else {
					row.marker_reliable_ref_breadth = row.marker_ref_breadth;
					row.marker_reliable_ref_mean_depth = row.marker_ref_mean_depth;
					row.marker_reliable_ref_hit_ctx = row.marker_raw_xny_ctx;
					row.marker_reliable_ref_hit_mean_depth =
						row.marker_ref_hit_mean_depth;
					row.marker_reliable_ref_hit_median_depth = 0.0;
					row.marker_reliable_ref_hit_depth_variance =
						row.marker_ref_depth_variance;
					row.marker_reliable_ref_zip_af = row.marker_ref_zip_af;
				}
				row.marker_effective_depth =
					ani_readwise_effective_abundance_depth_values(
						row.marker_reliable_ref_hit_median_depth,
						row.marker_reliable_ref_mean_depth,
						row.marker_reliable_ref_zip_af);
			}
		}
		const uint32_t unique_overlap_u32 = ani_clamp_u64_to_u32(unique_overlap);
		row.mash_dist = get_mashD(Bitslen.ctx / 2, ref_ctx_total[rn],
								  qry_ctx_for_dist, unique_overlap_u32);
		row.aaf_dist = get_aafD(Bitslen.ctx / 2, ref_ctx_total[rn],
								qry_ctx_for_dist, unique_overlap_u32);
		finalize_row_selected_metric(&row, ani_opt);
		if (row.selected_ani < ani_opt->anicut)
			continue;
		kv_push(ani_row_t, survivors, row);
	}
	if (auto_readwise_abundance_report)
	{
		keep_readwise_default_calls(&survivors);
		if (kv_size(survivors))
			qsort(&kv_A(survivors, 0), kv_size(survivors),
				  sizeof(ani_row_t), cmp_readwise_default_call_desc);
	}
	else if (kv_size(survivors))
	{
		qsort(&kv_A(survivors, 0), kv_size(survivors), sizeof(ani_row_t), cmp_ani_desc);
	}
	size_t out_n = kv_size(survivors);
	if (ani_opt->ntop > 0 && (size_t)ani_opt->ntop < out_n)
		out_n = (size_t)ani_opt->ntop;
	if (ani_opt->abundance_model != ANI_ABUNDANCE_NONE)
		normalize_readwise_abundance_depth(&survivors, out_n);
	ani_write_cami_profile(ani_opt, query_path, refname, refanno, &survivors, out_n);
	for (size_t i = 0; i < out_n; ++i) {
		const ani_row_t *r = &kv_A(survivors, i);
		print_unified_detail_row(outfp, ani_opt, query_path, refname[r->rn],
								 r, annotation_at(refanno, r->rn));
	}
	ani_write_readwise_unique_sidecar(
		ani_opt, query_path, refname, refanno, ref_n, ref_uses_ctxobj96,
		index, index96, index_n, ref_ctx_total,
		unique_acc, unique_ref_ctx_cov,
		total_reads, total_density_blocks,
		density_block_mode, density_block_ctx,
		ref_infile_meta, ref_domain_profiles,
		auto_readwise_abundance_report);
	ani_write_readwise_exact_split_sidecar(
		ani_opt, query_path, refname, refanno, ref_n, ref_uses_ctxobj96,
		index, index96, index_n, ref_ctx_total,
		exact_split_acc, exact_split_ref_ctx_cov, exact_split_ref_ctx_hit_weight,
		total_reads, exact_split_total_density_blocks,
		ref_infile_meta, ref_domain_profiles,
		auto_readwise_abundance_report);
	kv_destroy(survivors);

	ani_opt->readwise_query = old_readwise;
	ani_opt->unassembled = old_unassembled;
	ani_opt->v = old_v;
	ani_opt->ctxcut = old_ctxcut;
	ani_opt->afcut = old_afcut;
	ani_opt->anicut = old_anicut;
	get_generic_dist_from_features = old_dist_fn;
	if (outfp != stdout)
		fclose(outfp);
	free(abundance_stats);
	free(reliable_abundance_stats);
	free(raw_feature_acc);
	free(filtered_feature_acc);
	free(marker_ref_ctx_total);
	free(marker_abundance_stats);
	free(marker_raw_feature_acc);
	free(marker_filtered_feature_acc);
	free(marker_reliable_abundance_stats);
	free(ctx_filter_stats);
	ani_u64_set_destroy(&qry_ctx_seen);
	ani_u64_set_destroy(&qry_ref_ctx_seen);
	free(acc);
	free(ref_hit_bits);
	free(ref_ctx_cov);
	free(ref_ctx_mindiff_cov);
	free(ref_ctx_hit_weight);
	free(unique_acc);
	free(unique_ref_hit_bits);
	free(unique_ref_ctx_cov);
	free(exact_split_acc);
	free(exact_split_ref_hit_bits);
	free(exact_split_ref_ctx_cov);
	free(exact_split_ref_ctx_hit_weight);
	free(fence);
	free(ref_ctx_total);
	if (refanno)
		free_read_from_file(refanno, (size_t)ref_n * PATHLEN);
	if (ref_infile_meta)
		free_read_from_file(ref_infile_meta, (size_t)ref_n * sizeof(ref_infile_meta[0]));
	if (ref_domain_profiles)
		free_read_from_file(ref_domain_profiles,
							(size_t)ref_n * sizeof(ref_domain_profiles[0]));
	free(ref_ctxmeta);
	free_reference_sorted_index((ctxgidobj_t *)index_mem, index_bytes, index_is_mmap);
	free_read_from_file(ref_stat, ref_stat_size);
	return 0;
}



/* comb_manysmall_sortedsketch64Xcomb_fewlarge_sortedsketch64_filter_and_sort_survivors using hash table  */
//helpers: ctx→(beg,end) hash table for query a

typedef struct {
    uint64_t *keys;   // store (ctx+1); 0 means empty
    uint64_t *vals;   // pack(beg,end)
    uint32_t  mask;   // cap-1 (cap power of 2)
    uint32_t  n_runs;
} ctxrun_ht_t;

static inline uint32_t next_pow2_u32(uint32_t x)
{
    if (x <= 1) return 1;
    --x;
    x |= x >> 1; x |= x >> 2; x |= x >> 4; x |= x >> 8; x |= x >> 16;
    return x + 1;
}

static inline uint64_t pack_be(uint32_t beg, uint32_t end)
{
    return ((uint64_t)beg << 32) | (uint64_t)end;
}
static inline uint32_t unpack_beg(uint64_t be) { return (uint32_t)(be >> 32); }
static inline uint32_t unpack_end(uint64_t be) { return (uint32_t)(be); }

static inline void ctxrun_ht_init_from_n(ctxrun_ht_t *ht, uint32_t n)
{
    // Over-allocate from n (len_qry). Very fast lookups; simple build.
    uint32_t cap = next_pow2_u32(n * 2u);         // ~0.5 load factor if nrun~n
    ht->keys = (uint64_t*)calloc(cap, sizeof(uint64_t));
    ht->vals = (uint64_t*)malloc((size_t)cap * sizeof(uint64_t));
    if (!ht->keys || !ht->vals) err(EXIT_FAILURE, "malloc ctxrun_ht");
    ht->mask = cap - 1;
    ht->n_runs = 0;
}

static inline void ctxrun_ht_destroy(ctxrun_ht_t *ht)
{
    free(ht->keys); free(ht->vals);
    ht->keys = NULL; ht->vals = NULL; ht->mask = 0; ht->n_runs = 0;
}

static inline void ctxrun_ht_put(ctxrun_ht_t *ht, uint64_t ctx, uint64_t val_be)
{
    const uint64_t k = ctx + 1; // reserve 0
    uint32_t pos = (uint32_t)mix64(k) & ht->mask;
    while (ht->keys[pos] && ht->keys[pos] != k) pos = (pos + 1) & ht->mask;
    ht->keys[pos] = k;
    ht->vals[pos] = val_be;
}

static inline int ctxrun_ht_get(const ctxrun_ht_t *ht, uint64_t ctx, uint64_t *out_be)
{
    const uint64_t k = ctx + 1;
    uint32_t pos = (uint32_t)mix64(k) & ht->mask;
    for (;;) {
        const uint64_t cur = ht->keys[pos];
        if (!cur) return 0;
        if (cur == k) { *out_be = ht->vals[pos]; return 1; }
        pos = (pos + 1) & ht->mask;
    }
}

// Build ctx->(beg,end) from sorted a[] in ONE pass.
// Must be called once per query qn; afterwards ht is read-only and thread-safe.
static inline void build_ctxrun_ht_from_sorted_query(const uint64_t *a, uint32_t n, uint8_t nobjbits, ctxrun_ht_t *ht)
{
    ctxrun_ht_init_from_n(ht, n);

    for (uint32_t i = 0; i < n; ) {
        const uint64_t ctx = a[i] >> nobjbits;
        const uint32_t beg = i;
        do { ++i; } while (i < n && (a[i] >> nobjbits) == ctx);
        const uint32_t end = i;
        ctxrun_ht_put(ht, ctx, pack_be(beg, end));
        ht->n_runs++;
    }
}

static inline uint32_t count_ctx_runs_sorted_ctxobj64(const uint64_t *a, uint32_t n, uint8_t nobjbits,
													  bool ignoreconflict)
{
    uint32_t runs = 0;
    for (uint32_t i = 0; i < n; ) {
        const uint64_t ctx = a[i] >> nobjbits;
        const uint32_t begin = i;
        do { ++i; } while (i < n && (a[i] >> nobjbits) == ctx);
        if (!ignoreconflict || i - begin == 1)
            runs++;
    }
    return runs;
}

// Return 1 if enough shared contexts remain for the report AF cutoff, else 0.
// Fills features regardless; caller can compute dist/ani.
static inline int get_features_scan_b_hash_a(
    const uint64_t *a, uint32_t n,
    const uint64_t *b, uint32_t m,
    const ctxrun_ht_t *ht,                    // built from a
	uint32_t need_X,                          // ceil(afcut * min(query ctx count, reference ctx count))
    bool ignore_ref_conflict,
    ani_features_t *f)
{
    const uint8_t  nobjbits = Bitslen.obj;
    const uint64_t objmask  = ((nobjbits == 64) ? ~0ULL : ((1ULL << nobjbits) - 1ULL));

    memset(f, 0, sizeof *f);

    for (uint32_t j = 0; j < m; ) {
        // Early fail: even if all remaining b entries hit, cannot reach need_X
        if ((uint32_t)f->XnY_ctx + (m - j) < need_X)
            return 0;

        const uint64_t ctxB = b[j] >> nobjbits;
        const uint32_t b_begin = j;
        do { ++j; } while (j < m && (b[j] >> nobjbits) == ctxB);
        const uint32_t b_end = j;

        if (ignore_ref_conflict && b_end - b_begin > 1)
            continue;

        uint64_t be;
        if (!ctxrun_ht_get(ht, ctxB, &be))
            continue;

        const uint32_t beg  = unpack_beg(be);
        const uint32_t end  = unpack_end(be);

        f->XnY_ctx++;
        const int min_diff = min_diff_sections_ctxobj64_runs(a, beg, end, b, b_begin, b_end, objmask);

        if (min_diff > 0) {
            f->N_diff_obj++;
            f->N_diff_obj_section += min_diff;
            if (min_diff > 1) f->N_mut2_ctx++;
        }
    }

    return ((uint32_t)f->XnY_ctx >= need_X);
}


void comb_manysmall_sortedsketch64Xcomb_fewlarge_sortedsketch64_filter_and_sort_survivors(ani_opt_t *ani_opt)
{
    unify_sketch_t *qry = generic_sketch_parse(ani_opt->qrydir, ani_query_parse_flags(ani_opt));
    unify_sketch_t *ref = generic_sketch_parse(ani_opt->refdir, ani_ref_parse_flags(ani_opt));
    pairwise_check_compatible(ref, qry);
    pairwise_prepare_minco_model(ref);
    load_infile_meta_for_best_guard(qry, ani_opt->qrydir, ani_opt);
    load_infile_meta_for_best_guard(ref, ani_opt->refdir, ani_opt);
    ani_ctxmeta_rec_t *qry_ctxmeta =
        read_optional_ani_ctxmeta_stats(ani_opt->qrydir, qry->infile_num);
    ani_ctxmeta_rec_t *ref_ctxmeta =
        read_optional_ani_ctxmeta_stats(ani_opt->refdir, ref->infile_num);

    const uint32_t Q = qry->infile_num;
    const uint32_t R = ref->infile_num;

    FILE *outfp = (ani_opt->outf[0] == '\0') ? stdout : fopen(ani_opt->outf, "w");
    if (!outfp) err(errno, "%s", ani_opt->outf);
    print_ani_detail_header(outfp, ani_opt, false);

    const int P = (ani_opt->p > 0) ? ani_opt->p : 1;

    // store formatted output per qn for ordered printing
    kstring_t *ks_arr = (kstring_t*)calloc(Q, sizeof(kstring_t));
    if (!ks_arr) err(EXIT_FAILURE, "calloc ks_arr");

    const uint8_t  nobjbits = Bitslen.obj;

    for (uint32_t qn = 0; qn < Q; ++qn) {

        // ---- query a (large) ----
        const uint64_t *a = qry->comb_sketch + qry->sketch_index[qn];
        const uint32_t  n = (uint32_t)(qry->sketch_index[qn + 1] - qry->sketch_index[qn]);

        // Build ctx->(beg,end) once for this qn (key speedup)
        ctxrun_ht_t ht;
        build_ctxrun_ht_from_sorted_query(a, n, nobjbits, &ht);

        // Thread-local survivor buffers
        kv_ani_row_t *tls = (kv_ani_row_t*)calloc((size_t)P, sizeof(*tls));
        if (!tls) err(EXIT_FAILURE, "calloc tls");
        for (int t = 0; t < P; ++t) kv_init(tls[t]);

        #pragma omp parallel num_threads(P)
        {
            const int tid = omp_get_thread_num();
            ani_row_t row;

            #pragma omp for schedule(dynamic, 256)
            for (uint32_t rn = 0; rn < R; ++rn) {

                // ---- ref b (small, conflict-free) ----
                const uint64_t *b = ref->comb_sketch + ref->sketch_index[rn];
                const uint32_t  m = (uint32_t)(ref->sketch_index[rn + 1] - ref->sketch_index[rn]);
                const uint32_t  m_ctx = ref->conflict ? count_ctx_runs_sorted_ctxobj64(b, m, nobjbits, ani_opt->ignoreconflict) : m;
                const uint32_t  n_ctx = qry->conflict ? ht.n_runs : n;
                if (m_ctx == 0 || n_ctx == 0)
                    continue;

                const ani_ctxmeta_rec_t *qmeta = ani_ctxmeta_at(qry_ctxmeta, qn);
                const ani_ctxmeta_rec_t *rmeta = ani_ctxmeta_at(ref_ctxmeta, rn);
				const uint32_t need_X = ani_density_af_needed_ctx(ani_opt, n_ctx, m_ctx,
                                                                  qmeta, rmeta);

                ani_features_t f;
                get_features_scan_b_hash_a(a, n, b, m, &ht, need_X, ani_opt->ignoreconflict && ref->conflict, &f);
                if (f.XnY_ctx < need_X) continue;

                const double af_q = (double)f.XnY_ctx / (double)n_ctx;
                const double af_r = (double)f.XnY_ctx / (double)m_ctx;
                const ani_density_af_t density_af = ani_estimate_density_af(
                    qmeta, rmeta, (uint32_t)f.XnY_ctx, af_q, af_r);

				if (!ani_report_af_pass(ani_opt, density_af.qry, density_af.ref)) continue;

	                ani_features_t tmp = f;
	                tmp.X_ctx = n_ctx;
	                const double blastn_af_q = lm3ways_af_ANIb_from_features(&tmp);
	                tmp = f;
	                tmp.X_ctx = m_ctx;
	                const double blastn_af_r = lm3ways_af_ANIb_from_features(&tmp);

	                row = make_selected_output_row(rn, &f, ani_opt, n_ctx, m_ctx,
	                                               af_q, blastn_af_q, af_r, blastn_af_r,
	                                               density_af,
	                                               infile_meta_at(qry, qn), infile_meta_at(ref, rn));
	                if (row.selected_ani <= ani_opt->anicut) continue;

	                kv_push(ani_row_t, tls[tid], row);
	            }
        } // omp parallel

        // merge tls to surv
        kv_ani_row_t surv; kv_init(surv);
        for (int t = 0; t < P; ++t) {
            if (kv_size(tls[t])) kv_append_rows(&surv, &tls[t]);
            kv_destroy(tls[t]);
        }
        free(tls);

        ctxrun_ht_destroy(&ht);

        // sort + format
        if (kv_size(surv))
            qsort(&kv_A(surv,0), kv_size(surv), sizeof(ani_row_t), cmp_ani_desc);

        kstring_t ks_out = (kstring_t){0,0,0};
        format_rows_to_kstr(qry, ref, qn, &surv, ani_opt, &ks_out);
        ks_arr[qn] = ks_out;

        kv_destroy(surv);
    }

    // ordered print
    for (uint32_t qn = 0; qn < Q; ++qn) {
        if (ks_arr[qn].l) fwrite(ks_arr[qn].s, 1, ks_arr[qn].l, outfp);
        free(ks_arr[qn].s);
    }
    free(ks_arr);

    if (outfp != stdout) fclose(outfp);
    free(qry_ctxmeta);
    free(ref_ctxmeta);
}
