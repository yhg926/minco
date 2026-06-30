/* Hidden sketch-set maintenance helpers for minco context-object sketches. */
#include "global_basic.h"
#include "command_ani.h"
#include "command_operate.h"
#include "command_sketch_wrapper.h"
#include "pairwise_graph.h"
#include "sketch_inspect.h"
#include "../klib/khash.h"
#include <fcntl.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

const char minco_pan_prefix[] = "lpan"; // uint64_t pan
const char minco_uniq_pan_prefix[] = "luniq_pan";
const char minco_pan96_prefix[] = "lpan.ctxobj96";
const char minco_uniq_pan96_prefix[] = "luniq_pan.ctxobj96";
// common vars
static size_t file_size;
static minco_sketch_stat_t minco_stat_readin, minco_stat_pan, minco_stat_origin;
static struct stat s;
static char outfpath[PATHLEN + 20];
static int ret;
extern const char sorted_comb_ctxgid64obj32[];

static bool operate_stat_uses_ctxobj96_payload(const minco_sketch_stat_t *stat)
{
	if (!stat)
		return false;
	const int ctx_bits = stat->coden_len > 0 ? 4 * stat->coden_len : 4 * stat->hclen;
	const int obj_bits = 2 * stat->klen - ctx_bits;
	return ctx_bits > 64 || obj_bits > 32 || ctx_bits + obj_bits > 64 ||
		   ctx_bits + GID_NBITS > 64;
}

typedef struct minco_payload_layout_spec
{
	bool use_ctxobj96;
	const char *combined_suffix;
	const char *idx_suffix;
	const char *pan_suffix;
	const char *uniq_pan_suffix;
	size_t record_size;
} minco_payload_layout_spec_t;

static minco_payload_layout_spec_t minco_payload_layout_for_stat(
	const minco_sketch_stat_t *stat)
{
	const bool use_ctxobj96 = operate_stat_uses_ctxobj96_payload(stat);
	minco_payload_layout_spec_t spec = {
		.use_ctxobj96 = use_ctxobj96,
		.combined_suffix = use_ctxobj96 ? combined_sketch96_suffix
										: combined_sketch_suffix,
		.idx_suffix = use_ctxobj96 ? idx_sketch96_suffix : idx_sketch_suffix,
		.pan_suffix = use_ctxobj96 ? minco_pan96_prefix : minco_pan_prefix,
		.uniq_pan_suffix = use_ctxobj96 ? minco_uniq_pan96_prefix
										: minco_uniq_pan_prefix,
		.record_size = use_ctxobj96 ? sizeof(ctxobj96_t) : sizeof(uint64_t),
	};
	return spec;
}

static void copy_minco_sketch_annotations(const char *indir, const char *outdir, int infile_num)
{
	if (infile_num <= 0 || !file_exists_in_folder(indir, sketch_anno_stat))
		return;

	size_t anno_file_size = 0;
	char *anno_path = test_get_fullpath(indir, sketch_anno_stat);
	char (*annotations)[PATHLEN] = read_from_file(anno_path, &anno_file_size);
	free(anno_path);

	const size_t expected_size = (size_t)infile_num * PATHLEN;
	if (anno_file_size != expected_size)
		err(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
			__func__, indir, sketch_anno_stat, anno_file_size, expected_size);

	write_to_file(test_create_fullpath(outdir, sketch_anno_stat), annotations, anno_file_size);
	free_read_from_file(annotations, anno_file_size);
}

static void copy_minco_domain_profiles(const char *indir, const char *outdir, int infile_num)
{
	if (infile_num <= 0 || !file_exists_in_folder(indir, sketch_domain_stat))
		return;

	size_t domain_file_size = 0;
	char *domain_path = test_get_fullpath(indir, sketch_domain_stat);
	uint8_t *profiles = read_from_file(domain_path, &domain_file_size);
	free(domain_path);

	const size_t expected_size = (size_t)infile_num * sizeof(profiles[0]);
	if (domain_file_size != expected_size)
		err(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
			__func__, indir, sketch_domain_stat, domain_file_size, expected_size);

	write_to_file(test_create_fullpath(outdir, sketch_domain_stat), profiles, domain_file_size);
	free_read_from_file(profiles, domain_file_size);
}

void sketch_inspect_print_samples(const char *sketch_path)
{
	void *mem_stat = read_from_file(test_get_fullpath(sketch_path, sketch_stat), &file_size);
	memcpy(&minco_stat_readin, mem_stat, sizeof(minco_stat_readin));
	char (*tmpname)[PATHLEN] = mem_stat + sizeof(minco_sketch_stat_t);
	const char *idx_suffix = operate_stat_uses_ctxobj96_payload(&minco_stat_readin)
								 ? idx_sketch96_suffix
								 : idx_sketch_suffix;
	uint64_t *mem_index = (uint64_t *)read_from_file(test_get_fullpath(sketch_path, idx_suffix), &file_size);
	for (int i = 0; i < minco_stat_readin.infile_num; i++)
		printf("%lu\t%s\n", mem_index[i + 1] - mem_index[i], tmpname[i]);
	free_all(mem_stat, mem_index, NULL);
}

static const char *inspect_ctxmeta_mode_name(uint8_t mode)
{
	switch (mode)
	{
	case MINCO_CTXMETA_PRECONFLICT:
		return "preconflict";
	case MINCO_CTXMETA_POSTCONFLICT:
		return "postconflict";
	case MINCO_CTXMETA_BOTH:
		return "both";
	case MINCO_CTXMETA_NONE:
	default:
		return "none";
	}
}

static long double inspect_threshold_density(uint64_t threshold, uint32_t hash_bits)
{
	if (hash_bits == 0 || hash_bits > 64 || threshold == UINT64_MAX)
		return 0.0L;
	return ((long double)threshold + 1.0L) / ldexpl(1.0L, (int)hash_bits);
}

void sketch_inspect_print_ctxmeta(const char *sketch_path)
{
	size_t stat_size = 0;
	char *stat_path = test_get_fullpath(sketch_path, sketch_stat);
	void *mem_stat = read_from_file(stat_path, &stat_size);
	free(stat_path);
	if (!minco_stat_decode_mem(mem_stat, stat_size, &minco_stat_readin, NULL))
		errx(EINVAL, "%s(): malformed %s/%s", __func__, sketch_path, sketch_stat);
	char (*names)[PATHLEN] = minco_stat_names_from_mem(mem_stat, stat_size);

	if (!file_exists_in_folder(sketch_path, minco_ctxmeta_bin_stat))
		errx(EXIT_FAILURE, "%s/%s does not exist; build the sketch with --ctxmeta",
			 sketch_path, minco_ctxmeta_bin_stat);
	char *ctxmeta_path = test_get_fullpath(sketch_path, minco_ctxmeta_bin_stat);
	size_t ctxmeta_size = 0;
	minco_ctxmeta_record_t *records = read_from_file(ctxmeta_path, &ctxmeta_size);
	free(ctxmeta_path);
	const size_t expected_size =
		(size_t)minco_stat_readin.infile_num * sizeof(records[0]);
	if (ctxmeta_size != expected_size)
		errx(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
			 __func__, sketch_path, minco_ctxmeta_bin_stat,
			 ctxmeta_size, expected_size);

	puts("sample_id\tsample_path\tmode\tvalid\thash_bits\tthreshold\tsketch_entries"
		 "\tselected_observed_ctx\tselected_estimated_unique_ctx"
		 "\tpreconflict_observed_ctx\tpreconflict_estimated_unique_ctx"
		 "\tpostconflict_observed_ctx\tpostconflict_estimated_unique_ctx");
	for (int i = 0; i < minco_stat_readin.infile_num; ++i)
	{
		const minco_ctxmeta_record_t *r = &records[i];
		printf("%d\t%s\t%s\t%u\t%u\t%" PRIu64 "\t%" PRIu64
			   "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64
			   "\t%" PRIu64 "\t%" PRIu64 "\n",
			   i,
			   names[i],
			   inspect_ctxmeta_mode_name(r->mode),
			   (unsigned)r->valid,
			   r->hash_bits,
			   r->threshold,
			   r->sketch_entries,
			   r->selected_observed_ctx,
			   r->selected_estimated_unique_ctx,
			   r->preconflict_observed_ctx,
			   r->preconflict_estimated_unique_ctx,
			   r->postconflict_observed_ctx,
			   r->postconflict_estimated_unique_ctx);
	}
	free_read_from_file(records, ctxmeta_size);
	free_read_from_file(mem_stat, stat_size);
}

static const char *inspect_density_policy_name(uint32_t policy)
{
	switch (policy)
	{
	case MINCO_STAT_DENSITY_POLICY_SINGLE_SAMPLE:
		return "single_sample_density";
	case MINCO_STAT_DENSITY_POLICY_LARGEST_SAMPLE:
		return "largest_sample_density";
	case MINCO_STAT_DENSITY_POLICY_EXPLICIT_THRESHOLD:
		return "explicit_threshold";
	case MINCO_STAT_DENSITY_POLICY_NONE:
	default:
		return "none";
	}
}

static const char *inspect_sample_name_or_empty(char (*names)[PATHLEN], int infile_num,
												uint32_t sample_id)
{
	if (!names || sample_id == MINCO_STAT_DENSITY_SAMPLE_ID_NONE ||
		sample_id >= (uint32_t)infile_num)
		return "";
	return names[sample_id];
}

void sketch_inspect_print_ctxsetmeta(const char *sketch_path)
{
	size_t stat_size = 0;
	char *stat_path = test_get_fullpath(sketch_path, sketch_stat);
	void *mem_stat = read_from_file(stat_path, &stat_size);
	free(stat_path);
	minco_sketch_info_t info = {0};
	if (!minco_stat_decode_mem(mem_stat, stat_size, &minco_stat_readin, &info))
		errx(EINVAL, "%s(): malformed %s/%s", __func__, sketch_path, sketch_stat);
	if (!info.has_minco_ext || info.density_valid_sample_count == 0)
		errx(EXIT_FAILURE, "%s has no density summary; build the sketch with --ctxmeta",
			 sketch_path);
	char (*names)[PATHLEN] = minco_stat_names_from_mem(mem_stat, stat_size);

	puts("key\tvalue");
	printf("meta_version\t1\n");
	printf("sample_count\t%d\n", minco_stat_readin.infile_num);
	printf("valid_sample_count\t%u\n", info.density_valid_sample_count);
	printf("hash_bits\t%u\n", info.density_hash_bits);
	printf("hash_bits_mixed\t%u\n",
		   (info.density_flags & MINCO_STAT_DENSITY_FLAG_MIXED_HASH_BITS) ? 1u : 0u);
	printf("universal_density_policy\t%s\n",
		   inspect_density_policy_name(info.density_universal_policy));
	printf("universal_sample_id\t%u\n", info.density_universal_sample_id);
	printf("universal_sample_path\t%s\n",
		   inspect_sample_name_or_empty(names, minco_stat_readin.infile_num,
										info.density_universal_sample_id));
	printf("universal_hash_bits\t%u\n", info.density_hash_bits);
	printf("universal_threshold\t%" PRIu64 "\n", info.density_universal_threshold);
	printf("universal_density\t%.18Le\n",
		   inspect_threshold_density(info.density_universal_threshold,
									 info.density_hash_bits));
	printf("min_sample_id\t%u\n", info.density_min_sample_id);
	printf("min_sample_path\t%s\n",
		   inspect_sample_name_or_empty(names, minco_stat_readin.infile_num,
										info.density_min_sample_id));
	printf("min_sample_hash_bits\t%u\n", info.density_hash_bits);
	printf("min_sample_threshold\t%" PRIu64 "\n", info.density_min_threshold);
	printf("min_sample_density\t%.18Le\n",
		   inspect_threshold_density(info.density_min_threshold,
									 info.density_hash_bits));
	printf("largest_sample_id\t%u\n", info.density_max_sample_id);
	printf("largest_sample_path\t%s\n",
		   inspect_sample_name_or_empty(names, minco_stat_readin.infile_num,
										info.density_max_sample_id));
	printf("largest_sample_hash_bits\t%u\n", info.density_hash_bits);
	printf("largest_sample_threshold\t%" PRIu64 "\n", info.density_max_threshold);
	printf("largest_sample_density\t%.18Le\n",
		   inspect_threshold_density(info.density_max_threshold,
									 info.density_hash_bits));
	free_read_from_file(mem_stat, stat_size);
}

void print_minco_sample_names(set_opt_t *set_opt)
{
	sketch_inspect_print_samples(set_opt->insketchpath);
}

void sketch_inspect_print_content(const char *sketch_path, int show_mode)
{
	FILE *fp = NULL, *fab = NULL, *fpos = NULL;
	uint64_t kmer;
	uint64_t pos;
	uint32_t abun;
	ctxgidobj_t ctxgidobj;
	bool abundance = 0;
	void *stat_mem = read_from_file(test_get_fullpath(sketch_path, sketch_stat), &file_size);
	memcpy(&minco_stat_readin, stat_mem, sizeof(minco_stat_readin));
	const bool use_ctxobj96 = operate_stat_uses_ctxobj96_payload(&minco_stat_readin);
	free(stat_mem);
	if (show_mode == 1)
	{
		const char *idx_suffix = use_ctxobj96 ? idx_sketch96_suffix : idx_sketch_suffix;
		uint64_t *sketch_index = read_from_file(test_get_fullpath(sketch_path, idx_suffix), &file_size);
		int infile_num = file_size / sizeof(uint64_t) - 1;
		if (use_ctxobj96)
		{
			if ((fp = fopen(test_get_fullpath(sketch_path, combined_sketch96_suffix), "rb")) == NULL)
				err(EXIT_FAILURE, "%s(): Failed to open file '%s/%s'", __func__, sketch_path, combined_sketch96_suffix);
			ctxobj96_t rec;
			for (int i = 0; i < infile_num; i++)
			{
				for (uint64_t j = sketch_index[i]; j < sketch_index[i + 1]; j++)
				{
					if (fread(&rec, sizeof(rec), 1, fp) != 1)
						err(EXIT_FAILURE, "%s(): Failed to read %s", __func__, combined_sketch96_suffix);
					printf("%d\t%016" PRIx64 "\t%08x\n", i, rec.ctx, rec.obj);
				}
			}
			free(sketch_index);
			goto cleanup;
		}
		if ((fp = fopen(test_get_fullpath(sketch_path, combined_sketch_suffix), "rb")) == NULL)
			err(EXIT_FAILURE, "%s(): Failed to open file '%s/%s'", __func__, sketch_path, combined_sketch_suffix);

		if (file_exists_in_folder(sketch_path, combined_ab_suffix))
		{
			if ((fab = fopen(test_get_fullpath(sketch_path, combined_ab_suffix), "rb")) == NULL)
				err(EXIT_FAILURE, "%s(): Failed to open file '%s/%s'", __func__, sketch_path, combined_ab_suffix);
			abundance = 1;
		}

		for (int i = 0; i < infile_num; i++)
		{
			for (uint64_t j = sketch_index[i]; j < sketch_index[i + 1]; j++)
			{
				fread(&kmer, sizeof(kmer), 1, fp);
				if (abundance)
				{
					fread(&abun, sizeof(abun), 1, fab);
					printf("%d\t%lx\t%u\n", i, kmer, abun);
				}
				else
					printf("%d\t%lx\n", i, kmer);
			}
		}
		free(sketch_index);
	}
	else if (show_mode == 2)
	{
		if (use_ctxobj96)
		{
			if ((fp = fopen(test_get_fullpath(sketch_path, sorted_comb_ctx64gid32obj32), "rb")) == NULL)
				err(EXIT_FAILURE, "%s(): Failed to open file '%s/%s'", __func__, sketch_path, sorted_comb_ctx64gid32obj32);
			ctxgidobj128_t rec;
			while (fread(&rec, sizeof(rec), 1, fp) == 1)
				printf("%016" PRIx64 "\t%x\t%08x\n", rec.ctx, rec.gid, rec.obj);
			goto cleanup;
		}
		if ((fp = fopen(test_get_fullpath(sketch_path, sorted_comb_ctxgid64obj32), "rb")) == NULL)
			err(EXIT_FAILURE, "%s(): Failed to open file '%s/%s'", __func__, sketch_path, sorted_comb_ctxgid64obj32);
		while (!feof(fp))
		{
			fread(&ctxgidobj, sizeof(ctxgidobj), 1, fp);
			printf("%lx\t%x\n", ctxgidobj.ctxgid, ctxgidobj.obj);
		}
	}
	else if (show_mode == 3)
	{
		if (use_ctxobj96)
			errx(EXIT_FAILURE, "%s(): ctxobj96 sketches do not support legacy ctxobj64 positions", __func__);
		uint64_t *sketch_index = read_from_file(test_get_fullpath(sketch_path, idx_sketch_suffix), &file_size);
		int infile_num = file_size / sizeof(uint64_t) - 1;
		const uint64_t sketch_entries = sketch_index[infile_num];
		if ((fp = fopen(test_get_fullpath(sketch_path, combined_sketch_suffix), "rb")) == NULL)
			err(EXIT_FAILURE, "%s(): Failed to open file '%s/%s'", __func__, sketch_path, combined_sketch_suffix);
		if ((fpos = fopen(test_get_fullpath(sketch_path, sketch_position_suffix), "rb")) == NULL)
			err(EXIT_FAILURE, "%s(): Failed to open file '%s/%s'; build this sketch with --position",
				__func__, sketch_path, sketch_position_suffix);

		struct stat comb_st;
		struct stat pos_st;
		if (fstat(fileno(fp), &comb_st) != 0)
			err(errno, "%s(): fstat %s/%s", __func__, sketch_path, combined_sketch_suffix);
		if (fstat(fileno(fpos), &pos_st) != 0)
			err(errno, "%s(): fstat %s/%s", __func__, sketch_path, sketch_position_suffix);
		const off_t expected_size = (off_t)(sketch_entries * sizeof(uint64_t));
		if (comb_st.st_size != expected_size)
			err(EINVAL, "%s(): %s/%s has %ld bytes, expected %ld",
				__func__, sketch_path, combined_sketch_suffix,
				(long)comb_st.st_size, (long)expected_size);
		if (pos_st.st_size != expected_size)
			err(EINVAL, "%s(): %s/%s has %ld bytes, expected %ld",
				__func__, sketch_path, sketch_position_suffix,
				(long)pos_st.st_size, (long)expected_size);

		for (int i = 0; i < infile_num; i++)
		{
			for (uint64_t j = sketch_index[i]; j < sketch_index[i + 1]; j++)
			{
				if (fread(&kmer, sizeof(kmer), 1, fp) != 1)
					err(EXIT_FAILURE, "%s(): Failed to read %s", __func__, combined_sketch_suffix);
				if (fread(&pos, sizeof(pos), 1, fpos) != 1)
					err(EXIT_FAILURE, "%s(): Failed to read %s", __func__, sketch_position_suffix);
				printf("%d\t%lx\t%lu\n", i, kmer, pos);
			}
		}
		free(sketch_index);
	}
	else
		err(EXIT_FAILURE, "%s(): only show modes 1, 2 and 3 are supported, show_mode =%d", __func__, show_mode);
cleanup:
	if (fp)
		fclose(fp);
	if (abundance)
		fclose(fab);
	if (fpos)
		fclose(fpos);
}

void show_content(set_opt_t *set_opt)
{
	sketch_inspect_print_content(set_opt->insketchpath, set_opt->show);
}

KHASH_SET_INIT_INT64(kmer_set)

static inline void u64_swap(uint64_t *a, uint64_t *b)
{
	uint64_t tmp = *a;
	*a = *b;
	*b = tmp;
}

static void u64_insertion_sort(uint64_t *arr, size_t n)
{
	for (size_t i = 1; i < n; ++i)
	{
		uint64_t key = arr[i];
		size_t j = i;
		while (j > 0 && key < arr[j - 1])
		{
			arr[j] = arr[j - 1];
			--j;
		}
		arr[j] = key;
	}
}

static size_t u64_partition(uint64_t *arr, size_t n)
{
	uint64_t *left = arr;
	uint64_t *mid = arr + n / 2;
	uint64_t *right = arr + n - 1;
	if (*mid < *left)
		u64_swap(mid, left);
	if (*right < *left)
		u64_swap(right, left);
	if (*mid < *right)
		u64_swap(mid, right);
	const uint64_t pivot = *right;

	size_t i = 0;
	size_t j = n - 1;
	while (1)
	{
		while (arr[i] < pivot)
			++i;
		while (pivot < arr[j])
			--j;
		if (i >= j)
			return j;
		u64_swap(&arr[i], &arr[j]);
		++i;
		--j;
	}
}

static void u64_parallel_quicksort_task(uint64_t *arr, size_t n)
{
	const size_t insertion_threshold = 32;
	const size_t parallel_threshold = (size_t)1 << 20;
	if (n <= insertion_threshold)
	{
		u64_insertion_sort(arr, n);
		return;
	}

	size_t p = u64_partition(arr, n);
#ifdef _OPENMP
#pragma omp task shared(arr) if (p + 1 > parallel_threshold)
#endif
	u64_parallel_quicksort_task(arr, p + 1);
#ifdef _OPENMP
#pragma omp task shared(arr) if (n - p - 1 > parallel_threshold)
#endif
	u64_parallel_quicksort_task(arr + p + 1, n - p - 1);
#ifdef _OPENMP
#pragma omp taskwait
#endif
}

static void sort_uint64_values(uint64_t *arr, size_t n, int threads)
{
#ifdef _OPENMP
	if (threads < 1)
		threads = 1;
#pragma omp parallel num_threads(threads)
	{
#pragma omp single nowait
		u64_parallel_quicksort_task(arr, n);
	}
#else
	u64_parallel_quicksort_task(arr, n);
#endif
}

static size_t compact_union_or_unique_uint64(uint64_t *sorted, size_t n, bool unique_only)
{
	size_t out = 0;
	for (size_t i = 0; i < n;)
	{
		size_t j = i + 1;
		while (j < n && sorted[j] == sorted[i])
			++j;
		if (!unique_only || j == i + 1)
			sorted[out++] = sorted[i];
		i = j;
	}
	return out;
}

static void kh_resize_for_size(khash_t(kmer_set) *h, size_t n)
{
	if (n == 0)
		return;
	const size_t requested = n + n / 3 + 1024;
	if (requested <= (size_t)UINT32_MAX)
		kh_resize(kmer_set, h, (khint_t)requested);
}

static size_t *build_u64_highbit_bucket_offsets(const uint64_t *sorted, size_t n, unsigned bits)
{
	const size_t bucket_count = (size_t)1 << bits;
	size_t *offsets = calloc(bucket_count + 1, sizeof(offsets[0]));
	if (!offsets)
		err(errno, "%s(): OOM markerdb lookup buckets", __func__);

	size_t bucket = 0;
	for (size_t i = 0; i < n; ++i)
	{
		const size_t value_bucket = sorted[i] >> (64 - bits);
		while (bucket < value_bucket)
			offsets[++bucket] = i;
	}
	while (bucket < bucket_count)
		offsets[++bucket] = n;
	return offsets;
}

static inline bool u64_bucketed_contains(
	const uint64_t *sorted,
	const size_t *bucket_offsets,
	unsigned bits,
	uint64_t value)
{
	const size_t bucket = value >> (64 - bits);
	size_t lo = bucket_offsets[bucket];
	size_t hi = bucket_offsets[bucket + 1];
	while (lo < hi)
	{
		const size_t mid = lo + (hi - lo) / 2;
		if (sorted[mid] < value)
			lo = mid + 1;
		else
			hi = mid;
	}
	return lo < bucket_offsets[bucket + 1] && sorted[lo] == value;
}

static inline int ctxobj96_compare_value(ctxobj96_t a, ctxobj96_t b)
{
	if (a.ctx != b.ctx)
		return (a.ctx > b.ctx) - (a.ctx < b.ctx);
	return (a.obj > b.obj) - (a.obj < b.obj);
}

static size_t compact_union_or_unique_ctxobj96(ctxobj96_t *sorted, size_t n,
											   bool unique_only)
{
	size_t out = 0;
	for (size_t i = 0; i < n;)
	{
		size_t j = i + 1;
		while (j < n && ctxobj96_compare_value(sorted[j], sorted[i]) == 0)
			++j;
		if (!unique_only || j == i + 1)
			sorted[out++] = sorted[i];
		i = j;
	}
	return out;
}

static inline size_t ctxobj96_bucket_id(uint64_t ctx, unsigned bits)
{
	if (bits == 0)
		return 0;
	return (size_t)(ctx >> (64u - bits));
}

static size_t *build_ctxobj96_bucket_offsets(const ctxobj96_t *sorted,
											 size_t n,
											 unsigned bits)
{
	const size_t bucket_count = (size_t)1 << bits;
	size_t *offsets = calloc(bucket_count + 1, sizeof(offsets[0]));
	if (!offsets)
		err(errno, "%s(): OOM ctxobj96 lookup buckets", __func__);

	size_t bucket = 0;
	for (size_t i = 0; i < n; ++i)
	{
		const size_t value_bucket = ctxobj96_bucket_id(sorted[i].ctx, bits);
		while (bucket < value_bucket)
			offsets[++bucket] = i;
	}
	while (bucket < bucket_count)
		offsets[++bucket] = n;
	return offsets;
}

static inline bool ctxobj96_bucketed_contains(
	const ctxobj96_t *sorted,
	const size_t *bucket_offsets,
	unsigned bits,
	ctxobj96_t value)
{
	const size_t bucket = ctxobj96_bucket_id(value.ctx, bits);
	size_t lo = bucket_offsets[bucket];
	size_t hi = bucket_offsets[bucket + 1];
	while (lo < hi)
	{
		const size_t mid = lo + (hi - lo) / 2;
		if (ctxobj96_compare_value(sorted[mid], value) < 0)
			lo = mid + 1;
		else
			hi = mid;
	}
	return lo < bucket_offsets[bucket + 1] &&
		   ctxobj96_compare_value(sorted[lo], value) == 0;
}

static uint64_t ctxobj96_count_ctx_runs(const ctxobj96_t *arr, size_t n)
{
	uint64_t count = 0;
	for (size_t i = 0; i < n;)
	{
		const uint64_t ctx = arr[i].ctx;
		do
		{
			++i;
		} while (i < n && arr[i].ctx == ctx);
		++count;
	}
	return count;
}

static uint32_t ctxobj96_count_shared_contexts(const ctxobj96_t *a, size_t n,
											   const ctxobj96_t *b, size_t m)
{
	size_t i = 0;
	size_t j = 0;
	uint64_t count = 0;
	while (i < n && j < m)
	{
		const uint64_t actx = a[i].ctx;
		const uint64_t bctx = b[j].ctx;
		size_t an = i + 1;
		while (an < n && a[an].ctx == actx)
			an++;
		size_t bn = j + 1;
		while (bn < m && b[bn].ctx == bctx)
			bn++;
		if (actx < bctx)
		{
			i = an;
			continue;
		}
		if (bctx < actx)
		{
			j = bn;
			continue;
		}
		if (count < UINT32_MAX)
			++count;
		i = an;
		j = bn;
	}
	return count > UINT32_MAX ? UINT32_MAX : (uint32_t)count;
}

typedef struct ctxobj96_hash_record
{
	uint64_t h;
	ctxobj96_t rec;
} ctxobj96_hash_record_t;

static int ctxobj96_hash_record_cmp(const void *pa, const void *pb)
{
	const ctxobj96_hash_record_t *a = pa;
	const ctxobj96_hash_record_t *b = pb;
	if (a->h != b->h)
		return (a->h > b->h) - (a->h < b->h);
	return ctxobj96_compare_value(a->rec, b->rec);
}

static uint64_t ctxobj96_hash_value(ctxobj96_t rec)
{
	return mix64(rec.ctx ^ (uint64_t)MINCO_SEED);
}

static void ctxobj96_select_bottom_by_hash(ctxobj96_t *arr, size_t n, size_t keep)
{
	if (keep >= n)
	{
		ctxobj96_sort_array(arr, n);
		return;
	}
	ctxobj96_hash_record_t *tmp = malloc(n * sizeof(tmp[0]));
	if (!tmp)
		err(errno, "%s(): OOM ctxobj96 hash buffer", __func__);
	for (size_t i = 0; i < n; ++i)
	{
		tmp[i].h = ctxobj96_hash_value(arr[i]);
		tmp[i].rec = arr[i];
	}
	qsort(tmp, n, sizeof(tmp[0]), ctxobj96_hash_record_cmp);
	for (size_t i = 0; i < keep; ++i)
		arr[i] = tmp[i].rec;
	free(tmp);
	ctxobj96_sort_array(arr, keep);
}

static size_t ctx_bucket_id(uint64_t ctx, unsigned bucket_bits, unsigned ctx_bits)
{
	if (bucket_bits == 0)
		return 0;
	if (ctx_bits > bucket_bits)
		return (size_t)(ctx >> (ctx_bits - bucket_bits));
	if (ctx_bits < bucket_bits)
		return (size_t)(ctx << (bucket_bits - ctx_bits));
	return (size_t)ctx;
}

static size_t *build_ctx_bucket_offsets(const uint64_t *sorted_ctx, size_t n,
										unsigned bucket_bits, unsigned ctx_bits)
{
	const size_t bucket_count = (size_t)1 << bucket_bits;
	size_t *offsets = calloc(bucket_count + 1, sizeof(offsets[0]));
	if (!offsets)
		err(errno, "%s(): OOM context markerdb lookup buckets", __func__);

	size_t bucket = 0;
	for (size_t i = 0; i < n; ++i)
	{
		const size_t value_bucket = ctx_bucket_id(sorted_ctx[i], bucket_bits, ctx_bits);
		while (bucket < value_bucket)
			offsets[++bucket] = i;
	}
	while (bucket < bucket_count)
		offsets[++bucket] = n;
	return offsets;
}

static inline bool ctx_bucketed_contains(
	const uint64_t *sorted_ctx,
	const size_t *bucket_offsets,
	unsigned bucket_bits,
	unsigned ctx_bits,
	uint64_t ctx)
{
	const size_t bucket = ctx_bucket_id(ctx, bucket_bits, ctx_bits);
	size_t lo = bucket_offsets[bucket];
	size_t hi = bucket_offsets[bucket + 1];
	while (lo < hi)
	{
		const size_t mid = lo + (hi - lo) / 2;
		if (sorted_ctx[mid] < ctx)
			lo = mid + 1;
		else
			hi = mid;
	}
	return lo < bucket_offsets[bucket + 1] && sorted_ctx[lo] == ctx;
}

static void pwrite_all(int fd, const void *data, size_t bytes, off_t offset, const char *path)
{
	const char *cursor = data;
	while (bytes > 0)
	{
		const ssize_t written = pwrite(fd, cursor, bytes, offset);
		if (written < 0)
		{
			if (errno == EINTR)
				continue;
			err(errno, "%s(): write %s", __func__, path);
		}
		if (written == 0)
			err(EIO, "%s(): short write %s", __func__, path);
		cursor += written;
		offset += written;
		bytes -= (size_t)written;
	}
}

static void minco_warn_markerdb_small_refs(const set_opt_t *set_opt,
										   const uint64_t *post_fco_pos,
										   const void *mem_stat,
										   size_t stat_size,
										   const char *mode_label);

typedef struct markerdb_ctx_pair_edge_collect
{
	pairwise_edge_list_t edges;
} markerdb_ctx_pair_edge_collect_t;

static void markerdb_ctx_pair_edge_cb(void *ctx, int qry, int ref,
									  const pairwise_eval_t *eval)
{
	(void)eval;
	markerdb_ctx_pair_edge_collect_t *collector = ctx;
	pairwise_edge_list_add(&collector->edges, qry, ref);
}

static char *minco_build_temp_ctxgid_index(const uint64_t *fco_pos,
										   const uint64_t *mem_ctxobj,
										   uint32_t infile_num,
										   size_t in_entry_ct)
{
	char tmpl[] = "/tmp/minco_markerdb_ctxpair_index.XXXXXX";
	int fd = mkstemp(tmpl);
	if (fd < 0)
		err(errno, "%s(): create temporary sorted context index", __func__);

	ctxgidobj_t *ctxgid = ctxobj64_2ctxgidobj(
		(uint64_t *)fco_pos, (uint64_t *)mem_ctxobj, (int)infile_num,
		(uint32_t)in_entry_ct);
	ctxgidobj_sort_array(ctxgid, in_entry_ct);

	const uint64_t bytes_u64 = (uint64_t)in_entry_ct * (uint64_t)sizeof(ctxgid[0]);
	if (in_entry_ct != 0 && bytes_u64 / in_entry_ct != sizeof(ctxgid[0]))
		errx(EINVAL, "%s(): temporary sorted index byte size overflow", __func__);
	if (ftruncate(fd, (off_t)bytes_u64) != 0)
		err(errno, "%s(): resize %s", __func__, tmpl);
	pwrite_all(fd, ctxgid, (size_t)bytes_u64, 0, tmpl);
	if (close(fd) != 0)
		err(errno, "%s(): close %s", __func__, tmpl);
	free(ctxgid);

	char *path = strdup(tmpl);
	if (!path)
		err(errno, "%s(): OOM temporary index path", __func__);
	return path;
}

static void minco_mark_pair_shared_contexts(uint8_t *drop_entry,
											const uint64_t *mem_ctxobj,
											const uint64_t *fco_pos,
											uint32_t a,
											uint32_t b)
{
	uint64_t ai = fco_pos[a];
	const uint64_t ae = fco_pos[a + 1];
	uint64_t bi = fco_pos[b];
	const uint64_t be = fco_pos[b + 1];
	while (ai < ae && bi < be)
	{
		const uint64_t actx = mem_ctxobj[ai] >> Bitslen.obj;
		const uint64_t bctx = mem_ctxobj[bi] >> Bitslen.obj;
		uint64_t an = ai + 1;
		while (an < ae && (mem_ctxobj[an] >> Bitslen.obj) == actx)
			an++;
		uint64_t bn = bi + 1;
		while (bn < be && (mem_ctxobj[bn] >> Bitslen.obj) == bctx)
			bn++;

		if (actx < bctx)
		{
			ai = an;
			continue;
		}
		if (bctx < actx)
		{
			bi = bn;
			continue;
		}
		for (uint64_t p = ai; p < an; ++p)
			drop_entry[p] = 1;
		for (uint64_t p = bi; p < bn; ++p)
			drop_entry[p] = 1;
		ai = an;
		bi = bn;
	}
}

static void minco_mark_pair_shared_contexts96(uint8_t *drop_entry,
											 const ctxobj96_t *mem_ctxobj,
											 const uint64_t *fco_pos,
											 uint32_t a,
											 uint32_t b)
{
	uint64_t ai = fco_pos[a];
	const uint64_t ae = fco_pos[a + 1];
	uint64_t bi = fco_pos[b];
	const uint64_t be = fco_pos[b + 1];
	while (ai < ae && bi < be)
	{
		const uint64_t actx = mem_ctxobj[ai].ctx;
		const uint64_t bctx = mem_ctxobj[bi].ctx;
		uint64_t an = ai + 1;
		while (an < ae && mem_ctxobj[an].ctx == actx)
			an++;
		uint64_t bn = bi + 1;
		while (bn < be && mem_ctxobj[bn].ctx == bctx)
			bn++;

		if (actx < bctx)
		{
			ai = an;
			continue;
		}
		if (bctx < actx)
		{
			bi = bn;
			continue;
		}
		for (uint64_t p = ai; p < an; ++p)
			drop_entry[p] = 1;
		for (uint64_t p = bi; p < bn; ++p)
			drop_entry[p] = 1;
		ai = an;
		bi = bn;
	}
}

static size_t ctxgidobj128_lower_bound_ctx(const ctxgidobj128_t *arr,
										   size_t n,
										   uint64_t ctx)
{
	size_t lo = 0;
	size_t hi = n;
	while (lo < hi)
	{
		const size_t mid = lo + (hi - lo) / 2;
		if (arr[mid].ctx < ctx)
			lo = mid + 1;
		else
			hi = mid;
	}
	return lo;
}

static size_t ctxgidobj128_upper_bound_ctx(const ctxgidobj128_t *arr,
										   size_t n,
										   uint64_t ctx,
										   size_t lo)
{
	size_t hi = n;
	while (lo < hi)
	{
		const size_t mid = lo + (hi - lo) / 2;
		if (arr[mid].ctx <= ctx)
			lo = mid + 1;
		else
			hi = mid;
	}
	return lo;
}

static void minco_collect_pairwise_context_edges96(
	const set_opt_t *set_opt,
	const uint64_t *fco_pos,
	const ctxobj96_t *mem_ctxobj,
	uint32_t n,
	size_t in_entry_ct,
	markerdb_ctx_pair_edge_collect_t *edge_collect,
	pairwise_index_scan_stats_t *scan_stats)
{
	if (!edge_collect || !scan_stats)
		errx(EINVAL, "%s(): invalid output arguments", __func__);

	ctxgidobj128_t *index = ctxobj96_2ctxgidobj128(
		(uint64_t *)fco_pos, (ctxobj96_t *)mem_ctxobj, (int)n,
		(uint32_t)in_entry_ct);
	ctxgidobj128_sort_array(index, in_entry_ct);

	uint32_t *ctx_counts = calloc((size_t)n, sizeof(ctx_counts[0]));
	if (!ctx_counts)
		err(errno, "%s(): OOM context count table", __func__);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(set_opt->p)
#endif
	for (uint32_t i = 0; i < n; ++i)
	{
		const size_t len = (size_t)(fco_pos[i + 1] - fco_pos[i]);
		const uint64_t count = ctxobj96_count_ctx_runs(mem_ctxobj + fco_pos[i], len);
		ctx_counts[i] = count > UINT32_MAX ? UINT32_MAX : (uint32_t)count;
	}

	const uint32_t index_min_votes =
		set_opt->markerdb_ctx_index_min_votes > 0
			? set_opt->markerdb_ctx_index_min_votes
			: (set_opt->markerdb_ctx_min_xny > 0
				   ? set_opt->markerdb_ctx_min_xny
				   : 1);
	const uint32_t sample_step =
		set_opt->markerdb_ctx_index_sample_step > 0
			? set_opt->markerdb_ctx_index_sample_step
			: 1;
	const int threads = set_opt->p > 0 ? set_opt->p : 1;

	uint64_t candidate_pairs = 0;
	uint64_t exact_pairs = 0;
	uint64_t accepted_edges = 0;
	uint64_t ctx_rejects = 0;
	uint64_t max_af_rejects = 0;

#ifdef _OPENMP
#pragma omp parallel num_threads(threads) reduction(+:candidate_pairs,exact_pairs,accepted_edges,ctx_rejects,max_af_rejects)
#endif
	{
		uint32_t *votes = calloc((size_t)n, sizeof(votes[0]));
		uint32_t *touched = malloc((size_t)n * sizeof(touched[0]));
		if (!votes || !touched)
			err(errno, "%s(): OOM pairwise ctxobj96 scan buffers", __func__);

#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
		for (uint32_t qn = 0; qn < n; ++qn)
		{
			size_t touched_n = 0;
			const uint64_t q_begin = fco_pos[qn];
			const uint64_t q_end = fco_pos[qn + 1];
			for (uint64_t qi = q_begin; qi < q_end; qi += sample_step)
			{
				const uint64_t ctx = mem_ctxobj[qi].ctx;
				const size_t first = ctxgidobj128_lower_bound_ctx(index, in_entry_ct, ctx);
				if (first == in_entry_ct || index[first].ctx != ctx)
					continue;
				const size_t ctx_end =
					ctxgidobj128_upper_bound_ctx(index, in_entry_ct, ctx, first);
				if (set_opt->markerdb_ctx_index_max_ctx_freq > 0 &&
					ctx_end - first > set_opt->markerdb_ctx_index_max_ctx_freq)
					continue;
				for (size_t d = first; d < ctx_end; ++d)
				{
					const uint32_t gid = index[d].gid;
					if (gid >= qn)
						break;
					if (votes[gid] == 0)
						touched[touched_n++] = gid;
					if (votes[gid] < UINT32_MAX)
						++votes[gid];
				}
			}

			for (size_t ti = 0; ti < touched_n; ++ti)
			{
				const uint32_t rn = touched[ti];
				++candidate_pairs;
				if (votes[rn] < index_min_votes)
				{
					votes[rn] = 0;
					continue;
				}
				++exact_pairs;
				const ctxobj96_t *qry = mem_ctxobj + fco_pos[qn];
				const size_t qry_n = (size_t)(fco_pos[qn + 1] - fco_pos[qn]);
				const ctxobj96_t *ref = mem_ctxobj + fco_pos[rn];
				const size_t ref_n = (size_t)(fco_pos[rn + 1] - fco_pos[rn]);
				const uint32_t xny =
					ctxobj96_count_shared_contexts(qry, qry_n, ref, ref_n);
				if (xny < set_opt->markerdb_ctx_min_xny)
				{
					++ctx_rejects;
					votes[rn] = 0;
					continue;
				}
				const double af_qry = ctx_counts[qn] > 0
										  ? (double)xny / (double)ctx_counts[qn]
										  : 0.0;
				const double af_ref = ctx_counts[rn] > 0
										  ? (double)xny / (double)ctx_counts[rn]
										  : 0.0;
				const double max_af = af_qry > af_ref ? af_qry : af_ref;
				if (max_af < set_opt->markerdb_ctx_min_af)
				{
					++max_af_rejects;
					votes[rn] = 0;
					continue;
				}
				++accepted_edges;
				pairwise_eval_t eval = {
					.distance = 0.0,
					.similarity = 1.0,
					.xny_ctx = xny,
					.af_qry = af_qry,
					.af_ref = af_ref,
					.max_af = max_af,
					.valid = true,
				};
#ifdef _OPENMP
#pragma omp critical(markerdb_ctxobj96_edge_collect)
#endif
				markerdb_ctx_pair_edge_cb(edge_collect, (int)qn, (int)rn, &eval);
				votes[rn] = 0;
			}
		}

		free(votes);
		free(touched);
	}

	scan_stats->candidate_pairs = candidate_pairs;
	scan_stats->exact_pairs = exact_pairs;
	scan_stats->accepted_edges = accepted_edges;
	scan_stats->distance_edges = accepted_edges;
	scan_stats->ctx_rejects = ctx_rejects;
	scan_stats->max_af_rejects = max_af_rejects;

	free(ctx_counts);
	free(index);
}

static void minco_write_pairwise_context_markerdb96(set_opt_t *set_opt,
													const void *mem_stat,
													size_t stat_file_size)
{
	const_comask_init(&minco_stat_readin);
	if (Bitslen.ctx == 0 || Bitslen.ctx > 64)
		errx(EXIT_FAILURE, "%s(): unsupported context bit width: %u",
			 __func__, (unsigned)Bitslen.ctx);
	if (minco_stat_readin.infile_num < 0)
		errx(EXIT_FAILURE, "%s(): invalid sample count %d",
			 __func__, minco_stat_readin.infile_num);
	if (minco_stat_readin.conflict)
		errx(EXIT_FAILURE, "%s(): --markerdb-ctx-pairwise currently requires non-conflict reference sketches",
			 __func__);

	size_t idx_file_size = 0;
	uint64_t *fco_pos = (uint64_t *)read_from_file(
		test_get_fullpath(set_opt->insketchpath, idx_sketch96_suffix),
		&idx_file_size);
	const uint32_t n = (uint32_t)minco_stat_readin.infile_num;
	const size_t expected_idx_size = ((size_t)n + 1) * sizeof(fco_pos[0]);
	if (idx_file_size != expected_idx_size)
		errx(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
			 __func__, set_opt->insketchpath, idx_sketch96_suffix,
			 idx_file_size, expected_idx_size);

	size_t comb_file_size = 0;
	ctxobj96_t *mem_ctxobj = (ctxobj96_t *)read_from_file(
		test_get_fullpath(set_opt->insketchpath, combined_sketch96_suffix),
		&comb_file_size);
	if (comb_file_size % sizeof(mem_ctxobj[0]) != 0)
		errx(EINVAL, "%s(): %s/%s size %zu is not a multiple of %zu",
			 __func__, set_opt->insketchpath, combined_sketch96_suffix,
			 comb_file_size, sizeof(mem_ctxobj[0]));
	const size_t in_entry_ct = comb_file_size / sizeof(mem_ctxobj[0]);
	if (in_entry_ct > UINT32_MAX)
		errx(EXIT_FAILURE, "%s(): ctxobj96 pairwise markerdb currently supports up to %u entries",
			 __func__, UINT32_MAX);
	if (in_entry_ct != fco_pos[n])
		errx(EINVAL, "%s(): %s/%s size does not match %s",
			 __func__, set_opt->insketchpath, combined_sketch96_suffix,
			 idx_sketch96_suffix);

	markerdb_ctx_pair_edge_collect_t edge_collect = {0};
	pairwise_index_scan_stats_t scan_stats = {0};
	minco_collect_pairwise_context_edges96(set_opt, fco_pos, mem_ctxobj, n,
										   in_entry_ct, &edge_collect,
										   &scan_stats);

	uint8_t *drop_entry = calloc(in_entry_ct ? in_entry_ct : 1, sizeof(drop_entry[0]));
	if (!drop_entry)
		err(errno, "%s(): OOM pairwise marker bitmap", __func__);
	for (size_t e = 0; e < edge_collect.edges.n; ++e)
	{
		const pairwise_edge_pair_t edge = edge_collect.edges.edges[e];
		minco_mark_pair_shared_contexts96(drop_entry, mem_ctxobj, fco_pos,
										  (uint32_t)edge.a, (uint32_t)edge.b);
	}

	uint64_t *post_fco_pos = calloc((size_t)n + 1, sizeof(post_fco_pos[0]));
	if (!post_fco_pos)
		err(errno, "%s(): OOM output index", __func__);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(set_opt->p)
#endif
	for (uint32_t i = 0; i < n; i++)
	{
		uint64_t count = 0;
		for (uint64_t p = fco_pos[i]; p < fco_pos[i + 1]; ++p)
			if (!drop_entry[p])
				count++;
		post_fco_pos[i + 1] = count;
	}
	for (uint32_t i = 0; i < n; i++)
		post_fco_pos[i + 1] += post_fco_pos[i];
	const uint64_t out_entry_ct = post_fco_pos[n];

	char *out_comb_path = test_create_fullpath(set_opt->outdir,
											   combined_sketch96_suffix);
	int outfd = open(out_comb_path, O_CREAT | O_TRUNC | O_WRONLY, 0644);
	if (outfd < 0)
		err(errno, "%s(): cannot create %s/%s", __func__,
			set_opt->outdir, combined_sketch96_suffix);
	const uint64_t out_bytes_u64 = out_entry_ct * (uint64_t)sizeof(mem_ctxobj[0]);
	if (out_entry_ct != 0 && out_bytes_u64 / out_entry_ct != sizeof(mem_ctxobj[0]))
		errx(EINVAL, "%s(): markerdb output byte size overflow", __func__);
	if (ftruncate(outfd, (off_t)out_bytes_u64) != 0)
		err(errno, "%s(): resize %s", __func__, out_comb_path);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(set_opt->p)
#endif
	for (uint32_t i = 0; i < n; i++)
	{
		const uint64_t sample_count = post_fco_pos[i + 1] - post_fco_pos[i];
		if (sample_count == 0)
			continue;
		ctxobj96_t *sample_markers =
			malloc(sample_count * sizeof(sample_markers[0]));
		if (!sample_markers)
			err(errno, "%s(): OOM sample context marker buffer", __func__);
		uint64_t out = 0;
		for (uint64_t p = fco_pos[i]; p < fco_pos[i + 1]; ++p)
			if (!drop_entry[p])
				sample_markers[out++] = mem_ctxobj[p];
		if (out != sample_count)
			err(EINVAL, "%s(): pairwise context marker count changed for sample %u",
				__func__, i);
		pwrite_all(
			outfd,
			sample_markers,
			sample_count * sizeof(sample_markers[0]),
			(off_t)(post_fco_pos[i] * sizeof(sample_markers[0])),
			out_comb_path);
		free(sample_markers);
	}
	if (close(outfd) != 0)
		err(errno, "%s(): close %s", __func__, out_comb_path);

	write_to_file(format_string("%s/%s", set_opt->outdir, idx_sketch96_suffix),
				  post_fco_pos, ((size_t)n + 1) * sizeof(post_fco_pos[0]));
	write_to_file(test_create_fullpath(set_opt->outdir, sketch_stat),
				  mem_stat, stat_file_size);
	copy_minco_sketch_annotations(set_opt->insketchpath, set_opt->outdir,
								  minco_stat_readin.infile_num);
	copy_minco_domain_profiles(set_opt->insketchpath, set_opt->outdir,
							   minco_stat_readin.infile_num);

	fprintf(stderr,
			"minco set: pairwise coden15 context markerdb kept %" PRIu64
			" ctxobj entries from %zu input entries; min_xny=%u min_af=%.12g "
			"accepted_pairs=%zu candidate_pairs=%" PRIu64
			" exact_pairs=%" PRIu64 " index_max_ctx_freq=%u index_min_votes=%u\n",
			out_entry_ct, in_entry_ct, set_opt->markerdb_ctx_min_xny,
			set_opt->markerdb_ctx_min_af, edge_collect.edges.n,
			scan_stats.candidate_pairs, scan_stats.exact_pairs,
			set_opt->markerdb_ctx_index_max_ctx_freq,
			set_opt->markerdb_ctx_index_min_votes > 0
				? set_opt->markerdb_ctx_index_min_votes
				: (set_opt->markerdb_ctx_min_xny > 0
					   ? set_opt->markerdb_ctx_min_xny
					   : 1));
	minco_warn_markerdb_small_refs(set_opt, post_fco_pos, mem_stat, stat_file_size,
								   "pairwise coden15 context-markerdb");

	free(out_comb_path);
	free(drop_entry);
	pairwise_edge_list_free(&edge_collect.edges);
	free_all(fco_pos, post_fco_pos, NULL);
	free_read_from_file(mem_ctxobj, comb_file_size);
}

static void minco_write_pairwise_context_markerdb(set_opt_t *set_opt,
												  const void *mem_stat,
												  size_t stat_file_size)
{
	const_comask_init(&minco_stat_readin);
	if (operate_stat_uses_ctxobj96_payload(&minco_stat_readin))
	{
		minco_write_pairwise_context_markerdb96(set_opt, mem_stat, stat_file_size);
		return;
	}
	if (Bitslen.ctx == 0 || Bitslen.ctx > 63)
		errx(EXIT_FAILURE, "%s(): unsupported context bit width: %u",
			 __func__, (unsigned)Bitslen.ctx);
	if (minco_stat_readin.infile_num >= (1 << GID_NBITS))
		errx(EXIT_FAILURE, "%s(): genome number %u exceeds maximum %u",
			 __func__, minco_stat_readin.infile_num, 1 << GID_NBITS);
	if (Bitslen.ctx + GID_NBITS > 64)
		errx(EXIT_FAILURE, "%s(): context_bits_len(%u)+gid_bits_len(%d) exceed 64",
			 __func__, (unsigned)Bitslen.ctx, GID_NBITS);
	if (minco_stat_readin.conflict)
		errx(EXIT_FAILURE, "%s(): --markerdb-ctx-pairwise currently requires non-conflict reference sketches",
			 __func__);

	size_t idx_file_size = 0;
	uint64_t *fco_pos = (uint64_t *)read_from_file(
		test_get_fullpath(set_opt->insketchpath, idx_sketch_suffix),
		&idx_file_size);
	const size_t expected_idx_size =
		((size_t)minco_stat_readin.infile_num + 1) * sizeof(fco_pos[0]);
	if (idx_file_size != expected_idx_size)
		errx(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
			 __func__, set_opt->insketchpath, idx_sketch_suffix,
			 idx_file_size, expected_idx_size);

	size_t comb_file_size = 0;
	uint64_t *mem_ctxobj = (uint64_t *)read_from_file(
		test_get_fullpath(set_opt->insketchpath, combined_sketch_suffix),
		&comb_file_size);
	const size_t in_entry_ct = comb_file_size / sizeof(mem_ctxobj[0]);
	if (comb_file_size != fco_pos[minco_stat_readin.infile_num] * sizeof(mem_ctxobj[0]))
		errx(EINVAL, "%s(): %s/%s size does not match %s",
			 __func__, set_opt->insketchpath, combined_sketch_suffix,
			 idx_sketch_suffix);
	if (in_entry_ct > UINT32_MAX)
		errx(EXIT_FAILURE, "%s(): pairwise context markerdb currently supports up to %u entries",
			 __func__, UINT32_MAX);

	minco_sketch_info_t sketch_info = {0};
	if (!minco_stat_decode_mem(mem_stat, stat_file_size, &minco_stat_readin, &sketch_info))
		errx(EINVAL, "%s(): malformed minco stat buffer", __func__);
	const uint32_t n = (uint32_t)minco_stat_readin.infile_num;
	unify_sketch_t sketch = {
		.stat_type = 2,
		.mem_stat = (void *)mem_stat,
		.gname = minco_stat_names_from_mem((void *)mem_stat, stat_file_size),
		.comb_sketch = mem_ctxobj,
		.sketch_index = fco_pos,
		.infile_num = minco_stat_readin.infile_num,
		.kmerlen = minco_stat_readin.klen,
		.hash_id = sketch_info.has_minco_ext ? sketch_info.sketch_id
											  : minco_stat_readin.hash_id,
		.conflict = minco_stat_readin.conflict,
		.minco_info = sketch_info,
	};
	sketch.stats.minco_stat = minco_stat_readin;

	bool temp_index = false;
	char *index_path = NULL;
	if (file_exists_in_folder(set_opt->insketchpath, (char *)sorted_comb_ctxgid64obj32))
	{
		index_path = test_get_fullpath(set_opt->insketchpath, sorted_comb_ctxgid64obj32);
	}
	else
	{
		fprintf(stderr,
				"minco set: --markerdb-ctx-pairwise building temporary sorted context index\n");
		index_path = minco_build_temp_ctxgid_index(
			fco_pos, mem_ctxobj, n, in_entry_ct);
		temp_index = true;
	}

	markerdb_ctx_pair_edge_collect_t edge_collect = {0};
	pairwise_metric_expr_t metric = pairwise_metric_expr_single(PAIRWISE_METRIC_AAF);
	pairwise_index_scan_options_t scan_opt = {
		.metric = &metric,
		.cut = 1e100,
		.ctxcut = set_opt->markerdb_ctx_min_xny,
		.max_afcut = set_opt->markerdb_ctx_min_af,
		.index_max_ctx_freq = set_opt->markerdb_ctx_index_max_ctx_freq,
		.index_min_votes = set_opt->markerdb_ctx_index_min_votes > 0
							   ? set_opt->markerdb_ctx_index_min_votes
							   : (set_opt->markerdb_ctx_min_xny > 0
									  ? set_opt->markerdb_ctx_min_xny
									  : 1),
		.index_sample_step = set_opt->markerdb_ctx_index_sample_step,
		.threads = set_opt->p,
	};
	pairwise_index_scan_stats_t scan_stats = {0};
	pairwise_indexed_self_scan(&sketch, index_path, &scan_opt,
							   markerdb_ctx_pair_edge_cb, &edge_collect,
							   NULL, NULL, &scan_stats);

	if (temp_index && unlink(index_path) != 0)
		fprintf(stderr, "minco set: warning: failed to remove temporary index %s: %s\n",
				index_path, strerror(errno));
	free(index_path);

	uint8_t *drop_entry = calloc(in_entry_ct ? in_entry_ct : 1, sizeof(drop_entry[0]));
	if (!drop_entry)
		err(errno, "%s(): OOM pairwise marker bitmap", __func__);

	for (size_t e = 0; e < edge_collect.edges.n; ++e)
	{
		const pairwise_edge_pair_t edge = edge_collect.edges.edges[e];
		minco_mark_pair_shared_contexts(drop_entry, mem_ctxobj, fco_pos,
										(uint32_t)edge.a, (uint32_t)edge.b);
	}

	uint64_t *post_fco_pos = calloc((size_t)n + 1, sizeof(post_fco_pos[0]));
	if (!post_fco_pos)
		err(errno, "%s(): OOM output index", __func__);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(set_opt->p)
#endif
	for (uint32_t i = 0; i < n; i++)
	{
		uint64_t count = 0;
		for (uint64_t p = fco_pos[i]; p < fco_pos[i + 1]; ++p)
			if (!drop_entry[p])
				count++;
		post_fco_pos[i + 1] = count;
	}
	for (uint32_t i = 0; i < n; i++)
		post_fco_pos[i + 1] += post_fco_pos[i];
	const uint64_t out_entry_ct = post_fco_pos[n];

	char *out_comb_path = test_create_fullpath(set_opt->outdir, combined_sketch_suffix);
	int outfd = open(out_comb_path, O_CREAT | O_TRUNC | O_WRONLY, 0644);
	if (outfd < 0)
		err(errno, "%s(): cannot create %s/%s", __func__,
			set_opt->outdir, combined_sketch_suffix);
	const uint64_t out_bytes_u64 = out_entry_ct * (uint64_t)sizeof(mem_ctxobj[0]);
	if (out_entry_ct != 0 && out_bytes_u64 / out_entry_ct != sizeof(mem_ctxobj[0]))
		errx(EINVAL, "%s(): markerdb output byte size overflow", __func__);
	if (ftruncate(outfd, (off_t)out_bytes_u64) != 0)
		err(errno, "%s(): resize %s", __func__, out_comb_path);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(set_opt->p)
#endif
	for (uint32_t i = 0; i < n; i++)
	{
		const uint64_t sample_count = post_fco_pos[i + 1] - post_fco_pos[i];
		if (sample_count == 0)
			continue;
		uint64_t *sample_markers = malloc(sample_count * sizeof(sample_markers[0]));
		if (!sample_markers)
			err(errno, "%s(): OOM sample context marker buffer", __func__);
		uint64_t out = 0;
		for (uint64_t p = fco_pos[i]; p < fco_pos[i + 1]; ++p)
			if (!drop_entry[p])
				sample_markers[out++] = mem_ctxobj[p];
		if (out != sample_count)
			err(EINVAL, "%s(): pairwise context marker count changed for sample %u",
				__func__, i);
		pwrite_all(
			outfd,
			sample_markers,
			sample_count * sizeof(sample_markers[0]),
			(off_t)(post_fco_pos[i] * sizeof(sample_markers[0])),
			out_comb_path);
		free(sample_markers);
	}
	if (close(outfd) != 0)
		err(errno, "%s(): close %s", __func__, out_comb_path);

	write_to_file(format_string("%s/%s", set_opt->outdir, idx_sketch_suffix),
				  post_fco_pos, ((size_t)n + 1) * sizeof(post_fco_pos[0]));
	write_to_file(test_create_fullpath(set_opt->outdir, sketch_stat),
				  mem_stat, stat_file_size);
	copy_minco_sketch_annotations(set_opt->insketchpath, set_opt->outdir,
								  minco_stat_readin.infile_num);
	copy_minco_domain_profiles(set_opt->insketchpath, set_opt->outdir,
							   minco_stat_readin.infile_num);

	fprintf(stderr,
			"minco set: pairwise context markerdb kept %" PRIu64
			" ctxobj entries from %zu input entries; min_xny=%u min_af=%.12g "
			"accepted_pairs=%zu candidate_pairs=%" PRIu64
			" exact_pairs=%" PRIu64 " index_max_ctx_freq=%u index_min_votes=%u\n",
			out_entry_ct, in_entry_ct, set_opt->markerdb_ctx_min_xny,
			set_opt->markerdb_ctx_min_af, edge_collect.edges.n,
			scan_stats.candidate_pairs, scan_stats.exact_pairs,
			set_opt->markerdb_ctx_index_max_ctx_freq,
			scan_opt.index_min_votes);
	minco_warn_markerdb_small_refs(set_opt, post_fco_pos, mem_stat, stat_file_size,
								   "pairwise context-markerdb");

	free(out_comb_path);
	free(drop_entry);
	pairwise_edge_list_free(&edge_collect.edges);
	free_all(fco_pos, post_fco_pos, NULL);
	free_read_from_file(mem_ctxobj, comb_file_size);
}

static void minco_warn_markerdb_small_refs(const set_opt_t *set_opt,
										   const uint64_t *post_fco_pos,
										   const void *mem_stat,
										   size_t stat_size,
										   const char *mode_label)
{
	if (!set_opt || set_opt->markerdb_warn_threshold == 0 || !post_fco_pos)
		return;

	minco_sketch_stat_t stat = {0};
	if (!minco_stat_decode_mem(mem_stat, stat_size, &stat, NULL))
		errx(EINVAL, "%s(): malformed minco stat buffer", __func__);
	if (stat.infile_num <= 0)
		return;

	const uint64_t threshold = set_opt->markerdb_warn_threshold;
	uint32_t below_ct = 0;
	for (int i = 0; i < stat.infile_num; ++i)
	{
		const uint64_t retained = post_fco_pos[i + 1] - post_fco_pos[i];
		if (retained < threshold)
			below_ct++;
	}
	if (below_ct == 0)
		return;

	const char (*names)[PATHLEN] = minco_stat_const_names_from_mem(mem_stat, stat_size);
	fprintf(stderr,
			"minco set: WARNING: %u/%d refs have markerdb sketch size below %"
			PRIu64 " after %s\n",
			below_ct, stat.infile_num, threshold,
			mode_label ? mode_label : "markerdb");
	fprintf(stderr, "minco set: low-marker-ref\tsample_id\tsketch_entries\tsample\n");
	for (int i = 0; i < stat.infile_num; ++i)
	{
		const uint64_t retained = post_fco_pos[i + 1] - post_fco_pos[i];
		if (retained < threshold)
			fprintf(stderr, "minco set: low-marker-ref\t%d\t%" PRIu64 "\t%s\n",
					i, retained, names ? names[i] : "");
	}
}

static uint64_t *minco_unique_ref_contexts_from_ctxobj(
	const uint64_t *fco_pos,
	const uint64_t *mem_ctxobj,
	uint32_t infile_num,
	size_t in_entry_ct,
	size_t *unique_ctx_ct)
{
	if (in_entry_ct > UINT32_MAX)
		errx(EXIT_FAILURE, "%s(): sketch has %zu entries; context markerdb currently supports up to %u",
			 __func__, in_entry_ct, UINT32_MAX);
	ctxgidobj_t *ctxgid = ctxobj64_2ctxgidobj(
		(uint64_t *)fco_pos, (uint64_t *)mem_ctxobj, (int)infile_num,
		(uint32_t)in_entry_ct);
	ctxgidobj_sort_array(ctxgid, in_entry_ct);

	uint64_t *unique_ctx = malloc(in_entry_ct * sizeof(unique_ctx[0]));
	if (!unique_ctx && in_entry_ct)
		err(errno, "%s(): OOM unique context list", __func__);

	const uint64_t gidmask_local = (1ULL << GID_NBITS) - 1ULL;
	size_t out = 0;
	for (size_t i = 0; i < in_entry_ct; )
	{
		const uint64_t ctx = ctxgid[i].ctxgid >> GID_NBITS;
		uint32_t prev_gid = UINT32_MAX;
		uint32_t distinct_gids = 0;
		while (i < in_entry_ct && (ctxgid[i].ctxgid >> GID_NBITS) == ctx)
		{
			const uint32_t gid = (uint32_t)(ctxgid[i].ctxgid & gidmask_local);
			if (gid != prev_gid)
			{
				distinct_gids++;
				prev_gid = gid;
			}
			++i;
		}
		if (distinct_gids == 1)
			unique_ctx[out++] = ctx;
	}

	free(ctxgid);
	uint64_t *shrunk = realloc(unique_ctx, out * sizeof(unique_ctx[0]));
	if (shrunk || out == 0)
		unique_ctx = shrunk;
	*unique_ctx_ct = out;
	return unique_ctx;
}

static uint64_t *minco_unique_ref_contexts_from_ctxobj96(
	const uint64_t *fco_pos,
	const ctxobj96_t *mem_ctxobj,
	uint32_t infile_num,
	size_t in_entry_ct,
	size_t *unique_ctx_ct)
{
	if (in_entry_ct > UINT32_MAX)
		errx(EXIT_FAILURE, "%s(): sketch has %zu entries; ctxobj96 context markerdb currently supports up to %u",
			 __func__, in_entry_ct, UINT32_MAX);
	ctxgidobj128_t *ctxgid = ctxobj96_2ctxgidobj128(
		(uint64_t *)fco_pos, (ctxobj96_t *)mem_ctxobj, (int)infile_num,
		(uint32_t)in_entry_ct);
	ctxgidobj128_sort_array(ctxgid, in_entry_ct);

	uint64_t *unique_ctx = malloc(in_entry_ct * sizeof(unique_ctx[0]));
	if (!unique_ctx && in_entry_ct)
		err(errno, "%s(): OOM unique context list", __func__);

	size_t out = 0;
	for (size_t i = 0; i < in_entry_ct;)
	{
		const uint64_t ctx = ctxgid[i].ctx;
		uint32_t prev_gid = UINT32_MAX;
		uint32_t distinct_gids = 0;
		while (i < in_entry_ct && ctxgid[i].ctx == ctx)
		{
			const uint32_t gid = ctxgid[i].gid;
			if (gid != prev_gid)
			{
				distinct_gids++;
				prev_gid = gid;
			}
			++i;
		}
		if (distinct_gids == 1)
			unique_ctx[out++] = ctx;
	}

	free(ctxgid);
	uint64_t *shrunk = realloc(unique_ctx, out * sizeof(unique_ctx[0]));
	if (shrunk || out == 0)
		unique_ctx = shrunk;
	*unique_ctx_ct = out;
	return unique_ctx;
}

static void minco_write_context_markerdb96(set_opt_t *set_opt,
										   const void *mem_stat,
										   size_t stat_file_size)
{
	const_comask_init(&minco_stat_readin);
	if (Bitslen.ctx == 0 || Bitslen.ctx > 64)
		errx(EXIT_FAILURE, "%s(): unsupported context bit width: %u",
			 __func__, (unsigned)Bitslen.ctx);

	size_t idx_file_size = 0;
	uint64_t *fco_pos = (uint64_t *)read_from_file(
		test_get_fullpath(set_opt->insketchpath, idx_sketch96_suffix),
		&idx_file_size);
	const size_t expected_idx_size =
		((size_t)minco_stat_readin.infile_num + 1) * sizeof(fco_pos[0]);
	if (idx_file_size != expected_idx_size)
		errx(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
			 __func__, set_opt->insketchpath, idx_sketch96_suffix,
			 idx_file_size, expected_idx_size);

	size_t comb_file_size = 0;
	ctxobj96_t *mem_ctxobj = (ctxobj96_t *)read_from_file(
		test_get_fullpath(set_opt->insketchpath, combined_sketch96_suffix),
		&comb_file_size);
	if (comb_file_size % sizeof(mem_ctxobj[0]) != 0)
		errx(EINVAL, "%s(): %s/%s size %zu is not a multiple of %zu",
			 __func__, set_opt->insketchpath, combined_sketch96_suffix,
			 comb_file_size, sizeof(mem_ctxobj[0]));
	const size_t in_entry_ct = comb_file_size / sizeof(mem_ctxobj[0]);
	if (in_entry_ct != fco_pos[minco_stat_readin.infile_num])
		errx(EINVAL, "%s(): %s/%s size does not match %s",
			 __func__, set_opt->insketchpath, combined_sketch96_suffix,
			 idx_sketch96_suffix);

	size_t unique_ctx_ct = 0;
	uint64_t *unique_ctx = minco_unique_ref_contexts_from_ctxobj96(
		fco_pos, mem_ctxobj, minco_stat_readin.infile_num,
		in_entry_ct, &unique_ctx_ct);

	const unsigned bucket_bits = Bitslen.ctx < 24 ? Bitslen.ctx : 24;
	size_t *bucket_offsets = build_ctx_bucket_offsets(
		unique_ctx, unique_ctx_ct, bucket_bits, Bitslen.ctx);

	uint64_t *post_fco_pos = calloc((size_t)minco_stat_readin.infile_num + 1,
								   sizeof(post_fco_pos[0]));
	if (!post_fco_pos)
		err(errno, "%s(): OOM output index", __func__);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(set_opt->p)
#endif
	for (uint32_t i = 0; i < minco_stat_readin.infile_num; i++)
	{
		uint64_t count = 0;
		for (uint64_t n = fco_pos[i]; n < fco_pos[i + 1]; n++)
		{
			const uint64_t ctx = mem_ctxobj[n].ctx;
			if (ctx_bucketed_contains(unique_ctx, bucket_offsets,
									  bucket_bits, Bitslen.ctx, ctx))
				count++;
		}
		post_fco_pos[i + 1] = count;
	}
	for (uint32_t i = 0; i < minco_stat_readin.infile_num; i++)
		post_fco_pos[i + 1] += post_fco_pos[i];
	const uint64_t out_entry_ct = post_fco_pos[minco_stat_readin.infile_num];

	char *out_comb_path = test_create_fullpath(set_opt->outdir,
											   combined_sketch96_suffix);
	int outfd = open(out_comb_path, O_CREAT | O_TRUNC | O_WRONLY, 0644);
	if (outfd < 0)
		err(errno, "%s(): cannot create %s/%s", __func__,
			set_opt->outdir, combined_sketch96_suffix);
	const uint64_t out_bytes_u64 = out_entry_ct * (uint64_t)sizeof(mem_ctxobj[0]);
	if (out_entry_ct != 0 && out_bytes_u64 / out_entry_ct != sizeof(mem_ctxobj[0]))
		errx(EINVAL, "%s(): markerdb output byte size overflow", __func__);
	if (ftruncate(outfd, (off_t)out_bytes_u64) != 0)
		err(errno, "%s(): resize %s", __func__, out_comb_path);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(set_opt->p)
#endif
	for (uint32_t i = 0; i < minco_stat_readin.infile_num; i++)
	{
		const uint64_t sample_count = post_fco_pos[i + 1] - post_fco_pos[i];
		if (sample_count == 0)
			continue;
		ctxobj96_t *sample_markers =
			malloc(sample_count * sizeof(sample_markers[0]));
		if (!sample_markers)
			err(errno, "%s(): OOM sample context marker buffer", __func__);
		uint64_t out = 0;
		for (uint64_t n = fco_pos[i]; n < fco_pos[i + 1]; n++)
		{
			const uint64_t ctx = mem_ctxobj[n].ctx;
			if (ctx_bucketed_contains(unique_ctx, bucket_offsets,
									  bucket_bits, Bitslen.ctx, ctx))
				sample_markers[out++] = mem_ctxobj[n];
		}
		if (out != sample_count)
			err(EINVAL, "%s(): context marker count changed for sample %u",
				__func__, i);
		pwrite_all(
			outfd,
			sample_markers,
			sample_count * sizeof(sample_markers[0]),
			(off_t)(post_fco_pos[i] * sizeof(sample_markers[0])),
			out_comb_path);
		free(sample_markers);
	}
	if (close(outfd) != 0)
		err(errno, "%s(): close %s", __func__, out_comb_path);

	write_to_file(format_string("%s/%s", set_opt->outdir, idx_sketch96_suffix),
				  post_fco_pos,
				  ((size_t)minco_stat_readin.infile_num + 1) * sizeof(post_fco_pos[0]));
	write_to_file(test_create_fullpath(set_opt->outdir, sketch_stat),
				  mem_stat, stat_file_size);
	copy_minco_sketch_annotations(set_opt->insketchpath, set_opt->outdir,
								  minco_stat_readin.infile_num);
	copy_minco_domain_profiles(set_opt->insketchpath, set_opt->outdir,
							   minco_stat_readin.infile_num);

	fprintf(stderr,
			"minco set: coden15 context markerdb kept %zu ref-specific contexts, %" PRIu64
			" ctxobj entries from %zu input entries\n",
			unique_ctx_ct, out_entry_ct, in_entry_ct);
	minco_warn_markerdb_small_refs(set_opt, post_fco_pos, mem_stat, stat_file_size,
								   "coden15 context-markerdb");

	free(out_comb_path);
	free(bucket_offsets);
	free(unique_ctx);
	free_all(fco_pos, post_fco_pos, NULL);
	free_read_from_file(mem_ctxobj, comb_file_size);
}

static void minco_write_context_markerdb(set_opt_t *set_opt, const void *mem_stat,
										 size_t stat_file_size)
{
	if (set_opt->q2markerdb_context_pairwise)
	{
		minco_write_pairwise_context_markerdb(set_opt, mem_stat, stat_file_size);
		return;
	}

	const_comask_init(&minco_stat_readin);
	if (operate_stat_uses_ctxobj96_payload(&minco_stat_readin))
	{
		minco_write_context_markerdb96(set_opt, mem_stat, stat_file_size);
		return;
	}
	if (Bitslen.ctx == 0 || Bitslen.ctx > 63)
		errx(EXIT_FAILURE, "%s(): unsupported context bit width: %u",
			 __func__, (unsigned)Bitslen.ctx);
	if (minco_stat_readin.infile_num >= (1 << GID_NBITS))
		errx(EXIT_FAILURE, "%s(): genome number %u exceeds maximum %u",
			 __func__, minco_stat_readin.infile_num, 1 << GID_NBITS);
	if (Bitslen.ctx + GID_NBITS > 64)
		errx(EXIT_FAILURE, "%s(): context_bits_len(%u)+gid_bits_len(%d) exceed 64",
			 __func__, (unsigned)Bitslen.ctx, GID_NBITS);

	size_t idx_file_size = 0;
	uint64_t *fco_pos = (uint64_t *)read_from_file(
		test_get_fullpath(set_opt->insketchpath, idx_sketch_suffix),
		&idx_file_size);
	const size_t expected_idx_size =
		((size_t)minco_stat_readin.infile_num + 1) * sizeof(fco_pos[0]);
	if (idx_file_size != expected_idx_size)
		errx(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
			 __func__, set_opt->insketchpath, idx_sketch_suffix,
			 idx_file_size, expected_idx_size);

	size_t comb_file_size = 0;
	uint64_t *mem_ctxobj = (uint64_t *)read_from_file(
		test_get_fullpath(set_opt->insketchpath, combined_sketch_suffix),
		&comb_file_size);
	const size_t in_entry_ct = comb_file_size / sizeof(mem_ctxobj[0]);
	if (comb_file_size != fco_pos[minco_stat_readin.infile_num] * sizeof(mem_ctxobj[0]))
		errx(EINVAL, "%s(): %s/%s size does not match %s",
			 __func__, set_opt->insketchpath, combined_sketch_suffix,
			 idx_sketch_suffix);

	size_t unique_ctx_ct = 0;
	uint64_t *unique_ctx = minco_unique_ref_contexts_from_ctxobj(
		fco_pos, mem_ctxobj, minco_stat_readin.infile_num,
		in_entry_ct, &unique_ctx_ct);
	free_read_from_file(mem_ctxobj, comb_file_size);
	mem_ctxobj = NULL;

	const unsigned bucket_bits = Bitslen.ctx < 24 ? Bitslen.ctx : 24;
	size_t *bucket_offsets = build_ctx_bucket_offsets(
		unique_ctx, unique_ctx_ct, bucket_bits, Bitslen.ctx);

	mem_ctxobj = (uint64_t *)read_from_file(
		test_get_fullpath(set_opt->insketchpath, combined_sketch_suffix),
		&comb_file_size);
	uint64_t *post_fco_pos = calloc((size_t)minco_stat_readin.infile_num + 1,
								   sizeof(post_fco_pos[0]));
	if (!post_fco_pos)
		err(errno, "%s(): OOM output index", __func__);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(set_opt->p)
#endif
	for (uint32_t i = 0; i < minco_stat_readin.infile_num; i++)
	{
		uint64_t count = 0;
		for (uint64_t n = fco_pos[i]; n < fco_pos[i + 1]; n++)
		{
			const uint64_t ctx = mem_ctxobj[n] >> Bitslen.obj;
			if (ctx_bucketed_contains(unique_ctx, bucket_offsets,
									  bucket_bits, Bitslen.ctx, ctx))
				count++;
		}
		post_fco_pos[i + 1] = count;
	}

	for (uint32_t i = 0; i < minco_stat_readin.infile_num; i++)
		post_fco_pos[i + 1] += post_fco_pos[i];
	const uint64_t out_entry_ct = post_fco_pos[minco_stat_readin.infile_num];

	char *out_comb_path = test_create_fullpath(set_opt->outdir, combined_sketch_suffix);
	int outfd = open(out_comb_path, O_CREAT | O_TRUNC | O_WRONLY, 0644);
	if (outfd < 0)
		err(errno, "%s(): cannot create %s/%s", __func__,
			set_opt->outdir, combined_sketch_suffix);
	const uint64_t out_bytes_u64 = out_entry_ct * (uint64_t)sizeof(mem_ctxobj[0]);
	if (out_entry_ct != 0 && out_bytes_u64 / out_entry_ct != sizeof(mem_ctxobj[0]))
		errx(EINVAL, "%s(): markerdb output byte size overflow", __func__);
	if (ftruncate(outfd, (off_t)out_bytes_u64) != 0)
		err(errno, "%s(): resize %s", __func__, out_comb_path);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(set_opt->p)
#endif
	for (uint32_t i = 0; i < minco_stat_readin.infile_num; i++)
	{
		const uint64_t sample_count = post_fco_pos[i + 1] - post_fco_pos[i];
		if (sample_count == 0)
			continue;
		uint64_t *sample_markers = malloc(sample_count * sizeof(sample_markers[0]));
		if (!sample_markers)
			err(errno, "%s(): OOM sample context marker buffer", __func__);
		uint64_t out = 0;
		for (uint64_t n = fco_pos[i]; n < fco_pos[i + 1]; n++)
		{
			const uint64_t ctx = mem_ctxobj[n] >> Bitslen.obj;
			if (ctx_bucketed_contains(unique_ctx, bucket_offsets,
									  bucket_bits, Bitslen.ctx, ctx))
				sample_markers[out++] = mem_ctxobj[n];
		}
		if (out != sample_count)
			err(EINVAL, "%s(): context marker count changed for sample %u",
				__func__, i);
		pwrite_all(
			outfd,
			sample_markers,
			sample_count * sizeof(sample_markers[0]),
			(off_t)(post_fco_pos[i] * sizeof(sample_markers[0])),
			out_comb_path);
		free(sample_markers);
	}
	if (close(outfd) != 0)
		err(errno, "%s(): close %s", __func__, out_comb_path);

	write_to_file(format_string("%s/%s", set_opt->outdir, idx_sketch_suffix),
				  post_fco_pos,
				  ((size_t)minco_stat_readin.infile_num + 1) * sizeof(post_fco_pos[0]));
	write_to_file(test_create_fullpath(set_opt->outdir, sketch_stat),
				  mem_stat, stat_file_size);
	copy_minco_sketch_annotations(set_opt->insketchpath, set_opt->outdir,
								  minco_stat_readin.infile_num);
	copy_minco_domain_profiles(set_opt->insketchpath, set_opt->outdir,
							   minco_stat_readin.infile_num);

	fprintf(stderr,
			"minco set: context markerdb kept %zu ref-specific contexts, %" PRIu64
			" ctxobj entries from %zu input entries\n",
			unique_ctx_ct, out_entry_ct, in_entry_ct);
	minco_warn_markerdb_small_refs(set_opt, post_fco_pos, mem_stat, stat_file_size,
								   "context-markerdb");

	free(out_comb_path);
	free(bucket_offsets);
	free(unique_ctx);
	free_all(fco_pos, post_fco_pos, NULL);
	free_read_from_file(mem_ctxobj, comb_file_size);
}

static uint64_t minco_downsample_unique_ctx_count(const uint64_t *values, size_t n,
												  uint32_t n_obj_bits)
{
	if (!values || n == 0)
		return 0;
	uint64_t count = 1;
	uint64_t prev = n_obj_bits == 64 ? 0 : (values[0] >> n_obj_bits);
	for (size_t i = 1; i < n; ++i)
	{
		const uint64_t ctx = n_obj_bits == 64 ? 0 : (values[i] >> n_obj_bits);
		if (ctx != prev)
		{
			++count;
			prev = ctx;
		}
	}
	return count;
}

static uint64_t minco_downsample_density_estimate(uint64_t observed,
												  uint64_t threshold,
												  uint32_t hash_bits)
{
	if (observed == 0 || hash_bits == 0)
		return observed;
	const long double hash_space = ldexpl(1.0L, (int)hash_bits);
	const long double denom = (long double)threshold + 1.0L;
	if (denom <= 0.0L)
		return 0;
	const long double estimate =
		((long double)observed * hash_space / denom) + 0.5L;
	if (estimate >= (long double)UINT64_MAX)
		return UINT64_MAX;
	return (uint64_t)estimate;
}

static void minco_downsample_copy_optional_sidecar(const char *indir,
												   const char *outdir,
												   const char *name,
												   size_t expected_size)
{
	if (!file_exists_in_folder(indir, name))
		return;
	char *in_path = test_get_fullpath(indir, name);
	size_t sidecar_size = 0;
	void *data = read_from_file(in_path, &sidecar_size);
	free(in_path);
	if (sidecar_size != expected_size)
		errx(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
			 __func__, indir, name, sidecar_size, expected_size);
	char *out_path = test_create_fullpath(outdir, name);
	write_to_file(out_path, data, sidecar_size);
	free(out_path);
	free_read_from_file(data, sidecar_size);
}

static void minco_downsample_prepare_outdir(const char *outdir)
{
	if (!outdir || outdir[0] == '\0')
		errx(EINVAL, "%s(): missing output directory", __func__);
	if (mkdir(outdir, 0777) != 0 && errno != EEXIST)
		err(errno, "%s(): cannot create %s", __func__, outdir);
	const char *core_names[] = {
		sketch_stat,
		combined_sketch_suffix,
		combined_sketch96_suffix,
		idx_sketch_suffix,
		idx_sketch96_suffix,
		minco_ctxmeta_bin_stat,
		sketch_domain_stat,
		combined_ab_suffix,
		sketch_position_suffix,
		sorted_comb_ctxgid64obj32,
		sorted_comb_ctx64gid32obj32,
	};
	for (size_t i = 0; i < sizeof(core_names) / sizeof(core_names[0]); ++i)
		if (file_exists_in_folder(outdir, core_names[i]))
			errx(EINVAL, "%s already contains %s; choose an empty output directory",
				 outdir, core_names[i]);
}

static size_t minco_downsample_file_size(const char *dir, const char *name)
{
	char *path = test_get_fullpath(dir, name);
	struct stat st;
	if (stat(path, &st) != 0)
		err(errno, "%s(): stat %s", __func__, path);
	free(path);
	if (st.st_size < 0)
		errx(EINVAL, "%s(): negative file size for %s/%s", __func__, dir, name);
	return (size_t)st.st_size;
}

static FILE *minco_downsample_open_optional_entry_file(const char *dir,
													   const char *name,
													   uint64_t total_entries,
													   size_t elem_size)
{
	if (!file_exists_in_folder(dir, name))
		return NULL;
	const size_t expected_size = (size_t)total_entries * elem_size;
	const size_t observed_size = minco_downsample_file_size(dir, name);
	if (observed_size != expected_size)
		errx(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
			 __func__, dir, name, observed_size, expected_size);
	char *path = test_get_fullpath(dir, name);
	FILE *fp = fopen(path, "rb");
	if (!fp)
		err(errno, "%s(): fopen %s", __func__, path);
	free(path);
	return fp;
}

static void minco_downsample_seek(FILE *fp, uint64_t entry_offset,
								  size_t elem_size, const char *label)
{
	if (!fp)
		return;
	const uint64_t byte_offset_u64 = entry_offset * (uint64_t)elem_size;
	if (entry_offset != 0 && byte_offset_u64 / entry_offset != elem_size)
		errx(EINVAL, "%s(): byte offset overflow for %s", __func__, label);
	if (fseeko(fp, (off_t)byte_offset_u64, SEEK_SET) != 0)
		err(errno, "%s(): seek %s", __func__, label);
}

static void minco_downsample_read_exact(FILE *fp, void *buf, size_t elem_size,
										size_t count, const char *label)
{
	if (count == 0)
		return;
	if (fread(buf, elem_size, count, fp) != count)
		err(errno, "%s(): short read from %s", __func__, label);
}

static void minco_downsample_write_exact(FILE *fp, const void *buf, size_t elem_size,
										 size_t count, const char *label)
{
	if (count == 0)
		return;
	if (fwrite(buf, elem_size, count, fp) != count)
		err(errno, "%s(): write %s", __func__, label);
}

static int minco_sketch_downsample96(set_opt_t *set_opt,
									 void *mem_stat,
									 size_t stat_size,
									 const minco_sketch_info_t *info,
									 const minco_sketch_stat_t *origin,
									 char (*names)[PATHLEN])
{
	const int infile_num = origin->infile_num;
	char *in_idx_path = test_get_fullpath(set_opt->insketchpath,
										  idx_sketch96_suffix);
	uint64_t *in_idx = read_from_file(in_idx_path, &file_size);
	free(in_idx_path);
	const size_t expected_idx_size = (size_t)(infile_num + 1) * sizeof(in_idx[0]);
	if (file_size != expected_idx_size)
		errx(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
			 __func__, set_opt->insketchpath, idx_sketch96_suffix,
			 file_size, expected_idx_size);
	for (int i = 0; i < infile_num; ++i)
		if (in_idx[i + 1] < in_idx[i])
			errx(EINVAL, "%s(): offsets are not monotonic at sample %d",
				 __func__, i);
	const uint64_t total_entries = in_idx[infile_num];
	const size_t comb_size = minco_downsample_file_size(set_opt->insketchpath,
														combined_sketch96_suffix);
	if (comb_size != (size_t)total_entries * sizeof(ctxobj96_t))
		errx(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
			 __func__, set_opt->insketchpath, combined_sketch96_suffix,
			 comb_size, (size_t)total_entries * sizeof(ctxobj96_t));

	minco_downsample_prepare_outdir(set_opt->outdir);
	minco_downsample_copy_optional_sidecar(set_opt->insketchpath, set_opt->outdir,
										   sketch_anno_stat, (size_t)infile_num * PATHLEN);
	minco_downsample_copy_optional_sidecar(set_opt->insketchpath, set_opt->outdir,
										   sketch_infile_meta_stat,
										   (size_t)infile_num * sizeof(infile_meta_t));
	minco_downsample_copy_optional_sidecar(set_opt->insketchpath, set_opt->outdir,
										   sketch_qc_stat,
										   (size_t)infile_num * sizeof(minco_sketch_qc_stat_t));
	minco_downsample_copy_optional_sidecar(set_opt->insketchpath, set_opt->outdir,
										   sketch_domain_stat,
										   (size_t)infile_num * sizeof(uint8_t));

	char *in_comb_path = test_get_fullpath(set_opt->insketchpath,
										   combined_sketch96_suffix);
	FILE *in_comb = fopen(in_comb_path, "rb");
	if (!in_comb)
		err(errno, "%s(): open %s", __func__, in_comb_path);
	char *out_comb_path = test_create_fullpath(set_opt->outdir,
											   combined_sketch96_suffix);
	FILE *out_comb = fopen(out_comb_path, "wb");
	if (!out_comb)
		err(errno, "%s(): open %s", __func__, out_comb_path);

	uint64_t *out_idx = calloc((size_t)infile_num + 1, sizeof(out_idx[0]));
	minco_ctxmeta_record_t *ctxmeta =
		calloc((size_t)infile_num, sizeof(ctxmeta[0]));
	if (!out_idx || !ctxmeta)
		err(errno, "%s(): OOM downsample output metadata", __func__);

	minco_ctxmeta_record_t *old_ctxmeta = NULL;
	size_t old_ctxmeta_size = 0;
	if (file_exists_in_folder(set_opt->insketchpath, minco_ctxmeta_bin_stat))
	{
		char *old_ctxmeta_path =
			test_get_fullpath(set_opt->insketchpath, minco_ctxmeta_bin_stat);
		old_ctxmeta = read_from_file(old_ctxmeta_path, &old_ctxmeta_size);
		free(old_ctxmeta_path);
		const size_t expected_ctxmeta_size =
			(size_t)infile_num * sizeof(old_ctxmeta[0]);
		if (old_ctxmeta_size != expected_ctxmeta_size)
			errx(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
				 __func__, set_opt->insketchpath, minco_ctxmeta_bin_stat,
				 old_ctxmeta_size, expected_ctxmeta_size);
	}

	minco_stat_density_summary_t density = {
		.min_threshold = UINT64_MAX,
		.max_threshold = UINT64_MAX,
		.universal_threshold = UINT64_MAX,
		.min_sample_id = MINCO_STAT_DENSITY_SAMPLE_ID_NONE,
		.max_sample_id = MINCO_STAT_DENSITY_SAMPLE_ID_NONE,
		.universal_sample_id = MINCO_STAT_DENSITY_SAMPLE_ID_NONE,
		.valid_sample_count = 0,
		.hash_bits = 64,
		.universal_policy = MINCO_STAT_DENSITY_POLICY_NONE,
		.flags = 0,
	};

	ctxobj96_t *buf = NULL;
	size_t buf_cap = 0;
	fprintf(stderr, "minco set downsample: %s -> %s; samples=%d; S=%u; payload=ctxobj96\n",
			set_opt->insketchpath, set_opt->outdir, infile_num,
			set_opt->sketch_size);
	for (int i = 0; i < infile_num; ++i)
	{
		const uint64_t sample_entries = in_idx[i + 1] - in_idx[i];
		if (sample_entries > buf_cap)
		{
			ctxobj96_t *nbuf = realloc(buf, (size_t)sample_entries * sizeof(buf[0]));
			if (!nbuf && sample_entries)
				err(errno, "%s(): OOM coden15 downsample buffer", __func__);
			buf = nbuf;
			buf_cap = (size_t)sample_entries;
		}
		minco_downsample_seek(in_comb, in_idx[i], sizeof(buf[0]),
							  combined_sketch96_suffix);
		minco_downsample_read_exact(in_comb, buf, sizeof(buf[0]),
									(size_t)sample_entries,
									combined_sketch96_suffix);
		const uint64_t keep64 =
			sample_entries < (uint64_t)set_opt->sketch_size
				? sample_entries
				: (uint64_t)set_opt->sketch_size;
		const size_t keep = (size_t)keep64;
		ctxobj96_select_bottom_by_hash(buf, (size_t)sample_entries, keep);
		out_idx[i + 1] = out_idx[i] + keep64;
		minco_downsample_write_exact(out_comb, buf, sizeof(buf[0]), keep,
									 combined_sketch96_suffix);

		uint8_t mode = MINCO_CTXMETA_BOTH;
		if (old_ctxmeta)
			mode = old_ctxmeta[i].mode;
		ctxmeta[i].mode = mode;
		ctxmeta[i].hash_bits = 64;
		ctxmeta[i].sketch_entries = keep64;
		if (keep > 0)
		{
			uint64_t threshold = 0;
			for (size_t j = 0; j < keep; ++j)
			{
				const uint64_t h = ctxobj96_hash_value(buf[j]);
				if (h > threshold)
					threshold = h;
			}
			const uint64_t unique_ctx = ctxobj96_count_ctx_runs(buf, keep);
			const uint64_t estimate =
				minco_downsample_density_estimate(unique_ctx, threshold, 64);
			ctxmeta[i].valid = 1;
			ctxmeta[i].threshold = threshold;
			ctxmeta[i].selected_observed_ctx = unique_ctx;
			ctxmeta[i].selected_estimated_unique_ctx = estimate;
			ctxmeta[i].preconflict_observed_ctx = unique_ctx;
			ctxmeta[i].preconflict_estimated_unique_ctx = estimate;
			ctxmeta[i].postconflict_observed_ctx = unique_ctx;
			ctxmeta[i].postconflict_estimated_unique_ctx = estimate;
			density.valid_sample_count++;
			if (threshold < density.min_threshold)
			{
				density.min_threshold = threshold;
				density.min_sample_id = (uint32_t)i;
			}
			if (threshold >= density.max_threshold ||
				density.max_sample_id == MINCO_STAT_DENSITY_SAMPLE_ID_NONE)
			{
				density.max_threshold = threshold;
				density.max_sample_id = (uint32_t)i;
			}
		}
		if ((i + 1) % 5000 == 0 || i + 1 == infile_num)
			fprintf(stderr,
					"minco set downsample: %d/%d samples processed; entries=%" PRIu64 "\n",
					i + 1, infile_num, out_idx[i + 1]);
	}
	if (fclose(in_comb) != 0 || fclose(out_comb) != 0)
		err(errno, "%s(): close coden15 sketch files", __func__);
	free(in_comb_path);
	free(out_comb_path);

	if (density.valid_sample_count > 0)
	{
		density.universal_threshold = density.max_threshold;
		density.universal_sample_id = density.max_sample_id;
		density.universal_policy =
			density.valid_sample_count == 1
				? MINCO_STAT_DENSITY_POLICY_SINGLE_SAMPLE
				: MINCO_STAT_DENSITY_POLICY_LARGEST_SAMPLE;
	}
	char *out_idx_path = test_create_fullpath(set_opt->outdir, idx_sketch96_suffix);
	write_to_file(out_idx_path, out_idx,
				  ((size_t)infile_num + 1) * sizeof(out_idx[0]));
	free(out_idx_path);
	char *out_ctxmeta_path =
		test_create_fullpath(set_opt->outdir, minco_ctxmeta_bin_stat);
	write_to_file(out_ctxmeta_path, ctxmeta,
				  (size_t)infile_num * sizeof(ctxmeta[0]));
	free(out_ctxmeta_path);

	minco_sketch_stat_t out_stat = *origin;
	out_stat.koc = 0;
	char *out_stat_path = test_create_fullpath(set_opt->outdir, sketch_stat);
	minco_stat_write_path_with_density(out_stat_path, &out_stat,
									   (const char (*)[PATHLEN])names,
									   set_opt->sketch_size,
									   MINCO_STAT_SELECTION_BOTTOMK,
									   info ? info->flags : 0,
									   UINT64_MAX,
									   &density);
	free(out_stat_path);
	fprintf(stderr,
			"minco set downsample: complete; entries=%" PRIu64
			"; universal_threshold=%" PRIu64 "\n",
			out_idx[infile_num], density.universal_threshold);

	if (old_ctxmeta)
		free_read_from_file(old_ctxmeta, old_ctxmeta_size);
	free(buf);
	free(out_idx);
	free(ctxmeta);
	free_read_from_file(in_idx, expected_idx_size);
	free_read_from_file(mem_stat, stat_size);
	return 1;
}

int minco_sketch_downsample(set_opt_t *set_opt)
{
	if (!set_opt)
		errx(EINVAL, "%s(): missing set options", __func__);
	if (set_opt->sketch_size == 0)
		errx(EINVAL, "%s(): --sketch-size must be positive", __func__);
	if (strcmp(set_opt->insketchpath, set_opt->outdir) == 0)
		errx(EINVAL, "%s(): output directory must differ from input sketch", __func__);

	size_t stat_size = 0;
	char *stat_path = test_get_fullpath(set_opt->insketchpath, sketch_stat);
	void *mem_stat = read_from_file(stat_path, &stat_size);
	free(stat_path);
	minco_sketch_info_t info = {0};
	if (!minco_stat_decode_mem(mem_stat, stat_size, &minco_stat_origin, &info))
		errx(EINVAL, "%s(): malformed %s/%s", __func__,
			 set_opt->insketchpath, sketch_stat);
	if (info.target_sketch_size && set_opt->sketch_size > info.target_sketch_size)
		errx(EINVAL, "%s(): requested -S %u is larger than source sketch target %u",
			 __func__, set_opt->sketch_size, info.target_sketch_size);
	const int infile_num = minco_stat_origin.infile_num;
	if (infile_num <= 0)
		errx(EINVAL, "%s(): input sketch has no samples", __func__);
	const uint32_t hash_bits = info.hash_bits ? info.hash_bits
											  : minco_stat_hash_bits_from_dim(&minco_stat_origin);
	if (hash_bits > 64)
		errx(EINVAL, "%s(): invalid hash_bits=%u", __func__, hash_bits);
	const uint32_t n_obj_bits = 64u - hash_bits;
	char (*names)[PATHLEN] = minco_stat_names_from_mem(mem_stat, stat_size);
	if (operate_stat_uses_ctxobj96_payload(&minco_stat_origin))
		return minco_sketch_downsample96(set_opt, mem_stat, stat_size, &info,
										 &minco_stat_origin, names);

	char *in_idx_path = test_get_fullpath(set_opt->insketchpath, idx_sketch_suffix);
	uint64_t *in_idx = read_from_file(in_idx_path, &file_size);
	free(in_idx_path);
	const size_t expected_idx_size = (size_t)(infile_num + 1) * sizeof(in_idx[0]);
	if (file_size != expected_idx_size)
		errx(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
			 __func__, set_opt->insketchpath, idx_sketch_suffix,
			 file_size, expected_idx_size);
	for (int i = 0; i < infile_num; ++i)
		if (in_idx[i + 1] < in_idx[i])
			errx(EINVAL, "%s(): offsets are not monotonic at sample %d",
				 __func__, i);
	const uint64_t total_entries = in_idx[infile_num];
	const size_t comb_size = minco_downsample_file_size(set_opt->insketchpath,
														combined_sketch_suffix);
	if (comb_size != (size_t)total_entries * sizeof(uint64_t))
		errx(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
			 __func__, set_opt->insketchpath, combined_sketch_suffix, comb_size,
			 (size_t)total_entries * sizeof(uint64_t));

	minco_downsample_prepare_outdir(set_opt->outdir);
	minco_downsample_copy_optional_sidecar(set_opt->insketchpath, set_opt->outdir,
										   sketch_anno_stat, (size_t)infile_num * PATHLEN);
	minco_downsample_copy_optional_sidecar(set_opt->insketchpath, set_opt->outdir,
										   sketch_infile_meta_stat,
										   (size_t)infile_num * sizeof(infile_meta_t));
	minco_downsample_copy_optional_sidecar(set_opt->insketchpath, set_opt->outdir,
										   sketch_qc_stat,
										   (size_t)infile_num * sizeof(minco_sketch_qc_stat_t));
	minco_downsample_copy_optional_sidecar(set_opt->insketchpath, set_opt->outdir,
										   sketch_domain_stat,
										   (size_t)infile_num * sizeof(uint8_t));

	char *in_comb_path = test_get_fullpath(set_opt->insketchpath, combined_sketch_suffix);
	FILE *in_comb = fopen(in_comb_path, "rb");
	if (!in_comb)
		err(errno, "%s(): open %s", __func__, in_comb_path);
	char *out_comb_path = test_create_fullpath(set_opt->outdir, combined_sketch_suffix);
	FILE *out_comb = fopen(out_comb_path, "wb");
	if (!out_comb)
		err(errno, "%s(): open %s", __func__, out_comb_path);
	FILE *in_ab = minco_downsample_open_optional_entry_file(
		set_opt->insketchpath, combined_ab_suffix, total_entries, sizeof(uint32_t));
	FILE *out_ab = NULL;
	if (in_ab)
	{
		char *out_ab_path = test_create_fullpath(set_opt->outdir, combined_ab_suffix);
		out_ab = fopen(out_ab_path, "wb");
		if (!out_ab)
			err(errno, "%s(): open %s", __func__, out_ab_path);
		free(out_ab_path);
	}
	FILE *in_pos = minco_downsample_open_optional_entry_file(
		set_opt->insketchpath, sketch_position_suffix, total_entries, sizeof(uint64_t));
	FILE *out_pos = NULL;
	if (in_pos)
	{
		char *out_pos_path = test_create_fullpath(set_opt->outdir, sketch_position_suffix);
		out_pos = fopen(out_pos_path, "wb");
		if (!out_pos)
			err(errno, "%s(): open %s", __func__, out_pos_path);
		free(out_pos_path);
	}

	uint64_t *buf = malloc((size_t)set_opt->sketch_size * sizeof(buf[0]));
	uint32_t *abuf = in_ab ? malloc((size_t)set_opt->sketch_size * sizeof(abuf[0])) : NULL;
	uint64_t *pbuf = in_pos ? malloc((size_t)set_opt->sketch_size * sizeof(pbuf[0])) : NULL;
	uint64_t *out_idx = calloc((size_t)infile_num + 1, sizeof(out_idx[0]));
	minco_ctxmeta_record_t *ctxmeta =
		calloc((size_t)infile_num, sizeof(ctxmeta[0]));
	if (!buf || (in_ab && !abuf) || (in_pos && !pbuf) || !out_idx || !ctxmeta)
		err(errno, "%s(): OOM downsample buffers", __func__);

	minco_ctxmeta_record_t *old_ctxmeta = NULL;
	size_t old_ctxmeta_size = 0;
	if (file_exists_in_folder(set_opt->insketchpath, minco_ctxmeta_bin_stat))
	{
		char *old_ctxmeta_path =
			test_get_fullpath(set_opt->insketchpath, minco_ctxmeta_bin_stat);
		old_ctxmeta = read_from_file(old_ctxmeta_path, &old_ctxmeta_size);
		free(old_ctxmeta_path);
		const size_t expected_ctxmeta_size =
			(size_t)infile_num * sizeof(old_ctxmeta[0]);
		if (old_ctxmeta_size != expected_ctxmeta_size)
			errx(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
				 __func__, set_opt->insketchpath, minco_ctxmeta_bin_stat,
				 old_ctxmeta_size, expected_ctxmeta_size);
	}

	minco_stat_density_summary_t density = {
		.min_threshold = UINT64_MAX,
		.max_threshold = UINT64_MAX,
		.universal_threshold = UINT64_MAX,
		.min_sample_id = MINCO_STAT_DENSITY_SAMPLE_ID_NONE,
		.max_sample_id = MINCO_STAT_DENSITY_SAMPLE_ID_NONE,
		.universal_sample_id = MINCO_STAT_DENSITY_SAMPLE_ID_NONE,
		.valid_sample_count = 0,
		.hash_bits = hash_bits,
		.universal_policy = MINCO_STAT_DENSITY_POLICY_NONE,
		.flags = 0,
	};

	fprintf(stderr, "minco set downsample: %s -> %s; samples=%d; S=%u\n",
			set_opt->insketchpath, set_opt->outdir, infile_num,
			set_opt->sketch_size);
	for (int i = 0; i < infile_num; ++i)
	{
		const uint64_t sample_entries = in_idx[i + 1] - in_idx[i];
		const uint64_t keep64 =
			sample_entries < (uint64_t)set_opt->sketch_size
				? sample_entries
				: (uint64_t)set_opt->sketch_size;
		const size_t keep = (size_t)keep64;
		out_idx[i + 1] = out_idx[i] + keep64;
		minco_downsample_seek(in_comb, in_idx[i], sizeof(uint64_t),
							  combined_sketch_suffix);
		minco_downsample_read_exact(in_comb, buf, sizeof(buf[0]), keep,
									combined_sketch_suffix);
		minco_downsample_write_exact(out_comb, buf, sizeof(buf[0]), keep,
									 combined_sketch_suffix);
		if (in_ab)
		{
			minco_downsample_seek(in_ab, in_idx[i], sizeof(uint32_t),
								  combined_ab_suffix);
			minco_downsample_read_exact(in_ab, abuf, sizeof(abuf[0]), keep,
										combined_ab_suffix);
			minco_downsample_write_exact(out_ab, abuf, sizeof(abuf[0]), keep,
										 combined_ab_suffix);
		}
		if (in_pos)
		{
			minco_downsample_seek(in_pos, in_idx[i], sizeof(uint64_t),
								  sketch_position_suffix);
			minco_downsample_read_exact(in_pos, pbuf, sizeof(pbuf[0]), keep,
										sketch_position_suffix);
			minco_downsample_write_exact(out_pos, pbuf, sizeof(pbuf[0]), keep,
										 sketch_position_suffix);
		}

		uint8_t mode = MINCO_CTXMETA_BOTH;
		if (old_ctxmeta)
			mode = old_ctxmeta[i].mode;
		ctxmeta[i].mode = mode;
		ctxmeta[i].hash_bits = hash_bits;
		ctxmeta[i].sketch_entries = keep64;
		if (keep > 0)
		{
			const uint64_t threshold =
				n_obj_bits == 64 ? 0 : (buf[keep - 1] >> n_obj_bits);
			const uint64_t unique_ctx =
				minco_downsample_unique_ctx_count(buf, keep, n_obj_bits);
			const uint64_t estimate =
				minco_downsample_density_estimate(unique_ctx, threshold, hash_bits);
			ctxmeta[i].valid = 1;
			ctxmeta[i].threshold = threshold;
			ctxmeta[i].selected_observed_ctx = unique_ctx;
			ctxmeta[i].selected_estimated_unique_ctx = estimate;
			ctxmeta[i].preconflict_observed_ctx = unique_ctx;
			ctxmeta[i].preconflict_estimated_unique_ctx = estimate;
			ctxmeta[i].postconflict_observed_ctx = unique_ctx;
			ctxmeta[i].postconflict_estimated_unique_ctx = estimate;
			density.valid_sample_count++;
			if (threshold < density.min_threshold)
			{
				density.min_threshold = threshold;
				density.min_sample_id = (uint32_t)i;
			}
			if (threshold >= density.max_threshold || density.max_sample_id == MINCO_STAT_DENSITY_SAMPLE_ID_NONE)
			{
				density.max_threshold = threshold;
				density.max_sample_id = (uint32_t)i;
			}
		}
		if ((i + 1) % 5000 == 0 || i + 1 == infile_num)
			fprintf(stderr,
					"minco set downsample: %d/%d samples processed; entries=%" PRIu64 "\n",
					i + 1, infile_num, out_idx[i + 1]);
	}
	if (fclose(in_comb) != 0 || fclose(out_comb) != 0)
		err(errno, "%s(): close sketch files", __func__);
	free(in_comb_path);
	free(out_comb_path);
	if (in_ab && (fclose(in_ab) != 0 || fclose(out_ab) != 0))
		err(errno, "%s(): close abundance files", __func__);
	if (in_pos && (fclose(in_pos) != 0 || fclose(out_pos) != 0))
		err(errno, "%s(): close position files", __func__);

	if (density.valid_sample_count > 0)
	{
		density.universal_threshold = density.max_threshold;
		density.universal_sample_id = density.max_sample_id;
		density.universal_policy =
			density.valid_sample_count == 1
				? MINCO_STAT_DENSITY_POLICY_SINGLE_SAMPLE
				: MINCO_STAT_DENSITY_POLICY_LARGEST_SAMPLE;
	}
	char *out_idx_path = test_create_fullpath(set_opt->outdir, idx_sketch_suffix);
	write_to_file(out_idx_path, out_idx, ((size_t)infile_num + 1) * sizeof(out_idx[0]));
	free(out_idx_path);
	char *out_ctxmeta_path =
		test_create_fullpath(set_opt->outdir, minco_ctxmeta_bin_stat);
	write_to_file(out_ctxmeta_path, ctxmeta, (size_t)infile_num * sizeof(ctxmeta[0]));
	free(out_ctxmeta_path);

	minco_sketch_stat_t out_stat = minco_stat_origin;
	out_stat.koc = in_ab != NULL;
	char *out_stat_path = test_create_fullpath(set_opt->outdir, sketch_stat);
	minco_stat_write_path_with_density(out_stat_path, &out_stat,
									   (const char (*)[PATHLEN])names,
									   set_opt->sketch_size,
									   MINCO_STAT_SELECTION_BOTTOMK,
									   info.flags,
									   UINT64_MAX,
									   &density);
	free(out_stat_path);
	fprintf(stderr,
			"minco set downsample: complete; entries=%" PRIu64
			"; universal_threshold=%" PRIu64 "\n",
			out_idx[infile_num], density.universal_threshold);

	if (old_ctxmeta)
		free_read_from_file(old_ctxmeta, old_ctxmeta_size);
	free(buf);
	free(abuf);
	free(pbuf);
	free(out_idx);
	free(ctxmeta);
	free_read_from_file(in_idx, expected_idx_size);
	free_read_from_file(mem_stat, stat_size);
	return 1;
}

int minco_sketch_union(set_opt_t *set_opt)
{ // for both union and uniq union

	void *mem_stat = read_from_file(test_get_fullpath(set_opt->insketchpath, sketch_stat), &file_size);
	const size_t stat_file_size = file_size;
	memcpy(&minco_stat_readin, mem_stat, sizeof(minco_stat_readin));
	const minco_payload_layout_spec_t payload =
		minco_payload_layout_for_stat(&minco_stat_readin);
	if (set_opt->operation == 3 && set_opt->q2markerdb_context)
	{
		minco_write_context_markerdb(set_opt, mem_stat, stat_file_size);
		free_all(mem_stat, NULL);
		return 1;
	}
	if (minco_stat_readin.infile_num == 1 && !set_opt->q2markerdb)
	{ // no need create
		printf("only 1 sketch, use %s as pan-sketch?(Y/N)\n", set_opt->insketchpath);
		char inpbuff;
		scanf(" %c", &inpbuff);
		if ((inpbuff == 'Y') || (inpbuff == 'y'))
		{
			chdir(set_opt->insketchpath);
			if (rename(payload.combined_suffix, payload.pan_suffix) != 0)
				err(errno, "minco_sketch_union()");
			printf("the union directory: %s created successfully\n", set_opt->insketchpath);
			return 1;
		}
	}
	if (payload.use_ctxobj96)
	{
		ctxobj96_t *mem_ctxobj = (ctxobj96_t *)read_from_file(
			test_get_fullpath(set_opt->insketchpath, payload.combined_suffix),
			&file_size);
		if (file_size % sizeof(mem_ctxobj[0]) != 0)
			errx(EINVAL, "%s(): %s/%s size %zu is not a multiple of %zu",
				 __func__, set_opt->insketchpath, payload.combined_suffix,
				 file_size, sizeof(mem_ctxobj[0]));
		const size_t in_kmer_ct = file_size / sizeof(mem_ctxobj[0]);
		ctxobj96_t *sorted_kmers = (ctxobj96_t *)malloc(file_size ? file_size : 1);
		if (!sorted_kmers)
			err(errno, "%s(): OOM copying %s/%s for set operation",
				__func__, set_opt->insketchpath, payload.combined_suffix);
		memcpy(sorted_kmers, mem_ctxobj, file_size);
		free_read_from_file(mem_ctxobj, file_size);
		mem_ctxobj = NULL;

		ctxobj96_sort_array(sorted_kmers, in_kmer_ct);

		if (set_opt->operation == 2)
		{
			const size_t union_ct =
				compact_union_or_unique_ctxobj96(sorted_kmers, in_kmer_ct, false);
			write_to_file(test_create_fullpath(set_opt->outdir, payload.pan_suffix),
						  sorted_kmers, union_ct * sizeof(sorted_kmers[0]));
		}
		else if (set_opt->operation == 3)
		{
			const size_t unique_ct =
				compact_union_or_unique_ctxobj96(sorted_kmers, in_kmer_ct, true);
			if (!set_opt->q2markerdb)
			{
				write_to_file(test_create_fullpath(set_opt->outdir, payload.uniq_pan_suffix),
							  sorted_kmers, unique_ct * sizeof(sorted_kmers[0]));
			}
			else
			{
				const size_t unique_bytes = unique_ct * sizeof(sorted_kmers[0]);
				ctxobj96_t *shrunk = realloc(sorted_kmers, unique_bytes ? unique_bytes : 1);
				if (shrunk || unique_bytes == 0)
					sorted_kmers = shrunk;

				const unsigned bucket_bits = 24;
				size_t *bucket_offsets = build_ctxobj96_bucket_offsets(
					sorted_kmers, unique_ct, bucket_bits);

				mem_ctxobj = (ctxobj96_t *)read_from_file(
					test_get_fullpath(set_opt->insketchpath, payload.combined_suffix),
					&file_size);
				const size_t comb_file_size = file_size;
				uint64_t *fco_pos = (uint64_t *)read_from_file(
					test_get_fullpath(set_opt->insketchpath, payload.idx_suffix),
					&file_size);
				uint64_t *post_fco_pos = calloc((minco_stat_readin.infile_num + 1),
												sizeof(uint64_t));
				if (!post_fco_pos)
					err(errno, "%s(): OOM output index", __func__);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(set_opt->p)
#endif
				for (uint32_t i = 0; i < minco_stat_readin.infile_num; i++)
				{
					uint64_t count = 0;
					for (uint64_t n = fco_pos[i]; n < fco_pos[i + 1]; n++)
					{
						if (ctxobj96_bucketed_contains(sorted_kmers, bucket_offsets,
													   bucket_bits, mem_ctxobj[n]))
							count++;
					}
					post_fco_pos[i + 1] = count;
				}

				for (uint32_t i = 0; i < minco_stat_readin.infile_num; i++)
					post_fco_pos[i + 1] += post_fco_pos[i];
				if (post_fco_pos[minco_stat_readin.infile_num] != unique_ct)
					err(EINVAL, "%s(): markerdb output count %lu != unique count %zu",
						__func__, post_fco_pos[minco_stat_readin.infile_num],
						unique_ct);

				char *out_comb_path =
					test_create_fullpath(set_opt->outdir, payload.combined_suffix);
				int outfd = open(out_comb_path, O_CREAT | O_TRUNC | O_WRONLY, 0644);
				if (outfd < 0)
					err(errno, "%s(): cannot create %s/%s", __func__,
						set_opt->outdir, payload.combined_suffix);
				if (ftruncate(outfd, (off_t)unique_bytes) != 0)
					err(errno, "%s(): resize %s", __func__, out_comb_path);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(set_opt->p)
#endif
				for (uint32_t i = 0; i < minco_stat_readin.infile_num; i++)
				{
					const uint64_t sample_count =
						post_fco_pos[i + 1] - post_fco_pos[i];
					if (sample_count == 0)
						continue;
					ctxobj96_t *sample_markers =
						malloc(sample_count * sizeof(sample_markers[0]));
					if (!sample_markers)
						err(errno, "%s(): OOM sample marker buffer", __func__);
					uint64_t out = 0;
					for (uint64_t n = fco_pos[i]; n < fco_pos[i + 1]; n++)
					{
						if (ctxobj96_bucketed_contains(sorted_kmers, bucket_offsets,
													   bucket_bits, mem_ctxobj[n]))
							sample_markers[out++] = mem_ctxobj[n];
					}
					if (out != sample_count)
						err(EINVAL, "%s(): marker count changed for sample %u",
							__func__, i);
					pwrite_all(
						outfd,
						sample_markers,
						sample_count * sizeof(sample_markers[0]),
						(off_t)(post_fco_pos[i] * sizeof(sample_markers[0])),
						out_comb_path);
					free(sample_markers);
				}
				if (close(outfd) != 0)
					err(errno, "%s(): close %s", __func__, out_comb_path);
				write_to_file(format_string("%s/%s", set_opt->outdir, payload.idx_suffix),
							  post_fco_pos,
							  (minco_stat_readin.infile_num + 1) * sizeof(post_fco_pos[0]));
				minco_warn_markerdb_small_refs(set_opt, post_fco_pos, mem_stat,
											   stat_file_size, "coden15 markerdb");
				free(out_comb_path);
				free(bucket_offsets);
				free_read_from_file(mem_ctxobj, comb_file_size);
				mem_ctxobj = NULL;
				free_all(fco_pos, post_fco_pos, NULL);
			}
		}
		else
			err(EINVAL, "operation value %d neither 2 (-u: union) nor 3 (-q :uniq uion )", set_opt->operation);
		write_to_file(test_create_fullpath(set_opt->outdir, sketch_stat),
					  mem_stat, stat_file_size);
		copy_minco_sketch_annotations(set_opt->insketchpath, set_opt->outdir,
									  minco_stat_readin.infile_num);
		copy_minco_domain_profiles(set_opt->insketchpath, set_opt->outdir,
								   minco_stat_readin.infile_num);
		free(sorted_kmers);
		free_all(mem_stat, NULL);
		return 1;
	}
	// union operation
	uint64_t *mem_ctxobj = (uint64_t *)read_from_file(test_get_fullpath(set_opt->insketchpath, combined_sketch_suffix), &file_size);
	const size_t in_kmer_ct = file_size / sizeof(mem_ctxobj[0]);
	uint64_t *sorted_kmers = (uint64_t *)malloc(file_size);
	if (!sorted_kmers)
		err(errno, "%s(): OOM copying %s/%s for set operation",
			__func__, set_opt->insketchpath, combined_sketch_suffix);
	memcpy(sorted_kmers, mem_ctxobj, file_size);
	free_read_from_file(mem_ctxobj, file_size);
	mem_ctxobj = NULL;

	sort_uint64_values(sorted_kmers, in_kmer_ct, set_opt->p);

	if (set_opt->operation == 2)
	{ // -u: normal union mode
		const size_t union_ct = compact_union_or_unique_uint64(sorted_kmers, in_kmer_ct, false);
		write_to_file(test_create_fullpath(set_opt->outdir, minco_pan_prefix),
					  sorted_kmers, union_ct * sizeof(sorted_kmers[0]));
	}
	else if (set_opt->operation == 3)
	{ // -q: uniq union mode
		const size_t unique_ct = compact_union_or_unique_uint64(sorted_kmers, in_kmer_ct, true);
		if (!set_opt->q2markerdb)
		{
			write_to_file(test_create_fullpath(set_opt->outdir, minco_uniq_pan_prefix),
						  sorted_kmers, unique_ct * sizeof(sorted_kmers[0]));
		}
		else
		{ // genereate markerdb directly instead of uniq union
			const size_t unique_bytes = unique_ct * sizeof(sorted_kmers[0]);
			uint64_t *shrunk = realloc(sorted_kmers, unique_bytes);
			if (shrunk || unique_bytes == 0)
				sorted_kmers = shrunk;

			const unsigned bucket_bits = 24;
			size_t *bucket_offsets = build_u64_highbit_bucket_offsets(
				sorted_kmers,
				unique_ct,
				bucket_bits);

			mem_ctxobj = (uint64_t *)read_from_file(test_get_fullpath(set_opt->insketchpath, combined_sketch_suffix), &file_size);
			const size_t comb_file_size = file_size;
			uint64_t *fco_pos = (uint64_t *)read_from_file(test_get_fullpath(set_opt->insketchpath, idx_sketch_suffix), &file_size);
			uint64_t *post_fco_pos = calloc((minco_stat_readin.infile_num + 1), sizeof(uint64_t));
			if (!post_fco_pos)
				err(errno, "%s(): OOM output index", __func__);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(set_opt->p)
#endif
			for (uint32_t i = 0; i < minco_stat_readin.infile_num; i++)
			{
				uint64_t count = 0;
				for (uint64_t n = fco_pos[i]; n < fco_pos[i + 1]; n++)
				{
					if (u64_bucketed_contains(sorted_kmers, bucket_offsets, bucket_bits, mem_ctxobj[n]))
						count++;
				}
				post_fco_pos[i + 1] = count;
			}

			for (uint32_t i = 0; i < minco_stat_readin.infile_num; i++)
				post_fco_pos[i + 1] += post_fco_pos[i];
			if (post_fco_pos[minco_stat_readin.infile_num] != unique_ct)
				err(EINVAL, "%s(): markerdb output count %lu != unique count %zu",
					__func__, post_fco_pos[minco_stat_readin.infile_num], unique_ct);

			char *out_comb_path = test_create_fullpath(set_opt->outdir, combined_sketch_suffix);
			int outfd = open(out_comb_path, O_CREAT | O_TRUNC | O_WRONLY, 0644);
			if (outfd < 0)
				err(errno, "%s(): cannot create %s/%s", __func__, set_opt->outdir, combined_sketch_suffix);
			if (ftruncate(outfd, (off_t)unique_bytes) != 0)
				err(errno, "%s(): resize %s", __func__, out_comb_path);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(set_opt->p)
#endif
			for (uint32_t i = 0; i < minco_stat_readin.infile_num; i++)
			{
				const uint64_t sample_count = post_fco_pos[i + 1] - post_fco_pos[i];
				if (sample_count == 0)
					continue;
				uint64_t *sample_markers = malloc(sample_count * sizeof(sample_markers[0]));
				if (!sample_markers)
					err(errno, "%s(): OOM sample marker buffer", __func__);
				uint64_t out = 0;
				for (uint64_t n = fco_pos[i]; n < fco_pos[i + 1]; n++)
				{
					if (u64_bucketed_contains(sorted_kmers, bucket_offsets, bucket_bits, mem_ctxobj[n]))
						sample_markers[out++] = mem_ctxobj[n];
				}
				if (out != sample_count)
					err(EINVAL, "%s(): marker count changed for sample %u", __func__, i);
				pwrite_all(
					outfd,
					sample_markers,
					sample_count * sizeof(sample_markers[0]),
					(off_t)(post_fco_pos[i] * sizeof(sample_markers[0])),
					out_comb_path);
				free(sample_markers);
			}
			if (close(outfd) != 0)
				err(errno, "%s(): close %s", __func__, out_comb_path);
			write_to_file(format_string("%s/%s", set_opt->outdir, idx_sketch_suffix), post_fco_pos, (minco_stat_readin.infile_num + 1) * sizeof(post_fco_pos[0]));
			minco_warn_markerdb_small_refs(set_opt, post_fco_pos, mem_stat, stat_file_size,
										   "markerdb");
			free(out_comb_path);
			free(bucket_offsets);
			free_read_from_file(mem_ctxobj, comb_file_size);
			mem_ctxobj = NULL;
			free_all(fco_pos, post_fco_pos, NULL);
		}
	}
	else
		err(EINVAL, "operation value %d neither 2 (-u: union) nor 3 (-q :uniq uion )", set_opt->operation);
	write_to_file(test_create_fullpath(set_opt->outdir, sketch_stat), mem_stat, stat_file_size);
	copy_minco_sketch_annotations(set_opt->insketchpath, set_opt->outdir, minco_stat_readin.infile_num);
	copy_minco_domain_profiles(set_opt->insketchpath, set_opt->outdir, minco_stat_readin.infile_num);
	free(sorted_kmers);
	free_all(mem_stat, NULL);
	return 1;
}

int minco_sketch_operate(set_opt_t *set_opt)
{
	clock_t start_time = clock();
	void *mem_stat_pan = read_from_file(test_get_fullpath(set_opt->pansketchpath, sketch_stat), &file_size);
	memcpy(&minco_stat_pan, mem_stat_pan, sizeof(minco_stat_pan));
	void *mem_stat_minco = read_from_file(test_get_fullpath(set_opt->insketchpath, sketch_stat), &file_size);
	memcpy(&minco_stat_origin, mem_stat_minco, sizeof(minco_stat_origin));
	if (minco_stat_pan.hash_id != minco_stat_origin.hash_id)
		err(EXIT_FAILURE, "%s(): %s sketcing id %u != %s id %u", __func__, set_opt->pansketchpath, minco_stat_origin.hash_id, set_opt->insketchpath, minco_stat_pan.hash_id);
	const minco_payload_layout_spec_t pan_payload =
		minco_payload_layout_for_stat(&minco_stat_pan);
	const minco_payload_layout_spec_t origin_payload =
		minco_payload_layout_for_stat(&minco_stat_origin);
	if (pan_payload.use_ctxobj96 != origin_payload.use_ctxobj96)
		errx(EXIT_FAILURE, "%s(): cannot combine ctxobj64 and ctxobj96 set inputs",
			 __func__);
	if (minco_stat_origin.koc)
		printf("%s() Warning: k-mer abundances are dropped in this sketch operation\n ", __func__);
	// copy sketch stat file to result sketch
	write_to_file(test_create_fullpath(set_opt->outdir, sketch_stat), mem_stat_minco, file_size);
	copy_minco_sketch_annotations(set_opt->insketchpath, set_opt->outdir, minco_stat_origin.infile_num);
	copy_minco_domain_profiles(set_opt->insketchpath, set_opt->outdir, minco_stat_origin.infile_num);

	if (origin_payload.use_ctxobj96)
	{
		char *ctxobj_fpath;
		if (file_exists_in_folder(set_opt->pansketchpath, pan_payload.pan_suffix))
			ctxobj_fpath = test_get_fullpath(set_opt->pansketchpath,
											 pan_payload.pan_suffix);
		else if (file_exists_in_folder(set_opt->pansketchpath,
									   pan_payload.uniq_pan_suffix))
			ctxobj_fpath = test_get_fullpath(set_opt->pansketchpath,
											 pan_payload.uniq_pan_suffix);
		else
			err(EXIT_FAILURE, "%s(): cannot find %s or %s under %s",
				__func__, pan_payload.pan_suffix, pan_payload.uniq_pan_suffix,
				set_opt->pansketchpath);

		size_t pan_file_size = 0;
		ctxobj96_t *mem_pan = (ctxobj96_t *)read_from_file(ctxobj_fpath,
														   &pan_file_size);
		if (pan_file_size % sizeof(mem_pan[0]) != 0)
			errx(EINVAL, "%s(): pan file %s size %zu is not a multiple of %zu",
				 __func__, ctxobj_fpath, pan_file_size, sizeof(mem_pan[0]));
		const size_t pan_ct = pan_file_size / sizeof(mem_pan[0]);
		ctxobj96_sort_array(mem_pan, pan_ct);
		const unsigned bucket_bits = 24;
		size_t *bucket_offsets =
			build_ctxobj96_bucket_offsets(mem_pan, pan_ct, bucket_bits);

		size_t idx_file_size = 0;
		uint64_t *fco_pos = (uint64_t *)read_from_file(
			test_get_fullpath(set_opt->insketchpath, origin_payload.idx_suffix),
			&idx_file_size);
		uint64_t *post_fco_pos =
			calloc((minco_stat_origin.infile_num + 1), sizeof(uint64_t));
		size_t ctxobj_file_size = 0;
		ctxobj96_t *tmp_ctxobj_mem = (ctxobj96_t *)read_from_file(
			test_get_fullpath(set_opt->insketchpath, origin_payload.combined_suffix),
			&ctxobj_file_size);
		ctxobj96_t *post_ctxobj_mem = (ctxobj96_t *)malloc(ctxobj_file_size ? ctxobj_file_size : 1);
		if (!post_fco_pos || !post_ctxobj_mem)
			err(errno, "%s(): OOM ctxobj96 set operation buffers", __func__);
		uint64_t post_kmer_ct = 0;

		for (uint32_t i = 0; i < minco_stat_origin.infile_num; i++)
		{
			for (uint64_t n = fco_pos[i]; n < fco_pos[i + 1]; n++)
			{
				if (set_opt->operation ==
					(ctxobj96_bucketed_contains(mem_pan, bucket_offsets,
												 bucket_bits, tmp_ctxobj_mem[n]) != 0))
					post_ctxobj_mem[post_kmer_ct++] = tmp_ctxobj_mem[n];
			}
			post_fco_pos[i + 1] = post_kmer_ct;
		}

		sprintf(outfpath, "%s/%s", set_opt->outdir, origin_payload.combined_suffix);
		write_to_file(outfpath, post_ctxobj_mem,
					  post_kmer_ct * sizeof(post_ctxobj_mem[0]));
		sprintf(outfpath, "%s/%s", set_opt->outdir, origin_payload.idx_suffix);
		write_to_file(outfpath, post_fco_pos,
					  (minco_stat_origin.infile_num + 1) * sizeof(post_fco_pos[0]));

		free(ctxobj_fpath);
		free(bucket_offsets);
		free_read_from_file(mem_pan, pan_file_size);
		free_read_from_file(fco_pos, idx_file_size);
		free_read_from_file(tmp_ctxobj_mem, ctxobj_file_size);
		free(post_fco_pos);
		free(post_ctxobj_mem);
		free_all(mem_stat_pan, mem_stat_minco, NULL);
		return 1;
	}

	char *ctxobj_fpath;
	if (file_exists_in_folder(set_opt->pansketchpath, minco_pan_prefix))
		ctxobj_fpath = test_get_fullpath(set_opt->pansketchpath, minco_pan_prefix);
	else if (file_exists_in_folder(set_opt->pansketchpath, minco_uniq_pan_prefix))
		ctxobj_fpath = test_get_fullpath(set_opt->pansketchpath, minco_uniq_pan_prefix);
	else
		err(EXIT_FAILURE, "%s():cannot find %s or %s under %s ", __func__, minco_pan_prefix, minco_uniq_pan_prefix, set_opt->pansketchpath);

	uint64_t *mem_pan = (uint64_t *)read_from_file(ctxobj_fpath, &file_size);
	khash_t(kmer_set) *h = kh_init(kmer_set);
	uint32_t kmer_ct = file_size / sizeof(uint64_t);
	for (int i = 0; i < kmer_ct; i++)
		kh_put(kmer_set, h, mem_pan[i], &ret);
		// read sketch offset table to memory
	uint64_t *fco_pos = (uint64_t *)read_from_file(test_get_fullpath(set_opt->insketchpath, idx_sketch_suffix), &file_size);
	// post operation index	calloc
	uint64_t *post_fco_pos = calloc((minco_stat_origin.infile_num + 1), sizeof(uint64_t));
		// read context-object payload to memory
	uint64_t *tmp_ctxobj_mem = (uint64_t *)read_from_file(test_get_fullpath(set_opt->insketchpath, combined_sketch_suffix), &file_size);
	uint64_t *post_ctxobj_mem = (uint64_t *)malloc(file_size);
	uint32_t post_kmer_ct = 0;

	// sketch operation
	for (uint32_t i = 0; i < minco_stat_origin.infile_num; i++)
	{
		for (uint64_t n = fco_pos[i]; n < fco_pos[i + 1]; n++)
		{
			// make sure set_opt->operation == 0 if subtract, == 1 if intersect
			if (set_opt->operation == (kh_get(kmer_set, h, tmp_ctxobj_mem[n]) != kh_end(h)))
				post_ctxobj_mem[post_kmer_ct++] = tmp_ctxobj_mem[n];
		}
		post_fco_pos[i + 1] = post_kmer_ct;
	}
	kh_destroy(kmer_set, h);
		// write filtered context-object payload
	sprintf(outfpath, "%s/%s", set_opt->outdir, combined_sketch_suffix);
	write_to_file(outfpath, post_ctxobj_mem, post_kmer_ct * sizeof(post_ctxobj_mem[0]));
	// write index
	sprintf(outfpath, "%s/%s", set_opt->outdir, idx_sketch_suffix);
	write_to_file(outfpath, post_fco_pos, (minco_stat_origin.infile_num + 1) * sizeof(post_fco_pos[0]));

	free_all(mem_stat_pan, mem_stat_minco, ctxobj_fpath, mem_pan, fco_pos, post_fco_pos, tmp_ctxobj_mem, post_ctxobj_mem, NULL);
	return 1;
}

int minco_group_samples(set_opt_t *set_opt)
{ // for sorted minco context-object sketches

	void *mem_stat = read_from_file(test_get_fullpath(set_opt->insketchpath, sketch_stat), &file_size);
	memcpy(&minco_stat_readin, mem_stat, sizeof(minco_stat_readin));

	compan_t *subset = organize_taxf(set_opt->subsetf);
	if (minco_stat_readin.infile_num != subset->gn)
		err(EXIT_FAILURE, "%s(): %s's genome number %d != %s's line number %d", __func__, set_opt->insketchpath, minco_stat_readin.infile_num, set_opt->subsetf, subset->gn);

	const minco_payload_layout_spec_t payload =
		minco_payload_layout_for_stat(&minco_stat_readin);
	if (payload.use_ctxobj96)
	{
		uint64_t *tmp_idx = (uint64_t *)read_from_file(
			test_get_fullpath(set_opt->insketchpath, payload.idx_suffix),
			&file_size);
		const size_t idx_file_size = file_size;
		ctxobj96_t *mem_ctxobj = (ctxobj96_t *)read_from_file(
			test_get_fullpath(set_opt->insketchpath, payload.combined_suffix),
			&file_size);
		const size_t ctxobj_file_size = file_size;

		if (ctxobj_file_size % sizeof(mem_ctxobj[0]) != 0 ||
			tmp_idx[subset->gn] * sizeof(mem_ctxobj[0]) != ctxobj_file_size)
			err(EXIT_FAILURE, "%s(): %s and %s sizes are inconsistent",
				__func__, payload.idx_suffix, payload.combined_suffix);

		char *ctxobj_fpath = test_create_fullpath(set_opt->outdir,
												  payload.combined_suffix);
		int fd = open(ctxobj_fpath, O_CREAT | O_RDWR, 0644);
		if (fd < 0)
			err(errno, "%s(): open %s", __func__, ctxobj_fpath);
		if (ftruncate(fd, (off_t)ctxobj_file_size) == -1)
			err(EXIT_FAILURE, "%s(): ftruncate %s", __func__, ctxobj_fpath);
		ctxobj96_t *grouped_ctxobj =
			(ctxobj96_t *)mmap(NULL, ctxobj_file_size, PROT_READ | PROT_WRITE,
							   MAP_SHARED, fd, 0);
		if (grouped_ctxobj == MAP_FAILED)
			err(EXIT_FAILURE, "%s(): mmap %s", __func__, ctxobj_fpath);

		uint64_t *out_idx = calloc((subset->taxn + 1), sizeof(uint64_t));
		if (!out_idx)
			err(errno, "%s(): OOM grouped index", __func__);
		int outfn = 0;

		for (int t = 0; t < subset->taxn; t++)
		{
			if (subset->tax[t].taxid == 0)
				continue;
			outfn++;
			uint64_t start_offset = out_idx[outfn] = out_idx[outfn - 1];
			for (int n = 1; n <= subset->tax[t].gids[0]; n++)
			{
				int gid = subset->tax[t].gids[n];
				const uint64_t n_entries = tmp_idx[gid + 1] - tmp_idx[gid];
				memcpy(grouped_ctxobj + start_offset, mem_ctxobj + tmp_idx[gid],
					   n_entries * sizeof(grouped_ctxobj[0]));
				start_offset += n_entries;
			}

			if (subset->tax[t].gids[0] > 1)
			{
				ctxobj96_sort_array(grouped_ctxobj + out_idx[outfn],
									start_offset - out_idx[outfn]);
				size_t len = compact_union_or_unique_ctxobj96(
					grouped_ctxobj + out_idx[outfn],
					start_offset - out_idx[outfn], false);
				out_idx[outfn] += len;
			}
			else
				out_idx[outfn] = start_offset;
		}
		if (msync(grouped_ctxobj, out_idx[outfn] * sizeof(grouped_ctxobj[0]),
				  MS_SYNC) == -1)
			err(EXIT_FAILURE, "%s(): msync error", __func__);
		if (ftruncate(fd, (off_t)(out_idx[outfn] * sizeof(grouped_ctxobj[0]))) == -1)
			err(EXIT_FAILURE, "%s(): ftruncate resize failed", __func__);
		if (munmap(grouped_ctxobj, ctxobj_file_size) == -1)
			err(EXIT_FAILURE, "%s(): munmap", __func__);
		close(fd);

		write_to_file(test_create_fullpath(set_opt->outdir, payload.idx_suffix),
					  out_idx, (outfn + 1) * sizeof(out_idx[0]));

		minco_stat_readin.infile_num = outfn;
		minco_stat_readin.koc = 0;
		char (*input_anno)[PATHLEN] = NULL;
		size_t input_anno_size = 0;
		if (file_exists_in_folder(set_opt->insketchpath, sketch_anno_stat))
		{
			char *anno_path = test_get_fullpath(set_opt->insketchpath, sketch_anno_stat);
			input_anno = read_from_file(anno_path, &input_anno_size);
			free(anno_path);
			const size_t expected_size = (size_t)subset->gn * PATHLEN;
			if (input_anno_size != expected_size)
				err(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
					__func__, set_opt->insketchpath, sketch_anno_stat,
					input_anno_size, expected_size);
		}

		char (*tmpfname)[PATHLEN] = calloc(outfn, PATHLEN);
		char (*tmpanno)[PATHLEN] = input_anno ? calloc(outfn, PATHLEN) : NULL;
		if (!tmpfname || (input_anno && !tmpanno))
			err(errno, "%s(): OOM output names/annotations", __func__);
		int idx = 0;
		for (int t = 0; t < subset->taxn; t++)
		{
			if (subset->tax[t].taxid != 0)
			{
				if (subset->tax[t].taxname != NULL)
					snprintf(tmpfname[idx], PATHLEN, "%d_%s",
							 subset->tax[t].taxid, subset->tax[t].taxname);
				else
					snprintf(tmpfname[idx], PATHLEN, "%d", subset->tax[t].taxid);
				if (tmpanno)
				{
					if (subset->tax[t].gids[0] == 1)
						memcpy(tmpanno[idx], input_anno[subset->tax[t].gids[1]],
							   PATHLEN);
					else
						snprintf(tmpanno[idx], PATHLEN, "%s", tmpfname[idx]);
				}
				idx++;
			}
			free_all(subset->tax[t].gids, subset->tax[t].taxname, NULL);
		}
		concat_and_write_to_file(test_create_fullpath(set_opt->outdir, sketch_stat),
								 &minco_stat_readin, sizeof(minco_stat_readin),
								 tmpfname, PATHLEN * outfn);
		if (tmpanno)
			write_to_file(test_create_fullpath(set_opt->outdir, sketch_anno_stat),
						  tmpanno, (size_t)outfn * PATHLEN);

		if (input_anno)
			free_read_from_file(input_anno, input_anno_size);
		free(ctxobj_fpath);
		free_read_from_file(mem_ctxobj, ctxobj_file_size);
		free_read_from_file(tmp_idx, idx_file_size);
		free_all(mem_stat, out_idx, subset->tax, subset, tmpfname, tmpanno, NULL);
		return outfn;
	}

		// read offsets and context-object payload
	uint64_t *tmp_idx = (uint64_t *)read_from_file(test_get_fullpath(set_opt->insketchpath, idx_sketch_suffix), &file_size);
	uint64_t *mem_ctxobj = (uint64_t *)read_from_file(test_get_fullpath(set_opt->insketchpath, combined_sketch_suffix), &file_size);

	if (tmp_idx[subset->gn] * sizeof(uint64_t) != file_size)
		err(EXIT_FAILURE, "%s(): %s last(%u) index(%lu) * sizeof(uint64_t) != %s file size (%lu) ",
			__func__, idx_sketch_suffix, subset->gn, tmp_idx[subset->gn], combined_sketch_suffix, file_size);
		// output offsets and context-object payload
	char *ctxobj_fpath = test_create_fullpath(set_opt->outdir, combined_sketch_suffix);
	int fd = open(ctxobj_fpath, O_CREAT | O_RDWR, 0644);
	ftruncate(fd, file_size); // allowing much larger grouped_ctxobj than using malloc
	uint64_t *grouped_ctxobj = (uint64_t *)mmap(NULL, file_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);

	uint64_t *out_idx = calloc((subset->taxn + 1), sizeof(uint64_t));
	int outfn = 0; // uint64_t grouped_kmer_ct = 0;

	for (int t = 0; t < subset->taxn; t++)
	{
		if (subset->tax[t].taxid == 0)
			continue; // ignore taxid 0
		outfn++;
		uint64_t start_offset = out_idx[outfn] = out_idx[outfn - 1];
		for (int n = 1; n <= subset->tax[t].gids[0]; n++)
		{
			int gid = subset->tax[t].gids[n];
			memcpy(grouped_ctxobj + start_offset, mem_ctxobj + tmp_idx[gid], (tmp_idx[gid + 1] - tmp_idx[gid]) * sizeof(grouped_ctxobj[0]));
			start_offset += (tmp_idx[gid + 1] - tmp_idx[gid]);
		}

		if (subset->tax[t].gids[0] > 1)
		{
			qsort(grouped_ctxobj + out_idx[outfn], start_offset - out_idx[outfn], sizeof(grouped_ctxobj[0]), qsort_comparator_uint64);
			size_t len = dedup_sorted_uint64(grouped_ctxobj + out_idx[outfn], start_offset - out_idx[outfn]);
			out_idx[outfn] += len;
		}
		else
			out_idx[outfn] = start_offset;
	}
	// write grouped kmer and index to result
	// write_to_file(test_create_fullpath(set_opt->outdir,combined_sketch_suffix), grouped_ctxobj, grouped_kmer_ct*sizeof(grouped_ctxobj[0]));
	if (msync(grouped_ctxobj, out_idx[outfn] * sizeof(grouped_ctxobj[0]), MS_SYNC) == -1)
		err(EXIT_FAILURE, "%s(): msync error", __func__);
	if (ftruncate(fd, out_idx[outfn] * sizeof(grouped_ctxobj[0])) == -1)
		err(EXIT_FAILURE, "%s(): ftrucate resize failed", __func__);
	if (munmap(grouped_ctxobj, file_size) == -1)
		err(EXIT_FAILURE, "%s(): munmap", __func__);
	close(fd);

	write_to_file(test_create_fullpath(set_opt->outdir, idx_sketch_suffix), out_idx, (outfn + 1) * sizeof(out_idx[0]));
	// write stat file
	minco_stat_readin.infile_num = outfn;
	minco_stat_readin.koc = 0;
	char (*input_anno)[PATHLEN] = NULL;
	size_t input_anno_size = 0;
	if (file_exists_in_folder(set_opt->insketchpath, sketch_anno_stat))
	{
		char *anno_path = test_get_fullpath(set_opt->insketchpath, sketch_anno_stat);
		input_anno = read_from_file(anno_path, &input_anno_size);
		free(anno_path);
		const size_t expected_size = (size_t)subset->gn * PATHLEN;
		if (input_anno_size != expected_size)
			err(EINVAL, "%s(): %s/%s has %zu bytes, expected %zu",
				__func__, set_opt->insketchpath, sketch_anno_stat,
				input_anno_size, expected_size);
	}

	char (*tmpfname)[PATHLEN] = calloc(outfn, PATHLEN);
	char (*tmpanno)[PATHLEN] = input_anno ? calloc(outfn, PATHLEN) : NULL;
	if (!tmpfname || (input_anno && !tmpanno))
		err(errno, "%s(): OOM output names/annotations", __func__);
	int idx = 0;
	for (int t = 0; t < subset->taxn; t++)
	{
		if (subset->tax[t].taxid != 0)
		{
			if (subset->tax[t].taxname != NULL)
				snprintf(tmpfname[idx], PATHLEN, "%d_%s", subset->tax[t].taxid, subset->tax[t].taxname);
			else
				snprintf(tmpfname[idx], PATHLEN, "%d", subset->tax[t].taxid);
			if (tmpanno)
			{
				if (subset->tax[t].gids[0] == 1)
					memcpy(tmpanno[idx], input_anno[subset->tax[t].gids[1]], PATHLEN);
				else
					snprintf(tmpanno[idx], PATHLEN, "%s", tmpfname[idx]);
			}
			idx++;
		}
		free_all(subset->tax[t].gids, subset->tax[t].taxname, NULL);
	}
	concat_and_write_to_file(test_create_fullpath(set_opt->outdir, sketch_stat), &minco_stat_readin, sizeof(minco_stat_readin), tmpfname, PATHLEN * outfn);
	if (tmpanno)
		write_to_file(test_create_fullpath(set_opt->outdir, sketch_anno_stat), tmpanno, (size_t)outfn * PATHLEN);
	if (munmap(mem_ctxobj, file_size) == -1)
	{
		if (mem_ctxobj != NULL)
			free(mem_ctxobj);
		else
			err(EXIT_FAILURE, "%s(): munmap", __func__);
	}

	if (input_anno)
		free_read_from_file(input_anno, input_anno_size);
	free_all(mem_stat, tmp_idx, out_idx, subset->tax, subset, tmpfname, tmpanno, NULL);
	return outfn;
}
