/* Hidden sketch-set maintenance helpers for minco context-object sketches. */
#include "global_basic.h"
#include "command_ani.h"
#include "command_operate.h"
#include "command_sketch_wrapper.h"
#include "sketch_inspect.h"
#include "../klib/khash.h"
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

const char minco_pan_prefix[] = "lpan"; // uint64_t pan
const char minco_uniq_pan_prefix[] = "luniq_pan";
// common vars
static size_t file_size;
static minco_sketch_stat_t minco_stat_readin, minco_stat_pan, minco_stat_origin;
static struct stat s;
static char outfpath[PATHLEN + 20];
static int ret;
extern const char sorted_comb_ctxgid64obj32[];

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

void sketch_inspect_print_samples(const char *sketch_path)
{
	void *mem_stat = read_from_file(test_get_fullpath(sketch_path, sketch_stat), &file_size);
	memcpy(&minco_stat_readin, mem_stat, sizeof(minco_stat_readin));
	char (*tmpname)[PATHLEN] = mem_stat + sizeof(minco_sketch_stat_t);
	uint64_t *mem_index = (uint64_t *)read_from_file(test_get_fullpath(sketch_path, idx_sketch_suffix), &file_size);
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
	if (show_mode == 1)
	{
		uint64_t *sketch_index = read_from_file(test_get_fullpath(sketch_path, idx_sketch_suffix), &file_size);
		int infile_num = file_size / sizeof(uint64_t) - 1;
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
		idx_sketch_suffix,
		minco_ctxmeta_bin_stat,
		combined_ab_suffix,
		sketch_position_suffix,
		sorted_comb_ctxgid64obj32,
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
	if (minco_stat_readin.infile_num == 1)
	{ // no need create
		printf("only 1 sketch, use %s as pan-sketch?(Y/N)\n", set_opt->insketchpath);
		char inpbuff;
		scanf(" %c", &inpbuff);
		if ((inpbuff == 'Y') || (inpbuff == 'y'))
		{
			chdir(set_opt->insketchpath);
			if (rename(combined_sketch_suffix, minco_pan_prefix) != 0)
				err(errno, "minco_sketch_union()");
			printf("the union directory: %s created successfully\n", set_opt->insketchpath);
			return 1;
		}
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
	if (minco_stat_origin.koc)
		printf("%s() Warning: k-mer abundances are dropped in this sketch operation\n ", __func__);
	// copy sketch stat file to result sketch
	write_to_file(test_create_fullpath(set_opt->outdir, sketch_stat), mem_stat_minco, file_size);
	copy_minco_sketch_annotations(set_opt->insketchpath, set_opt->outdir, minco_stat_origin.infile_num);

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
