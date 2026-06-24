#ifndef GLOBAL_BASIC
#define GLOBAL_BASIC
#include <stdbool.h>
#include <stdint.h>
#include <limits.h>
/*CTX/OBJ/Kmer maskers*/
#define _64MASK 0xffffffffffffffffLLU
// #define TUPMASK (_64MASK >> (64-BITTL))
#define BIT1MASK 0x0000000000000001LLU
/*compile mode for different alphabet: 1:nt reduction,2:AA*/
#if ALPHABET == 1 // objs compatible mode
#define DEFAULT 15
#define OBJ_ALPH 16
#define OBJ_BITS 4
#define BIN_SZ 4096 // maximum sequences allowed in one bin
                    //	#include <stdbool.h>
extern const bool Objdist[16][16];
#define IS_DIFF(X, Y) (Objdist[(X)][(Y)])

#elif ALPHABET == 2 // AA seq
#define DEFAULT (-1)
#define OBJ_BITS 5
#define BIN_SZ 2048
#define IS_DIFF(X, Y) ((X) != (Y))
#else
#define DEFAULT (-1)
#define OBJ_ALPH 4
#define OBJ_BITS 2
#define BIN_SZ 65536 // 50 //for test //16384
#define IS_DIFF(X, Y) ((X) != (Y))
#endif
/**/
#define LMAX 4096      // characters limit per line for any documents might be analysis
#define PATHLEN 256    // limit for input file path length
#define MCO_BUF_S 4096 // buf size of reading/analysis mco
                       // COMPONENT_SZ control .co file divided how many components
#ifndef COMPONENT_SZ
#define COMPONENT_SZ 7 // 6 or 7,default unit component dimension size = 16^(COMPONENT_SZ) or 1<<4*(COMPONENT_SZ)
#endif

#ifndef CTX_SPC_USE_L
#define CTX_SPC_USE_L 8 // 8 for small mem.(< 1g), 4 for larger mem.  //ctx space occupy rate limit = 1/(1<<CTX_SPC_USE_L)
#endif

#define CTX_DR_LMT 100 // limit for least CTX after dimensionality reduction
#define LD_FCTR 0.6    // hash function load factor

//*****core data struct *****//
// gid array linked list: .mco genome id+obj array size
#define GID_ARR_SZ 16       // 32 //short arr[GID_ARR_SZ]
#define BBILLION 1073741824 // binary billion

/*type define */
typedef unsigned long long int llong;
/*argp wrapper*/
#include <argp.h>
/** Local Prototypes **/
struct arg_global
{
  int verbosity;
};

void log_printf(struct arg_global *g, int level, const char *fmt, ...);
#define ARGP_KEY_INVALID 16777219

// const char co_dstat[] = "cofiles.stat"; //stat file name in co dir
// const char mco_dstat[] = "mcofiles.stat";

//*****basic function prototypes ****/

FILE *fpathopen(const char *dpath, const char *fname, const char *mode);
double get_sys_mmry(void);
/*completeReverse*/
#define SWAP2 0x3333333333333333ULL
#define SWAP4 0x0F0F0F0F0F0F0F0FULL
#define SWAP8 0x00FF00FF00FF00FFULL
#define SWAP16 0x0000FFFF0000FFFFULL
#define SWAP32 0x00000000FFFFFFFFULL
static inline llong crvs64bits(llong n)
{

  n = ((n >> 2) & SWAP2) | ((n & SWAP2) << 2);
  n = ((n >> 4) & SWAP4) | ((n & SWAP4) << 4);
  n = ((n >> 8) & SWAP8) | ((n & SWAP8) << 8);
  n = ((n >> 16) & SWAP16) | ((n & SWAP16) << 16);
  n = ((n >> 32) & SWAP32) | ((n & SWAP32) << 32);
  return ~n;
}

/*basemap for different alphabet */
extern const int Basemap[128];
extern const char Mapbase[];
extern const unsigned int primer[25];
llong find_lgst_primer_2pow(int w);
#include <stdint.h>
uint32_t nextPrime(uint32_t);
/* orgnized infile table*/

typedef struct infile_entry
{
  size_t fsize;
  char *fpath;
} infile_entry_t;

typedef struct infile_tab
{
  int infile_num;
  infile_entry_t *organized_infile_tab; // infile_entry_t arr
} infile_tab_t;

#define BASENAME_LEN 128 // largest allowed bytes for input file basename
typedef struct bin_stat
{
  // estimated kmer count sum accross files in the bin (dimension reduction rate not comsidered here)
  llong est_kmc_bf_dr;
  char (*seqfilebasename)[BASENAME_LEN];
  llong exact_co_mem; // exact co file size in memory, estimated after stage I
} bin_stat_t;

infile_tab_t *organize_infile_list(char *list_path, int fmt_ck);
infile_tab_t *organize_infile_frm_arg(int num_remaining_args, char **remaining_args, int fmt_ck);
// per bin get file basename and estimate .co files (before dimension reduction) memory usage.
bin_stat_t *get_bin_basename_stat(infile_entry_t *organized_infile_tab, int *shuffle_arr, int binsz);

/*legency combco stat file type*/
typedef struct co_dirstat
{
  unsigned int shuf_id;
  bool koc;       // kmer occurence or not
  int kmerlen;    // 2*k
  int dim_rd_len; // 2*drlevel
  int comp_num;   // components number
  int infile_num; // .co file num, .co file with many components only count as 1
  // int comp_sz[comp_num * infile_num]; // do this if storage all .co per mco bin in one file
  llong all_ctx_ct;
} co_dstat_t;

/* minco sketch stat file type */
typedef struct minco_sketch_stat
{
  uint32_t hash_id; // sketching type coding
  bool koc;
  bool conflict; // with conflict objs ? 
  int coden_len; // Num of coden in the coden ctxobj pattern
  int klen;      // full length of kmer, 8..31
  int hclen;     // half context length,
  int holen;     // half outer object length, 1..64
	  int compat_filter_shift;    // compatibility slot for the old hash-filter shift; minco writes 0
  int infile_num;
} minco_sketch_stat_t;

#define MINCO_STAT_EXT_MAGIC 0x4d434f53u /* "MCOS" */
#define MINCO_STAT_EXT_VERSION 1u
#define MINCO_STAT_HASH_FUNCTION_SPLITMIX64 1u
#define MINCO_STAT_SKETCH_MODEL_CTX_BOTTOMK 1u
#define MINCO_STAT_SELECTION_BOTTOMK 1u
#define MINCO_STAT_SELECTION_DENSITY_THRESHOLD 2u
#define MINCO_STAT_SELECTION_READWISE_DENSITY 3u
#define MINCO_STAT_DENSITY_POLICY_NONE 0u
#define MINCO_STAT_DENSITY_POLICY_SINGLE_SAMPLE 1u
#define MINCO_STAT_DENSITY_POLICY_LARGEST_SAMPLE 2u
#define MINCO_STAT_DENSITY_POLICY_EXPLICIT_THRESHOLD 3u
#define MINCO_STAT_DENSITY_SAMPLE_ID_NONE UINT32_MAX
#define MINCO_STAT_DENSITY_FLAG_MIXED_HASH_BITS 0x01u
#define MINCO_STAT_FLAG_MINCO_HASH_BOTTOMK 0x01u
#define MINCO_STAT_FLAG_STREAM_BOTTOMK 0x02u
#define MINCO_STAT_FLAG_KEEP_SOURCE_FILTER 0x04u
#define MINCO_STAT_FLAG_SPARSE_CTX_HASH 0x08u

#ifndef MINCO_SEED
#define MINCO_SEED 0x9e3779b97f4a7c15ULL
#endif

typedef struct minco_stat_density_summary
{
  uint64_t min_threshold;
  uint64_t max_threshold;
  uint64_t universal_threshold;
  uint32_t min_sample_id;
  uint32_t max_sample_id;
  uint32_t universal_sample_id;
  uint32_t valid_sample_count;
  uint32_t hash_bits;
  uint32_t universal_policy;
  uint32_t flags;
} minco_stat_density_summary_t;

/* Authoritative minco metadata appended after the legacy-compatible stat
 * prefix and fixed-width sample-name table:
 *   minco_sketch_stat_t prefix
 *   char sample_names[infile_num][PATHLEN]
 *   minco_stat_ext_v1_t extension
 *
 * The prefix keeps current payload readers working during this remodel; the
 * extension carries minco-specific identity such as -S and selection mode.
 */
typedef struct minco_stat_ext_v1
{
  uint32_t magic;
  uint16_t version;
  uint16_t struct_size;
  uint32_t target_sketch_size;
  uint32_t feature_id;
  uint32_t sketch_id;
  uint32_t hash_function;
  uint32_t hash_bits;
  uint64_t hash_seed;
  uint32_t sketch_model;
  uint32_t selection_mode;
  uint32_t flags;
  uint32_t reserved;
  uint64_t density_threshold;
  uint64_t density_min_threshold;
  uint64_t density_max_threshold;
  uint64_t density_universal_threshold;
  uint32_t density_min_sample_id;
  uint32_t density_max_sample_id;
  uint32_t density_universal_sample_id;
  uint32_t density_valid_sample_count;
  uint32_t density_hash_bits;
  uint32_t density_universal_policy;
  uint32_t density_flags;
  uint32_t density_reserved;
} minco_stat_ext_v1_t;

typedef struct minco_sketch_info
{
  bool has_minco_ext;
  uint32_t stat_version;
  uint32_t target_sketch_size;
  uint32_t feature_id;
  uint32_t sketch_id;
  uint32_t hash_function;
  uint32_t hash_bits;
  uint64_t hash_seed;
  uint32_t sketch_model;
  uint32_t selection_mode;
  uint32_t flags;
  uint64_t density_threshold;
  uint64_t density_min_threshold;
  uint64_t density_max_threshold;
  uint64_t density_universal_threshold;
  uint32_t density_min_sample_id;
  uint32_t density_max_sample_id;
  uint32_t density_universal_sample_id;
  uint32_t density_valid_sample_count;
  uint32_t density_hash_bits;
  uint32_t density_universal_policy;
  uint32_t density_flags;
} minco_sketch_info_t;

#define MINCO_PAYLOAD_CTXOBJ64 1u
#define MINCO_PAYLOAD_CTXOBJ96 2u

#ifndef MINCO_CTXOBJ_EXT_TYPES
#define MINCO_CTXOBJ_EXT_TYPES
typedef struct __attribute__((packed)) { uint64_t ctx; uint32_t obj; } ctxobj96_t;
typedef struct { uint64_t ctx; uint32_t gid; uint32_t obj; } ctxgidobj128_t;
_Static_assert(sizeof(ctxobj96_t) == 12, "ctxobj96_t must be a 12-byte record");
_Static_assert(sizeof(ctxgidobj128_t) == 16, "ctxgidobj128_t must be a 16-byte record");
#endif

typedef struct minco_ctxmeta_record
{
  uint8_t mode;
  uint8_t valid;
  uint16_t reserved16;
  uint32_t hash_bits;
  uint64_t threshold;
  uint64_t sketch_entries;
  uint64_t selected_observed_ctx;
  uint64_t selected_estimated_unique_ctx;
  uint64_t preconflict_observed_ctx;
  uint64_t preconflict_estimated_unique_ctx;
  uint64_t postconflict_observed_ctx;
  uint64_t postconflict_estimated_unique_ctx;
} minco_ctxmeta_record_t;
_Static_assert(sizeof(minco_ctxmeta_record_t) == 72,
               "minco_ctxmeta_record_t must remain a 72-byte v1 binary record");

#define DIM_SKETCH_QC_RANGE_VALID 1u
#define DIM_SKETCH_QC_RANGE_APPLIED 2u

/* Per-sample auxiliary sketch QC statistics. One record per sample, in the
 * same order as minco.stat names and minco.ctxobj64.offsets entries. */
typedef struct minco_sketch_qc_stat
{
  uint32_t flags;
  uint32_t reads_qc_lower;
  uint32_t reads_qc_upper;
  uint32_t reads_qc_mode;
} minco_sketch_qc_stat_t;

#define MINCO_INFILE_META_VERSION 1u
#define MINCO_INFILE_FMT_UNKNOWN 0
#define MINCO_INFILE_FMT_FASTA 1
#define MINCO_INFILE_FMT_FASTQ 2
#define MINCO_CREATE_NORMAL 0
#define MINCO_CREATE_ASONE 1
#define MINCO_CREATE_SPLITMFA 2
#define MINCO_CREATE_PIPECMD 3
#define MINCO_INFILE_FLAG_PIPECMD 0x01u
#define MINCO_INFILE_FLAG_STDIN 0x02u
#define MINCO_INFILE_FLAG_COMPRESSED 0x04u
#define MINCO_INFILE_FLAG_MEDIAN_APPROX 0x08u
#define MINCO_INFILE_FLAG_MIXED_FORMAT 0x10u

/* Per-sample input metadata. One record per sample, in the same order as
 * minco.stat names and minco.ctxobj64.offsets entries. meta_fmt_version == 0 means
 * metadata is unavailable/invalid for this sample.
 */
typedef struct infile_meta
{
  uint64_t total_length_bp;
  uint32_t record_count;
  uint32_t median_length_bp;
  float asm_level;
  float length_cv;
  uint8_t meta_fmt_version;
  int8_t infile_fmt;
  int8_t create_type;
  uint8_t infile_flags;
} infile_meta_t;

//********** input file formats test ********************//

/*input file formats array size */
#define ACPT_FMT_SZ 7
#define FAS_FMT_SZ 4
#define FQ_FMT_SZ 2
#define CO_FMT_SZ 1
#define MCO_FMT_SZ 1
#define CMPRESS_FMT_SZ 2

/* input file formats */
extern const char
    *acpt_infile_fmt[ACPT_FMT_SZ],
    *fasta_fmt[FAS_FMT_SZ],
    *fastq_fmt[FQ_FMT_SZ],
    *co_fmt[CO_FMT_SZ],
    *mco_fmt[MCO_FMT_SZ],
    *compress_fmt[CMPRESS_FMT_SZ];

/* format test fun.*/
#include <string.h>
static inline int isCompressfile(char *fname)
{
  int ret = 0;
  for (int i = 0; i < CMPRESS_FMT_SZ; i++)
  {
    int basename_len = strlen(fname) - strlen(compress_fmt[i]);
    if (strcmp((fname + basename_len), compress_fmt[i]) == 0)
      return 1;
  }
  return ret;
}

static inline int isOK_fmt_infile(char *fname, const char *test_fmt[], int test_fmt_arr_sz)
{
  int ret = 0;
  char suftmp[10];

  for (int i = 0; i < CMPRESS_FMT_SZ; i++)
  {
    int basename_len = strlen(fname) - strlen(compress_fmt[i]);
    // if infile suffix with .gz or other compress format
    if (strcmp((fname + basename_len), compress_fmt[i]) == 0)
    {
      char cp_fname[PATHLEN];
      strcpy(cp_fname, fname);
      *(cp_fname + basename_len) = '\0';
      fname = cp_fname;
      break;
    };
  };

  for (int i = 0; i < test_fmt_arr_sz; i++)
  {
    sprintf(suftmp, ".%s", test_fmt[i]);
    if (strcmp((char *)(fname + strlen(fname) - strlen(suftmp)), suftmp) == 0)
    {
      return 1;
    }
  };
  return ret;
};


// char check functions
static inline int is_dangerous(char c) {
    switch (c) {
        case '\0': case '\n': case '\r': case '\t': case '\\':
        case '"':  case '\'': case '|': case '>': case '<':
        case '$':  case '&':  case ';': case '*': case '?':
        case '#':  case '!':  case '%': case '@':
            return 1;
        default:
            return 0;
    }
}

static inline int is_bad_for_tokenizing(char c) {
    return (c == ' '  || c == '\t' || c == '\n' ||
            c == '\r' || c == '|'  || c == ','  ||
            c == '"'  || c == '\\' || (unsigned char)c < 32);
}

static inline int is_allowed_char(char c) {
    if (c >= 'A' && c <= 'Z') return 1;
    if (c >= 'a' && c <= 'z') return 1;
    if (c >= '0' && c <= '9') return 1;
    if (c == '_' || c == '-' || c == '.') return 1;
    return 0;
}

static inline int is_special(char c) {
    return !is_allowed_char(c);
}


#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>  // open function
#include <unistd.h> //close function
#include <stdint.h>

static inline void check(int test, const char *message, ...)
{
  if (test)
  {
    va_list args;
    va_start(args, message);
    vfprintf(stderr, message, args);
    va_end(args);
    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);
  }
};

typedef struct mmpco // member type match writeco2file()
{
  size_t fsize;
  unsigned int *mmpco;
} mmp_uint_t;

static inline mmp_uint_t mmp_uint_arr(char *cofname)
{
  mmp_uint_t cofilemmp; // char* fpath; llong fsize
  int fd;
  struct stat s;
  fd = open(cofname, O_RDONLY);
  check(fd < 0, "open %s failed: %s", cofname, strerror(errno));
  fstat(fd, &s);
  cofilemmp.fsize = s.st_size;
  cofilemmp.mmpco = mmap(NULL, s.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  check(cofilemmp.mmpco == MAP_FAILED, "mmap %s failed: %s", cofname, strerror(errno));
  close(fd);
  return cofilemmp;
};

typedef struct mmpany
{
  size_t fsize;
  void *mmp;
} mmp_any_t;

static inline mmp_any_t mmp_any(char *fname)
{
  mmp_any_t filemmp;
  int fd;
  struct stat s;
  fd = open(fname, O_RDONLY);
  check(fd < 0, "open %s failed: %s", fname, strerror(errno));
  fstat(fd, &s);
  filemmp.fsize = s.st_size;
  filemmp.mmp = mmap(NULL, s.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  check(filemmp.mmp == MAP_FAILED, "mmap %s failed: %s", fname, strerror(errno));
  close(fd);
  return filemmp;
};

void replaceChar(char *str, char oldChar, char newChar);

int str_suffix_match(char *str, const char *suf);
const char *get_pathname(const char *fullpath, const char *suf);
char *test_get_fullpath(const char *parent_path, const char *dstat_f);
char *test_create_fullpath(const char *parent_path, const char *dstat_f);
char *sketch_existing_fullpath(const char *parent_path, const char *dstat_f);
uint32_t minco_stat_hash_bits_from_dim(const minco_sketch_stat_t *stat);
uint32_t minco_stat_feature_id_from_dim(const minco_sketch_stat_t *stat);
uint32_t minco_stat_sketch_id_from_dim(const minco_sketch_stat_t *stat,
                                       uint32_t target_sketch_size,
                                       uint32_t selection_mode,
                                       uint32_t flags);
minco_stat_ext_v1_t minco_stat_make_ext(const minco_sketch_stat_t *stat,
                                        uint32_t target_sketch_size,
                                        uint32_t selection_mode,
                                        uint32_t flags,
                                        uint64_t density_threshold);
minco_stat_ext_v1_t minco_stat_make_ext_with_density(const minco_sketch_stat_t *stat,
                                                     uint32_t target_sketch_size,
                                                     uint32_t selection_mode,
                                                     uint32_t flags,
                                                     uint64_t density_threshold,
                                                     const minco_stat_density_summary_t *density_summary);
bool minco_stat_decode_mem(const void *mem, size_t stat_size,
                           minco_sketch_stat_t *legacy_out,
                           minco_sketch_info_t *info_out);
size_t minco_stat_base_size(int infile_num);
size_t minco_stat_full_size(int infile_num);
char (*minco_stat_names_from_mem(void *mem, size_t stat_size))[PATHLEN];
const char (*minco_stat_const_names_from_mem(const void *mem, size_t stat_size))[PATHLEN];
void minco_stat_write_path(const char *path, const minco_sketch_stat_t *stat,
                           const char (*names)[PATHLEN],
                           uint32_t target_sketch_size,
                           uint32_t selection_mode,
                           uint32_t flags,
                           uint64_t density_threshold);
void minco_stat_write_path_with_density(const char *path, const minco_sketch_stat_t *stat,
                                        const char (*names)[PATHLEN],
                                        uint32_t target_sketch_size,
                                        uint32_t selection_mode,
                                        uint32_t flags,
                                        uint64_t density_threshold,
                                        const minco_stat_density_summary_t *density_summary);
char *format_string(const char *format, ...);
int file_exists_in_folder(const char *folder, const char *filename);
void *read_from_file(const char *file_path, size_t *file_size);
void *read_file_mode(const char *file_path, size_t *file_size, const char *mode);
char **read_lines_from_file(const char *file_path, int *line_count);
void write_to_file(const char *file_path, const void *data, size_t data_size);
void concat_and_write_to_file(const char *file_path, const void *block1, size_t size1, const void *block2, size_t size2);
void replace_special_chars_with_underscore(char *str);
// infile fmt count struct
typedef struct
{
  int fasta;
  int fastq;
  int co;
  int mco;
} infile_fmt_count_t;
// infile fmt count function
infile_fmt_count_t *infile_fmt_count(infile_tab_t *infile_tab);

// uint64_t sketch
extern const char sketch_stat[];
extern const char sketch_qc_stat[];
extern const char sketch_anno_stat[];
extern const char sketch_infile_meta_stat[];
extern const char minco_ctxmeta_bin_stat[];
extern const char minco_ctxsetmeta_legacy_tsv_stat[];
extern const char sketch_position_suffix[];
extern const char combined_sketch_suffix[];
extern const char combined_sketch96_suffix[];
extern const char combined_ab_suffix[];
extern const char idx_sketch_suffix[];
extern const char idx_sketch96_suffix[];
extern const char sorted_comb_ctx64gid32obj32[];
extern const char minco_pan_prefix[];
extern const char minco_uniq_pan_prefix[];
// legency uint32_t sketch
extern const char co_dstat[];
extern const char skch_prefix[];
extern const char idx_prefix[];
extern const char pan_prefix[];
extern const char uniq_pan_prefix[];

extern const char mco_dstat[];
extern const char mco_gids_prefix[];
extern const char mco_idx_prefix[];

typedef unsigned int ctx_obj_ct_t;

int mkdir_p(const char *path);
void free_all(void *first, ...);

#define H1(K, HASH_SZ) ((K) % (HASH_SZ))
#define H2(K, HASH_SZ) (1 + (K) % ((HASH_SZ) - 1))
#define HASH(K, I, HASH_SZ) ((H1(K, HASH_SZ) + I * H2(K, HASH_SZ)) % HASH_SZ)
#define LOG2(X) ((unsigned)(8 * sizeof(unsigned long long) - __builtin_clzll((X)) - 1))

// general vector type and methods
typedef struct
{
  void *data;          // Pointer to the array data
  size_t element_size; // Size of each element
  size_t size;         // Current number of elements
  size_t capacity;     // Allocated capacity
} Vector;

void vector_init(Vector *vec, size_t element_size);
void vector_free(Vector *vec);
void vector_push(Vector *vec, const void *element);
void *vector_get(Vector *vec, size_t index);
void vector_reserve(Vector *v, size_t new_capacity);

//u64vec vector 
typedef struct{uint64_t *a;size_t n, cap; uint64_t threshold; size_t prune_at;} u64vec;

static inline void v_init(u64vec *v, size_t cap)
{
    v->a = cap ? (uint64_t *)malloc(cap * sizeof(uint64_t)) : NULL;
    v->n = 0;
    v->cap = cap;
    v->threshold = UINT64_MAX;
    v->prune_at = 0;
}
static inline void v_free(u64vec *v)
{
    free(v->a);
    v->a = NULL;
    v->n = v->cap = 0;
    v->threshold = UINT64_MAX;
    v->prune_at = 0;
}
static inline void v_reserve(u64vec *v, size_t need)
{
    if (need > v->cap)
    {
        size_t nc = v->cap ? v->cap : 8192;
        while (nc < need)
            nc <<= 1;
        v->a = (uint64_t *)realloc(v->a, nc * sizeof(uint64_t));
        v->cap = nc;
    }
}
static inline void v_push(u64vec *v, uint64_t x)
{
    if (v->n == v->cap)
        v_reserve(v, v->cap ? (v->cap << 1) : 8192);
    v->a[v->n++] = x;
}

// union type for in-memory sketch views
typedef struct
{
  int stat_type;          // 1 = 32-bit sketch, 2 = 64-bit sketch
  void *mem_stat;         // Pointer to memory-mapped or read data
  char (*gname)[PATHLEN]; // Query names
  uint64_t *comb_sketch;  // Combined k-mer counts
  ctxobj96_t *comb_sketch96;
  uint64_t *positions;    // Optional positions aligned with comb_sketch
  uint64_t *sketch_index; // Combined k-mer index
  uint32_t *abundance;
  minco_sketch_qc_stat_t *sample_qc;
  infile_meta_t *infile_meta;
  char (*annotation)[PATHLEN];
  int infile_num;   // Number of input files
  int kmerlen;      // K-mer length
  uint32_t hash_id; // hash_id or shuf_id
  uint32_t payload_layout;
  size_t payload_item_size;
  bool conflict; //if keep conflict obj
  minco_sketch_info_t minco_info;
  union
  {
    co_dstat_t co_stat_val;         // 32-bit sketch statistics
    minco_sketch_stat_t minco_stat; // 64-bit sketch statistics
  } stats;                          // Union of the two possible types
} unify_sketch_t;

typedef struct
{
  uint64_t *comb_sketch; // Combined k-mer counts
  uint64_t *sketch_index;
  int infile_num;   // Number of input files
} simple_sketch_t;
// ctx--diff obj count type
typedef struct
{
  uint32_t ctx_ct;
  uint32_t diff_obj;
} co_distance_t;

typedef enum sketch_parse_flags
{
  SKETCH_PARSE_NONE = 0,
  SKETCH_PARSE_POSITIONS = 1u << 0,
  SKETCH_PARSE_ABUNDANCE = 1u << 1,
  SKETCH_PARSE_SAMPLE_QC = 1u << 2,
  SKETCH_PARSE_ANNOTATION = 1u << 3,
  SKETCH_PARSE_INFILE_META = 1u << 4,
  SKETCH_PARSE_SAMPLE_SIDECARS =
      SKETCH_PARSE_SAMPLE_QC | SKETCH_PARSE_ANNOTATION | SKETCH_PARSE_INFILE_META,
  SKETCH_PARSE_ALL =
      SKETCH_PARSE_POSITIONS | SKETCH_PARSE_ABUNDANCE | SKETCH_PARSE_SAMPLE_SIDECARS,
  SKETCH_PARSE_ALL_EXCEPT_POSITIONS =
      SKETCH_PARSE_ABUNDANCE | SKETCH_PARSE_SAMPLE_SIDECARS
} sketch_parse_flags_t;

unify_sketch_t *generic_sketch_parse(const char *sketch_dir, unsigned flags);
void free_unify_sketch(unify_sketch_t *result);
void free_read_from_file(void *buffer, size_t file_size);

uint64_t GetAvailableMemory();
#endif
