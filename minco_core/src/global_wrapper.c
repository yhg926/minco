/*use by global cmd only*/
#include "global_wrapper.h"
#include "global_basic.h"
#include "command_sketch_wrapper.h"
#include "command_composite.h"
#include "command_matrix.h"
#include "command_ani.h"
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <limits.h>
#include <argp.h>
#include <argz.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif


char *domain;
char *long_domain;
/***===========global wraper=========== ***/

static struct argp_option opt_global[] = {
	{"usage",'u',0,OPTION_NO_USAGE,0},
	{"help",'?',0,OPTION_NO_USAGE,0},
	{"license",'l',0, OPTION_NO_USAGE,"license and copyright information.",0},
	{ 0 }
};

static char doc_license[] = 
"\n"
		"  Copyright 2019 Huiguang Yi. All Rights Reservered.\n\n"

  "Licensed under the Apache License, Version 2.0 (the \"License\");\n"
  "you may not use this file except in compliance with the License.\n"
  "You may obtain a copy of the License at\n\n"

  "	http://www.apache.org/licenses/LICENSE-2.0\n\n"

  "Unless required by applicable law or agreed to in writing, software\n"
  "distributed under the License is distributed on an \"AS IS\" BASIS,\n"
  "WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.\n"
  "See the License for the specific language governing permissions and\n"
  "limitations under the License.\n"
"\v"
;

static char doc_global[] =
"\n"
      "minco: pure C fixed-size context-minhash sketches for ANI analysis."
"\v"
      "The public sketch keeps the bottom 10,000 context hashes per sample by default; use --sketch-size to change it.\n"
      "Use `minco examples` for common workflows and `<subcommand> --help` for options.\n"
      "\n"
      "Public subcommands:\n"
"\n"
      "  sketch   \tCreate, inspect, index, filter, append, and deduplicate sketches.\n"
"\n"
      "  ani      \tEstimate ANI from sketches or direct FASTA/FASTQ inputs.\n"
"\n"
      "  matrix   \tReport context-distance matrices, sparse edges, clusters, and dedup plans.\n"

"\n"
      "  examples\tPrint common command workflows.\n"

"\n"
      "  doctor  \tCheck build/runtime environment basics.\n"

"\n"
;

static void print_global_examples(const char *prog)
{
  if (prog == NULL || prog[0] == '\0')
    prog = "minco";

  printf("\nCommon minco workflows\n\n");
  printf("Assembled genomes, reference/query ANI:\n");
  printf("  %s sketch -p8 --ctxmeta both -o ref_sketches refs/*.fasta.gz\n", prog);
  printf("  %s sketch -p8 --ctxmeta both -o qry_sketches queries/*.fasta.gz\n", prog);
  printf("  %s sketch -i ref_sketches\n", prog);
  printf("  %s ani -r ref_sketches -q qry_sketches -m0 -p8 -o ani.tsv\n\n", prog);

  printf("Many assembled genomes, fast all-vs-all context distance:\n");
  printf("  %s sketch -p8 --sketch-size 10000 -l genomes.list -o genomes.minco\n", prog);
  printf("  %s sketch -i genomes.minco\n", prog);
  printf("  %s matrix --format triangle -q genomes.minco -d -p8 -o dist.triangle.tsv\n\n", prog);

  printf("Many assembled genomes, all-vs-all ANI triangle:\n");
  printf("  %s ani -q genomes.minco -m2 -s -1 -d -p8 -o ani.triangle.tsv\n\n", prog);

  printf("Direct ANI from FASTA/FASTQ or sketch inputs:\n");
  printf("  %s ani -f0 -n0 -t0 -o pair_ani.tsv ref.fasta query1.fasta query2.fasta\n\n", prog);

  printf("Raw-read query ANI against assemblies:\n");
  printf("  %s sketch -p8 -o ref_sketches refs/*.fasta\n", prog);
  printf("  %s sketch --conflict --readsQC -p8 -o read_sketches reads/*.fastq.gz\n", prog);
  printf("  %s sketch -i ref_sketches\n", prog);
  printf("  %s ani -r ref_sketches --qraw read_sketches -m0 -p8 -o read_vs_ref.tsv\n\n", prog);
  printf("  %s ani -r ref_sketches --qraw reads.fastq.gz --readsQC -m0 -p8 -o one_read_file.tsv\n\n", prog);

  printf("Reads-to-reads ANI:\n");
  printf("  %s sketch --conflict --readsQC -p8 -o ref_read_sketches ref_reads/*.fastq.gz\n", prog);
  printf("  %s sketch --conflict --readsQC -p8 -o qry_read_sketches query_reads/*.fastq.gz\n", prog);
  printf("  %s sketch -i ref_read_sketches\n", prog);
  printf("  %s ani -r ref_read_sketches --qraw qry_read_sketches -m0 -p8 -o reads_to_reads.tsv\n\n", prog);

  printf("Streamed or tool-converted input:\n");
  printf("  samtools fastq reads.bam | %s sketch --conflict --readsQC -o reads_sketch -\n", prog);
  printf("  %s ani --pipecmd 'samtools fastq {}' \\\n", prog);
  printf("    -r ref_sketches --qraw reads.bam --readsQC -o raw_read_ani.tsv\n\n");

  printf("Notes:\n");
  printf("  - Sketch directories contain comblco, lcofiles.stat, and optional sidecars.\n");
  printf("  - Build an index with `sketch -i` before large reference or all-vs-all runs.\n");
  printf("  - For noisy reads, use `--conflict --readsQC` before `ani --qraw`.\n");
  printf("  - The default public sketch is fixed-size bottom-k: 10,000 contexts/sample.\n");
  printf("  - Change sketch size with `sketch --sketch-size N` or `ani -S N` for direct inputs.\n");
  printf("  - Use %s <subcommand> --help for full options.\n", prog);
}

static int path_has_entry(const char *path, const char *entry)
{
  if (path == NULL || entry == NULL || entry[0] == '\0')
    return 0;

  size_t entry_len = strlen(entry);
  const char *cursor = path;
  while (cursor[0] != '\0')
  {
    const char *end = strchr(cursor, ':');
    size_t len = end ? (size_t)(end - cursor) : strlen(cursor);
    if (len == entry_len && strncmp(cursor, entry, len) == 0)
      return 1;
    if (end == NULL)
      break;
    cursor = end + 1;
  }
  return 0;
}

static int path_has_executable_dir(const char *path, const char *prog)
{
  const char *slash = prog ? strrchr(prog, '/') : NULL;
  if (slash == NULL)
    return 0;

  char dir[PATH_MAX];
  size_t dir_len = (size_t)(slash - prog);
  if (dir_len == 0)
  {
    snprintf(dir, sizeof(dir), "/");
  }
  else
  {
    if (dir_len >= sizeof(dir))
      return 0;
    memcpy(dir, prog, dir_len);
    dir[dir_len] = '\0';
  }

  if (path_has_entry(path, dir))
    return 1;

  char resolved_dir[PATH_MAX];
  if (realpath(dir, resolved_dir) != NULL && path_has_entry(path, resolved_dir))
    return 1;

  return 0;
}

static void run_global_doctor(const char *prog)
{
  if (prog == NULL || prog[0] == '\0')
    prog = "minco";

  printf("\nminco doctor\n\n");
  printf("Version: %s\n", argp_program_version ? argp_program_version : "unknown");
  printf("Invocation: %s\n", prog);

#ifdef _OPENMP
  printf("OpenMP: enabled");
  printf(" (_OPENMP=%d, max_threads=%d)", _OPENMP, omp_get_max_threads());
  printf("\n");
#else
  printf("OpenMP: not enabled; -p/--threads will not speed up this build\n");
#endif

  if (access(".", W_OK) == 0)
    printf("Current directory writable: yes\n");
  else
    printf("Current directory writable: no\n");

  const char *path = getenv("PATH");
  if (path == NULL || path[0] == '\0')
    printf("PATH: not set\n");
  else if (path_has_executable_dir(path, prog))
    printf("PATH: includes the executable directory\n");
  else if (strchr(prog, '/') == NULL)
    printf("PATH: minco was invoked by command name\n");
  else
    printf("PATH: run `make install_env` for the recommended export line\n");

  printf("\nTry `%s examples` for common workflows.\n", prog);
}

static error_t parse_global(int key, char* arg, struct argp_state* state)
{
 // struct arg_global* global = state->input;
  if(key == '?' || key=='u')
  	state->name = long_domain;
	else if (key == 'l'){
		printf("%s\n",doc_license);
		return EINVAL; //EINVAL will end the parse loop, so the key never be assigned ARGP_KEY_NO_ARGS 
	}
  else if(key==ARGP_KEY_ARG){
      assert( arg );
      state->name = domain;
      if(strcmp(arg, "examples") == 0 || strcmp(arg, "example") == 0) {
         print_global_examples(domain);
         exit(0);
      }
      else if(strcmp(arg, "doctor") == 0) {
         run_global_doctor(domain);
         exit(0);
      }
      else if(strcmp(arg, "help") == 0) {
         state->name = long_domain;
         argp_state_help(state, stdout, ARGP_HELP_STD_HELP);
         exit(0);
      }
      if(strcmp(arg, "sketch") == 0)
      {
        cmd_sketch(state); //cmd_align
      }
			else if(strcmp(arg, "composite") == 0)	
				cmd_composite(state);		
			else if(strcmp(arg, "matrix") == 0)
				cmd_matrix(state);
	    else if(strcmp(arg, "ani") == 0)
        cmd_ani(state);
			else if(strcmp(arg, "primer") == 0)
					for(int i = 8;i<52;i++ )
				 		printf("%llu\n",find_lgst_primer_2pow(i));
			else {
        argp_error(state, "%s is not a valid command", arg);
      }
  }
  else if(key == ARGP_KEY_NO_ARGS){
        state->name = long_domain;
				printf("\n%s\n\n",argp_program_version);
				printf("Type '%s --license' for license and copyright information.\n\n", domain);
				printf("Type '%s examples' for common workflows.\n\n", domain);
        argp_state_help(state,stdout,ARGP_HELP_SHORT_USAGE);
        printf("\n");
        argp_state_help(state,stdout,ARGP_HELP_POST_DOC);
				//printf("\v");
        return EINVAL;
  }
  else if(key == ARGP_KEY_INVALID){
      state->name = state->argv[0] = domain;
      return ARGP_ERR_UNKNOWN;
  }
  else return ARGP_ERR_UNKNOWN;
  return 0;
}

static struct argp argp =
{
  opt_global,
  parse_global,
  "[arguments ...]",
  doc_global,
  0,
  0,
  0
};

void cmd_global(int argc, char**argv)
{
  struct arg_global global = {  };
  argp_parse(&argp, argc, argv, ARGP_IN_ORDER, &argc, &global);
}
