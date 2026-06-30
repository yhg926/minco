#ifndef COMMAND_SET_WRAPPER_H
#define COMMAND_SET_WRAPPER_H

#include "global_basic.h"
#include "command_dist.h"
#include <sys/stat.h>
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

typedef struct set_opt
{
  int operation;
  bool q2markerdb;
  bool q2markerdb_context;
  bool q2markerdb_context_pairwise;
  uint32_t markerdb_warn_threshold;
  uint32_t markerdb_ctx_min_xny;
  double markerdb_ctx_min_af;
  uint32_t markerdb_ctx_index_max_ctx_freq;
  uint32_t markerdb_ctx_index_min_votes;
  uint32_t markerdb_ctx_index_sample_step;
  uint8_t domain_profile;
  int p;
  int P;
	int show;
  uint32_t sketch_size;
  int num_remaining_args;
  char ** remaining_args;
  char insketchpath[PATHLEN];
  char pansketchpath[PATHLEN];
  char subsetf[PATHLEN]; //subset infile
  char outdir[PATHLEN];
} set_opt_t ;




int cmd_set(struct argp_state* state);
#endif 
