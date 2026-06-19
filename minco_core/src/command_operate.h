#ifndef COMMAND_OPERATE_H
#define COMMAND_OPERATE_H

#include "command_set.h"

//extern const char sketch_stat[];
//extern const char idx_sketch_suffix[];
void print_minco_sample_names(set_opt_t* set_opt);
void show_content(set_opt_t* set_opt);
int minco_sketch_union(set_opt_t* set_opt);
int minco_sketch_operate(set_opt_t* set_opt);
int minco_sketch_downsample(set_opt_t* set_opt);
int minco_group_samples(set_opt_t* set_opt);
#endif
