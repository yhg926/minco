#ifndef COMMAND_PROGRESS_H
#define COMMAND_PROGRESS_H

#include <stdbool.h>
#include <stdint.h>
#include <time.h>

typedef struct minco_progress
{
    bool enabled;
    const char *prefix;
    const char *label;
    const char *unit;
    uint64_t total;
    time_t started_at;
    time_t last_at;
} minco_progress_t;

minco_progress_t minco_progress_start(bool enabled,
                                    const char *prefix,
                                    const char *label,
                                    const char *unit,
                                    uint64_t total);
void minco_progress_update(minco_progress_t *progress, uint64_t done, bool force);
void minco_progress_done(minco_progress_t *progress);
void minco_progress_cb(void *ctx, uint64_t done, uint64_t total);

#endif
