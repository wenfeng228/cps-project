#ifndef PPI_TIMER_H
#define PPI_TIMER_H

#include "common.h"
#include <time.h>
#include <sys/time.h>

// Return a monotonic wall time in seconds
static inline double now_sec(void) {
#if defined(CLOCK_MONOTONIC)
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0) {
        return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
    }
#endif
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
}

static inline u64 now_us(void) {
    return (u64)(now_sec() * 1e6);
}

// Simple scope timer helper
typedef struct {
    const char* label;
    double t0;
} scope_timer_t;

static inline scope_timer_t timer_start(const char* label) {
    scope_timer_t t = { label, now_sec() };
    return t;
}
static inline double timer_stop(scope_timer_t* t) {
    return now_sec() - t->t0;
}

#endif // PPI_TIMER_H
