#ifndef PPI_LOG_H
#define PPI_LOG_H

#include "common.h"
#include <stdarg.h>

static inline int* _ppi_log_verbosity_slot(void) {
    static int v = 2;
    return &v;
}
static inline void log_set_verbosity(int v) { *_ppi_log_verbosity_slot() = v; }
static inline int  log_get_verbosity(void)  { return *_ppi_log_verbosity_slot(); }

#define _CLR_RESET  "\x1b[0m"
#define _CLR_RED    "\x1b[31m"
#define _CLR_YEL    "\x1b[33m"
#define _CLR_CYN    "\x1b[36m"
#define _CLR_MAG    "\x1b[35m"

static inline void _log_print(int level_gate,
                              const char* color,
                              const char* tag,
                              const char* fmt, va_list ap) {
    if (log_get_verbosity() < level_gate) return;
    fprintf(stderr, "%s[%s] %s", color, tag, _CLR_RESET);
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
}

static inline void log_err(const char* fmt, ...) {
    va_list ap; va_start(ap, fmt);
    _log_print(0, _CLR_RED, "ERR", fmt, ap);
    va_end(ap);
}
static inline void log_warn(const char* fmt, ...) {
    va_list ap; va_start(ap, fmt);
    _log_print(1, _CLR_YEL, "WRN", fmt, ap);
    va_end(ap);
}
static inline void log_info(const char* fmt, ...) {
    va_list ap; va_start(ap, fmt);
    _log_print(2, _CLR_CYN, "INF", fmt, ap);
    va_end(ap);
}
static inline void log_debug(const char* fmt, ...) {
    va_list ap; va_start(ap, fmt);
    _log_print(3, _CLR_MAG, "DBG", fmt, ap);
    va_end(ap);
}

#endif // PPI_LOG_H
