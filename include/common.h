#ifndef PPI_COMMON_H
#define PPI_COMMON_H

// C11 headers
#include <assert.h>
#include <errno.h>
#include <float.h>
#include <inttypes.h>
#include <limits.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h> 

// Short aliases
typedef uint64_t u64;
typedef uint32_t u32;
typedef uint16_t u16;
typedef uint8_t  u8;
typedef int64_t  i64;
typedef int32_t  i32;

// Generic helpers
#define ARRAY_LEN(a)    (sizeof(a) / sizeof((a)[0]))
#define MIN(a,b)        ((a) < (b) ? (a) : (b))
#define MAX(a,b)        ((a) > (b) ? (a) : (b))
#define CLAMP(x,l,h)    (MIN((h), MAX((l), (x))))
#define UNUSED(x)       ((void)(x))

// Error handling - abort on fatal
static inline void _ppi_fatalf(const char* tag,
                               const char* file, int line,
                               const char* fmt, ...) {
    fprintf(stderr, "[%s] %s:%d: ", tag, file, line);
    va_list ap; va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
    abort();
}

// Usage:
//   DIEF("message %d", code);
//   CHECK(ptr != NULL, "alloc failed: %zu bytes", n);
#define DIEF(...)  _ppi_fatalf("FATAL", __FILE__, __LINE__, __VA_ARGS__)
#define CHECK(cond, ...) do { if (!(cond)) _ppi_fatalf("ERROR", __FILE__, __LINE__, __VA_ARGS__); } while (0)

#endif // PPI_COMMON_H
