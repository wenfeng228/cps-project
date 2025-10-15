#include "io.h"
#include "log.h"
#include <ctype.h>

static void trim(char* s) {
    if (!s) return;
    size_t n = strlen(s);
    size_t i = 0, j = n ? n - 1 : 0;
    while (i < n && isspace((unsigned char)s[i])) i++;
    while (j + 1 > i && isspace((unsigned char)s[j])) j--;
    if (i > 0) memmove(s, s + i, j + 1 - i + 1);
    s[j + 1 - i] = '\0';
}

static int is_blank_or_comment(const char* s) {
    while (*s && isspace((unsigned char)*s)) ++s;
    return (!*s) || *s == '#' || *s == ';';
}

static int starts_with_number(const char* s) {
    while (*s && isspace((unsigned char)*s)) ++s;
    if (*s == '+' || *s == '-') ++s;
    return isdigit((unsigned char)*s) ? 1 : 0;
}

int io_read_edges_csv(const char* path, Edge** out_edges, size_t* out_m, int* out_max_node) {
    *out_edges = NULL; *out_m = 0; if (out_max_node) *out_max_node = -1;

    FILE* f = fopen(path, "r");
    if (!f) { log_err("io: cannot open edges: %s", path); return -1; }

    size_t cap = 1024, m = 0;
    Edge* arr = (Edge*)malloc(cap * sizeof(Edge));
    CHECK(arr, "edges alloc");

    char line[1024];
    size_t bad = 0, header_skipped = 0;

    while (fgets(line, sizeof(line), f)) {
        trim(line);
        if (is_blank_or_comment(line)) continue;

        for (char* q = line; *q; ++q) if (*q == ';') *q = ',';

        // header detection
        if (!starts_with_number(line)) {
            header_skipped++;
            continue;
        }

        int u, v; double w;
        int nmatch = sscanf(line, " %d , %d , %lf ", &u, &v, &w);
        if (nmatch < 2) { bad++; continue; }
        if (nmatch == 2) w = 1.0;

        if (m == cap) {
            cap *= 2;
            arr = (Edge*)realloc(arr, cap * sizeof(Edge));
            CHECK(arr, "edges realloc");
        }
        arr[m++] = (Edge){ .u = u, .v = v, .w = w };
        if (out_max_node) {
            if (u > *out_max_node) *out_max_node = u;
            if (v > *out_max_node) *out_max_node = v;
        }
    }
    fclose(f);

    *out_edges = arr;
    *out_m = m;

    if (header_skipped) log_debug("io: edges header lines skipped: %zu", header_skipped);
    if (bad) log_warn("io: edges bad lines skipped: %zu", bad);
    log_info("io: read %zu edges from %s", m, path);
    return 0;
}

int io_read_vector_csv(const char* path, double** out_x, size_t* out_n) {
    *out_x = NULL; *out_n = 0;

    FILE* f = fopen(path, "r");
    if (!f) { log_err("io: cannot open vector: %s", path); return -1; }

    size_t cap = 256, n = 0;
    double* x = (double*)malloc(cap * sizeof(double));
    CHECK(x, "vector alloc");

    char line[1024];
    size_t bad = 0, header_skipped = 0;

    while (fgets(line, sizeof(line), f)) {
        trim(line);
        if (is_blank_or_comment(line)) continue;

        for (char* q = line; *q; ++q) if (*q == ';') *q = ',';

        double val; int idx_dummy;
        int parsed = 0;

        if (strchr(line, ',')) {
            if (starts_with_number(line) &&
                sscanf(line, " %d , %lf ", &idx_dummy, &val) == 2) {
                parsed = 1;
            } else if (!starts_with_number(line)) {
                header_skipped++;
                continue;
            }
            if (!parsed) {
                bad++;
                continue;
            }
        } else {
            if (starts_with_number(line) &&
                sscanf(line, " %lf ", &val) == 1) {
                parsed = 1;
            } else if (!starts_with_number(line)) {
                header_skipped++;
                continue;
            } else {
                bad++;
                continue;
            }
        }

        if (n == cap) {
            cap *= 2;
            x = (double*)realloc(x, cap * sizeof(double));
            CHECK(x, "vector realloc");
        }
        x[n++] = val;
    }
    fclose(f);

    *out_x = x;
    *out_n = n;

    if (header_skipped) log_debug("io: vector header lines skipped: %zu", header_skipped);
    if (bad) log_warn("io: vector bad lines skipped: %zu", bad);
    log_info("io: read %zu values from %s", n, path);
    return 0;
}
