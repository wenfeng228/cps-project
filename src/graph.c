#include "graph.h"
#include "log.h"

struct graph {
    size_t n; // nodes
    size_t m; // edges kept
    size_t* off; // len n+1
    int* to; // len m
    double* ww; // len m
};

static int cmp_edge_uv(const void* a, const void* b) {
    const Edge* ea = (const Edge*)a;
    const Edge* eb = (const Edge*)b;
    if (ea->u != eb->u) return (ea->u < eb->u) ? -1 : 1;
    if (ea->v != eb->v) return (ea->v < eb->v) ? -1 : 1;
    if (ea->w == eb->w) return 0;
    return (ea->w < eb->w) ? -1 : 1;
}

graph_t* graph_from_edges(size_t n, const Edge* edges, size_t m) {
    graph_t* g = (graph_t*)calloc(1, sizeof(*g));
    CHECK(g, "graph alloc failed");
    g->n = n;

    // copy & filter edges into a temp buffer
    Edge* buf = NULL;
    if (m) {
        buf = (Edge*)malloc(m * sizeof(Edge));
        CHECK(buf, "edge buf alloc failed");
    }
    size_t keep = 0;
    for (size_t i = 0; i < m; ++i) {
        int u = edges[i].u, v = edges[i].v;
        if (u < 0 || v < 0 || (size_t)u >= n || (size_t)v >= n) {
            log_warn("drop out-of-range edge (%d -> %d)", u, v);
            continue;
        }
        buf[keep++] = edges[i];
    }

    // sort by (u,v) to fill CSR
    qsort(buf, keep, sizeof(Edge), cmp_edge_uv);

    g->m  = keep;
    g->off = (size_t*)calloc(n + 1, sizeof(size_t));
    g->to  = (int*)   malloc(keep * sizeof(int));
    g->ww  = (double*)malloc(keep * sizeof(double));
    CHECK(g->off && (keep == 0 || (g->to && g->ww)), "CSR alloc failed");

    // degree count
    for (size_t i = 0; i < keep; ++i) g->off[buf[i].u]++;

    // prefix sum
    size_t acc = 0;
    for (size_t u = 0; u < n; ++u) {
        size_t du = g->off[u];
        g->off[u] = acc;
        acc += du;
    }
    g->off[n] = acc; // = m

    // fill neighbors, stable by sort
    size_t* cur = (size_t*)malloc((n + 1) * sizeof(size_t));
    CHECK(cur, "cur alloc failed");
    memcpy(cur, g->off, (n + 1) * sizeof(size_t));

    for (size_t i = 0; i < keep; ++i) {
        int u = buf[i].u;
        size_t k = cur[u]++;
        g->to[k] = buf[i].v;
        g->ww[k] = buf[i].w;
    }

    free(cur);
    free(buf);
    log_info("graph: n=%zu, m=%zu (kept)", g->n, g->m);
    return g;
}

void graph_free(graph_t* g) {
    if (!g) return;
    free(g->off);
    free(g->to);
    free(g->ww);
    free(g);
}

size_t graph_num_nodes(const graph_t* g) { return g ? g->n : 0; }
size_t graph_num_edges(const graph_t* g) { return g ? g->m : 0; }

size_t graph_out_begin(const graph_t* g, int u, size_t* deg) {
    if (!g || u < 0 || (size_t)u >= g->n) { if (deg) *deg = 0; return 0; }
    size_t k0 = g->off[u];
    size_t k1 = g->off[u + 1];
    if (deg) *deg = k1 - k0;
    return k0;
}
int graph_out_v(const graph_t* g, size_t k) { return g->to[k]; }
double graph_out_w(const graph_t* g, size_t k) { return g->ww[k]; }
