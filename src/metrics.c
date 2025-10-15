#include "metrics.h"

void metrics_degrees(const graph_t* g, int* kout, int* kin) {
    size_t n = graph_num_nodes(g);
    size_t m = graph_num_edges(g);
    for (size_t i = 0; i < n; ++i) { kout[i] = 0; kin[i] = 0; }

    for (size_t u = 0; u < n; ++u) {
        size_t deg = 0, k0 = graph_out_begin(g, (int)u, &deg);
        kout[u] = (int)deg;
        for (size_t t = 0; t < deg; ++t) {
            int v = graph_out_v(g, k0 + t);
            (void)m; 
            kin[v] += 1;
        }
    }
}

void metrics_out_weight_sum(const graph_t* g, double* s_out) {
    size_t n = graph_num_nodes(g);
    for (size_t i = 0; i < n; ++i) s_out[i] = 0.0;

    for (size_t u = 0; u < n; ++u) {
        size_t deg = 0, k0 = graph_out_begin(g, (int)u, &deg);
        double acc = 0.0;
        for (size_t t = 0; t < deg; ++t) {
            acc += graph_out_w(g, k0 + t);
        }
        s_out[u] = acc;
    }
}
