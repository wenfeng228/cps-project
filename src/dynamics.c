#include "dynamics.h"
#include "log.h"

int step_update(const graph_t* g,
                const double* x, // I_t
                const double* b, // P_t
                const double* target, // T
                double gamma, // γ
                int clamp01,            
                double* x_next) { // I_{t+1}
    size_t n = graph_num_nodes(g);

    // 1) inflow[i] = Σ_j a_{ji} * P_j
    for (size_t i = 0; i < n; ++i) x_next[i] = 0.0;
    for (size_t u = 0; u < n; ++u) {
        size_t deg = 0, k0 = graph_out_begin(g, (int)u, &deg);
        double Pu = (b ? b[u] : 0.0);
        if (Pu == 0.0) continue;
        for (size_t t = 0; t < deg; ++t) {
            int v = graph_out_v(g, k0 + t);
            double w = graph_out_w(g, k0 + t);
            x_next[v] += w * Pu;
        }
    }

    // 2) I_{i,t+1} = I_{i,t} + γ * (T_i - I_{i,t}) * ( P_i + inflow_i )
    for (size_t i = 0; i < n; ++i) {
        double Ii = x[i];
        double gap = target ? (target[i] - Ii) : 1.0;
        if (gap < 0.0) gap = 0.0;
        double eff = (b ? b[i] : 0.0) + x_next[i];
        double val = Ii + gamma * gap * eff;
        // optional
        if (clamp01) { 
            if (val < 0.0) val = 0.0;
            if (val > 1.0) val = 1.0;
        }
        x_next[i] = val;
    }
    return 0;
}
