#include "alloc_rule.h"
#include "log.h"

static void allocate_uniform(size_t n, double budget, double* b) {
    double per = (n > 0) ? (budget / (double)n) : 0.0;
    for (size_t i = 0; i < n; ++i) b[i] = per;
}

int alloc_compute(const char* rule, size_t n,
                  const double* x, const double* target,
                  const int* kout, double budget,
                  double* b) {
    (void)x;
    if (!rule) rule = "uniform";

    if (n == 0) return 0;
    if (budget < 0.0) budget = 0.0;

    if (strcmp(rule, "uniform") == 0) {
        allocate_uniform(n, budget, b);
        return 0;
    }

    if (strcmp(rule, "degree") == 0) {
        double sum = 0.0;
        if (!kout) { allocate_uniform(n, budget, b); return 0; }
        for (size_t i = 0; i < n; ++i) sum += (double)(kout[i] > 0 ? kout[i] : 0);
        if (sum <= 0.0) { allocate_uniform(n, budget, b); return 0; }
        for (size_t i = 0; i < n; ++i) {
            double w = (double)(kout[i] > 0 ? kout[i] : 0);
            b[i] = budget * (w / sum);
        }
        return 0;
    }

    if (strcmp(rule, "gap") == 0) {
        if (!target || !x) { allocate_uniform(n, budget, b); return 0; }
        double sum = 0.0;
        for (size_t i = 0; i < n; ++i) {
            double g = target[i] - x[i];
            if (g < 0.0) g = 0.0;
            double k = (kout ? (double)kout[i] : 0.0);
            b[i] = g * (k + 1.0);   // q_i = (T_i - x_i) * (K_out_i + 1)
            sum += b[i];
        }
        if (sum <= 0.0) { allocate_uniform(n, budget, b); return 0; }
        for (size_t i = 0; i < n; ++i) b[i] = budget * (b[i] / sum); // P_i
        return 0;
    }

    // unknown rule -> uniform
    log_warn("alloc: unknown rule '%s', fallback to uniform", rule);
    allocate_uniform(n, budget, b);
    return 0;
}
