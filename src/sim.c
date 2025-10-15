#include "sim.h"
#include "io.h"
#include "graph.h"
#include "metrics.h"
#include "alloc_rule.h"
#include "dynamics.h"
#include "log.h"
#include "timer.h"
#include <sys/stat.h>
#include <math.h>

static int ensure_dir(const char* dir) {
#if defined(_WIN32)
    int rc = _mkdir(dir);
#else
    int rc = mkdir(dir, 0777);
#endif
    (void)rc; return 0;
}

// write per-step efficiency E_t (T rows)
static int write_efficiency(const char* out_path,
                            size_t T,
                            const double* eff_hist) {
    FILE* f = fopen(out_path, "w");
    if (!f) return -1;
    fprintf(f, "t,efficiency\n");
    for (size_t t = 0; t < T; ++t) {
        double e = eff_hist ? eff_hist[t] : 0.0;
        fprintf(f, "%zu,%.10g\n", t, e);
    }
    fclose(f);
    return 0;
}

static int write_trajectory(const char* out_path,
                            size_t T, size_t n,
                            const double* x_hist) {
    FILE* f = fopen(out_path, "w");
    if (!f) return -1;
    fprintf(f, "t,node,value\n");
    for (size_t t = 0; t <= T; ++t) {
        const double* xt = x_hist + t * n;
        for (size_t i = 0; i < n; ++i) {
            fprintf(f, "%zu,%zu,%.10g\n", t, i, xt[i]);
        }
    }
    fclose(f);
    return 0;
}

// write allocations history (T steps)
static int write_allocations(const char* out_path,
                             size_t T, size_t n,
                             const double* b_hist) {
    FILE* f = fopen(out_path, "w");
    if (!f) return -1;
    fprintf(f, "t,node,alloc\n");
    for (size_t t = 0; t < T; ++t) {
        const double* bt = b_hist + t * n;
        for (size_t i = 0; i < n; ++i) {
            fprintf(f, "%zu,%zu,%.10g\n", t, i, bt[i]);
        }
    }
    fclose(f);
    return 0;
}

// write per-step metrics (T+1 rows for t=0..T)
static int write_metrics(const char* out_path,
                         size_t T,
                         const double* sum_x_hist,
                         const double* dist_hist,
                         const double* rmse_hist, 
                         const double* gap_hist) { 
    FILE* f = fopen(out_path, "w");
    if (!f) return -1;
    fprintf(f, "t,sum_x,dist_to_target,rmse,gap_mean\n");
    for (size_t t = 0; t <= T; ++t) {
        double sx = sum_x_hist ? sum_x_hist[t] : 0.0;
        double d  = dist_hist   ? dist_hist[t]   : 0.0;
        double r   = rmse_hist   ? rmse_hist[t]  : 0.0;
        double gap = gap_hist    ? gap_hist[t]   : 0.0;
        fprintf(f, "%zu,%.10g,%.10g,%.10g,%.10g\n", t, sx, d, r, gap);
    }
    fclose(f);
    return 0;
}

static int write_vector_pair(const char* out_path,
                             size_t n,
                             const char* h1, const char* h2,
                             const double* a, const double* b) {
    FILE* f = fopen(out_path, "w");
    if (!f) return -1;
    fprintf(f, "node,%s,%s\n", h1, h2);
    for (size_t i = 0; i < n; ++i) {
        fprintf(f, "%zu,%.10g,%.10g\n", i, a ? a[i] : 0.0, b ? b[i] : 0.0);
    }
    fclose(f);
    return 0;
}

int sim_run(const config_t* cfg) {
    // 1) load inputs
    Edge* E = NULL; size_t M = 0; int max_node = -1;
    if (io_read_edges_csv(cfg->edges_file, &E, &M, &max_node) != 0) return -1;

    double* x0 = NULL; size_t N0 = 0;
    if (io_read_vector_csv(cfg->init_file, &x0, &N0) != 0) { free(E); return -1; }

    double* target = NULL; size_t Ny = 0;
    if (io_read_vector_csv(cfg->targets_file, &target, &Ny) != 0) { free(E); free(x0); return -1; }

    size_t n = N0 ? N0 : (size_t)(max_node + 1);
    if (n == 0) { log_err("sim: cannot infer n"); free(E); free(x0); free(target); return -1; }
    if (Ny && Ny < n) n = Ny;

    // 2) build graph
    graph_t* G = graph_from_edges(n, E, M);
    free(E);

    // 3) precompute metrics
    int* kout = (int*)malloc(n * sizeof(int));
    int* kin  = (int*)malloc(n * sizeof(int));
    double* s_out = (double*)malloc(n * sizeof(double));
    CHECK(kout && kin && s_out, "sim: metrics alloc failed");
    metrics_degrees(G, kout, kin);
    metrics_out_weight_sum(G, s_out);

    // 4) buffers
    size_t T = (size_t)cfg->T;
    double* x      = (double*)malloc(n * sizeof(double));
    double* x_next = (double*)malloc(n * sizeof(double));
    double* b      = (double*)malloc(n * sizeof(double));
    double* x_hist = (double*)malloc((T + 1) * n * sizeof(double));

    // 5) histories
    double* b_hist      = (double*)malloc(T * n * sizeof(double));
    double* sum_x_hist  = (double*)malloc((T + 1) * sizeof(double));
    double* dist_hist   = (double*)malloc((T + 1) * sizeof(double));

    CHECK(x && x_next && b && x_hist && b_hist && sum_x_hist && dist_hist,
          "sim: state alloc failed");

    double* rmse_hist = (double*)malloc((T + 1) * sizeof(double));
    double* gap_hist  = (double*)malloc((T + 1) * sizeof(double));
    double* eff_hist  = (double*)malloc(T * sizeof(double)); 
    double se_total   = 0.0; 

    CHECK(rmse_hist && gap_hist && eff_hist, "sim: extra metrics alloc failed");

    for (size_t i = 0; i < n; ++i) x[i] = (i < N0 ? x0[i] : 0.0);
    memcpy(x_hist + 0, x, n * sizeof(double));

    // 6) metrics at t=0
    double sx0 = 0.0, d0 = 0.0;
    for (size_t i = 0; i < n; ++i) {
        sx0 += x[i];
        if (target) {
            double diff = x[i] - target[i];
            d0 += diff * diff; 
        }
    }
    sum_x_hist[0] = sx0;
    dist_hist[0]  = target ? sqrt(d0) : 0.0;

    if (target) {
        rmse_hist[0] = dist_hist[0] / sqrt((double)n);
        double gap0 = 0.0;
        for (size_t i = 0; i < n; ++i) gap0 += (target[i] - x[i]);
        gap_hist[0] = gap0 / (double)n;
        se_total += d0; // d0 is ∑(x - T)^2 (no sqrt)
    } else {
        rmse_hist[0] = 0.0;
        gap_hist[0]  = 0.0;
    }

    // 7) loop
    scope_timer_t tt = timer_start("sim-steps");
    for (size_t t = 0; t < T; ++t) {
        // compute allocation and record it
        alloc_compute(cfg->alloc_rule, n, x, target, kout, cfg->budget, b);
        memcpy(b_hist + t * n, b, n * sizeof(double));

        // step update
        step_update(G, x, b, target, cfg->gamma, cfg->clamp01, x_next);

        memcpy(x, x_next, n * sizeof(double));
        memcpy(x_hist + (t + 1) * n, x, n * sizeof(double));

        // metrics at t+1
        double sx = 0.0, d = 0.0;
        for (size_t i = 0; i < n; ++i) {
            sx += x[i];
            if (target) { double diff = x[i] - target[i]; d += diff * diff; }
        }
        sum_x_hist[t + 1] = sx;
        dist_hist[t + 1]  = target ? sqrt(d) : 0.0;
        // per-step RMSE and mean gap at t+1
        if (target) {
            rmse_hist[t + 1] = dist_hist[t + 1] / sqrt((double)n);
            double gap = 0.0;
            for (size_t i = 0; i < n; ++i) gap += (target[i] - x[i]);
            gap_hist[t + 1] = gap / (double)n;
            se_total += d; // accumulate squared error for overall RMSE
        } else {
            rmse_hist[t + 1] = 0.0;
            gap_hist[t + 1]  = 0.0;
        }

        // efficiency Et = (sum_x[t+1] - sum_x[t]) / budget
        double B = cfg->budget;
        double sx_prev = sum_x_hist[t];
        eff_hist[t] = (B > 0.0) ? ((sx - sx_prev) / B) : 0.0;

    }
    double elapsed = timer_stop(&tt);
    log_info("sim: ran %zu steps in %.6f s", T, elapsed);

    // overall & steady RMSE
    double rmse_overall = 0.0, rmse_steady = 0.0;
    if (target && n > 0) {
        rmse_overall = sqrt(se_total / ((double)n * (double)(T + 1)));
        rmse_steady  = rmse_hist[T];
    }
    log_info("metrics: RMSE_overall=%.6g, RMSE_steady=%.6g",
            rmse_overall, rmse_steady);


    // 8) outputs
    ensure_dir(cfg->results_dir);

    char path[512];

    snprintf(path, sizeof(path), "%s/trajectory.csv", cfg->results_dir);
    if (write_trajectory(path, T, n, x_hist) == 0) log_info("sim: wrote %s", path);
    else log_warn("sim: cannot write %s", path);

    snprintf(path, sizeof(path), "%s/allocations.csv", cfg->results_dir);
    if (write_allocations(path, T, n, b_hist) == 0) log_info("sim: wrote %s", path);
    else log_warn("sim: cannot write %s", path);

    snprintf(path, sizeof(path), "%s/metrics.csv", cfg->results_dir);
    if (write_metrics(path, T, sum_x_hist, dist_hist, rmse_hist, gap_hist) == 0) log_info("sim: wrote %s", path);
    else log_warn("sim: cannot write %s", path);

    snprintf(path, sizeof(path), "%s/efficiency.csv", cfg->results_dir);
    if (write_efficiency(path, T, eff_hist) == 0)
        log_info("sim: wrote %s", path);
    else
        log_warn("sim: cannot write %s", path);

    snprintf(path, sizeof(path), "%s/final_vs_target.csv", cfg->results_dir);
    if (write_vector_pair(path, n, "final", "target", x, target) == 0) log_info("sim: wrote %s", path);
    else log_warn("sim: cannot write %s", path);

    snprintf(path, sizeof(path), "%s/degrees.csv", cfg->results_dir);
    double* kout_d = (double*)malloc(n * sizeof(double));
    double* kin_d  = (double*)malloc(n * sizeof(double));
    if (kout_d && kin_d) {
        for (size_t i = 0; i < n; ++i) { kout_d[i] = (double)kout[i]; kin_d[i] = (double)kin[i]; }
        if (write_vector_pair(path, n, "kout", "kin", kout_d, kin_d) == 0) log_info("sim: wrote %s", path);
        else log_warn("sim: cannot write %s", path);
    }
    free(kout_d); free(kin_d);

    // 9) cleanup
    free(kout); free(kin); free(s_out);
    free(x); free(x_next); free(b); free(x_hist);
    free(b_hist); free(sum_x_hist); free(dist_hist);

    free(rmse_hist);
    free(gap_hist);
    free(eff_hist);

    graph_free(G);
    free(x0); free(target);

    return 0;
}
