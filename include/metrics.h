#ifndef PPI_METRICS_H
#define PPI_METRICS_H

#include "common.h"
#include "graph.h"

// Compute out-degree and in-degree as counts
// kout/kin are length n
void metrics_degrees(const graph_t* g, int* kout, int* kin);

// Sum of outgoing weights per node
void metrics_out_weight_sum(const graph_t* g, double* s_out);

#endif // PPI_METRICS_H
