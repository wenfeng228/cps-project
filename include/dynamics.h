#ifndef PPI_DYNAMICS_H
#define PPI_DYNAMICS_H

#include "common.h"
#include "graph.h"

int step_update(const graph_t* g,
                const double* x,
                const double* b,
                const double* target,
                double gamma,
                int clamp01, 
                double* x_next);


#endif // PPI_DYNAMICS_H
