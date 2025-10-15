#ifndef PPI_ALLOC_RULE_H
#define PPI_ALLOC_RULE_H

#include "common.h"
#include "graph.h"

int alloc_compute(const char* rule, size_t n,
                  const double* x, const double* target,
                  const int* kout, double budget,
                  double* b);

#endif // PPI_ALLOC_RULE_H
