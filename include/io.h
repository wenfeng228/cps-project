#ifndef PPI_IO_H
#define PPI_IO_H

#include "common.h"
#include "graph.h"

// Read edges "u,v,w" csv 
int io_read_edges_csv(const char* path, Edge** out_edges, size_t* out_m, int* out_max_node);

// Read vector csv: either "value" per line or "idx,value"
int io_read_vector_csv(const char* path, double** out_x, size_t* out_n);

#endif // PPI_IO_H
