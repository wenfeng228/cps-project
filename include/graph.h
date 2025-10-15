#ifndef PPI_GRAPH_H
#define PPI_GRAPH_H

#include "common.h"

// Edge triple
typedef struct {
    int     u;   // source (0..n-1)
    int     v;   // target (0..n-1)
    double  w;   // weight
} Edge;

// Opaque CSR graph
typedef struct graph graph_t;

// Build CSR from edge list
graph_t* graph_from_edges(size_t n, const Edge* edges, size_t m);

// Free
void graph_free(graph_t* g);

// Basic info
size_t graph_num_nodes(const graph_t* g);
size_t graph_num_edges(const graph_t* g);

// Out-neighbor iteration - CSR:
size_t graph_out_begin(const graph_t* g, int u, size_t* deg);
int graph_out_v(const graph_t* g, size_t k);
double graph_out_w(const graph_t* g, size_t k);

#endif // PPI_GRAPH_H
