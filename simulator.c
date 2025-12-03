/**
* @file simulator.c
* @brief A simulator of network using OpenMP for parallel implementation.
*
* The main funtion contains spillover computation and iteration computation.
* It supports two graph representations: adjacency list and compressed sparse row.
* It uses three allocation strategies: uniform, local spillover, and global spillover.
*
* @note
* - Compile: mpicc -fopenmp -fopenmp simulator.c -o simulator -lm
* - Usage: srun -n #rank -c #thread ./simulator [gamma] [lambda] [budget] [output_mode] [nodes_file] [edges_file] ...
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <omp.h>

#define MAX_LINE_LEN 4096   /**< maximum length for line buffer reading */
#define DIST_INF 1e300  /**< representation of infinity for distance calculation */

/**********
min-heap definitions
**********/

/**
* @struct HeapNode
* @brief represents a single element in the priority queue.
*/
typedef struct {
    int32_t node; /**< node id */
    double dist;  /**< the priority value: distance/cost */
} HeapNode;

/**
* @struct MinHeap
* @brief a binary min-heap structure for efficient priority queue operations for Dijkstra's algorithm during global spillover calculation.
*/
typedef struct {
    HeapNode *data; /**< array of heap nodes */
    int32_t size;   /**< current number of elements in the heap */
    int32_t capacity;   /**< maximum capacity of the heap */
} MinHeap;

/** @brief returns the parent index of a heap node. */
static inline int32_t heap_parent(int32_t i) { return (i - 1) / 2; }
/** @brief returns the left child index of a heap node. */
static inline int32_t heap_left(int32_t i)   { return 2*i + 1; }
/** @brief returns the right child index of a heap node. */
static inline int32_t heap_right(int32_t i)  { return 2*i + 2; }

/**
* @brief allocates and initializes a new MinHeap.
* @param capacity The maximum number of elements the heap can hold.
* @return pointer to the created MinHeap.
*/
MinHeap* heap_create(int32_t capacity) {
    MinHeap* h = (MinHeap*) malloc(sizeof(MinHeap));
    if (!h) { perror("malloc heap"); exit(1); }
    h->data = (HeapNode*) malloc(capacity * sizeof(HeapNode));
    if (!h->data) { perror("malloc heap data"); exit(1); }
    h->size = 0;
    h->capacity = capacity;
    return h;
}

/**
* @brief frees the memory associated with a MinHeap.
* @param h pointer to the heap to free.
*/
void heap_free(MinHeap* h) {
    if (!h) return;
    free(h->data);
    free(h);
}

/**
* @brief pushes a new node into the min-heap.
* @param h pointer to the heap.
* @param node the node identifier.
* @param dist the priority/distance value.
*/
void heap_push(MinHeap* h, int32_t node, double dist) {
    if (h->size >= h->capacity) return;
    int i = h->size++;
    while (i > 0) {
        int p = heap_parent(i);
        if (h->data[p].dist <= dist) break;
        h->data[i] = h->data[p];
        i = p;
    }
    h->data[i].node = node;
    h->data[i].dist = dist;
}

/**
* @brief pops the node with the smallest distance from the heap.
* @param h pointer to the heap.
* @return the HeapNode with the minimum distance.
*/
HeapNode heap_pop(MinHeap* h) {
    HeapNode ret = h->data[0];
    int last = --h->size;
    HeapNode tmp = h->data[last];
    int i = 0;
    int child = heap_left(i);
    while (child < h->size) {
        if (child + 1 < h->size && h->data[child + 1].dist < h->data[child].dist) child++;
        if (tmp.dist <= h->data[child].dist) break;
        h->data[i] = h->data[child];
        i = child;
        child = heap_left(i);
    }
    h->data[i] = tmp;
    return ret;
}

/**********
adjacency list
**********/

/**
* @struct EdgeNode
* @brief a node in the adjacency list linked list.
*/
typedef struct EdgeNode {
    int32_t dst;    /**< destination node id */
    double w;   /**< weighted value of the edge */
    double uw;  /**< unweighted value used for calculations */
    struct EdgeNode *next;  /**< pointer to the next edge in the list */
} EdgeNode;

/**
* @struct AdjGraph
* @brief graph representation using an adjacency list.
*/
typedef struct {
    int32_t n_nodes;    /**< num of nodes in the graph */
    int32_t n_edges;    /**< num of edges in the graph */
    EdgeNode **head;    /**< array of pointers to the head of the linked lists for each node */
} AdjGraph;

/**
* @brief creates an empty adjacency graph.
* @param n_nodes number of nodes.
* @return pointer to the new AdjGraph.
*/
AdjGraph* create_adj_graph(int32_t n_nodes) {
    AdjGraph *g = (AdjGraph*) malloc(sizeof(AdjGraph));
    g->n_nodes = n_nodes;
    g->n_edges = 0;
    g->head = (EdgeNode**) calloc(n_nodes, sizeof(EdgeNode*));
    return g;
}

/**
* @brief adds an edge to the adjacency graph.
* @param g pointer to the graph.
* @param src source node id.
* @param dst destination node id.
* @param w edge weight.
* @param uw unweighted value.
*/
void add_edge_adj(AdjGraph *g, int src, int dst, double w, double uw) {
    EdgeNode *e = (EdgeNode*) malloc(sizeof(EdgeNode));
    e->dst = dst;
    e->w = w;
    e->uw = uw;
    e->next = g->head[src];
    g->head[src] = e;
    g->n_edges++;
}

/**
* @brief frees all memory associated with an adjacency graph.
* @param g pointer to the graph.
*/
void free_adj(AdjGraph *g) {
    if (!g) return;
    for (int i = 0; i < g->n_nodes; i++) {
        EdgeNode *curr = g->head[i];
        while (curr) {
            EdgeNode *tmp = curr;
            curr = curr->next;
            free(tmp);
        }
    }
    free(g->head);
    free(g);
}

/**********
csr
**********/

/**
* @struct CSRGraph
* @brief graph representation using compressed sparse row/CSR.
*/
typedef struct {
    int32_t n_nodes;    /**< num of nodes */
    int32_t n_edges;    /**< num of edges */
    int32_t *row_ptr;   /**< row pointers (size n_nodes + 1) */
    int32_t *col_ind;   /**< column indices (size n_edges) */
    double *val;    /**< edge weights (size n_edges) */
    double *val_unsigned;   /**< unweighted (size n_edges) */
} CSRGraph;

/**
* @brief frees all memory associated with a csr graph.
* @param g pointer to the graph.
*/
void free_csr(CSRGraph *g) {
    if (!g) return;
    free(g->row_ptr);
    free(g->col_ind);
    free(g->val);
    free(g->val_unsigned);
    free(g);
}

/**********
node
**********/

/**
* @struct Node
* @brief represents a single entity in the simulation network.
*/
typedef struct {
    int32_t id; /**< unique node id */
    double initial; /**< initial value */
    double target;  /**< target value */
    
    // simulation state
    double indicator;   /**< current state indicator */
    double alloc;   /**< budget allocated to this node */
    double gap; /**< difference between target and indicator */
    
    // pre-computed factor
    double spill_weight; /**< weight based on uniform, local, or global strategy */
} Node;

/**
* @brief resets node states for a new simulation run.
* @param nodes array of nodes.
* @param n num of nodes.
*/
void reset_nodes(Node *nodes, int n) {
    for (int i = 0; i < n; i++) {
        nodes[i].indicator = nodes[i].initial;
        nodes[i].alloc = 0.0;
        nodes[i].gap = 0.0;
        nodes[i].spill_weight = 1.0; // default to uniform
    }
}

/**
* @struct RankPair
* @brief helper structure for sorting nodes by spillover weight.
*/
typedef struct {
    int index;  /**< original index of the node */
    double value;   /**< value to sort by spillover weight */
} RankPair;

/**
* @brief comparator function for qsort to sort RankPairs descending.
*/
int compare_rank(const void *a, const void *b) {
    const RankPair *ra = (const RankPair *)a;
    const RankPair *rb = (const RankPair *)b;
    // Descending order
    if (rb->value > ra->value) return 1;
    if (rb->value < ra->value) return -1;
    return 0;
}

/**********
loaders
**********/

/**
* @brief loads node data from a csv file.
* @param path path to the csv file.
* @param[out] out_nodes pointer to the array of nodes.
* @param[out] out_n pointer to the integer storing node count.
* @return 0 on success, -1 on failure.
*/
int load_nodes(const char* path, Node **out_nodes, int *out_n) {
    FILE* f = fopen(path, "r");
    if (!f) return -1;
    char line[MAX_LINE_LEN];
    int n = 0;
    if (!fgets(line, sizeof(line), f)) { fclose(f); return -1; }
    while (fgets(line, sizeof(line), f)) n++;
    rewind(f);
    fgets(line, sizeof(line), f); // skip header

    Node *nodes = (Node*) malloc(n * sizeof(Node));
    if (!nodes) { fclose(f); return -1; }

    int i = 0;
    while (fgets(line, sizeof(line), f) && i < n) {
        int id; double init, tgt;
        if (sscanf(line, "%d,%lf,%lf", &id, &init, &tgt) >= 3) {
            nodes[i].id = id;
            nodes[i].initial = init;
            nodes[i].target = tgt;
            i++;
        }
    }
    fclose(f);
    *out_nodes = nodes;
    *out_n = i;
    return 0;
}

/**
* @struct EdgeRaw
* @brief temporary structure for reading edges before graph construction.
*/
typedef struct { int src; int dst; double w; double uw; } EdgeRaw;

/**
* @brief loads graph data and builds both csr and adjacency list representations.
* @param path path to the edges csv file.
* @param n_nodes total num of nodes.
* @param[out] csr_f pointer to forward csv graph.
* @param[out] csr_b pointer to backward csv graph.
* @param[out] adj_f pointer to forward adjacency graph.
* @param[out] adj_b pointer to backward adjacency graph.
* @return 0 on success, -1 on failure.
*/
int load_graphs(const char* path, int n_nodes, 
                CSRGraph** csr_f, CSRGraph** csr_b,
                AdjGraph** adj_f, AdjGraph** adj_b) {
    FILE* f = fopen(path, "r");
    if (!f) return -1;
    char line[MAX_LINE_LEN];
    fgets(line, sizeof(line), f); // skip header

    // read all edges first
    int cap = 1024, count = 0;
    EdgeRaw *edges = (EdgeRaw*) malloc(sizeof(EdgeRaw) * cap);
    
    while (fgets(line, sizeof(line), f)) {
        if (count == cap) {
            cap *= 2;
            edges = (EdgeRaw*) realloc(edges, sizeof(EdgeRaw) * cap);
        }
        int s, d; double w, uw;
        if (sscanf(line, "%d,%d,%lf,%lf", &s, &d, &w, &uw) >= 3) {
            edges[count].src = s;
            edges[count].dst = d;
            edges[count].w = w;
            edges[count].uw = uw;
            count++;
        }
    }
    fclose(f);

    // build adj graphs
    AdjGraph *af = create_adj_graph(n_nodes);
    AdjGraph *ab = create_adj_graph(n_nodes);
    
    // build csr graphs
    CSRGraph *cf = (CSRGraph*) malloc(sizeof(CSRGraph));
    CSRGraph *cb = (CSRGraph*) malloc(sizeof(CSRGraph));
    cf->n_nodes = n_nodes; cf->n_edges = count;
    cb->n_nodes = n_nodes; cb->n_edges = count;
    
    cf->row_ptr = (int32_t*) calloc(n_nodes + 1, sizeof(int32_t));
    cf->col_ind = (int32_t*) malloc(count * sizeof(int32_t));
    cf->val = (double*) malloc(count * sizeof(double));
    cf->val_unsigned = (double*) malloc(count * sizeof(double));

    cb->row_ptr = (int32_t*) calloc(n_nodes + 1, sizeof(int32_t));
    cb->col_ind = (int32_t*) malloc(count * sizeof(int32_t));
    cb->val = (double*) malloc(count * sizeof(double));
    cb->val_unsigned = (double*) malloc(count * sizeof(double));

    int32_t *deg_out = (int32_t*) calloc(n_nodes, sizeof(int32_t));
    int32_t *deg_in  = (int32_t*) calloc(n_nodes, sizeof(int32_t));

    // populate data
    for (int i = 0; i < count; i++) {
        int s = edges[i].src;
        int d = edges[i].dst;
        if (s < 0 || s >= n_nodes || d < 0 || d >= n_nodes) continue;

        // adj list
        add_edge_adj(af, s, d, edges[i].w, edges[i].uw);
        add_edge_adj(ab, d, s, edges[i].w, edges[i].uw);

        // csr deg count
        deg_out[s]++;
        deg_in[d]++;
    }

    // csr row ptrs
    for (int i = 0; i < n_nodes; i++) {
        cf->row_ptr[i+1] = cf->row_ptr[i] + deg_out[i];
        cb->row_ptr[i+1] = cb->row_ptr[i] + deg_in[i];
    }

    // csr fill
    int32_t *cur_out = (int32_t*) calloc(n_nodes, sizeof(int32_t));
    int32_t *cur_in  = (int32_t*) calloc(n_nodes, sizeof(int32_t));
    for (int i = 0; i < n_nodes; i++) { cur_out[i] = cf->row_ptr[i]; cur_in[i] = cb->row_ptr[i]; }

    for (int i = 0; i < count; i++) {
        int s = edges[i].src; 
        int d = edges[i].dst;
        if (s < 0 || s >= n_nodes || d < 0 || d >= n_nodes) continue;
        
        int p_f = cur_out[s]++;
        cf->col_ind[p_f] = d;
        cf->val[p_f] = edges[i].w;
        cf->val_unsigned[p_f] = edges[i].uw;

        int p_b = cur_in[d]++;
        cb->col_ind[p_b] = s;
        cb->val[p_b] = edges[i].w;
        cb->val_unsigned[p_b] = edges[i].uw;
    }

    free(edges);
    free(deg_out); free(deg_in);
    free(cur_out); free(cur_in);
    
    *csr_f = cf; *csr_b = cb;
    *adj_f = af; *adj_b = ab;
    return 0;
}

/**********
local spillover computation (csr)
**********/

/**
* @brief computes local spillover weights using csr representation.
* @param g the csr graph.
* @param nodes the array of nodes to update.
* @param n num of nodes.
*/
void compute_local_spillover_csr(const CSRGraph *g, Node *nodes, int n) {
    #pragma omp parallel for schedule(dynamic, 4)
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = g->row_ptr[i]; j < g->row_ptr[i+1]; j++) {
            sum += g->val[j];
        }
        nodes[i].spill_weight = sum;
    }
}

/**
* @brief computes global spillover weights using CSR representation.
* @param g the csr graph.
* @param nodes the array of nodes to update.
* @param n num of nodes.
*/
void compute_global_spillover_csr(const CSRGraph *g, Node *nodes, int n) {
    #pragma omp parallel
    {
        double *dist = (double*) malloc(n * sizeof(double));
        MinHeap *pq = heap_create(g->n_edges + n);

        #pragma omp for schedule(static, 1)
        for (int src = 0; src < n; src++) {
            for (int v = 0; v < n; v++) dist[v] = DIST_INF;
            dist[src] = 0.0;
            pq->size = 0;
            heap_push(pq, src, 0.0);

            while (pq->size > 0) {
                HeapNode top = heap_pop(pq);
                int u = top.node;
                if (top.dist > dist[u]) continue;

                for (int idx = g->row_ptr[u]; idx < g->row_ptr[u+1]; idx++) {
                    int v = g->col_ind[idx];
                    double w = g->val_unsigned[idx];
                    if (dist[u] + w < dist[v]) {
                        dist[v] = dist[u] + w;
                        heap_push(pq, v, dist[v]);
                    }
                }
            }
            double sum_dist = 0.0;
            for (int v = 0; v < n; v++) {
                if (v != src && dist[v] < DIST_INF) sum_dist += dist[v];
            }
            nodes[src].spill_weight = sum_dist;
        }
        free(dist);
        heap_free(pq);
    }
}

/**********
local spillover computation (adj)
**********/

/**
* @brief computes local spillover weights using adj list representation.
* @param g the adjacency graph.
* @param nodes the array of nodes to update.
* @param n num of nodes.
*/
void compute_local_spillover_adj(const AdjGraph *g, Node *nodes, int n) {
    #pragma omp parallel for schedule(dynamic, 4)
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        EdgeNode *curr = g->head[i];
        while (curr) {
            sum += curr->w;
            curr = curr->next;
        }
        nodes[i].spill_weight = sum;
    }
}

/**
* @brief computes global spillover weights using adj list representation.
* @param g the adj graph.
* @param nodes the array of nodes to update.
* @param n num of nodes.
*/
void compute_global_spillover_adj(const AdjGraph *g, Node *nodes, int n) {
    #pragma omp parallel
    {
        double *dist = (double*) malloc(n * sizeof(double));
        MinHeap *pq = heap_create(g->n_edges + n);

        #pragma omp for schedule(static, 1)
        for (int src = 0; src < n; src++) {
            for (int v = 0; v < n; v++) dist[v] = DIST_INF;
            dist[src] = 0.0;
            pq->size = 0;
            heap_push(pq, src, 0.0);

            while (pq->size > 0) {
                HeapNode top = heap_pop(pq);
                int u = top.node;
                if (top.dist > dist[u]) continue;

                EdgeNode *curr = g->head[u];
                while (curr) {
                    int v = curr->dst;
                    double w = curr->uw;
                    if (dist[u] + w < dist[v]) {
                        dist[v] = dist[u] + w;
                        heap_push(pq, v, dist[v]);
                    }
                    curr = curr->next;
                }
            }
            double sum_dist = 0.0;
            for (int v = 0; v < n; v++) {
                if (v != src && dist[v] < DIST_INF) sum_dist += dist[v];
            }
            nodes[src].spill_weight = sum_dist;
        }
        free(dist);
        heap_free(pq);
    }
}

/**********
simulation computation (adj, csr)
**********/

/**
* @brief computes the l2/Euclidean dist of the gap.
* @param nodes array of nodes.
* @param n num of nodes.
* @return the square root of the sum of squared gaps.
*/
double compute_l2_gap(Node *nodes, int n) {
    double sum_sq = 0.0;
    #pragma omp parallel for reduction(+:sum_sq)
    for (int i = 0; i < n; i++) {
        double g = nodes[i].target - nodes[i].indicator;
        if (g < 0) g = 0;
        sum_sq += g * g;
    }
    return sqrt(sum_sq);
}

/**
* @brief executes the simulation loop for budget allocation.
* @param rep_type the graph representation (0 for AdjList, 1 for CSR).
* @param rep_name string name of representation for logging.
* @param fac_name string name of the spillover factor strategy.
* @param bwd_graph void pointer to the backward graph (cast based on rep_type).
* @param nodes array of Node structures containing state.
* @param n total number of nodes.
* @param budget total budget to distribute per iteration.
* @param lambda impact factor for spillover weights.
* @param gamma step size for indicator updates.
* @param max_iter max num of iterations allowed.
* @param epsilon convergence threshold or tolerance for l2 gap.
* @param output_mode if 1, writes detailed csv logs; if 0, runs silently.
* @param output_dir directory path for output files.
* @param data_val the dataset scale identifier.
* @param num_candidates num of top candidates to track in logs.
* @param num_snapshots num of snapshots to record in the log file.
* @param[out] out_steps pointer to store the final iteration count.
* @param[out] out_l2 pointer to store the final l2 error.
*/
void simulation_iteration(
    int rep_type, // 0=adj, 1=csr
    const char *rep_name,
    const char *fac_name,
    const void *bwd_graph, 
    Node *nodes, int n, 
    double budget, double lambda, double gamma,
    int max_iter, double epsilon,
    int output_mode,
    const char *output_dir,
    int data_val,
    int num_candidates,
    int num_snapshots,
    int *out_steps, double *out_l2
) {
    double *new_ind = (double*) malloc(n * sizeof(double));
    int step = 0;
    double l2 = 1e9;

    // file pointer for logging
    FILE *log_file = NULL;
    int *cand_idx = NULL;
    int actual_candidates = (n < num_candidates) ? n : num_candidates;
    
    // calculate snapshot interval
    int snapshot_interval = max_iter / num_snapshots;
    if (snapshot_interval < 1) snapshot_interval = 1;

    // only open file if output_mode is 1 as it impacts runtime measure
    if (output_mode == 1) {
        char filename[1024]; 
        snprintf(filename, sizeof(filename), "%s/%s_%s_trace_%d.csv", output_dir, rep_name, fac_name, data_val);

        log_file = fopen(filename, "w");
        if (log_file) { 
            cand_idx = (int*) malloc(sizeof(int) * actual_candidates);

            if (strcmp(fac_name, "uniform") == 0) {
                // Take IDs 0 to K-1
                for (int i = 0; i < actual_candidates; i++) {
                    cand_idx[i] = i; 
                }
            } else {
                // rank by spill_weight (descending)
                RankPair *ranks = (RankPair*) malloc(sizeof(RankPair) * n);
                for (int i = 0; i < n; i++) {
                    ranks[i].index = i;
                    ranks[i].value = nodes[i].spill_weight;
                }
                qsort(ranks, n, sizeof(RankPair), compare_rank);

                for (int i = 0; i < actual_candidates; i++) {
                    cand_idx[i] = ranks[i].index;
                }
                free(ranks);
            }

            // write Header for Candidates
            fprintf(log_file, "step,l2");
            for(int k = 0; k < actual_candidates; k++) {
                int nid = nodes[cand_idx[k]].id;
                // output indicator and alloc for each candidate
                fprintf(log_file, ",n%d_ind,n%d_alloc", nid, nid);
            }
            fprintf(log_file, "\n");
        } else {
            fprintf(stderr, "Error opening file: %s\n", filename);
        }
    }

    while (step < max_iter) {
        // 1. calc gap & factor
        double sum_factor = 0.0;
        #pragma omp parallel for reduction(+:sum_factor)
        for (int i = 0; i < n; i++) {
            double g = nodes[i].target - nodes[i].indicator;
            if (g < 0) g = 0;
            nodes[i].gap = g;
            double factor = g * (lambda * nodes[i].spill_weight + 1.0);
            nodes[i].alloc = factor; 
            sum_factor += factor;
        }

        // 2. distribute budget
        if (sum_factor > 1e-9) {
            #pragma omp parallel for
            for (int i = 0; i < n; i++) {
                nodes[i].alloc = budget * (nodes[i].alloc / sum_factor);
            }
        } else {
            #pragma omp parallel for
            for (int i = 0; i < n; i++) nodes[i].alloc = 0.0;
        }

        // 3. update indicators
        if (rep_type == 0) { // adj list
            const AdjGraph *bg = (const AdjGraph*) bwd_graph;
            #pragma omp parallel for schedule(dynamic, 4)
            for (int i = 0; i < n; i++) {
                double incoming = 0.0;
                EdgeNode *e = bg->head[i];
                while (e) {
                    incoming += e->w * nodes[e->dst].alloc;
                    e = e->next;
                }
                new_ind[i] = nodes[i].indicator + gamma * nodes[i].gap * (nodes[i].alloc + incoming);
            }
        } else { // csr
            const CSRGraph *bg = (const CSRGraph*) bwd_graph;
            #pragma omp parallel for schedule(dynamic, 4)
            for (int i = 0; i < n; i++) {
                double incoming = 0.0;
                for (int idx = bg->row_ptr[i]; idx < bg->row_ptr[i+1]; idx++) {
                    int j = bg->col_ind[idx];
                    incoming += bg->val[idx] * nodes[j].alloc;
                }
                new_ind[i] = nodes[i].indicator + gamma * nodes[i].gap * (nodes[i].alloc + incoming);
            }
        }

        #pragma omp parallel for
        for (int i = 0; i < n; i++) nodes[i].indicator = new_ind[i];

        l2 = compute_l2_gap(nodes, n);

        if (log_file && (step % snapshot_interval == 0)) {
            fprintf(log_file, "%d,%.6f", step, l2);
            for (int k = 0; k < actual_candidates; k++) {
                int idx = cand_idx[k];
                fprintf(log_file, ",%.6f,%.6f", nodes[idx].indicator, nodes[idx].alloc);
            }
            fprintf(log_file, "\n");
        }

        if (l2 < epsilon) break;
        step++;
    }        
    
    if (log_file) {
        // ensure final step is logged if it wasn't caught by the modulo
        if (step % snapshot_interval != 0) {
            fprintf(log_file, "%d,%.6f", step, l2);
            for (int k = 0; k < actual_candidates; k++) {
                int idx = cand_idx[k];
                fprintf(log_file, ",%.6f,%.6f", nodes[idx].indicator, nodes[idx].alloc);
            }
            fprintf(log_file, "\n");
        }
        fclose(log_file);
        free(cand_idx);
    }
    
    free(new_ind);
    *out_steps = step;
    *out_l2 = l2;
}

/**********
main
**********/
/**
* @brief entry point of the simulation.
* @param argc argument count.
* @param argv argument vector.
* @return 0 on success, 1 on error.
*/
int main(int argc, char** argv) {
    // defaults
    double gamma = 0.1;
    double lambda = 0.25;
    double budget = 100.0;
    int max_iter = 5000;
    int output_mode = 0; 

    const char *nodes_path = "../data/nodes.csv";
    const char *graph_path = "../data/edges.csv";
    const char *output_dir = "../output";

    int n_candidates = 20;
    int n_snapshots = 100;
    
    // Parse CLI args
    if (argc >= 6) {
        gamma = atof(argv[1]);
        lambda = atof(argv[2]);
        budget = atof(argv[3]);
        max_iter = atoi(argv[4]);
        output_mode = atoi(argv[5]);
    }

    if (argc >= 8) {
        nodes_path = argv[6];
        graph_path = argv[7];
    }

    if (argc >= 9) {
        output_dir = argv[8];
    }

    if (argc >= 10) n_candidates = atoi(argv[9]);
    if (argc >= 11) n_snapshots = atoi(argv[10]);

    double epsilon = 1e-3;

    int n_threads = omp_get_max_threads();

    // load data
    Node *nodes = NULL;
    int n = 0;
    if (load_nodes(nodes_path, &nodes, &n) != 0) {
        printf("error: failed to load nodes\n");
        return 1;
    }

    CSRGraph *csr_fwd = NULL, *csr_bwd = NULL;
    AdjGraph *adj_fwd = NULL, *adj_bwd = NULL;
    if (load_graphs(graph_path, n, &csr_fwd, &csr_bwd, &adj_fwd, &adj_bwd) != 0) {
        printf("error: failed to load graphs\n");
        return 1;
    }

    long long num_nodes = 0;
    long long num_edges = 0;
    int data_val = 0;
    // find "nodes_" inside the string and read integer
    const char *p_m = strstr(nodes_path, "nodes_");
    if (p_m) {
        if (sscanf(p_m, "nodes_%d.csv", &data_val) == 1) {
            num_nodes = 69LL * (1LL << data_val);
            num_edges = 215LL * (1LL << data_val);
        }
    }

    // experiment matrix
    // representation: 0 = adj list, 1 = csr
    // factor: 0 = uniform, 1 = local, 2 = global
    const char* rep_names[] = {"adjacency_list", "compressed_sparse_row"};
    const char* fac_names[] = {"uniform", "local_spillover", "global_spillover"};

    for (int r = 0; r < 2; r++) {     
        for (int f = 0; f < 3; f++) { 
            
            reset_nodes(nodes, n);

            double t_spill_start = omp_get_wtime();
            double t_spill_end = t_spill_start;

            // 1. compute allocation factor
            if (f == 0) {
                // uniform
            } 
            else if (f == 1) { // local
                if (r == 0) compute_local_spillover_adj(adj_fwd, nodes, n);
                else        compute_local_spillover_csr(csr_fwd, nodes, n);
                t_spill_end = omp_get_wtime();
            } 
            else if (f == 2) { // global
                if (r == 0) compute_global_spillover_adj(adj_fwd, nodes, n);
                else        compute_global_spillover_csr(csr_fwd, nodes, n);
                t_spill_end = omp_get_wtime();
            }

            double time_spillover = t_spill_end - t_spill_start;
            
            // 2. run iteration
            double t_iter_start = omp_get_wtime();
            int final_step = 0;
            double final_l2 = 0.0;
            
            void *bwd_ptr = (r == 0) ? (void*)adj_bwd : (void*)csr_bwd;
            
            simulation_iteration(r, rep_names[r], fac_names[f], bwd_ptr, nodes, n, 
                                budget, lambda, gamma, max_iter, epsilon, 
                                output_mode, output_dir,
                                data_val,
                                n_candidates, n_snapshots,
                                &final_step, &final_l2);

            double t_iter_end = omp_get_wtime();
            double time_iteration = t_iter_end - t_iter_start;
            double time_total = time_spillover + time_iteration;

            // 3. print output
            if (f == 0) {
                printf("%d, %s, %s, nan, %.6f, %.6f, %d, %.6f, %.4f, %.4f, %.4f, %d, %.6f, %d, %lld, %lld\n",
                    (r*3 + f + 1), rep_names[r], fac_names[f], 
                    time_iteration, time_total, final_step, final_l2,
                    gamma, lambda, budget, max_iter, epsilon, n_threads,
                    num_nodes, num_edges);
            } else {
                printf("%d, %s, %s, %.6f, %.6f, %.6f, %d, %.6f, %.4f, %.4f, %.4f, %d, %.6f, %d, %lld, %lld\n",
                    (r*3 + f + 1), rep_names[r], fac_names[f], 
                    time_spillover, time_iteration, time_total, final_step, final_l2,
                    gamma, lambda, budget, max_iter, epsilon, n_threads,
                    num_nodes, num_edges);
            }
        }
    }

    // free memory
    free(nodes);
    free_csr(csr_fwd); free_csr(csr_bwd);
    free_adj(adj_fwd); free_adj(adj_bwd);
    
    return 0;
}