#ifndef PPI_CONFIG_H
#define PPI_CONFIG_H

#include "common.h"

typedef struct {
    // simulation
    int     T;
    double  budget;
    char    alloc_rule[32]; // "gap", "degree", "uniform"
    double  gamma;
    int     clamp01;
    // paths
    char    results_dir[256];
    char    edges_file[256];
    char    init_file[256];
    char    targets_file[256];
    // misc
    int     signed_edges;
    int     verbosity;
    // internal
    char    ini_path[256];
} config_t;

void  config_set_defaults(config_t* c);
int   config_load_ini(config_t* c, const char* path);
int   config_parse_argv(config_t* c, int argc, char** argv);
int   config_validate(const config_t* c, char* why, size_t why_cap);

void  config_print(const config_t* c, FILE* out);

#endif // PPI_CONFIG_H
