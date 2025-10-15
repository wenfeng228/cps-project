#include "config.h"
#include "log.h"
#include <ctype.h>
#include <sys/stat.h>

// helpers
static int file_exists(const char* p) {
    struct stat st;
    return (p && p[0] && stat(p, &st) == 0 && S_ISREG(st.st_mode));
}
static int dir_exists_or_mkdir(const char* p) {
    struct stat st;
    if (!p || !p[0]) return 0;
    if (stat(p, &st) == 0) return S_ISDIR(st.st_mode);
#if defined(_WIN32)
    int rc = _mkdir(p);
#else
    int rc = mkdir(p, 0777);
#endif
    if (rc == 0) return 1;
    return 0;
}
static void trim(char* s) {
    if (!s) return;
    size_t n = strlen(s);
    size_t i = 0, j = n ? n - 1 : 0;
    while (i < n && isspace((unsigned char)s[i])) i++;
    while (j + 1 > i && isspace((unsigned char)s[j])) j--;
    if (i > 0) memmove(s, s + i, j + 1 - i + 1);
    s[j + 1 - i] = '\0';
}
static void set_str(char dst[], size_t cap, const char* src) {
    if (!src) { dst[0] = '\0'; return; }
    snprintf(dst, cap, "%s", src);
}

// defaults
void config_set_defaults(config_t* c) {
    memset(c, 0, sizeof(*c));
    c->T = 20;
    c->budget = 0.000001;
    c->gamma = 0.1;
    c->clamp01 = 1;
    set_str(c->alloc_rule, sizeof(c->alloc_rule), "gap");

    set_str(c->results_dir,  sizeof(c->results_dir), "results/run1");
    set_str(c->edges_file,   sizeof(c->edges_file),  "data/edges_signed.csv");
    set_str(c->init_file,    sizeof(c->init_file),   "data/init.csv");
    set_str(c->targets_file, sizeof(c->targets_file),"data/targets.csv");

    c->signed_edges = 1;
    c->verbosity = 2;
    c->ini_path[0] = '\0';
}

// INI loader
static void apply_kv(config_t* c, const char* k, const char* v) {
    if (strcmp(k, "T") == 0)               { c->T = (int)strtol(v, NULL, 10); return; }
    if (strcmp(k, "budget") == 0)          { c->budget = strtod(v, NULL); return; }
    if (strcmp(k, "alloc_rule") == 0)      { set_str(c->alloc_rule, sizeof(c->alloc_rule), v); return; }
    if (strcmp(k, "results_dir") == 0)     { set_str(c->results_dir, sizeof(c->results_dir), v); return; }
    if (strcmp(k, "edges_file") == 0)      { set_str(c->edges_file, sizeof(c->edges_file), v); return; }
    if (strcmp(k, "init_file") == 0)       { set_str(c->init_file, sizeof(c->init_file), v); return; }
    if (strcmp(k, "targets_file") == 0)    { set_str(c->targets_file, sizeof(c->targets_file), v); return; }
    if (strcmp(k, "signed_edges") == 0)    { c->signed_edges = (int)strtol(v, NULL, 10); return; }
    if (strcmp(k, "verbosity") == 0)       { c->verbosity = (int)strtol(v, NULL, 10); return; }
    if (strcmp(k, "gamma") == 0)           { c->gamma = strtod(v, NULL); return; }
    if (strcmp(k, "clamp01") == 0)         { c->clamp01 = (int)strtol(v, NULL, 10); return; }
    // unknown keys are ignored by design
}

int config_load_ini(config_t* c, const char* path) {
    FILE* f = fopen(path, "r");
    if (!f) { log_warn("config: cannot open ini '%s'", path); return -1; }
    set_str(c->ini_path, sizeof(c->ini_path), path);

    char line[1024];
    while (fgets(line, sizeof(line), f)) {
        // strip comments: # or ;
        char* p = line;
        for (; *p; ++p) {
            if (*p == '#' || *p == ';') { *p = '\0'; break; }
        }
        trim(line);
        if (!line[0]) continue;

        // section ignored but allowed
        if (line[0] == '[') continue;

        // key = value
        char* eq = strchr(line, '=');
        if (!eq) continue;
        *eq = '\0';
        char* k = line;
        char* v = eq + 1;
        trim(k); trim(v);
        if (!k[0]) continue;

        apply_kv(c, k, v);
    }

    fclose(f);
    return 0;
}

// CLI parser
static void print_help(const char* prog) {
    fprintf(stderr,
        "Usage: %s [--config PATH] [--T N] [--budget X] [--alloc_rule NAME]\n"
        "          [--results_dir DIR] [--edges_file PATH] [--init_file PATH]\n"
        "          [--targets_file PATH] [--signed 0|1] [--verbosity 0..3]\n"
        "          [--gamma X] [--clamp01 0|1]\n"
        "          [--help]\n", prog);
}

int config_parse_argv(config_t* c, int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        if (strcmp(a, "--help") == 0 || strcmp(a, "-h") == 0) {
            print_help(argv[0]);
            return 1; // signal: asked for help
        } else if (strcmp(a, "--config") == 0 && i + 1 < argc) {
            const char* p = argv[++i];
            if (config_load_ini(c, p) == 0) {
                log_info("Loaded ini: %s", p);
            }
        } else if (strcmp(a, "--T") == 0 && i + 1 < argc) {
            c->T = (int)strtol(argv[++i], NULL, 10);
        } else if (strcmp(a, "--budget") == 0 && i + 1 < argc) {
            c->budget = strtod(argv[++i], NULL);
        } else if (strcmp(a, "--alloc_rule") == 0 && i + 1 < argc) {
            set_str(c->alloc_rule, sizeof(c->alloc_rule), argv[++i]);
        } else if (strcmp(a, "--results_dir") == 0 && i + 1 < argc) {
            set_str(c->results_dir, sizeof(c->results_dir), argv[++i]);
        } else if (strcmp(a, "--edges_file") == 0 && i + 1 < argc) {
            set_str(c->edges_file, sizeof(c->edges_file), argv[++i]);
        } else if (strcmp(a, "--init_file") == 0 && i + 1 < argc) {
            set_str(c->init_file, sizeof(c->init_file), argv[++i]);
        } else if (strcmp(a, "--targets_file") == 0 && i + 1 < argc) {
            set_str(c->targets_file, sizeof(c->targets_file), argv[++i]);
        } else if (strcmp(a, "--signed") == 0 && i + 1 < argc) {
            c->signed_edges = (int)strtol(argv[++i], NULL, 10);
        } else if (strcmp(a, "--verbosity") == 0 && i + 1 < argc) {
            c->verbosity = (int)strtol(argv[++i], NULL, 10);
        } else if (strcmp(a, "--gamma") == 0 && i + 1 < argc) {
            c->gamma = strtod(argv[++i], NULL);

        } else if (strcmp(a, "--clamp01") == 0 && i + 1 < argc) {
            c->clamp01 = (int)strtol(argv[++i], NULL, 10);
        } else {
            log_warn("Unknown arg: %s (use --help)", a);
        }
    }
    return 0;
}

int config_validate(const config_t* c, char* why, size_t cap) {
    // basic
    if (c->T <= 0) { snprintf(why, cap, "T must be > 0"); return -1; }
    if (c->budget < 0) { snprintf(why, cap, "budget must be >= 0"); return -1; }
    if (!(c->gamma > 0.0 && c->gamma <= 1.0)) {
        snprintf(why, cap, "gamma must be in (0,1]"); 
        return -1;
    }
    if (!(c->clamp01 == 0 || c->clamp01 == 1)) {  
        snprintf(why, cap, "clamp01 must be 0 or 1");
        return -1;
    }

    if (!dir_exists_or_mkdir(c->results_dir)) {
        snprintf(why, cap, "results_dir not exist (and mkdir failed): %s", c->results_dir);
        return -1;
    }

    if (!file_exists(c->init_file)) {
        snprintf(why, cap, "init_file not found: %s", c->init_file);
        return -1;
    }
    if (!file_exists(c->targets_file)) {
        snprintf(why, cap, "targets_file not found: %s", c->targets_file);
        return -1;
    }
    if (!file_exists(c->edges_file)) {
        snprintf(why, cap, "edges_file not found: %s", c->edges_file);
        return -1;
    }
    return 0;
}

void config_print(const config_t* c, FILE* out) {
    fprintf(out,
        "--- CONFIG ---\n"
        "T            = %d\n"
        "budget       = %.6f\n"
        "alloc_rule   = %s\n"
        "gamma        = %.6f\n"      
        "clamp01      = %d\n"       
        "results_dir  = %s\n"
        "edges_file   = %s\n"
        "init_file    = %s\n"
        "targets_file = %s\n"
        "signed_edges = %d\n"
        "verbosity    = %d\n"
        "ini_path     = %s\n",
        c->T, c->budget, c->alloc_rule,
        c->gamma, c->clamp01,   
        c->results_dir, c->edges_file, c->init_file, c->targets_file,
        c->signed_edges, c->verbosity,
        (c->ini_path[0] ? c->ini_path : "(none)"));
}
