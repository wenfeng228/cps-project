#include "common.h"
#include "log.h"
#include "config.h"
#include "sim.h"

int main(int argc, char** argv) {
    config_t cfg;
    config_set_defaults(&cfg);

    // load default ini if present, then override by argv.
    (void)config_load_ini(&cfg, "config/default.ini");
    int rc = config_parse_argv(&cfg, argc, argv);
    if (rc == 1) return 0; // --help

    log_set_verbosity(cfg.verbosity);

    char why[256];
    if (config_validate(&cfg, why, sizeof(why)) != 0) {
        log_err("Config invalid: %s", why);
        return 1;
    }

    log_info("PPI booting with config:");
    config_print(&cfg, stderr);

    if (sim_run(&cfg) != 0) {
        log_err("Simulation failed.");
        return 1;
    }
    log_info("Simulation OK.");
    return 0;
}
