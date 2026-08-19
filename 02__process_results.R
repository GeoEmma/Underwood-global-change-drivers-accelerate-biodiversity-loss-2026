# # ---
# title: "02__process_results"
# manuscript: "Title: Climate change, defaunation and disease interact to accelerate biodiversity loss in a tropical hotspot"
# corresponding_author: "Emma L Underwood"
# coauthors: "Kerry A Brown, Rebekka Allgayer, Mark Mulligan, Nigel Walford, Jette Wolff"
# last update: "2026-08-19"
# output: Aggregates the raw RangeShiftR Range outputs across replicates,
#         writes summary tables, and produces the per-simulation trajectory
#         plots. The per-simulation supplementary plots (Fig S2 population
#         dynamics, Fig S3 stage composition) are generated during 01 by
#         create_all_plots_with_rasters(); this script re-summarises the
#         Range files and produces the trajectory figures.
# prerequisite: 01__simulations.R
# next script: 03__occupancy_probability.R, then 04__publication_figures.R
# ---

#==============================================================================#
#                           0. Workspace set up ----
#==============================================================================#

source("config.R")
source(file.path("R", "functions_utils.R"))
source(file.path("R", "functions_results.R"))

## Paths ----
masterfolder <- here::here("outputs", "simulations", "Master")
results_dir  <- here::here("outputs", "analysis")
range_dir    <- file.path(masterfolder, "Outputs")
analysis_dir <- file.path(results_dir, "Analysis")   # used inside plotting functions
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

#==============================================================================#
#                           1. Load scenario mapping ----
#==============================================================================#

all_scenarios <- fread(
  here::here("outputs", "simulations", "Master_simulation_mapping.csv")
)
scenario_mapping <- all_scenarios
cat("Loaded", nrow(all_scenarios), "simulation scenarios\n")

#==============================================================================#
#                           2. Aggregate range data across replicates ----
#==============================================================================#

# process_range_data() reads every *_Range.txt, merges scenario metadata, and
# writes mean / SD summaries (per batch x sim x scenario x year) to analysis_dir.
range_data <- process_range_data(range_dir, scenario_mapping)

#==============================================================================#
#                           3. Per-simulation trajectory plots ----
#==============================================================================#

# Global extent for fixed-scale centroid plots (max extent across final ranges).
global_extent <- with(range_data$mean, list(
  min_x = min(min_X_mean, na.rm = TRUE), max_x = max(max_X_mean, na.rm = TRUE),
  min_y = min(min_Y_mean, na.rm = TRUE), max_y = max(max_Y_mean, na.rm = TRUE)
))

for (i in seq_len(nrow(all_scenarios))) {
  plot_scenario_trajectory(
    range_data_mean = range_data$mean,
    batch_num       = all_scenarios$batch_num[i],
    sim_num         = all_scenarios$sim_id[i],
    analysis_dir    = analysis_dir,
    save_plots      = TRUE,
    global_extent   = global_extent
  )
}

# Multi-simulation centroid comparisons (Batch 1 gradients)
plot_multi_centroid(range_data$mean, batch_num = 1, sim_ids = c(1, 2, 4))   # climate gradient (low defaunation)
plot_multi_centroid(range_data$mean, batch_num = 1, sim_ids = c(1, 5, 9))   # defaunation gradient (current climate)


## Clean up ----
gc()

#==============================================================================#
#                        ----  End of script ----
#==============================================================================#
