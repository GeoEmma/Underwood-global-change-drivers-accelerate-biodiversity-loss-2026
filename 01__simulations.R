# # ---
# title: "01__simulations"
# manuscript: "Title: Climate change, defaunation and disease interact to accelerate biodiversity loss in a tropical hotspot"
# corresponding_author: "Emma L Underwood"
# coauthors: "Kerry A Brown, Rebekka Allgayer, Mark Mulligan, Nigel Walford, Jette Wolff"
# last update: "2026-08-19"
# output: Runs all 30 scenario simulations (climate x defaunation x pathogen
#         timing), 10 replicates each, with memory-managed parallelism.
# prerequisite: 00__data_preparation.R
# next script: 01b__simulation_validation.R
# ---

#==============================================================================#
#                           0. Workspace set up ----
#==============================================================================#

# Working directory, library paths and packages are set in config.R. On a HPC,
# set any cluster-specific library path in the job script.
source("config.R")
source(file.path("R", "functions_utils.R"))
source(file.path("R", "functions_simulation.R"))
source(file.path("R", "functions_results.R"))   # per-sim plots run inside run_single_simulation()

# RangeShiftR v3.0.0 is installed via config.R (pinned to the released tag).
library(RangeShiftR)

## Paths ----
basefolder   <- here::here()
masterfolder <- here::here("outputs", "simulations", "Master")

## Create the RangeShiftR output directory tree ----
# FIX: previously masterfolder was assigned the (logical) return value of
# dir.create(); it is now a path, and the directory is created separately.
dir.create(masterfolder, recursive = TRUE, showWarnings = FALSE)
for (d in c("Inputs", "Outputs", "Output_Maps", "Plots")) {
  dir.create(file.path(masterfolder, d), recursive = TRUE, showWarnings = FALSE)
}

# =============================================================================#
#              1. Config unique batch/sim scenarios  ----
# =============================================================================#

### Set batch and sim numbers ----
ssp <- c('current','ssp1', 'ssp5')
emig_prob <- c(0.2, 0.1, 0.05)

# Batch 1: All 9 combinations (3 climate × 3 lemur)
batch1_scenarios <- expand.grid(batch_num = 1, ssp = ssp, emig_prob = emig_prob, stringsAsFactors = FALSE)
# Batch 2: All 6 combinations (2 climate × 3 lemur)
batch2_scenarios <- expand.grid(batch_num = 2, ssp = ssp, emig_prob = emig_prob, stringsAsFactors = FALSE)
# Batch 3: Only SSP scenarios (2 climate × 3 lemur) - no current climate repeats
batch3_scenarios <- expand.grid(batch_num = 3, ssp = c('ssp1', 'ssp5'), emig_prob = emig_prob, stringsAsFactors = FALSE)
# Batch 4: Only SSP scenarios (2 climate × 3 lemur) - no current climate repeats
batch4_scenarios <- expand.grid(batch_num = 4, ssp = c('ssp1', 'ssp5'), emig_prob = emig_prob, stringsAsFactors = FALSE)

# Combine all scenarios
all_scenarios <- rbind(batch1_scenarios, batch2_scenarios, batch3_scenarios, batch4_scenarios)
all_scenarios$sim_id <- 1:nrow(all_scenarios) # Add simulation IDs

# Create and save mapping file for the unique sims
mapping_file <- data.frame(
  sim_id = all_scenarios$sim_id,
  batch_num = all_scenarios$batch_num,
  ssp = all_scenarios$ssp,
  emig_prob = all_scenarios$emig_prob,
  description = paste0("Batch", all_scenarios$batch_num, "_", all_scenarios$ssp, "_emig", all_scenarios$emig_prob),
  scenario_code = paste0(
    # Emigration coding: 0.2=1, 0.1=2, 0.05=3
    ifelse(all_scenarios$emig_prob == 0.2, "1", 
           ifelse(all_scenarios$emig_prob == 0.1, "2", "3")),
    # SSP coding: current=a, ssp1=b, ssp5=c  
    ifelse(all_scenarios$ssp == "current", "a",
           ifelse(all_scenarios$ssp == "ssp1", "b",
                  ifelse(all_scenarios$ssp == "ssp5", "c")))
  ),
  pathogen_description = case_when(
    all_scenarios$batch_num == 1 ~ "No pathogens (baseline)",
    all_scenarios$batch_num == 2 & all_scenarios$ssp == "current" ~ "Pathogens at year 100 - baseline climate",
    all_scenarios$batch_num == 3 & all_scenarios$ssp == "current" ~ "Pathogens at year 100 - baseline climate",
    all_scenarios$batch_num == 4 & all_scenarios$ssp == "current" ~ "Pathogens at year 100 - baseline climate",
    all_scenarios$batch_num == 2 ~ "Pathogens at year 101 - dynamic landscape",
    all_scenarios$batch_num == 3 ~ "Pathogens at year 131 - dynamic landscape", 
    all_scenarios$batch_num == 4 ~ "Pathogens at year 161 - dynamic landscape"
  )
)
# write to file
write.csv(mapping_file, file = paste0(basefolder, "/01__simulations/Master_simulation_mapping.csv"), row.names = FALSE)

# Clean up
cleanup_objects(batch1_scenarios, batch2_scenarios, batch3_scenarios, batch4_scenarios, 
                mapping_file, batch_num, ssp, emig_prob, required_dirs, dir)

# Path to landscape, pathogen scaling and species distribution data 
data_files <- list(
  pathogen_file = file.path(basefolder, "01__simulations", "landscapes", "pathogen_scaling_data.RDS"),
  landscape_file = file.path(basefolder, "01__simulations", "landscapes", "landscapefile.RDS"),
  spdist_file = file.path(basefolder, "01__simulations", "landscapes", "spdistfile_binary_40.RDS")
)

# Temp superseding to add more memory management to function (11.08.2025)
# NOTE: extract_scenario_data() moved to R/functions_simulation.R




#==============================================================================#
#                         2. Main simulation set up ----
#==============================================================================#

### Start timing
start_time <- Sys.time()
check_memory("Start")


# Function to process simulations with periodic worker refresh
# NOTE: run_simulations_with_memory_management() moved to R/functions_simulation.R



# Function to run a single simulation in full
# NOTE: run_single_simulation() moved to R/functions_simulation.R



#==============================================================================#
#               3. Run model scenarios with memory management ----
#==============================================================================#


# Use function with mem management
simulation_results <- run_simulations_with_memory_management(
  sim_indices = which(all_scenarios$sim_id %in% c(8:12)), # 8-12, 13-21 on UP HPC, 22-31, 32-42 on ecoc9
  all_scenarios = all_scenarios,
  data_files = data_files,
  masterfolder = masterfolder,
  ncores = 4,  # 4 - slightly fewer cores (than 7 on previous runs) for stability and mem sharing
  sims_per_worker_refresh = 2  # Refresh workers every 2 sims
)

# Clean up
cleanup_objects(all_scenarios, cl, ncores)

# Calculate successful runs
successful_sims <- sapply(simulation_results, function(x) !is.null(x) && !inherits(x, "try-error") && 
                            (is.list(x) && x$status == "completed"))

# Final summary
end_time <- Sys.time()
total_duration <- end_time - start_time

cat("\n=== Final Summary ===\n")
cat("Total wallclock time:", format(total_duration), "\n")
cat("Completed successfully:", sum(successful_sims), "out of", length(simulation_results), "scenarios\n")

check_memory("Final")

cat("\n Simulation batch completed!\n")

# # Clean up
gc()

# =============================================================================#
#                                 End of script  ----
# =============================================================================#