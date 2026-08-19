# # ---
# title: "03__occupancy_probability"
# manuscript: "Title: Climate change, defaunation and disease interact to accelerate biodiversity loss in a tropical hotspot"
# corresponding_author: "Emma L Underwood"
# coauthors: "Kerry A Brown, Rebekka Allgayer, Mark Mulligan, Nigel Walford, Jette Wolff"
# last update: "2026-08-19"
# output: Occupancy probability rasters (GeoTiff) and point layers (GeoPackage),
#         summarised across replicates for selected simulation years.
# how:
#   - Reads the raw RangeShiftR Pop files produced by 01__simulations.R
#   - For each selected year, computes the proportion of replicates in which
#     each 100 m cell was occupied (mean occupancy probability)
#   - Writes one GeoTiff per simulation-year and one GeoPackage per simulation
#     (one layer per year) into 02__plotting_results for mapping
# prerequisite: 01__simulations.R must have produced the Batch*_Sim*_Land1_Pop.txt
#               and _Range.txt files. Run after 02 or independently.
# next script: 04__publication_figures.R (reads the OccRasters/OccPoints below)
# ---

#==============================================================================#
#                           0. Workspace set up ----
#==============================================================================#

## Load project configuration and functions ----
# config.R sets the single project root (via here::here) and loads packages
source("config.R")
source(file.path("R", "functions_occupancy.R"))

## Paths ----
masterfolder <- here::here("outputs", "simulations", "Master")
results_dir  <- here::here("outputs", "occupancy")
template_dir <- here::here("data", "StudyArea")

raster_dir <- file.path(results_dir, "OccRasters")   # GeoTiff outputs
points_dir <- file.path(results_dir, "OccPoints")    # GeoPackage outputs

dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(points_dir, recursive = TRUE, showWarnings = FALSE)

#==============================================================================#
#                           1. Load inputs ----
#==============================================================================#

## Scenario mapping ----
all_scenarios <- fread(
  here::here("outputs", "simulations", "Master_simulation_mapping.csv")
)
cat("Loaded", nrow(all_scenarios), "simulation scenarios\n")

## Template raster ----
# Must be the 100 m template so it matches the Pop file resolution and provides
# the projected extent/CRS (ESRI:102022) for the false-origin coordinate fix.
template_r <- terra::rast(file.path(template_dir, "calo_current_CAZ_100m.tif"))

#==============================================================================#
#                           2. Run configuration ----
#==============================================================================#

## Which simulations and years to summarise ----
# target_years are nominal; get_key_years() substitutes the last viable year
# for any scenario that goes extinct before a target (see functions_occupancy.R).
sim_ids      <- all_scenarios$sim_id      # all simulations
target_years <- c(110L, 140L, 170L, 200L) # post-spin-up snapshots

#==============================================================================#
#                           3. Generate rasters and points ----
#==============================================================================#

occ_results <- run_occ_rasters_batch(
  sim_ids       = sim_ids,
  all_scenarios = all_scenarios,
  masterfolder  = masterfolder,
  template_r    = template_r,
  raster_dir    = raster_dir,
  points_dir    = points_dir,
  points_format = "gpkg",       # one GeoPackage per sim, one layer per year
  target_years  = target_years,
  fix_coords    = TRUE,         # convert RangeShiftR grid coords to ESRI:102022
  occ_col       = "NInd",       # any individual present = occupied
  n_reps        = 10L
)

#==============================================================================#
#                           4. Outputs ----
#==============================================================================#

# >>> QGIS: load the following for the spatial occupancy figures
#     (manuscript Figure 4 and Supplementary Figure S4). These are also read
#     automatically by 04__publication_figures.R.
#       Rasters: outputs/occupancy/OccRasters/OccProb_Sim<XX>_Year<YYY>.tif
#       Points : outputs/occupancy/OccPoints/OccPoints_Sim<XX>.gpkg
#                (each GeoPackage holds one point layer per year:
#                 Sim<XX>_Year110, Sim<XX>_Year140, ...)
#

cat("\nOccupancy stage complete.\n")
cat("Rasters:", raster_dir, "\n")
cat("Points :", points_dir, "\n")

## Clean up ----
gc()

#==============================================================================#
#                           ----  End of script ----
#==============================================================================#
