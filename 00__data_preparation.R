# # ---
# title: "00__data_preparation"
# manuscript: "Title: Climate change, defaunation and disease interact to accelerate biodiversity loss in a tropical hotspot"
# corresponding_author: "Emma L Underwood"
# coauthors: "Kerry A Brown, Rebekka Allgayer, Mark Mulligan, Nigel Walford, Jette Wolff"
# last update: "2026-08-19"
# output: Cropped and aligned SDM data converted into the matrices required by
#         the RangeShiftR workflow. Next script: 01__simulations.R
# ---

# Aim
# - Copy the fungal-wilt SDM outputs into the project folder structure
#   (Zenodo 10.5281/zenodo.21160109)
# - Convert habitat suitability to 0-100 integer matrices on the study grid
# - Build the species-distribution, landscape and pathogen-scaling input files

#==============================================================================#
#                           0. Workspace set up ----
#==============================================================================#

# Paths, working directory and packages are handled centrally in config.R.
source("config.R")
source(file.path("R", "functions_utils.R"))
source(file.path("R", "functions_data_prep.R"))

# NOTE: path variables below are inherited from the original workflow. Confirm
# they resolve under the repository's data/ and outputs/ layout (see config.R)
# before running on a new machine.

### Base folders for 00__data preparation ----
basefolder <- "./Underwood-global-change-drivers-accelerate-biodiversity-loss-2026/" # Change as required

### Create data preparation folder structure ----
dir.create(file.path(basefolder, "data", "00__data_preparation"), recursive = TRUE, showWarnings = FALSE)
data_prep <- file.path(basefolder, "data", "00__data_preparation")

data_prep_dirs <- c(
  "StudyArea", "HabitatMaps", "SpDist"
)

for(dir in data_prep_dirs) {
  dir.create(file.path(data_prep, dir), recursive = TRUE, showWarnings = FALSE)
}

### Input folders ----
folder_habmaps <- file.path(data_prep, "HabitatMaps")
folder_spdis <- file.path(data_prep, "SpDist")
folder_AOI <- file.path(data_prep, "StudyArea")

### Download Underwood et al. (2026) ensemble SDMs ---- 
# habitat suitability raster maps from Zenodo (10.5281/zenodo.21160109)
# Users must save to file prior to moving past this point
# Folders from  SDM  - output habitat suitability to be converted to matrices
folder_SDMs <- file.path("./fungal-pathogen-SDM-ensemble-outputs/")

# Current (1981-2010 projection across MDG - 1 per spp)
folder_current <- file.path(folder_SDMs, "current_projection")
current_calo <- list.files(file.path(folder_current), pattern="calo", full.names = T) # 1 x current SDM ensemble prediction for C. paniculatum (endemic tree)
current_ <- list.files(file.path(folder_current), pattern="vertlept", full.names = T) # 1 x current SDM ensemble prediction for L. calophylli (fungal wilt)

# Future (2011-2040, 2041-2070, 2071-2100 for SSPs 126, 370 and 585 across MDG - 9 per spp)
folder_future <- file.path(folder_SDMs, "future_prediction")
future_calo <- list.files(file.path(folder_future), pattern="calo", full.names = T)
future_ <- list.files(file.path(folder_future), pattern="vertlept", full.names = T)


#==============================================================================#
#                       1. Load, crop, prepare maps ----
#==============================================================================#

# List of habitat maps
habmaps <- list()

# Calo current
calocurrent <- terra::rast(current_calo)
# Add to list with base name
habmaps[["calocurrent"]] <- calocurrent

# Calo future
for (future in future_calo) {
  # Extract base name without extension
  base_name <- tools::file_path_sans_ext(basename(future))
  # Load raster
  future_rast <- terra::rast(future)
  # Add to list using base name
  habmaps[[base_name]] <- future_rast
}

#  current
current <- terra::rast(current_)
# Add to list with base name
habmaps[["current"]] <- current

#  future
for (future in future_) {
  # Extract base name without extension
  base_name <- tools::file_path_sans_ext(basename(future))
  # Load raster
  future_rast <- terra::rast(future)
  # Add to list using base name
  habmaps[[base_name]] <- future_rast
}


#==============================================================================#
##               Rasterize CAZ corridor (study area/AOI) ----
#==============================================================================#


# Load vector of CAZ corridor created in GIS (Protected Areas of CAZ and RNP joined via Convex Hull)- made available within this repository
CAZ_AOI <- terra::vect(file.path(folder_AOI, "CAZ_AOI.shp"))

# Example climate output (with 1km resolution to match/be a multiple of cell)
calocurrent <- terra::rast(current_calo)

# Identify a single layer from the spatraster to use
calo_curr <- calocurrent$ensemble_occ_weighted # layer name within the raster

# Reproject CAZ_AOI to same projection as climate raster
CAZ_AOI_ <- terra::project(CAZ_AOI, calo_curr)

# Crop climate data to CAZ corridor
calo_cur_CAZ <- terra::crop(calo_curr, CAZ_AOI_, snap="near")

# Apply mask to set values outside the polygon to NA
calo_cur_CAZ <- terra::mask(calo_cur_CAZ, CAZ_AOI_)

# Set values to 0 so it is just a template (remove calo current predictions!)
terra::values(calo_cur_CAZ) <- 0

# Disagg to 100m 
caz_100 <- terra::disagg(calo_cur_CAZ, fact=10)

# Write 1km template of CAZ corridor to file to be used in future
terra::writeRaster(calo_cur_CAZ, filename = file.path(folder_AOI, "CAZ_1km.tif"), overwrite=T)

# Write 100m template of CAZ corridor to file
terra::writeRaster(caz_100, filename = file.path(folder_AOI, "CAZ_100m.tif"), overwrite=T)

## Load CAZ corridors
caz_100 <- terra::rast(file.path(folder_AOI, "CAZ_100m.tif"))
caz_1k <- terra::rast(file.path(folder_AOI, "CAZ_1km.tif"))


#==============================================================================#
##               Convert full hab maps to 1km matrix ----
#==============================================================================#

# This section leaves habitat suitability maps at their original extent (MDG) and resolution (1km)
# Convert all habitat maps to matrices and write to files
for (map_name in names(habmaps)) {
  # Get the current raster from the list
  current_rast <- habmaps[[map_name]]
  
  # Extract only the first layer (predictions have 5 layers from SDM output including error and uncertainty layers)
  if (terra::nlyr(current_rast) > 1) {
    current_rast <- current_rast[[1]]  # Extract just the first layer
    cat("Extracted first layer from", map_name, "\n") 
  }
  
  # Create a reference 1km template raster for all of Madagascar with aligned cells
  # Use first map's extent but ensure consistent cell alignment
  if (map_name == names(habmaps)[1]) {
    # Create template from first raster
    mdg_template_1km <- current_rast
    # Reset values (optional, just to make it clear it's a template)
    terra::values(mdg_template_1km) <- NA
    # Save the template for future use
    terra::writeRaster(mdg_template_1km, file.path(folder_AOI, "MDG_template_1km.tif"), overwrite=TRUE)
  }
  
  # Align to calo_cur_CAZ first (resample to match resolution and extent if needed)
  current_rast_aligned <- terra::resample(current_rast, mdg_template_1km, method="near")
  
  # Convert raster to matrix
  mat_data <- terra::as.matrix(current_rast_aligned)
  
  # Create file path for output
  output_file <- file.path(folder_habmaps, paste0(map_name, "_matrix_MDG_1km.txt"))
  
  # Write matrix to file
  write.table(mat_data, file = output_file, row.names = FALSE)
}

gc()


#==============================================================================#
##               ----  Resample/disagg full maps to 100m maps ----
#==============================================================================#


# Create 100m template for Madagascar by disaggregating the 1km template
mdg_template_100m <- terra::disagg(mdg_template_1km, fact=10)
terra::writeRaster(mdg_template_100m, file.path(folder_AOI, "MDG_template_100m.tif"), overwrite=TRUE)

# These predictions are full across all of MDG, but disaggregated to a finer scale resolution:

# Convert all habitat maps to matrices and write to files
for (map_name in names(habmaps)) {
  # Get the current raster from the list
  current_rast <- habmaps[[map_name]]
  
  # Extract only the first layer (predictions have 5 layers from SDM output including error and uncertainty layers)
  if (terra::nlyr(current_rast) > 1) {
    current_rast <- current_rast[[1]]  # Extract just the first layer
    cat("Extracted first layer from", map_name, "\n") 
  }
  
  # First ensure alignment at 1km
  current_rast_aligned <- terra::resample(current_rast, mdg_template_1km, method="near")
  
  # Then disaggregate to 100m - each 1km cell will be split into exactly 10x10 cells
  current_rast_100m <- terra::disagg(current_rast_aligned, fact=10)
  
  # Final alignment check with 100m template
  current_rast_100m <- terra::resample(current_rast_100m, mdg_template_100m, method="near")
  
  # Convert to matrix
  mat_data <- terra::as.matrix(current_rast_100m)
  
  # Create file path for output
  output_file <- file.path(folder_habmaps, paste0(map_name, "_matrix_MDG_100m.txt"))
  
  # Write matrix to file
  write.table(mat_data, file = output_file, row.names = FALSE)
  
  # Write raster to file
  terra::writeRaster(current_rast_100m, file=file.path(folder_habmaps, paste0(map_name, "_MDG_100m.tif")))
  
  cat("Processed aligned full extent 100m for:", map_name, "\n")
}

gc() # clean up 


#==============================================================================#
##               ----  Crop 100m full maps to CAZ + %----
#==============================================================================#

# This section also converts decimal SDM outputs (decimals, 0-1) to percentages of suitability between 0-100 (whole numbers)
# These predictions are cropped to CAZ corridor only, and are already disaggregated to 100m (see above):
# 
# # Load existing cropped raster files (assuming they end with "_CAZ_100m.tif")
# # These files are saved intermediate outputs from earlier in script - they are disagged but not yet cropped
existing_rasters <- list()

# Get list of all existing cropped raster files
raster_files <- list.files(folder_habmaps, pattern = "_MDG_100m\\.tif$", full.names = TRUE)


# Crop to CAZ and convert all 100m habitat maps to matrices and write to file
for (raster_file in raster_files) {
  # Extract base name (remove "_CAZ_100m.tif" suffix)
  base_name <- gsub("_MDG_100m\\.tif$", "", basename(raster_file))
  
  # Load the raster
  current_rast <- terra::rast(raster_file)
  
  # Add to existing raster list (semi-processed, saved to chp6 habmaps)
  existing_rasters[[base_name]] <- current_rast
}

for (map_name in names(existing_rasters)) {
  # Get the current raster
  current_rast <- existing_rasters[[map_name]]
  
  # Crop from MDG ext to CAZ
  current_rast <- terra::crop(current_rast, CAZ_AOI_)
  current_rast_cropped_100m <- terra::mask(current_rast, CAZ_AOI_)
  
  # Final alignment with CAZ 100m template
  current_rast_cropped_100m <- terra::resample(current_rast_cropped_100m, caz_100, method="near")
  current_rast_cropped_100m <- terra::mask(current_rast_cropped_100m, caz_100)
  
  # Convert from 0-1 decimal values to 0-100 integer percentages
  current_rast_perc <- current_rast_cropped_100m * 100
  current_rast_perc <- round(current_rast_perc)
  
  # Ensure values are constrained to 0-100 range (RangeShiftR ImportedLandscape() req)
  current_rast_perc <- terra::clamp(current_rast_perc, lower=0, upper=100)
  
  # Convert to matrix
  mat_data <- terra::as.matrix(current_rast_perc)
  
  # Create file paths for percentage outputs
  output_raster_file <- file.path(folder_habmaps, paste0(map_name, "_CAZ_100m_perc.tif"))
  output_matrix_file <- file.path(folder_habmaps, paste0(map_name, "_matrix_CAZ_100m_perc.txt"))
  
  # Write percentage raster to file
  terra::writeRaster(current_rast_perc, file = output_raster_file, overwrite = TRUE)
  
  # Write percentage matrix to file
  write.table(mat_data, file = output_matrix_file, row.names = FALSE)
  
}

gc() # clean up


#==============================================================================#
# 2. Create binary presence/absence file from habitat quality threshold ----
#==============================================================================#

# NOTE: create_species_distribution_file() moved to R/functions_data_prep.R


#==============================================================================#
# 3. Create landscape files from multiple input landscapes ----
#==============================================================================#

# Creates habitat quality matrix lists (and back up of other file formats)
# Ready for RangeShiftR ImportedLandscape() function
# NOTE: create_landscape_files() moved to R/functions_data_prep.R



#==============================================================================#
#          4. Create pathogen spatial demographic scaling layers ----
#==============================================================================#

# ## Original demographic parameters from field data 
# ## # article: https://www.cambridge.org/core/journals/journal-of-tropical-ecology/article/demography-phenology-and-sex-of-calophyllum-brasiliense-clusiaceae-trees-in-the-atlantic-forest/19321CB50A079DE83187C242BC6A6486
# # Extracted from Table 1: Annual mortality and transition rates
# # Stage mapping: 0=Seed, 1=Seedling, 2=Juvenile, 3=sub-adult, 4=Adult
# 
# ## Original survival probabilities (σ) = 1 - mortality rate ----
# original_survival <- c(
#   1.000,  # Seed (0) - assumed perfect survival to next stage
#   0.493,  # Seedling (1) - from 1 - 0.507 mortality
#   0.802,  # Juvenile (2) - from 1 - 0.198 mortality  
#   0.936,  # sub-adult (3) - from 1 - 0.064 mortality
#   1.000   # Adult (4) - from 1 - 0 mortality (adjusted to 0.99 in matrix)
# )
# 
# ## Original development probabilities (γ) ----
# # Calculated as: transition rate / survival probability
# original_development <- c(
#   1.000,  # Seed→Seedling (0→1) - all seeds develop in same year
#   0.200,  # Seedling→Juvenile (1→2) - 0.0986 / 0.493 = 0.200
#   0.156,  # Juvenile→sub-adult (2→3) - 0.125 / 0.802 = 0.156
#   0.118,  # sub-adult→Adult (3→4) - 0.11 / 0.936 = 0.118
#   0.000   # Adult - no further development
# )
# 
# ## Pathogen survival multipliers by stage ----
# pathogen_survival_multipliers <- c(
#   1.000000000,  # Seed (0) - no effect, seedlings have full survival
#   0.756578947,  # Seedling (1)
#   0.553990610,  # Juvenile (2)
#   0.388559952,  # sub-adult (3)
#   0.251415797   # Adult (4) # adults are observed to have highest mortality
# )

# NOTE: create_pathogen_scaling_files() moved to R/functions_data_prep.R


#==============================================================================#
#                 Run 00__data_preparation script in full----
#==============================================================================#

#  to create species distribution file
spdist_result <- create_species_distribution_file(presence_threshold = 40)

#  to create input landscape matrix list
landscape_result <- create_landscape_files()

#  to create pathogen scaling files (multiple files, pathogen suitability + 3D arrays):
pathogen_result <- create_pathogen_scaling_files()

# # Clean up
gc()

#==============================================================================#
#                        ----  End of script ----
#==============================================================================#