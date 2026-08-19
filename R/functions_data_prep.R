# # ---
# title: "functions_data_prep.R"
# manuscript: "Title: Climate change, defaunation and disease interact to accelerate biodiversity loss in a tropical hotspot"
# corresponding_author: "Emma L Underwood"
# coauthors: "Kerry A Brown, Rebekka Allgayer, Mark Mulligan, Nigel Walford, Jette Wolff"
# last update: "2026-08-19"
# output: Landscape, species-distribution and pathogen-scaling layer builders
#         used by 00__data_preparation.R.
# requires: terra, data.table  (loaded centrally in config.R)
# ---

# Functions are grouped by pipeline stage. Bodies are unchanged from the
# original working code; only their location, headers and section markers
# have been standardised for the public repository.


#==============================================================================#
#   1. Filename parsing ----
#==============================================================================#

parse_filename_to_matrix_name <- function(filename, ssp_mapping, year_mapping) {
  
  # Handle current landscape
  if(grepl("current", filename)) {
    return("current_land0")
  }
  
  # Parse SSP scenario
  ssp_match <- regexpr("ssp[0-9]+", filename)
  if(ssp_match != -1) {
    ssp_extracted <- regmatches(filename, ssp_match)
    ssp_code <- ssp_mapping[[ssp_extracted]]
    if(is.null(ssp_code)) {
      stop("Unknown SSP scenario: ", ssp_extracted)
    }
  } else {
    stop("Could not parse SSP scenario from filename: ", filename)
  }
  
  # Parse year period
  year_match <- regexpr("[0-9]{4}-[0-9]{4}", filename)
  if(year_match != -1) {
    year_extracted <- regmatches(filename, year_match)
    year_code <- year_mapping[[year_extracted]]
    if(is.null(year_code)) {
      stop("Unknown year period: ", year_extracted)
    }
  } else {
    stop("Could not parse year period from filename: ", filename)
  }
  
  # Create matrix name
  matrix_name <- paste0(ssp_code, "_land", year_code) # for calophyllum landscape
  return(matrix_name)
}  # End of function parse_filename_to_matrix_name


parse_pathogen_filename_to_matrix_name <- function(filename, ssp_mapping, year_mapping) {
  
  # Handle current landscape
  if(grepl("current", filename)) {
    return("wilt_current_land0")
  }
  
  # Parse SSP scenario
  ssp_match <- regexpr("ssp[0-9]+", filename)
  if(ssp_match != -1) {
    ssp_extracted <- regmatches(filename, ssp_match)
    ssp_code <- ssp_mapping[[ssp_extracted]]
    if(is.null(ssp_code)) {
      stop("Unknown SSP scenario: ", ssp_extracted)
    }
  } else {
    stop("Could not parse SSP scenario from filename: ", filename)
  }
  
  # Parse year period
  year_match <- regexpr("[0-9]{4}-[0-9]{4}", filename)
  if(year_match != -1) {
    year_extracted <- regmatches(filename, year_match)
    year_code <- year_mapping[[year_extracted]]
    if(is.null(year_code)) {
      stop("Unknown year period: ", year_extracted)
    }
  } else {
    stop("Could not parse year period from filename: ", filename)
  }
  
  # Create matrix name with "wilt_" prefix
  matrix_name <- paste0("wilt_", ssp_code, "_land", year_code)
  return(matrix_name)
}  # End of function parse_pathogen_filename_to_matrix_name


#==============================================================================#
#   2. Pathogen scaling arrays ----
#==============================================================================#

convert_pathogen_to_scaling <- function(pathogen_matrix, stage_multiplier) {
  # pathogen_matrix: 0-1 scale (0 = no pathogen, 1 = max pathogen)
  # stage_multiplier: minimum survival fraction (0.251-1.000)
  # Returns: 0-100 percentage scale for RangeShiftR
  
  scaling_percentage <- (stage_multiplier * 100) + 
    (1 - pathogen_matrix) * (100 - stage_multiplier * 100)
  
  return(scaling_percentage)
}  # End of function convert_pathogen_to_scaling


create_no_pathogen_array <- function(template_matrix, n_stages = 5) {
  dims <- c(nrow(template_matrix), ncol(template_matrix), n_stages)
  scaling_array <- array(100, dim = dims)  # 100% survival = no pathogen effect
  
  # Copy NA pattern from template to all stages 
  # # (this is so the study landscape and pathogen layers are identical across arrays)
  na_mask <- is.na(template_matrix)
  for (stage in 1:n_stages) {
    scaling_array[,,stage][na_mask] <- NA
  }
  
  return(scaling_array)
}  # End of function create_no_pathogen_array


create_pathogen_array <- function(pathogen_matrix, pathogen_multipliers, n_stages = 5) {
  dims <- c(nrow(pathogen_matrix), ncol(pathogen_matrix), n_stages)
  scaling_array <- array(1, dim = dims)
  
  # Fill with stage-specific scaling
  for (stage in 1:n_stages) {
    scaling_array[,,stage] <- convert_pathogen_to_scaling(
      pathogen_matrix, 
      pathogen_multipliers[stage]
    )
  }
  
  na_mask <- is.na(pathogen_matrix) # copy the NA pattern (outside of study area shape) to all stages
  for (stage in 1:n_stages) {
    scaling_array[,,stage][na_mask] <- NA
  }
  return(scaling_array)
}  # End of function create_pathogen_array


#==============================================================================#
#   3. Input-file builders (called by 00__data_preparation.R) ----
#==============================================================================#

create_species_distribution_file <- function(landscape_path = NULL,
                                             presence_threshold = 40) {
  # Default path if not provided
  if(is.null(landscape_path)) {
    landscape_path <- file.path(folder_habmaps, "calo_current_CAZ_100m_perc.tif")
  }
  cat("Input landscape:", landscape_path, "\n")
  cat("Presence threshold: >=", presence_threshold, "\n\n")
  
  # Check if file exists
  if(!file.exists(landscape_path)) {
    stop("Landscape file not found: ", landscape_path)
  }

  #============================================================================#
  ## Load landscape ----
  #============================================================================#
  
  # Load landscape for checking purposes
  landscape <- terra::rast(landscape_path)
  all_values <- terra::values(landscape)
  
  # Convert values outside 0-100 range to NA
  all_values[all_values < 0 | all_values > 100] <- NA
  
  # Calculate habitat quality after conversion
  valid_values <- all_values[!is.na(all_values)]
  
  cat("AFTER CONVERSION:\n")
  cat("- Total cells:", format(length(all_values), big.mark = ","), "\n")
  cat("- Valid cells (0-100):", format(length(valid_values), big.mark = ","), "\n")
  cat("- NA cells:", format(sum(is.na(all_values)), big.mark = ","), "\n")
  cat("- Quality range:", min(valid_values), "-", max(valid_values), "\n")
  cat("- Mean quality:", round(mean(valid_values), 1), "\n\n")
  
  #============================================================================#
  ## Apply threshold to create binary pres/abs ----
  #============================================================================#
  
  # Apply threshold to create binary presence
  # Create binary values
  binary_values <- rep(NA, length(all_values))  # Start with NA
  binary_values[!is.na(all_values) & all_values >= presence_threshold] <- 1  # Presence
  binary_values[!is.na(all_values) & all_values < presence_threshold] <- 0   # Absence
  # NA values remain NA
  
  # Calculate summary
  presence_cells <- sum(binary_values == 1, na.rm = TRUE)
  absence_cells <- sum(binary_values == 0, na.rm = TRUE)
  na_cells <- sum(is.na(binary_values))
  
  # Summary of cells
  cat("- Presence cells (1):", format(presence_cells, big.mark = ","), 
      "(", round(100 * presence_cells / length(valid_values), 1), "%)\n")
  cat("- Absence cells (0):", format(absence_cells, big.mark = ","), 
      "(", round(100 * absence_cells / length(valid_values), 1), "%)\n")
  cat("- NA cells:", format(na_cells, big.mark = ","), "\n\n")
  
  #============================================================================#
  ## Convert to matrix and create file paths ----
  #============================================================================#
  
  # Convert to matrix format
  # Create binary raster
  binary_raster <- landscape
  terra::values(binary_raster) <- binary_values
  
  # Convert to matrix
  binary_matrix <- as.matrix(binary_raster, wide = TRUE)
  
  # Create list with 1 matrix (RangeShiftR version requirement)
  spdist_list <- list(binary_matrix)
  
  # Create output paths
  input_dir <- dirname(landscape_path)
  output_dir <- file.path(basefolder,"01__simulations","landscapes")
  input_name <- "spdistfile"
  base_name <- paste0(input_name, "_binary_", presence_threshold)
  
  rds_path <- file.path(output_dir, paste0(base_name, ".RDS"))
  asc_path <- file.path(output_dir, paste0(base_name, ".asc"))
  tif_path <- file.path(output_dir, paste0(base_name, ".tif"))
  txt_path <- file.path(output_dir, paste0(base_name, ".txt"))
  
  #============================================================================#
  ## Save files ----
  #============================================================================#
  
  # 1. Save as .RDS (primary file)
  cat("- Saving .RDS file:", rds_path, "\n")
  saveRDS(spdist_list, rds_path)
  
  # 2. Save as .asc (ASCII grid)
  cat("- Saving .asc file:", asc_path, "\n")
  terra::writeRaster(binary_raster, asc_path,
                     filetype = "AAIGrid",
                     overwrite = TRUE,
                     datatype = "INT2S",
                     NAflag = -9999)  # Explicit NA flag
  
  # 3. Save as .tif
  cat("- Saving .tif file:", tif_path, "\n")
  terra::writeRaster(binary_raster, tif_path,
                     filetype = "GTiff",
                     overwrite = TRUE,
                     datatype = "INT2S",
                     NAflag = -9999)  # Explicit NA flag
  
  # 4. Save as .txt (same as .asc but with .txt extension)
  cat("- Saving .txt file:", txt_path, "\n")
  terra::writeRaster(binary_raster, txt_path,
                     filetype = "AAIGrid",
                     overwrite = TRUE,
                     datatype = "INT2S",
                     NAflag = -9999)  # Explicit NA flag
  
  #============================================================================#
  ## Verify files ----
  #============================================================================#
  
  # Check .RDS file
  if(file.exists(rds_path)) {
    test_load <- readRDS(rds_path)
    cat("✓ .RDS file: List with", length(test_load), "matrix(es),", 
        "dimensions:", paste(dim(test_load[[1]]), collapse = " x "), "\n")
  }
  
  # Check .asc file
  if(file.exists(asc_path)) {
    header_lines <- readLines(asc_path, n = 6)
    cat("✓ .asc file created, header:\n")
    for(line in header_lines) {
      cat("    ", line, "\n")
    }
  }
  
  # Check other files
  if(file.exists(tif_path)) cat(".tif file created\n")
  if(file.exists(txt_path)) cat(".txt file created\n")
  
  cat("\nCOMPLETE! Binary distribution files created.\n")
  cat("Primary file: spdistfile.RDS (list with 1 matrix)\n")
  cat("Backups: .asc, .tif, .txt\n")
  cat("Values: 1 = presence (>=", presence_threshold, "), 0 = absence, NA = no data\n\n")
  
  return(list(
    spdist_list = spdist_list,
    binary_raster = binary_raster,
    rds_path = rds_path,
    asc_path = asc_path,
    tif_path = tif_path,
    txt_path = txt_path,
    presence_cells = presence_cells,
    absence_cells = absence_cells,
    na_cells = na_cells,
    threshold = presence_threshold
  ))
}  # End of function create_species_distribution_file


create_landscape_files <- function(landscape_dir = NULL) {
  
  cat("CREATING LANDSCAPE FILES\n")
  cat("========================\n")
  
  # Default directory if not provided
  if(is.null(landscape_dir)) {
    landscape_dir <- folder_habmaps
  }
  
  cat("Input directory:", landscape_dir, "\n\n")
  
  # Define SSP and year mappings
  ssps <- c("ssp5", "ssp1")
  years <- c(0, 101, 131, 161)  # 0=current, 101=2011-2040, 131=2041-2070, 161=2071-2100
  
  # Create mapping dictionaries for each year/SSP combination
  ssp_mapping <- list( # SSPs
    "ssp585" = "ssp5",
    "ssp126" = "ssp1",
    "current" = "current"
  )
  
  year_mapping <- list( # years
    "current" = 0,
    "2011-2040" = 101,
    "2041-2070" = 131,
    "2071-2100" = 161
  )
  
  # Find all landscape files
  landscape_files <- list.files(landscape_dir, 
                                pattern = "calo_.*_CAZ_100m_perc\\.tif$",  # calo only
                                full.names = TRUE)
  
  cat("Found", length(landscape_files), "landscape files to process:\n")
  for(i in seq_along(landscape_files)) {
    cat("  ", i, ":", basename(landscape_files[i]), "\n")
  }
  cat("\n")
  
  # Initialize storage for all matrices
  all_matrices <- list()
  matrix_names <- character()
  
  # Create output directory
  output_dir <- file.path(basefolder, "01__simulations", "landscapes")
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  #============================================================================#
  ## Process each landscape file ----
  #============================================================================#
  
  for(i in seq_along(landscape_files)) { # adds checks and skips files if needed
    
    landscape_path <- landscape_files[i]
    filename <- basename(landscape_path)
    
    cat("Processing file", i, "of", length(landscape_files), ":", filename, "\n")
    cat("===============================================\n")
    
    # Parse filename to extract SSP and year information
    matrix_name <- parse_filename_to_matrix_name(filename, ssp_mapping, year_mapping)
    
    cat("Matrix name:", matrix_name, "\n")
    
    # Check if file exists
    if(!file.exists(landscape_path)) {
      cat("WARNING: File not found, skipping:", landscape_path, "\n\n")
      next
    }
    
    #==========================================================================#
    ## Load and process landscape ----
    #==========================================================================#
    
    # Load landscape
    landscape <- terra::rast(landscape_path)
    all_values <- terra::values(landscape)
    
    # Convert values outside 0-100 range to NA
    all_values[all_values < 0 | all_values > 100] <- NA
    
    # Calculate habitat quality after conversion
    valid_values <- all_values[!is.na(all_values)]
    
    cat("LANDSCAPE ANALYSIS:\n")
    cat("- Total cells:", format(length(all_values), big.mark = ","), "\n")
    cat("- Valid cells (0-100):", format(length(valid_values), big.mark = ","), "\n")
    cat("- NA cells:", format(sum(is.na(all_values)), big.mark = ","), "\n")
    cat("- Quality range:", min(valid_values), "-", max(valid_values), "\n")
    cat("- Mean quality:", round(mean(valid_values), 1), "\n\n")
    
    #==========================================================================#
    ## Convert to matrix and store ----
    #==========================================================================#
    
    # Ready for next stage of modelling - version of RangeShifter requires matrix tables
    # Update raster values
    terra::values(landscape) <- all_values
    
    # Convert to matrix
    landscape_matrix <- as.matrix(landscape, wide = TRUE)
    
    # Store matrix and name
    all_matrices[[matrix_name]] <- landscape_matrix
    matrix_names <- c(matrix_names, matrix_name)
    
    #==========================================================================#
    ## Save additional backup files with extra formats ----
    #==========================================================================#
    
    # Create individual file paths
    asc_path <- file.path(output_dir, paste0(matrix_name, ".asc"))
    tif_path <- file.path(output_dir, paste0(matrix_name, ".tif"))
    txt_path <- file.path(output_dir, paste0(matrix_name, ".txt"))
    
    # Save individual backup files
    cat("Saving individual backup files:\n")
    
    # Save as .asc
    cat("- Saving .asc file:", basename(asc_path), "\n")
    terra::writeRaster(landscape, asc_path,
                       filetype = "AAIGrid",
                       overwrite = TRUE,
                       datatype = "FLT4S",
                       NAflag = -9999)
    
    # Save as .tif
    cat("- Saving .tif file:", basename(tif_path), "\n")
    terra::writeRaster(landscape, tif_path,
                       filetype = "GTiff",
                       overwrite = TRUE,
                       datatype = "FLT4S",
                       NAflag = -9999)
    
    # Save as .txt
    cat("- Saving .txt file:", basename(txt_path), "\n")
    terra::writeRaster(landscape, txt_path,
                       filetype = "AAIGrid",
                       overwrite = TRUE,
                       datatype = "FLT4S",
                       NAflag = -9999)
    
    cat("Completed processing:", matrix_name, "\n")
    cat("===========================================\n\n")
  }
  
  
  #==========================================================================#
  ### Create dummy "current_land100" landscape for pathogen scenarios  ----
  #==========================================================================#
  
  # Create a duplicate "current_land100" which is a copy of "current_land0".
  # This is required to exactly match the number of dynamic land years and demographic
  # scaling for batches 2-4 where pathogen is present.
  
  if("current_land0" %in% names(all_matrices)) {
    cat("Creating current_land100 as duplicate of current_land0 (for pathogen synchronization)...\n")
    all_matrices[["current_land100"]] <- all_matrices[["current_land0"]]
    matrix_names <- c(matrix_names, "current_land100")
    
    # Verify they're identical
    identical_check <- identical(all_matrices[["current_land0"]], all_matrices[["current_land100"]])
    cat("✓ current_land100 created - identical to current_land0:", identical_check, "\n")
    
    # Save individual backup files for current_land100
    current_land100_matrix <- all_matrices[["current_land100"]]
    
    # Convert back to raster for saving backups
    if(exists("landscape")) {  # Use the last processed landscape as template
      template_raster <- landscape
      terra::values(template_raster) <- as.vector(current_land100_matrix)
      
      # Create file paths
      asc_path <- file.path(output_dir, "current_land100.asc")
      tif_path <- file.path(output_dir, "current_land100.tif") 
      txt_path <- file.path(output_dir, "current_land100.txt")
      
      cat("Saving backup files for current_land100:\n")
      
      # Save as .asc
      cat("- Saving .asc file: current_land100.asc\n")
      terra::writeRaster(template_raster, asc_path,
                         filetype = "AAIGrid", 
                         overwrite = TRUE,
                         datatype = "FLT4S",
                         NAflag = -9999)
      
      # Save as .tif  
      cat("- Saving .tif file: current_land100.tif\n")
      terra::writeRaster(template_raster, tif_path,
                         filetype = "GTiff",
                         overwrite = TRUE, 
                         datatype = "FLT4S",
                         NAflag = -9999)
      
      # Save as .txt
      cat("- Saving .txt file: current_land100.txt\n")
      terra::writeRaster(template_raster, txt_path,
                         filetype = "AAIGrid",
                         overwrite = TRUE,
                         datatype = "FLT4S", 
                         NAflag = -9999)
    }
    
  } else {
    cat("WARNING: current_land0 not found - cannot create current_land100\n")
  }
  
  cat("✓ Duplicate landscape key creation completed\n\n")
  
  #============================================================================#
  ## Save combined RDS file ----
  #============================================================================#
  
  cat("SAVING COMBINED LANDSCAPE FILES\n")
  cat("===============================\n")
  
  # Create final RDS path
  rds_path <- file.path(output_dir, "landscapefile.RDS")
  
  # Save combined list as .RDS (primary file)
  cat("- Saving combined .RDS file:", rds_path, "\n")
  cat("- Contains", length(all_matrices), "matrices:", paste(names(all_matrices), collapse = ", "), "\n")
  saveRDS(all_matrices, rds_path)
  
  
  #============================================================================#
  ## Verify files ----
  #============================================================================#
  
  cat("\nVERIFYING LANDSCAPE FILES\n")
  cat("=========================\n")
  
  # Check .RDS file
  if(file.exists(rds_path)) {
    test_load <- readRDS(rds_path)
    cat("✓ .RDS file: List with", length(test_load), "matrix(es)\n")
    for(name in names(test_load)) {
      cat("  -", name, ":", paste(dim(test_load[[name]]), collapse = " x "), "\n")
    }
  }
  
  # Check backup files
  backup_files <- list.files(output_dir, pattern = "_landscape\\.(asc|tif|txt)$")
  cat("✓ Backup files created:", length(backup_files), "files\n")
  
  cat("\nCOMPLETE! Landscape files created.\n")
  cat("Primary file: landscapefile.RDS (list with", length(all_matrices), "matrices)\n")
  cat("Individual backups: .asc, .tif, .txt for each matrix\n")
  cat("Values: 0-100 = habitat quality, NA = no data\n")
  cat("NODATA_value = -9999 in backup files\n\n")
  
  return(list(
    all_matrices = all_matrices,
    matrix_names = matrix_names,
    rds_path = rds_path,
    output_dir = output_dir,
    files_processed = length(landscape_files)
  ))
}  # End of function create_landscape_files

## 4a. Create pathogen scaling files ----
# Creates pathogen spatial demographic scaling layers

create_pathogen_scaling_files <- function(pathogen_dir = NULL) {
  
  cat("CREATING PATHOGEN SCALING FILES\n")
  cat("===============================\n")
  
  # Default directory if not provided
  if(is.null(pathogen_dir)) {
    pathogen_dir <- file.path(basefolder, "00__data_preparation", "HabitatMaps", "PathogenMaps")
  }
  
  cat("Input directory:", pathogen_dir, "\n\n")
  
  # Define SSP and year mappings
  ssp_mapping <- list(
    "ssp585" = "ssp5",
    "ssp126" = "ssp1",
    "current" = "current"
  )
  
  year_mapping <- list(
    "current" = 0,
    "2011-2040" = 101,
    "2041-2070" = 131,
    "2071-2100" = 161
  )
  
  # Find all pathogen files
  pathogen_files <- list.files(pathogen_dir, 
                               pattern = "vertlept_.*_CAZ_100m_perc\\.tif$", 
                               full.names = TRUE)
  
  cat("Found", length(pathogen_files), "pathogen files to process:\n")
  for(i in seq_along(pathogen_files)) {
    cat("  ", i, ":", basename(pathogen_files[i]), "\n")
  }
  cat("\n")
  
  # Initialize storage for all matrices
  pathogen_matrices <- list()
  matrix_names <- character()
  
  # Create output directory
  output_dir <- file.path(basefolder, "01__simulations", "landscapes")
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  #============================================================================#
  ### Process each pathogen file to create matrices ----
  #============================================================================#
  
  for(i in seq_along(pathogen_files)) {
    
    pathogen_path <- pathogen_files[i]
    filename <- basename(pathogen_path)
    
    cat("Processing file", i, "of", length(pathogen_files), ":", filename, "\n")
    cat("===============================================\n")
    
    # Parse filename to extract SSP and year information
    matrix_name <- parse_pathogen_filename_to_matrix_name(filename, ssp_mapping, year_mapping)
    
    cat("Matrix name:", matrix_name, "\n")
    
    # Check if file exists
    if(!file.exists(pathogen_path)) {
      cat("WARNING: File not found, skipping:", pathogen_path, "\n\n")
      next
    }
    
    #==========================================================================#
    ### Load and process pathogen file ----
    #==========================================================================#
    
    # Load pathogen raster
    pathogen_raster <- terra::rast(pathogen_path)
    all_values <- terra::values(pathogen_raster)
    
    # Convert values outside 0-100 range to NA
    all_values[all_values < 0 | all_values > 1] <- NA  # sdm output range should be 0-1 continuous decimal data
    
    # Calculate pathogen suitability after conversion
    valid_values <- all_values[!is.na(all_values)]
    
    cat("PATHOGEN ANALYSIS:\n")
    cat("- Total cells:", format(length(all_values), big.mark = ","), "\n")
    cat("- Valid cells (0-100):", format(length(valid_values), big.mark = ","), "\n")
    cat("- NA cells:", format(sum(is.na(all_values)), big.mark = ","), "\n")
    cat("- Suitability range:", min(valid_values), "-", max(valid_values), "\n")
    cat("- Mean suitability:", round(mean(valid_values), 1), "\n\n")
    
    #==========================================================================#
    ### Convert to matrix and store ----
    #==========================================================================#
    
    # Update raster values
    terra::values(pathogen_raster) <- all_values
    
    # Convert to matrix
    pathogen_matrix <- as.matrix(pathogen_raster, wide = TRUE)
    
    # Store matrix and name
    pathogen_matrices[[matrix_name]] <- pathogen_matrix
    matrix_names <- c(matrix_names, matrix_name)
    
    #==========================================================================#
    ### Save individual backup files ----
    #==========================================================================#
    
    # Create individual file paths
    asc_path <- file.path(output_dir, paste0(matrix_name, ".asc"))
    tif_path <- file.path(output_dir, paste0(matrix_name, ".tif"))
    txt_path <- file.path(output_dir, paste0(matrix_name, ".txt"))
    
    # Save individual backup files
    cat("Saving individual backup files:\n")
    
    # Save as .asc
    cat("- Saving .asc file:", basename(asc_path), "\n")
    terra::writeRaster(pathogen_raster, asc_path,
                       filetype = "AAIGrid",
                       overwrite = TRUE,
                       datatype = "FLT4S",
                       NAflag = -9999)
    
    # Save as .tif
    cat("- Saving .tif file:", basename(tif_path), "\n")
    terra::writeRaster(pathogen_raster, tif_path,
                       filetype = "GTiff",
                       overwrite = TRUE,
                       datatype = "FLT4S",
                       NAflag = -9999)
    
    # Save as .txt
    cat("- Saving .txt file:", basename(txt_path), "\n")
    terra::writeRaster(pathogen_raster, txt_path,
                       filetype = "AAIGrid",
                       overwrite = TRUE,
                       datatype = "FLT4S",
                       NAflag = -9999)
    
    cat("Completed processing:", matrix_name, "\n")
    cat("===========================================\n\n")
  }
  
  #============================================================================#
  ##            4b. Create demographic scaling layers ----
  #============================================================================#
  
  cat("CREATING DEMOGRAPHIC SCALING LAYERS\n")
  cat("===================================\n")
  
  # Experimental design: 30 unique scenarios require 6 pathogen scaling scenarios
  # Baseline climate: x 2 pathogen timing (years 0, 100) serves batches 2-4
  # SSP climate: 3 pathogen temporal timing options (starting in years 101, 131, 161) across batches 2-4
  # Year 0 pathogen is a dummy layer with 100% survival scaling, "current" acts from Year 100
  
  # Pathogen survival multipliers by stage
  pathogen_multipliers <- c(1.000, 0.757, 0.554, 0.389, 0.251)  # Stages 0-4 
  # (Computed at top of section, see  "4. Create pathogen spatial demographic scaling layers") 
  
  # Create scaling layers for each scenario combination
  scaling_data <- list()
  
  #============================================================================#
  #### Baseline climate pathogen scenario ----
  #============================================================================#
  
  # Single baseline scenario serves all baseline climate pathogen simulations
  if ("wilt_current_land0" %in% names(pathogen_matrices)) {
    current_matrix <- pathogen_matrices[["wilt_current_land0"]]
    
    # DynamicLandYears for baseline climate with pathogen
    baseline_years <- c(0, 100)
    
    # Create demogScaleLayers
    demog_scale_layers <- list()
    
    # Year 0: Dummy layer (100% scaling, no pathogen effects)
    demog_scale_layers[["0"]] <- create_no_pathogen_array(current_matrix)
    
    # Year 100: Pathogen introduction using baseline climate pathogen
    demog_scale_layers[["100"]] <- create_pathogen_array(current_matrix, pathogen_multipliers)
    
    # Store scaling data
    scaling_data[["current_pathogen"]] <- list(
      layers = demog_scale_layers,
      n_layers = 5,
      years = baseline_years,
      climate = "current",
      pathogen_start_year = 100,
      surv_layer_mapping = 1:5
    )
    
    # Validate array dimensions
    for (year in names(demog_scale_layers)) {
      arr_dims <- dim(demog_scale_layers[[year]])
      cat("  current_pathogen Year", year, "array dimensions:", paste(arr_dims, collapse = "x"), "\n")
    }
    
    cat("Created scaling layers for current_pathogen - Pathogen starts at year 100\n")
    cat("  DynamicLandYears:", paste(baseline_years, collapse = ", "), "\n")
    
  } else {
    cat("WARNING: Missing wilt_current_land0 matrix for baseline pathogen scenario\n")
  }
  
  #============================================================================#
  #### SSP climate pathogen scenarios ----
  #============================================================================#
  
  # SSP scenarios with different pathogen introduction timings
  ssp_scenarios <- c("ssp1", "ssp5")
  
  for (ssp in ssp_scenarios) {
    for (batch in 2:4) {
      # Determine pathogen start year based on batch
      pathogen_start_year <- if (batch == 2) 101 else if (batch == 3) 131 else 161
      
      # DynamicLandYears for SSP scenarios include all climate transition years
      dynamic_years <- c(0, 101, 131, 161)
      
      # Check if all required pathogen matrices exist for this SSP
      required_years <- c(101, 131, 161)
      pathogen_matrix_names <- paste0("wilt_", ssp, "_land", required_years)
      
      if (all(pathogen_matrix_names %in% names(pathogen_matrices))) {
        
        # Create demogScaleLayers for each time point
        demog_scale_layers <- list()
        
        for (year in dynamic_years) {
          if (year < pathogen_start_year) {
            # Before pathogen introduction: Dummy layer (100% scaling)
            template_matrix <- pathogen_matrices[[pathogen_matrix_names[1]]]
            demog_scale_layers[[as.character(year)]] <- create_no_pathogen_array(template_matrix)
            
          } else {
            # Pathogen introduction or continuation: Use appropriate pathogen matrix
            if (year <= 101) {
              pathogen_matrix <- pathogen_matrices[[paste0("wilt_", ssp, "_land101")]]
            } else if (year <= 131) {
              pathogen_matrix <- pathogen_matrices[[paste0("wilt_", ssp, "_land131")]]
            } else {
              pathogen_matrix <- pathogen_matrices[[paste0("wilt_", ssp, "_land161")]]
            }
            demog_scale_layers[[as.character(year)]] <- create_pathogen_array(pathogen_matrix, pathogen_multipliers)
          }
        }
        
        # Store scaling data
        key <- paste0(ssp, "_batch", batch)
        scaling_data[[key]] <- list(
          layers = demog_scale_layers,
          n_layers = 5,
          years = dynamic_years,
          batch = batch,
          climate = ssp,
          pathogen_start_year = pathogen_start_year,
          surv_layer_mapping = 1:5
        )
        
        # Validate array dimensions
        for (year in names(demog_scale_layers)) {
          arr_dims <- dim(demog_scale_layers[[year]])
          cat("  ", key, "Year", year, "array dimensions:", paste(arr_dims, collapse = "x"), "\n")
        }
        
        cat("Created scaling layers for", key, "- Pathogen starts at year", pathogen_start_year, "\n")
        cat("  DynamicLandYears:", paste(dynamic_years, collapse = ", "), "\n")
        
      } else {
        missing_matrices <- pathogen_matrix_names[!pathogen_matrix_names %in% names(pathogen_matrices)]
        cat("WARNING: Missing pathogen matrices for", ssp, "batch", batch, ":", paste(missing_matrices, collapse = ", "), "\n")
      }
    }
  }
  
  #============================================================================#
  #### Save combined RDS files ----
  #============================================================================#
  
  cat("\nSAVING COMBINED PATHOGEN FILES\n")
  cat("==============================\n")
  
  # Save pathogen matrices
  pathogen_rds_path <- file.path(output_dir, "pathogen_matrices.RDS")
  cat("- Saving pathogen matrices .RDS file:", pathogen_rds_path, "\n")
  saveRDS(pathogen_matrices, pathogen_rds_path)
  
  # Save scaling data
  scaling_rds_path <- file.path(output_dir, "pathogen_scaling_data.RDS")
  cat("- Saving pathogen scaling data .RDS file:", scaling_rds_path, "\n")
  saveRDS(scaling_data, scaling_rds_path)
  
  #============================================================================#
  #### Verify files ----
  #============================================================================#
  
  cat("\nVERIFYING PATHOGEN FILES\n")
  cat("========================\n")
  
  # Check pathogen matrices RDS file
  if(file.exists(pathogen_rds_path)) {
    test_load <- readRDS(pathogen_rds_path)
    cat("✓ Pathogen matrices .RDS file: List with", length(test_load), "matrix(es)\n")
    for(name in names(test_load)) {
      cat("  -", name, ":", paste(dim(test_load[[name]]), collapse = " x "), "\n")
    }
  }
  
  # Check scaling data RDS file
  if(file.exists(scaling_rds_path)) {
    test_load <- readRDS(scaling_rds_path)
    cat("✓ Pathogen scaling data .RDS file: List with", length(test_load), "scenario(s)\n")
    for(name in names(test_load)) {
      scenario_data <- test_load[[name]]
      cat("  -", name, ": Years", paste(scenario_data$years, collapse = ", "), 
          ", Pathogen starts:", scenario_data$pathogen_start_year, "\n")
    }
  }
  
  # Check backup files
  backup_files <- list.files(output_dir, pattern = "_pathogen\\.(asc|tif|txt)$")
  cat("✓ Backup files created:", length(backup_files), "files\n")
  
  cat("\nCOMPLETE! Pathogen files created.\n")
  cat("Primary files: pathogen_matrices.RDS, pathogen_scaling_data.RDS\n")
  cat("Individual backups: .asc, .tif, .txt for each matrix\n")
  cat("Values: 0-100 = pathogen suitability, NA = no data\n")
  cat("NODATA_value = -9999 in backup files\n")
  cat("Total pathogen scaling scenarios:", length(scaling_data), "\n")
  cat("Serves 42 unique simulation scenarios across 4 batches\n\n")
  
  # Print summary
  cat("PATHOGEN SCALING SCENARIOS:\n")
  cat("Baseline climate (serves batches 2-4):\n")
  if ("current_pathogen" %in% names(scaling_data)) {
    data <- scaling_data[["current_pathogen"]]
    cat(sprintf("  %-20s: %d layers, years: %s, pathogen starts: year %d\n", 
                "current_pathogen", data$n_layers, paste(data$years, collapse = ", "), 
                data$pathogen_start_year))
  }
  
  cat("SSP climate (batch-specific):\n")
  for (key in names(scaling_data)) {
    if (key != "current_pathogen") {
      data <- scaling_data[[key]]
      cat(sprintf("  %-20s: %d layers, years: %s, pathogen starts: year %d\n", 
                  key, data$n_layers, paste(data$years, collapse = ", "), 
                  data$pathogen_start_year))
    }
  }
  
  return(list(
    pathogen_matrices = pathogen_matrices,
    scaling_data = scaling_data,
    pathogen_rds_path = pathogen_rds_path,
    scaling_rds_path = scaling_rds_path,
    output_dir = output_dir,
    files_processed = length(pathogen_files)
  ))
}  # End of function create_pathogen_scaling_files


#==============================================================================#
#                        ----  End of functions_data_prep.R ----
#==============================================================================#
