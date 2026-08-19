# # ---
# title: "functions_simulation.R"
# manuscript: "Title: Climate change, defaunation and disease interact to accelerate biodiversity loss in a tropical hotspot"
# corresponding_author: "Emma L Underwood"
# coauthors: "Kerry A Brown, Rebekka Allgayer, Mark Mulligan, Nigel Walford, Jette Wolff"
# last update: "2026-08-19"
# output: RangeShiftR simulation setup, execution and memory management used by
#         01__simulations.R.
# requires: RangeShiftR, foreach, doParallel, data.table, terra  (loaded centrally in config.R)
# ---

# Functions are grouped by pipeline stage. Bodies are unchanged from the
# original working code; only their location, headers and section markers
# have been standardised for the public repository.


#==============================================================================#
#   1. Scenario data extraction ----
#==============================================================================#

extract_scenario_data <- function(sim_id, batch_num, ssp, data_files) {
  
  # Use local environment to ensure cleanup
  local({
    
    if (ssp == "current") {
      if (batch_num == 1) {
        required_years <- c(0)
      } else {
        required_years <- c(0, 100)
      }
    } else {
      required_years <- c(0, 101, 131, 161)
    }
    
    # Extract data with immediate cleanup
    result <- list()
    
    # Pathogen data - load, extract, immediately remove
    if (batch_num > 1) {
      pathogen_temp <- readRDS(data_files$pathogen_file)
      pathogen_key <- if (ssp == "current") "current_pathogen" else paste0(ssp, "_batch", batch_num)
      
      if (pathogen_key %in% names(pathogen_temp)) {
        result$pathogen_data <- pathogen_temp[[pathogen_key]]
        result$demog_layers_list <- if (ssp == "current" && batch_num > 1) {
          result$pathogen_data$layers[as.character(c(0, 100))]
        } else {
          result$pathogen_data$layers[as.character(required_years)]
        }
      }
      rm(pathogen_temp)
      for(i in 1:2) gc(verbose = FALSE, full = TRUE)
    }
    
    # Landscape data - process one at a time
    landscape_temp <- readRDS(data_files$landscape_file)
    result$landscape_data <- list()
    
    for (j in seq_along(required_years)) {
      year <- required_years[j]
      key <- if (ssp == "current") {
        if (year == 0) "current_land0" else paste0("current_land", year)
      } else {
        if (year == 0) "current_land0" else paste0(ssp, "_land", year)
      }
      
      if (key %in% names(landscape_temp)) {
        result$landscape_data[[j]] <- landscape_temp[[key]]
      }
    }
    rm(landscape_temp)
    for(i in 1:2) gc(verbose = FALSE, full = TRUE)
    
    # Species distribution
    spdist_temp <- readRDS(data_files$spdist_file)
    result$spdist_data <- list(spdist_temp[[1]])
    rm(spdist_temp)
    for(i in 1:2) gc(verbose = FALSE, full = TRUE)
    
    result$years <- required_years
    return(result)
  })
}  # End of function extract_scenario_data


#==============================================================================#
#   2. Single-simulation runner ----
#==============================================================================#

run_single_simulation <- function(sim_num, all_scenarios, data_files, masterfolder) {
  
  # Wrap entire function in error handling
  result <- tryCatch({
    
    ## 0. Extract single sim params ----
    sim_id <- all_scenarios[sim_num,]$sim_id
    batch_num <- all_scenarios[sim_num,]$batch_num
    ssp <- all_scenarios[sim_num,]$ssp
    emig_prob <- all_scenarios[sim_num,]$emig_prob
    
    cat("Worker", Sys.getpid(), "starting sim", sim_id, "- Batch:", batch_num, "SSP:", ssp, "Emig:", emig_prob, "\n")
    
    # Force garbage collection at start
    gc(verbose = FALSE)
    
    # Check available memory (Linux only - optional, won't break on other systems)
    if (Sys.info()["sysname"] == "Linux") {
      tryCatch({
        meminfo <- system("free -m", intern = TRUE)
        cat("Available memory before sim", sim_id, ":\n", meminfo[2], "\n")
      }, error = function(e) {
        # Silent fail if memory check doesn't work
      })
    }
    
    # Use extract_scenario_data() to get sim specific data
    extracted_data <- extract_scenario_data(sim_id, batch_num, ssp, data_files)
    landscape_data <- extracted_data$landscape_data
    spdist_data <- extracted_data$spdist_data
    pathogen_data <- extracted_data$pathogen_data
    demog_layers_list <- extracted_data$demog_layers_list
    years <- extracted_data$years
    
    # Clean up extraction and force GC
    rm(extracted_data)
    gc(verbose = FALSE)
    
    ## 1. Simulation object ----
    sim_obj <- safe_execute({
      RangeShiftR::Simulation(
        Simulation = as.integer(sim_id),
        Years = 201,
        Replicates = 10, # 10 reps proven to be suitable
        OutIntPop = 10, # 10 for size constraints
        OutIntRange = 1 # 1 as file is small
      )
    }, "Creating simulation object")
    
    if (is.null(sim_obj)) {
      stop(paste("Failed to create simulation object for sim", sim_id))
    }
    
    # Force GC after each major object creation
    gc(verbose = FALSE)
    
    ## 2. Landscape object ----
    # Determine dynamic years based on climate scenario and pathogen presence
    if (ssp == "current") {
      # Current climate scenarios
      if (batch_num == 1) {
        dynamic_years <- NULL  # Static baseline
      } else {
        dynamic_years <- c(0, 100)  # Dynamic due to pathogen introduction
      }
    } else {
      # SSP scenarios - always dynamic due to climate change
      dynamic_years <- c(0, 101, 131, 161)
    }
    
    # Create landscape object 
    # Updated April 2026 for v3.0.0 (replaced "LandscapeFile" with "LandscapeMatrix =" and added "OriginCoords")
    land_obj <- safe_execute({
      if (is.null(demog_layers_list)) { # no pathogen scaling
        if (is.null(dynamic_years)) { # Batch 1 - current climate only
          RangeShiftR::ImportedLandscape(
            LandscapeMatrix = landscape_data,
            Resolution = 100,
            OriginCoords = c(0,0),
            HabPercent = TRUE,
            K_or_DensDep = 6,
            SpDistMatrix = spdist_data,
            SpDistResolution = 100
          )
        } else {
          # dynamic landscape - include DynamicLandYears parameter
          RangeShiftR::ImportedLandscape(
            LandscapeMatrix = landscape_data,
            Resolution = 100,
            OriginCoords = c(0,0),
            HabPercent = TRUE,
            K_or_DensDep = 6,
            DynamicLandYears = dynamic_years,
            SpDistMatrix = spdist_data,
            SpDistResolution = 100
          )
        }
      } else {
        # With dynamic and pathogen scaling (Batch 2+ always dynamic)
        RangeShiftR::ImportedLandscape(
          LandscapeMatrix = landscape_data,
          Resolution = 100,
          OriginCoords = c(0,0),
          HabPercent = TRUE,
          K_or_DensDep = 6,
          DynamicLandYears = dynamic_years,
          SpDistMatrix = spdist_data,
          SpDistResolution = 100,
          demogScaleLayers = demog_layers_list,
          nrDemogScaleLayers = pathogen_data$n_layers
        )
      }
    }, paste("Creating landscape object for sim", sim_id))
    
    if (is.null(land_obj)) {
      stop(paste("Failed to create landscape object for sim", sim_id))
    }
    
    # Clean up landscape data immediately after landscape object creation
    cleanup_objects(landscape_data, spdist_data, demog_layers_list, dynamic_years)
    gc(verbose = FALSE)
    check_memory(paste("After landscape object - sim", sim_id))
    
    ## 3. Transition matrices ----
    # Taken from literature
    base_transition_matrix <- matrix(c(
      0, 0, 1.5, 18, 30,
      1, 0.3944, 0, 0, 0,
      0, 0.0986, 0.677, 0, 0,
      0, 0, 0.125, 0.826, 0,
      0, 0, 0, 0.11, 0.99
    ), nrow = 5, byrow = TRUE)
    
    # Developed during sensitivity testing
    tier4_stg_weights_modified <- matrix(c(
      0, 0.00, 0.00, 0.000, 0.000,
      0, 0.08, 0.06, 0.090, 0.110,
      0, 0.03, 0.22, 0.100, 0.120,
      0, 0.01, 0.06, 0.200, 0.110,
      0, 0.00, 0.03, 0.150, 0.090
    ), nrow = 5, byrow = TRUE)
    
    ## 4. Stage structure ----
    # Taken from genus-wide literature
    if (is.null(pathogen_data)) {
      stage_obj <- RangeShiftR::StageStructure(
        Stages = 5,
        TransMatrix = base_transition_matrix,
        MaxAge = 400,
        MinAge = c(0, 0, 2, 5, 10),
        RepSeasons = 1,
        SurvSched = 0,
        FecDensDep = FALSE,
        DevDensDep = TRUE,
        SurvDensDep = TRUE,
        DevDensCoeff = 1,
        SurvDensCoeff = 1,
        DevStageWtsMatrix = tier4_stg_weights_modified,
        SurvStageWtsMatrix = tier4_stg_weights_modified
      )
    } else { # Batches 2-4
      surv_layer_matrix <- matrix(pathogen_data$surv_layer_mapping, nrow = 5, ncol = 1)
      
      stage_obj <- safe_execute({RangeShiftR::StageStructure(
        Stages = 5,
        TransMatrix = base_transition_matrix,
        MaxAge = 400,
        MinAge = c(0, 0, 2, 5, 10),
        RepSeasons = 1,
        SurvSched = 0,
        FecDensDep = FALSE,
        DevDensDep = TRUE,
        SurvDensDep = TRUE,
        DevDensCoeff = 1,
        SurvDensCoeff = 1,
        DevStageWtsMatrix = tier4_stg_weights_modified,
        SurvStageWtsMatrix = tier4_stg_weights_modified,
        SurvLayer = surv_layer_matrix
      )
      }, "Creating stage structure")
      cleanup_objects(surv_layer_matrix, force_gc = FALSE)
    }
    
    if (is.null(stage_obj)) {
      stop(paste("Failed to create stage structure for sim", sim_id))
    }
    
    # Clean up
    cleanup_objects(base_transition_matrix, tier4_stg_weights_modified, pathogen_data)
    gc(verbose = FALSE)
    
    ## 5. Emigration and Demography ----
    emigration_matrix <- matrix(c(0, 1, 2, 3, 4,
                                  emig_prob, 0, 0, 0, 0), nrow = 5)
    
    demo_obj <- safe_execute({
      RangeShiftR::Demography(StageStruct = stage_obj, ReproductionType = 0)
    }, paste("Creating demography object for sim", sim_id))
    
    cleanup_objects(stage_obj)
    
    emigration_obj <- safe_execute({
      Emigration(
        EmigProb = emigration_matrix,
        SexDep = FALSE,
        StageDep = TRUE,
        DensDep = FALSE,
        IndVar = FALSE,
        UseFullKern = FALSE
      )
    }, paste("Creating emigration object for sim", sim_id))
    
    cleanup_objects(emigration_matrix)
    
    ## 6. Dispersal ----
    # mean dispersal distance taken from Chapter 4 Estimating Long-Distance Dispersal in Calophyllum: An Individual-Based Gut
    # Retention Model in: Underwood, E 2026, 
    # 'Climate change, deforestation and disease: modelling triple threats to an endemic tree in a tropical biodiversity hotspot', 
    # Doctor of Philosophy (PhD), Kingston University Higher Education Corporation, Kingston upon Thames, U.K..
    # https://researchinnovation.kingston.ac.uk/en/publications/climate-change-deforestation-and-disease-modelling-triple-threats/
    disp_obj <- safe_execute({
      RangeShiftR::Dispersal(
        Emigration = emigration_obj,
        Transfer = DispersalKernel(Distances = matrix(c(183)), DoubleKernel = FALSE, StageDep = FALSE),
        Settlement = Settlement() # default: Settle = 0, to "die when unsuitable" when using DispersalKernel()
      )}, paste("Creating dispersal object for sim", sim_id))
    
    cleanup_objects(emigration_obj)
    
    ## 7. Initialisation ----
    init_obj <- safe_execute({
      RangeShiftR::Initialise(
        InitType = 1,
        SpType = 0,
        InitDens = 0,
        PropStages = c(0, 0.97, 0.01, 0.01, 0.01),
        InitAge = 2
      )
    }, "Creating initialisation object")
    
    check_memory(paste("After all objects created - sim", sim_id))
    
    ## 8. Create sim object ----
    rs_sim <- safe_execute(RangeShiftR::RSsim(
      batchnum = as.integer(batch_num),
      simul = sim_obj,
      land = land_obj,
      demog = demo_obj,
      dispersal = disp_obj,
      init = init_obj,
      seed = 10040+sim_id
    ), "Creating complete simulation")
    
    # Clean up intermediate objects
    cleanup_objects(sim_obj, land_obj, demo_obj, disp_obj, init_obj)
    
    ## 9. Validate parameters ----
    validateRSparams(rs_sim)
    check_memory(paste("Before running simulation - sim", sim_id))
    
    ## 10. Run simulation ----
    sim_start_time <- Sys.time()
    
    # Simple timeout protection (optional - remove if it causes issues)
    sim_result <- tryCatch({
      RangeShiftR::RunRS(rs_sim, masterfolder)
      "success"
    }, error = function(e) {
      cat("Simulation", sim_id, "failed with error:", e$message, "\n")
      return("error")
    })
    
    if (sim_result != "success") {
      cleanup_objects(rs_sim)
      gc(verbose = FALSE)
      return(list(sim_id = sim_id, status = sim_result, error = "Simulation failed"))
    }
    
    sim_end_time <- Sys.time()
    sim_duration <- sim_end_time - sim_start_time
    
    ## 11. Create all single-sim plots ----
    plot_result <- tryCatch({
      create_all_plots_with_rasters(
        sim_id = sim_id,
        batch_num = batch_num,
        ssp = ssp,
        emig_prob = emig_prob,
        masterfolder = masterfolder,
        include_spatial = TRUE, # spatial plots adapted from JW code
        selected_years = c(100, 120, 140, 160, 180, 200),
        spatial_sample_size = NULL
      )
      "success"
    }, error = function(e) {
      cat("Plotting failed for sim", sim_id, ":", e$message, "\n")
      "plot_error"
    })
    
    # CRITICAL: Final cleanup
    cleanup_objects(rs_sim, years, sim_start_time, sim_end_time)
    
    # Aggressive garbage collection at end
    for (i in 1:2) {
      gc(verbose = FALSE)
    }
    
    cat("Completed sim", sim_id, "- Runtime:", format(sim_duration), "\n")
    
    return(list(
      sim_id = sim_id,
      batch_num = batch_num,
      ssp = ssp,
      emig_prob = emig_prob,
      runtime = sim_duration,
      status = "completed",
      plot_status = plot_result
    ))
    
  }, error = function(e) {
    # Global error handler
    cat("FATAL ERROR in sim", ifelse(exists("sim_id"), sim_id, sim_num), ":", e$message, "\n")
    
    # Emergency cleanup
    gc(verbose = FALSE)
    
    return(list(
      sim_id = ifelse(exists("sim_id"), sim_id, sim_num),
      status = "fatal_error",
      error = e$message,
      runtime = NA
    ))
  })
  
  return(result)
}  # End of function run_single_simulation


#==============================================================================#
#   3. Memory-managed batch runner ----
#==============================================================================#

run_simulations_with_memory_management <- function(sim_indices, all_scenarios, data_files, 
                                                   masterfolder, ncores = 4, 
                                                   sims_per_worker_refresh = 2) {
  
  all_results <- list()
  
  # Process in chunks to periodically refresh workers
  chunks <- split(sim_indices, ceiling(seq_along(sim_indices) / sims_per_worker_refresh))
  
  for (chunk_idx in seq_along(chunks)) {
    chunk_sims <- chunks[[chunk_idx]]
    
    cat("\n=== Processing chunk", chunk_idx, "of", length(chunks), "===\n")
    cat("Simulations:", paste(chunk_sims, collapse = ", "), "\n")
    
    # Create fresh cluster for each chunk
    cl <- makeCluster(ncores, outfile = file.path(masterfolder, 
                                                  paste0("cluster_debug_chunk", chunk_idx, ".txt")))
    
    # Configure workers with memory limits
    clusterEvalQ(cl, {
      options(warn = 1)
      gc(reset = TRUE)  # Reset gc stats
      NULL
    })
    
    registerDoParallel(cl)
    
    # Run simulations in this chunk
    chunk_results <- foreach(sim_num = chunk_sims,
                             .packages = c("RangeShiftR", "dplyr", "terra", 
                                           "data.table", "ggplot2"),
                             .export = c("all_scenarios", "data_files", "masterfolder", 
                                         "safe_execute", "check_memory", "cleanup_objects", 
                                         "extract_scenario_data", "run_single_simulation", 
                                         "create_all_plots_with_rasters", "calculate_extinction_years", 
                                         "plot_stage_dynamics_dual", "plot_occupancy_dual",
                                         "plot_stage_composition", "plot_spatial_distribution_by_stage"
                             ),
                             .errorhandling = "pass",
                             .verbose = FALSE) %dopar% {
                               
                               # Aggressive gc at start
                               for(i in 1:3) gc(verbose = FALSE, full = TRUE)
                               
                               # Run simulation
                               result <- run_single_simulation(sim_num, all_scenarios, data_files, masterfolder)
                               
                               # Aggressive gc at end
                               for(i in 1:3) gc(verbose = FALSE, full = TRUE)
                               
                               return(result)
                             }
    
    # Stop cluster to fully release memory
    stopCluster(cl)
    registerDoSEQ()  # Reset to sequential
    
    # Store results
    all_results <- c(all_results, chunk_results)
    
    # Force R session cleanup
    gc(verbose = TRUE, full = TRUE)
    
    # Brief pause for OS memory management
    Sys.sleep(3)
    
    cat("Chunk", chunk_idx, "completed. Memory freed.\n")
  }
  
  return(all_results)
}  # End of function run_simulations_with_memory_management


#==============================================================================#
#   4. Extinction timing (shared with results stage) ----
#==============================================================================#

# NOTE: calculate_extinction_years() moved to R/functions_utils.R (shared helper)


#==============================================================================#
#                        ----  End of functions_simulation.R ----
#==============================================================================#
