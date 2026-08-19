# # ---
# title: "functions_results.R"
# manuscript: "Title: Climate change, defaunation and disease interact to accelerate biodiversity loss in a tropical hotspot"
# corresponding_author: "Emma L Underwood"
# coauthors: "Kerry A Brown, Rebekka Allgayer, Mark Mulligan, Nigel Walford, Jette Wolff"
# last update: "2026-08-19"
# output: Range-data aggregation and per-simulation supplementary plots
#         (Supplementary Figures S2 and S3) used by 02__process_results.R.
# requires: data.table, dplyr, tidyr, ggplot2, patchwork, scales, viridis, terra, cowplot  (loaded centrally in config.R)
# ---

# Functions are grouped by pipeline stage. Bodies are unchanged from the
# original working code; only their location, headers and section markers
# have been standardised for the public repository.


#==============================================================================#
#   1. Templates ----
#==============================================================================#

prepare_aggregated_templates <- function(template_dir = "./00__data_preparation/StudyArea",
                                         source_dir = "./00__data_preparation/HabitatMaps",
                                         template_100m_filename = "calo_current_CAZ_100m.tif",
                                         output_100m_filename = "CAZ_100m.tif",  # output name for written copies of templates
                                         output_1km_filename = "CAZ_1km_aggregated.tif",
                                         aggregation_factor = 10,
                                         overwrite = FALSE) { 
  
  cat("\n=== PREPARING AGGREGATED TEMPLATES ===\n")
  
  # Define file paths
  source_100m_path <- file.path(source_dir, template_100m_filename)  # Source from HabitatMaps
  template_100m_path <- file.path(template_dir, output_100m_filename)  # Destination in StudyArea
  template_1km_path <- file.path(template_dir, output_1km_filename)
  
  # Check if both already exist
  if (file.exists(template_100m_path) && file.exists(template_1km_path) && !overwrite) {
    cat("  Both templates already exist in", template_dir, "\n")
    cat("  100m:", template_100m_path, "\n")
    cat("  1km:", template_1km_path, "\n")
    cat("  Set overwrite=TRUE to regenerate\n")
    
    # Load and return both templates
    template_100m <- terra::rast(template_100m_path)
    template_1km <- terra::rast(template_1km_path)
    
    return(list(
      template_100m = template_100m,
      template_1km = template_1km,
      paths = list(
        path_100m = template_100m_path,
        path_1km = template_1km_path
      )
    ))
  }
  
  # Load 100m template from source
  cat("Loading 100m template from source...\n")
  if (!file.exists(source_100m_path)) {
    stop("Source 100m template not found: ", source_100m_path)
  }
  
  template_100m <- terra::rast(source_100m_path)
  
  # Standardize to 1/NA without pulling values into R where possible
  template_100m_standardized <- terra::ifel(!is.na(template_100m), 1, NA)
  
  terra::writeRaster(template_100m_standardized, template_100m_path, overwrite = TRUE)
  
  template_1km <- terra::aggregate(template_100m_standardized, fact = aggregation_factor, fun = "max", na.rm = TRUE)
  template_1km <- terra::ifel(!is.na(template_1km), 1, NA)
  terra::writeRaster(template_1km, template_1km_path, overwrite = TRUE)
  
  crs_info <- tryCatch({
    list(
      wkt  = terra::crs(template_1km),
      epsg = tryCatch(terra::crs(template_1km, proj=TRUE), error = function(e) NA_character_)
    )
  }, error = function(e) list(wkt = NA_character_, epsg = NA_character_))
  
  # Create metadata file
  metadata <- list(
    created = Sys.time(),
    source_100m = source_100m_path,
    output_100m = template_100m_path,
    output_1km = template_1km_path,
    aggregation_factor = aggregation_factor,
    method = "max with na.rm=TRUE (binary presence)",
    dimensions_100m = dim(template_100m_standardized),
    dimensions_1km = dim(template_1km),
    resolution_100m = terra::res(template_100m_standardized),
    resolution_1km = terra::res(template_1km),
    extent = as.vector(terra::ext(template_1km)),
    crs_info = crs_info
  )
  
  metadata_path <- file.path(template_dir, "template_aggregation_metadata.rds")
  saveRDS(metadata, metadata_path)
  cat("✓ Metadata saved to:", metadata_path, "\n")
  
  return(list(
    template_100m = template_100m_standardized,
    template_1km = template_1km,
    paths = list(
      path_100m = template_100m_path,
      path_1km = template_1km_path
    ),
    metadata = metadata
  ))
}  # End of function prepare_aggregated_templates


#==============================================================================#
#   2. Range-data aggregation ----
#==============================================================================#

process_range_data <- function(range_dir, scenario_mapping) {
  
  # Get all range files
  range_files <- list.files(range_dir, pattern = "_Range.txt$", full.names = TRUE)
  cat(sprintf("Found %d range files\n", length(range_files)))
  
  if(length(range_files) == 0) {
    stop("No range files found. Check the directory path and file pattern.")
  }
  
  # Read and combine all range data with scenario info
  all_range_data <- lapply(range_files, function(file) {
    # Extract batch_num and sim from filename
    filename <- basename(file)
    batch_num <- as.numeric(gsub(".*Batch([0-9]+).*", "\\1", filename))
    sim_id <- as.numeric(gsub(".*Sim([0-9]+).*", "\\1", filename))
    
    batch_num <- as.numeric(gsub(".*[bB]atch([0-9]+).*", "\\1", filename))
    sim_id <- as.numeric(gsub(".*[sS]im([0-9]+).*", "\\1", filename))
    
    # Read data
    data <- read.table(file, header = TRUE, sep = "\t")
    data$batch_num <- batch_num
    data$sim_id <- sim_id
    
    # Merge with scenario mapping using correct column names
    data <- merge(data, scenario_mapping, by = c("batch_num", "sim_id"))
    return(data)
  }) %>% bind_rows()
  
  # Create aggregated versions (mean and sd across replicates)
  range_data_mean <- all_range_data %>%
    group_by(batch_num, sim_id, ssp, emig_prob, Year) %>%
    summarise(
      NInds_mean = mean(NInds, na.rm = TRUE),
      NInd_stage4_mean = mean(NInd_stage4, na.rm = TRUE),
      NOccupCells_mean = mean(NOccupCells, na.rm = TRUE),
      min_X_mean = mean(min_X, na.rm = TRUE),
      max_X_mean = mean(max_X, na.rm = TRUE),
      min_Y_mean = mean(min_Y, na.rm = TRUE),
      max_Y_mean = mean(max_Y, na.rm = TRUE),
      n_reps = n(),
      .groups = 'drop'
    )
  
  range_data_sd <- all_range_data %>%
    group_by(batch_num, sim_id, ssp, emig_prob, Year) %>%
    summarise(
      NInds_sd = sd(NInds, na.rm = TRUE),
      NInd_stage4_sd = sd(NInd_stage4, na.rm = TRUE),
      NOccupCells_sd = sd(NOccupCells, na.rm = TRUE),
      min_X_sd = sd(min_X, na.rm = TRUE),
      max_X_sd = sd(max_X, na.rm = TRUE),
      min_Y_sd = sd(min_Y, na.rm = TRUE),
      max_Y_sd = sd(max_Y, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Write summarised range_data to file (Masterfolder > Plots > Analysis)
  readr::write_rds(all_range_data, file = file.path(analysis_dir, "master_all_range_data.RDS"))
  readr::write_rds(range_data_mean, file = file.path(analysis_dir, "master_all_range_data_mean.RDS"))
  readr::write_rds(range_data_sd,   file = file.path(analysis_dir, "master_all_range_data_sd.RDS"))  # FIX: was range_data_mean
  
  return(list(
    raw = all_range_data,
    mean = range_data_mean,
    sd = range_data_sd
  ))
}  # End of function process_range_data


#==============================================================================#
#   3. Per-simulation trajectory / centroid plots ----
#==============================================================================#

plot_scenario_trajectory <- function(range_data_mean, batch_num, sim_num,
                                     analysis_dir = "./analysis_output",
                                     save_plots = TRUE, global_extent = NULL) {
  
  scenario_data <- range_data_mean %>%
    filter(batch_num == !!batch_num, sim_id == sim_num)
  
  if(nrow(scenario_data) == 0) {
    cat("No data found for specified scenario\n")
    return(NULL)
  }
  
  # Create sim-specific folder within main analysis folder
  sim_folder <- file.path(analysis_dir, paste0("sim_", sim_num))
  if(save_plots && !dir.exists(sim_folder)) {
    dir.create(sim_folder, recursive = TRUE)
  }
  
  # Extract scenario info
  scenario_info <- unique(scenario_data[c("ssp", "emig_prob")])
  
  # Calculate derived metrics
  plot_data <- scenario_data %>%
    mutate(
      range_area = (max_X_mean - min_X_mean) * (max_Y_mean - min_Y_mean),
      centroid_x = (max_X_mean + min_X_mean) / 2,
      centroid_y = (max_Y_mean + min_Y_mean) / 2,
      has_adults = NInd_stage4_mean > 0
    )
  
  # Define larger text theme
  large_text_theme <- theme_minimal() +
    theme(
      axis.title = element_text(size = 36),  # 3x default (12 -> 36)
      axis.text = element_text(size = 30),   # 3x default (10 -> 30)
      plot.title = element_text(size = 42),  # 3x default (14 -> 42)
      plot.subtitle = element_text(size = 33) # 3x default (11 -> 33)
    )
  
  # Plot (a) - Population size
  p1 <- ggplot(plot_data, aes(x = Year)) +
    geom_ribbon(aes(ymin = 0, 
                    ymax = ifelse(has_adults, NInds_mean, 0)), 
                fill = "forestgreen", alpha = 0.3) +
    geom_ribbon(aes(ymin = 0, 
                    ymax = ifelse(!has_adults, NInds_mean, 0)), 
                fill = "coral", alpha = 0.3) +
    geom_line(aes(y = NInds_mean), size = 1) +
    geom_vline(xintercept = 100, linetype = "dashed", alpha = 0.5) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = sprintf("Population Dynamics - Batch %d, Sim %d", batch_num, sim_num),
         subtitle = sprintf("Climate: %s | Emigration: %s", scenario_info$ssp, scenario_info$emig_prob),
         y = "Total Individuals", x = "Year") +
    large_text_theme
  
  # Plot (b) - Range area
  p2 <- ggplot(plot_data, aes(x = Year, y = range_area/1e6)) +
    geom_area(fill = "steelblue", alpha = 0.5) +
    geom_line(size = 1) +
    geom_vline(xintercept = 100, linetype = "dashed", alpha = 0.5) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = sprintf("Range area - Batch %d, Sim %d", batch_num, sim_num),
         subtitle = sprintf("Climate: %s | Defaunation scenario: %s", scenario_info$ssp, scenario_info$emig_prob),
         y = expression(Range~area~(km^2)), x = "Year") +
    large_text_theme
  
  # Plot (c) - Centroid movement
  # (c-1) dynamic scaling
  p3_dynamic <- ggplot(plot_data, aes(x = centroid_x/1000, y = centroid_y/1000, colour = Year)) +
    geom_path(size = 0.8, alpha = 0.8,
              arrow = arrow(type = "open", length = unit(0.15, "inches"))) +
    geom_point(data = filter(plot_data, Year %in% seq(0, max(Year), by = 20)),
               size = 2, alpha = 0.6) +
    geom_point(data = filter(plot_data, Year %in% c(100, max(Year))),
               size = 4, shape = 21, fill = "white", stroke = 1.5) +
    scale_colour_viridis_c() +
    labs(title = sprintf("Centroid movement (dynamic) - Batch %d, Sim %d", batch_num, sim_num),
         subtitle = sprintf("Climate scenario: %s | Defaunation scenario: %s", scenario_info$ssp, scenario_info$emig_prob),
         x = "X (km)", y = "Y (km)") +
    coord_equal() +
    scale_x_continuous(expand = expansion(mult = 0.05)) +
    large_text_theme +
    theme(legend.title = element_text(size = 30),
          legend.text = element_text(size = 27))
  
  # (c-2) fixed scaling (uses global_extent if supplied)
  p3_fixed <- p3_dynamic +
    labs(title = sprintf("Centroid movement (Fixed extent) - Batch %d, Sim %d", batch_num, sim_num)) +
    if(!is.null(global_extent)) {
      coord_equal(
        xlim = c(global_extent$min_x, global_extent$max_x) / 1000,
        ylim = c(global_extent$min_y, global_extent$max_y) / 1000
      )
    } else {
      coord_equal()
      + scale_x_continuous(expand = expansion(mult = 0.05))
    }
  
  # Plot (d) - Occupied cells
  p4 <- ggplot(plot_data, aes(x = Year, y = NOccupCells_mean)) +
    geom_area(fill = "coral", alpha = 0.5) +
    geom_line(size = 1) +
    geom_vline(xintercept = 100, linetype = "dashed", alpha = 0.5) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = sprintf("Occupied cells - Batch %d, Sim %d", batch_num, sim_num),
         subtitle = sprintf("Climate scenario: %s | Defaunation scenario: %s", scenario_info$ssp, scenario_info$emig_prob),
         y = "Number of cells", x = "Year") +
    large_text_theme
  
  # Save individual plots with larger dimensions to accommodate text
  if(save_plots) {
    ggsave(file.path(sim_folder, "population_dynamics.png"), p1, 
           width = 14, height = 8, dpi = 300)
    ggsave(file.path(sim_folder, "range_area.png"), p2, 
           width = 14, height = 8, dpi = 300)
    ggsave(file.path(sim_folder, "centroid_movement_dynamic_ext.png"), p3_dynamic, 
           width = 14, height = 12, dpi = 300)
    ggsave(file.path(sim_folder, "centroid_movement_fixed_ext.png"), p3_fixed, 
           width = 14, height = 12, dpi = 300)
    ggsave(file.path(sim_folder, "occupied_cells.png"), p4, 
           width = 14, height = 8, dpi = 300)
    
    cat(paste0("Plots saved in: ", sim_folder, "\n"))
  }
  
  return(list(population = p1, range_area = p2, centroid_fixed = p3_fixed, centroid_dynamic = p3_dynamic, occupied = p4))
}  # End of function plot_scenario_trajectory


plot_multi_centroid <- function(range_data_mean, batch_num, sim_ids) {
  
  library(RcolourBrewer)
  n_sims <- length(sim_ids)
  cb_colours <- colourRampPalette(c("#E69F00", "#56B4E9", "#009E73"))(n_sims)
  colour_map <- tibble(sim_id = sim_ids, colour = cb_colours)
  
  # filter to sims of interest
  scenario_data <- range_data_mean %>%
    filter(batch_num == !!batch_num, sim_id %in% sim_ids) %>% # !! ensures the filtering still works inside function
    mutate(
      range_area = (max_X_mean - min_X_mean) * (max_Y_mean - min_Y_mean),
      centroid_x = (max_X_mean + min_X_mean) / 2,
      centroid_y = (max_Y_mean + min_Y_mean) / 2
    )
  
  scenario_data <- scenario_data %>%
    filter(Year >= 100, Year <= 200) %>% # show only movement post-spin up for multi-sim plot
    arrange(sim_id, Year) %>%
    group_by(sim_id) %>%
    mutate(
      xend = lead(centroid_x)/1000,
      yend = lead(centroid_y)/1000,
      centroid_x_km = centroid_x/1000,
      centroid_y_km = centroid_y/1000
    ) %>%
    ungroup()
  
  # Assign line types for multi-sim plotting
  linetypes <- rep(c("solid","dashed","dotted","dotdash"), length.out = n_sims)
  scenario_data <- scenario_data %>% left_join(tibble(sim_id = sim_ids, linetype = linetypes), by="sim_id")
  
  # assign colour by sim_id (in the order supplied)
  colour_map <- tibble(sim_id = sim_ids,
                       colour = cb_colours[seq_along(sim_ids)])
  
  # Join sim_id
  scenario_data <- scenario_data %>%
    left_join(scenario_mapping %>% select(sim_id, ssp, emig_prob), by = "sim_id")
  
  # remove NAs from data
  scenario_data <- scenario_data %>%
    filter(!is.na(centroid_x_km) & !is.na(centroid_y_km))
  
  # plot
  p <- ggplot(scenario_data, aes(x = centroid_x_km, y = centroid_y_km,
                                 colour = factor(sim_id), linetype = factor(sim_id))) +
    geom_segment(
      data = filter(scenario_data, !is.na(xend) & !is.na(yend)),
      aes(x = centroid_x_km, y = centroid_y_km, xend = xend, yend = yend,
          colour = factor(sim_id), linetype = factor(sim_id)),
      arrow = arrow(type="open", length=unit(0.12, "inches")),
      size = 1, alpha = 0.7
    )+
    labs(title = sprintf("Centroid movement - Batch %d", batch_num),
         x = "X (km)", y = "Y (km)") +
    theme_minimal(base_size = 16)
  
  # for sim specific file name
  sim_str <- paste(sim_ids, collapse = "_") 
  
  ggsave(file.path(analysis_dir, paste0("Multi-sim_centroid_movement_",sim_str,".png")), p, 
         width = 10, height = 12, dpi = 300)
  
  return(p)
}  # End of function plot_multi_centroid


plot_geographic_comparison <- function(refugial_areas, 
                                       batch_to_show = 1,
                                       save_plot = TRUE,
                                       filename = file.path(analysis_dir,"_batch_", batch_num,"_max_final_geographic_ranges.png")) {
  
  ssp_colours <- c(
    "current" = "black",
    "ssp1" = "#FDB863",
    #"ssp3" = "#E66101", # removing SSP3 for thesis submission
    "ssp5" = "#E66101"  # choosing lighter value for submission "#8B0000"
  )
  
  range_rectangles <- refugial_areas %>%
    filter(batch_num == batch_to_show)
  
  p <- ggplot(range_rectangles) +
    geom_rect(aes(xmin = min_X_mean/1000, xmax = max_X_mean/1000,
                  ymin = min_Y_mean/1000, ymax = max_Y_mean/1000,
                  colour = ssp, linetype = factor(emig_prob)),
              fill = NA, size = 1.2, alpha = 0.8) +
    scale_colour_manual(values = ssp_colours, name = "Climate scenario") +
    scale_linetype_manual(values = c("0.05" = "dotted", 
                                     "0.1" = "dashed", 
                                     "0.2" = "solid"),
                          name = "Defaunation",
                          labels = c("High (dotted)",
                                     "Medium (dashed)", # medium defaunation
                                     "Low (solid)")) +
    coord_equal() +
    labs(title = sprintf("Geographic Range Comparison - Batch %d", batch_to_show),
         subtitle = "Final range extents for all scenarios",
         x = "X coordinate (km)", 
         y = "Y coordinate (km)") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 36),
      axis.text = element_text(size = 30),
      plot.title = element_text(size = 42),
      plot.subtitle = element_text(size = 33),
      legend.title = element_text(size = 33),
      legend.text = element_text(size = 30),
      legend.position = "right"
    )
  
  if(save_plot) {
    ggsave(filename, p, width = 14, height = 12, dpi = 300)
    cat(paste0("Geographic comparison saved as: ", filename, "\n"))
  }
  
  return(p)
}  # End of function plot_geographic_comparison


theme_publication <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "sans", size = 10),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 9),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      strip.text = element_text(size = 10, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, colour = "grey80", size = 0.5)
    )
}  # End of function theme_publication


#==============================================================================#
#   4. Stage dynamics and composition (Supplementary Figs S2, S3) ----
#==============================================================================#

plot_stage_dynamics_dual <- function(pop_data, sim_id, plot_dir, start_year = 0,
                                     formats = c("pdf", "png", "svg")) {
  
  suffix <- ifelse(start_year == 0, "full", "post100")
  cat(sprintf("  Creating stage dynamics plots (%s)...\n", suffix))
  
  stage_dir <- file.path(plot_dir, paste0("stage_dynamics_sim", sim_id))
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  
  base_filename <- file.path(stage_dir, paste0("stage_dynamics_sim", sim_id, "_", suffix))
  
  # Store extinction info to return later
  extinction_info <- list()
  
  # wrap plotting code  ---
  plotting_code <- function() {
    if (dev.cur() > 1) dev.off()  # Close any open devices / refresh
    
    # Set up layout with adjusted margins
    # Reduce bottom margin for panels 1-4, increase for panel 5
    # Left margin increased to accommodate y-axis label
    par(mfrow = c(5, 1), 
        mar = c(2, 6, 2, 2),  # Default margins for panels 1-4 (bottom, left, top, right)
        oma = c(5, 4, 2, 2))  # Outer margins for overall labels
    
    stage_cols <- c("NJuvs", "NInd_stage1", "NInd_stage2", "NInd_stage3", "NInd_stage4")
    stage_names <- c("Stage 0 (Seeds)", "Stage 1 (Seedlings)", "Stage 2 (Juveniles)", 
                     "Stage 3 (Sub-adults)", "Stage 4 (Adults)")
    
    # Define n_stages
    n_stages <- length(stage_cols)
    extinction_years <- list()
    
    if (is.data.table(pop_data)) {
      plot_data <- pop_data[Year >= start_year]
    } else {
      plot_data <- pop_data %>% filter(Year >= start_year)
    }
    
    for (i in seq_along(stage_cols)) {
      stage_col <- stage_cols[i]
      stage_label <- stage_names[i]
      
      # Adjust bottom margin for last plot
      if (i == n_stages) {
        par(mar = c(5, 6, 2, 2))  # Increase bottom margin for x-axis label
      }
      
      if (!stage_col %in% names(plot_data)) {
        plot.new()
        text(0.5, 0.5, paste(stage_col, "not found"), cex = 1.2)
        extinction_years[[stage_label]] <- NA
        next
      }
      
      if (is.data.table(plot_data)) {
        stage_summary <- plot_data[, .(stage_total = sum(get(stage_col), na.rm = TRUE)), 
                                   by = .(Year, Rep)][
                                     , .(mean_stage   = as.numeric(mean(stage_total, na.rm = TRUE)),
                                         median_stage = as.numeric(median(stage_total, na.rm = TRUE)),
                                         q25          = as.numeric(quantile(stage_total, 0.25, na.rm = TRUE)),
                                         q75          = as.numeric(quantile(stage_total, 0.75, na.rm = TRUE))), 
                                     by = Year]
      } else {
        stage_summary <- plot_data %>%
          group_by(Year, Rep) %>%
          summarise(stage_total = sum(.data[[stage_col]], na.rm = TRUE), .groups = 'drop') %>%
          group_by(Year) %>%
          summarise(
            mean_stage   = as.numeric(mean(stage_total, na.rm = TRUE)),
            median_stage = as.numeric(median(stage_total, na.rm = TRUE)),
            q25          = as.numeric(quantile(stage_total, 0.25, na.rm = TRUE)),
            q75          = as.numeric(quantile(stage_total, 0.75, na.rm = TRUE)),
            .groups = 'drop'
          )
      }
      
      post_spinup <- stage_summary[stage_summary$Year >= 100,]
      extinct_year <- if (nrow(post_spinup) > 0) {
        years_extinct <- post_spinup$Year[post_spinup$mean_stage < 0.5]
        if (length(years_extinct) > 0) min(years_extinct) else NA
      } else NA
      
      extinction_years[[stage_label]] <- extinct_year
      
      y_max <- max(stage_summary$q75, na.rm = TRUE) * 1.1
      if (!is.finite(y_max) || y_max <= 0) y_max <- 1
      
      # Show x-axis only for bottom plot
      show_x <- (i == n_stages)
      
      # Create the plot
      plot(stage_summary$Year, stage_summary$mean_stage,
           type = "l", lwd = 2, col = "darkblue",
           xlab = "", ylab = "",  # We'll add labels separately
           xaxt = ifelse(show_x, "s", "n"), yaxt = 'n',
           main = stage_label,
           xlim = c(start_year, 200), ylim = c(0, y_max),
           cex.main = 1.3, cex.axis = 1.2)
      
      # Custom Y axis with formatted labels
      y_ticks <- pretty(c(0, y_max))
      axis(2, at = y_ticks, 
           labels = format_axis_labels(y_ticks),
           cex.axis = 1.2, las = 1)  # las=1 for horizontal labels
      
      # Add interquartile range
      polygon(c(stage_summary$Year, rev(stage_summary$Year)),
              c(stage_summary$q25, rev(stage_summary$q75)),
              col = rgb(0, 0, 1, alpha = 0.3), border = NA)
      
      # Add median line
      lines(stage_summary$Year, stage_summary$median_stage, col = "orange", lty = 2)
      
      # Add extinction marker
      if (!is.na(extinct_year)) {
        abline(v = extinct_year, col = "darkred", lty = 2)
        text(extinct_year + 5, y_max * 0.7, paste("Extinct\nyear", extinct_year),
             col = "darkred", cex = 0.8, adj = 0)
      }
      
      # Add year 100 reference line
      if (start_year == 0) {
        abline(v = 100, col = "grey", lty = 3)
        if (i == 1) {  # Add label only on first plot
          text(100, y_max * 0.95, "Year 100", col = "grey40", cex = 0.8, adj = 0.5)
        }
      }
      
      # Add x-axis label only for bottom plot
      if (show_x) {
        mtext("Year", side = 1, line = 3, cex = 1.3)
      }
    }
    
    # Add overall y-axis label in outer margin (rotated)
    mtext("Total individuals (mean across all reps)", 
          side = 2, outer = TRUE, line = 2, cex = 1.3, las = 0)
    
    # Store extinction info in parent environment
    extinction_info <<- extinction_years
    
    invisible(NULL)  # Don't return anything from plotting function
  } # End plotting_code
  
  # Save with taller height
  save_plot_multi_format(plotting_code, base_filename,
                         width = 10, height = 16, dpi = 300, formats = formats)
  
  cat(sprintf("    ✓ Stage dynamics plot created (%s)\n", suffix))
  
  # Return extinction info
  return(extinction_info)
}  # End of function plot_stage_dynamics_dual


plot_stage_composition <- function(pop_data, sim_id, plot_dir, extinction_years = NULL, start_year = 0,
                                   formats = c("pdf", "png", "svg"),
                                   save_legend_separate = TRUE) {
  
  suffix <- ifelse(start_year == 0, "full", "post100")
  cat(sprintf("  Creating stage composition plot (%s)...\n", suffix))
  
  comp_dir <- file.path(plot_dir, paste0("stage_composition_sim", sim_id))
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Define base_filename for output
  base_filename <- file.path(comp_dir, paste0("stage_composition_sim", sim_id, "_", suffix))
  
  # Filter by start year
  if (is.data.table(pop_data)) {
    plot_data <- pop_data[Year >= start_year]
    years_to_plot <- seq(start_year, 200, by = 10)  # Every 10th year for readability
    composition <- plot_data[Year %in% years_to_plot, 
                             .(Stage0 = sum(NJuvs, na.rm = TRUE),
                               Stage1 = sum(NInd_stage1, na.rm = TRUE),
                               Stage2 = sum(NInd_stage2, na.rm = TRUE),
                               Stage3 = sum(NInd_stage3, na.rm = TRUE),
                               Stage4 = sum(NInd_stage4, na.rm = TRUE)), 
                             by = Year]
  } else {
    plot_data <- pop_data %>% filter(Year >= start_year)
    years_to_plot <- seq(start_year, 200, by = 10)
    composition <- plot_data %>%
      filter(Year %in% years_to_plot) %>%
      group_by(Year) %>%
      summarise(
        Stage0 = sum(NJuvs, na.rm = TRUE),
        Stage1 = sum(NInd_stage1, na.rm = TRUE),
        Stage2 = sum(NInd_stage2, na.rm = TRUE),
        Stage3 = sum(NInd_stage3, na.rm = TRUE),
        Stage4 = sum(NInd_stage4, na.rm = TRUE),
        .groups = 'drop'
      )
  }
  
  # wrap plotting code  ---
  plotting_code <- function() {
    # Increased left margin
    stage_mat <- as.matrix(composition[, -1])
    y_max <- max(colSums(t(stage_mat)), na.rm = TRUE)
    max_label_width <- max(nchar(format_axis_labels(pretty(c(0, y_max)))))
    left_margin <- max(12, max_label_width * 0.8)  # Dynamic left margin for handling different length values
    
    par(mfrow = c(1, 1), mar = c(5, left_margin, 3, 2),
        cex.lab = 2, cex.axis = 2)
    
    stage_mat <- as.matrix(composition[, -1])
    
    # Calculate y-axis range
    y_max <- max(colSums(t(stage_mat)), na.rm = TRUE)
    
    # Create barplot with fixed axes
    bp <- barplot(t(stage_mat), 
                  col = c("lightblue", "skyblue", "steelblue", "darkblue", "purple"),
                  main = paste("Stage composition over time - Sim", sim_id, 
                               ifelse(start_year > 0, "(Post-100)", "")),
                  xlab = "",  # Temporarily suppress
                  ylab = "",  # Temporarily suppress
                  names.arg = composition$Year,
                  las = 1,  # Horizontal labels
                  yaxt = 'n',  # Suppress default y-axis
                  ylim = c(0, y_max * 1.05))  # Add 5% headroom
    
    # Add custom y-axis without scientific notation
    y_ticks <- pretty(c(0, y_max))
    axis(2, at = y_ticks, 
         labels = format_axis_labels(y_ticks),
         las = 1, # Horizontal y-axis labels
         cex.axis = 2)  # tick label size
    
    # Supress for now
    # Add axis labels with proper positioning
    #mtext("Total Individuals", side = 2, line = 5, cex = 2) # change axis label sizes (default =1)
    #mtext("Year", side = 1, line = 3, cex = 2)
    
    # Legend being saved separately due to too fitting issues in batch 1
    # Add legend
    #legend("topright", 
    #      legend = c("Seeds", "Seedlings", "Juveniles", "Sub-adults", "Adults"),
    #     fill = c("lightblue", "skyblue", "steelblue", "darkblue", "purple"),
    #    cex = 0.8, bty = "n")
    
    # Add reference line at year 100 if applicable
    if (start_year == 0 && 100 %in% composition$Year) {
      year_100_pos <- which(composition$Year == 100)
      abline(v = bp[year_100_pos], col = "grey", lty = 3, lwd = 0.8)
    }
    
  } # End plotting_code
  # Save in multi-format ---
  save_plot_multi_format(plotting_code, base_filename,
                         width = 12, height = 8, dpi = 300, formats = formats)
  
  cat(sprintf("    ✓ Stage composition plot created (%s)\n", suffix))
}  # End of function plot_stage_composition


plot_stage_composition_mean <- function(pop_data, sim_id, plot_dir, extinction_years = NULL, start_year = 0,
                                        formats = c("pdf", "png", "svg"),
                                        save_legend_separate = TRUE,
                                        use_mean_across_reps = TRUE) {
  
  suffix <- ifelse(start_year == 0, "full", "post100")
  if (use_mean_across_reps) {
    suffix <- paste0(suffix, "_mean")
  }
  
  cat(sprintf("  Creating stage composition plot (%s)...\n", suffix))
  
  # Updated directory structure - plots go directly into the provided plot_dir
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Define base_filename for output
  base_filename <- file.path(plot_dir, paste0("stage_composition_sim", sim_id, "_", suffix))
  
  # Filter by start year
  if (is.data.table(pop_data)) {
    plot_data <- pop_data[Year >= start_year]
    
    years_to_plot <- seq(start_year, 200, by = 10)  # Every 10th year for readability
    
    if (use_mean_across_reps && "Rep" %in% names(plot_data)) {
      # Calculate mean across replicates for each year
      composition <- plot_data[Year %in% years_to_plot, 
                               .(Stage0_mean = mean(NJuvs, na.rm = TRUE),
                                 Stage1_mean = mean(NInd_stage1, na.rm = TRUE),
                                 Stage2_mean = mean(NInd_stage2, na.rm = TRUE),
                                 Stage3_mean = mean(NInd_stage3, na.rm = TRUE),
                                 Stage4_mean = mean(NInd_stage4, na.rm = TRUE),
                                 Stage0_sd = sd(NJuvs, na.rm = TRUE),
                                 Stage1_sd = sd(NInd_stage1, na.rm = TRUE),
                                 Stage2_sd = sd(NInd_stage2, na.rm = TRUE),
                                 Stage3_sd = sd(NInd_stage3, na.rm = TRUE),
                                 Stage4_sd = sd(NInd_stage4, na.rm = TRUE),
                                 n_reps = .N), 
                               by = Year]
      
      # Rename columns for consistent processing
      setnames(composition, 
               c("Stage0_mean", "Stage1_mean", "Stage2_mean", "Stage3_mean", "Stage4_mean"),
               c("Stage0", "Stage1", "Stage2", "Stage3", "Stage4"))
      
    } else {
      # Original approach - sum across all observations
      composition <- plot_data[Year %in% years_to_plot, 
                               .(Stage0 = sum(NJuvs, na.rm = TRUE),
                                 Stage1 = sum(NInd_stage1, na.rm = TRUE),
                                 Stage2 = sum(NInd_stage2, na.rm = TRUE),
                                 Stage3 = sum(NInd_stage3, na.rm = TRUE),
                                 Stage4 = sum(NInd_stage4, na.rm = TRUE)), 
                               by = Year]
    }
    
  } else {
    # dplyr approach
    plot_data <- pop_data %>% filter(Year >= start_year)
    
    years_to_plot <- seq(start_year, 200, by = 10)
    
    if (use_mean_across_reps && "Rep" %in% names(plot_data)) {
      # Calculate mean across replicates for each year
      composition <- plot_data %>%
        filter(Year %in% years_to_plot) %>%
        group_by(Year) %>%
        summarise(
          Stage0 = mean(NJuvs, na.rm = TRUE),
          Stage1 = mean(NInd_stage1, na.rm = TRUE),
          Stage2 = mean(NInd_stage2, na.rm = TRUE),
          Stage3 = mean(NInd_stage3, na.rm = TRUE),
          Stage4 = mean(NInd_stage4, na.rm = TRUE),
          Stage0_sd = sd(NJuvs, na.rm = TRUE),
          Stage1_sd = sd(NInd_stage1, na.rm = TRUE),
          Stage2_sd = sd(NInd_stage2, na.rm = TRUE),
          Stage3_sd = sd(NInd_stage3, na.rm = TRUE),
          Stage4_sd = sd(NInd_stage4, na.rm = TRUE),
          n_reps = n(),
          .groups = 'drop'
        )
    } else {
      # Original approach - sum across all observations
      composition <- plot_data %>%
        filter(Year %in% years_to_plot) %>%
        group_by(Year) %>%
        summarise(
          Stage0 = sum(NJuvs, na.rm = TRUE),
          Stage1 = sum(NInd_stage1, na.rm = TRUE),
          Stage2 = sum(NInd_stage2, na.rm = TRUE),
          Stage3 = sum(NInd_stage3, na.rm = TRUE),
          Stage4 = sum(NInd_stage4, na.rm = TRUE),
          .groups = 'drop'
        )
    }
  }
  
  # wrap plotting code  ---
  plotting_code <- function() {
    # Increased left margin
    stage_mat <- as.matrix(composition[, c("Stage0", "Stage1", "Stage2", "Stage3", "Stage4")])
    y_max <- max(rowSums(stage_mat), na.rm = TRUE)
    max_label_width <- max(nchar(format_axis_labels(pretty(c(0, y_max)))))
    left_margin <- max(12, max_label_width * 0.8)  # Dynamic left margin
    
    par(mfrow = c(1, 1), mar = c(5, left_margin, 3, 2),
        cex.lab = 2, cex.axis = 2)
    
    # Calculate y-axis range
    y_max <- max(rowSums(stage_mat), na.rm = TRUE)
    
    # Create title based on whether using means
    main_title <- paste("Stage composition over time - Sim", sim_id)
    if (use_mean_across_reps) {
      main_title <- paste(main_title, "(Mean across reps)")
    }
    if (start_year > 0) {
      main_title <- paste(main_title, "(Post-100)")
    }
    
    # Create barplot with fixed axes
    bp <- barplot(t(stage_mat), 
                  col = c("lightblue", "skyblue", "steelblue", "darkblue", "purple"),
                  main = main_title,
                  xlab = "",  # Temporarily suppress
                  ylab = "",  # Temporarily suppress
                  names.arg = composition$Year,
                  las = 1,  # Horizontal labels
                  yaxt = 'n',  # Suppress default y-axis
                  ylim = c(0, y_max * 1.05))  # Add 5% headroom
    
    # Add custom y-axis without scientific notation
    y_ticks <- pretty(c(0, y_max))
    axis(2, at = y_ticks, 
         labels = format_axis_labels(y_ticks),
         las = 1, # Horizontal y-axis labels
         cex.axis = 2)  # tick label size
    
    # Add reference line at year 100 if applicable
    if (start_year == 0 && 100 %in% composition$Year) {
      year_100_pos <- which(composition$Year == 100)
      abline(v = bp[year_100_pos], col = "grey", lty = 3, lwd = 0.8)
    }
    
  } # End plotting_code
  
  # Save in multi-format ---
  save_plot_multi_format(plotting_code, base_filename,
                         width = 12, height = 8, dpi = 300, formats = formats)
  
  # Print summary information
  if (use_mean_across_reps && "n_reps" %in% names(composition)) {
    cat(sprintf("    ✓ Stage composition plot created (%s) - averaged across %d replicates\n", 
                suffix, composition$n_reps[1]))
  } else {
    cat(sprintf("    ✓ Stage composition plot created (%s)\n", suffix))
  }
  
  # Optionally return the processed data for further analysis
  invisible(composition)
}  # End of function plot_stage_composition_mean


create_stg_comp_plots_all_yrs <- function(pop_data, sim_id, plot_dir, 
                                          formats = c("pdf", "png", "svg"),
                                          save_legend_separate = TRUE,
                                          use_mean_across_reps = TRUE,
                                          create_both_timeframes = TRUE) {
  
  cat(sprintf("\n=== Creating stage composition plots for Sim %d ===\n", sim_id))
  
  # Create full timeline plot (Year 0-200)
  plot_stage_composition_mean(pop_data = pop_data, 
                              sim_id = sim_id, 
                              plot_dir = plot_dir, 
                              start_year = 0,
                              formats = formats,
                              save_legend_separate = save_legend_separate,
                              use_mean_across_reps = use_mean_across_reps)
  
  # Create post-100 plot (Year 100-200) if requested
  if (create_both_timeframes) {
    plot_stage_composition_mean(pop_data = pop_data, 
                                sim_id = sim_id, 
                                plot_dir = plot_dir, 
                                start_year = 100,
                                formats = formats,
                                save_legend_separate = save_legend_separate,
                                use_mean_across_reps = use_mean_across_reps)
  }
  
  # Create legend if requested
  if (save_legend_separate) {
    plot_stage_composition_legend(plot_dir = plot_dir, 
                                  sim_id = sim_id, 
                                  formats = formats)
  }
  
  cat(sprintf("=== Completed stage composition plots for Sim %d ===\n\n", sim_id))
}  # End of function create_stg_comp_plots_all_yrs


process_stg_comp_mean_plots <- function(all_scenarios, masterfolder, 
                                        formats = c("pdf", "png", "svg"),
                                        save_legend_separate = FALSE,
                                        use_mean_across_reps = TRUE,
                                        create_both_timeframes = TRUE) {
  
  cat("Starting stage composition plot generation for all scenarios...\n")
  cat(sprintf("Total scenarios to process: %d\n\n", nrow(all_scenarios)))
  
  # Initialize progress tracking
  successful_plots <- 0
  failed_plots <- 0
  failed_sims <- c()
  
  # Loop through each scenario
  for (i in 1:nrow(all_scenarios)) {
    current_scenario <- all_scenarios[i, ]
    sim_id <- current_scenario$sim_id
    batch_num <- current_scenario$batch_num
    
    cat(sprintf("Processing %d/%d: Sim %d (Batch %d) - %s\n", 
                i, nrow(all_scenarios), sim_id, batch_num, current_scenario$description))
    
    # Construct file path for range data
    range_file <- file.path(masterfolder, "Outputs", 
                            paste0("Batch", batch_num, "_Sim", sim_id, "_Land1_Range.txt"))
    
    # Construct output directory
    plot_dir <- file.path(masterfolder, "Plots", paste0("sim_", sim_id), paste0("stage_composition_mean_sim", sim_id))
    
    # Check if file exists
    if (!file.exists(range_file)) {
      cat(sprintf("  WARNING: Range file not found: %s\n", range_file))
      failed_plots <- failed_plots + 1
      failed_sims <- c(failed_sims, sim_id)
      next
    }
    
    # Try to read and process the data
    tryCatch({
      # Read the range data
      range_data <- fread(range_file)
      
      # Check if required columns exist
      required_cols <- c("Year", "Rep", "NJuvs", "NInd_stage1", "NInd_stage2", "NInd_stage3", "NInd_stage4")
      missing_cols <- setdiff(required_cols, names(range_data))
      
      if (length(missing_cols) > 0) {
        cat(sprintf("  WARNING: Missing required columns: %s\n", paste(missing_cols, collapse = ", ")))
        failed_plots <- failed_plots + 1
        failed_sims <- c(failed_sims, sim_id)
        next
      }
      
      # Create the plots
      create_stg_comp_plots_all_yrs(pop_data = range_data,
                                    sim_id = sim_id,
                                    plot_dir = plot_dir,
                                    formats = formats,
                                    save_legend_separate = save_legend_separate,
                                    use_mean_across_reps = use_mean_across_reps,
                                    create_both_timeframes = create_both_timeframes)
      
      successful_plots <- successful_plots + 1
      cat(sprintf("  Successfully processed Sim %d\n", sim_id))
      
    }, error = function(e) {
      cat(sprintf("  ERROR processing Sim %d: %s\n", sim_id, e$message))
      failed_plots <- failed_plots + 1
      failed_sims <- c(failed_sims, sim_id)
    })
  }
  
  # Print summary
  cat("\n", "="*60, "\n")
  cat("PROCESSING SUMMARY\n")
  cat("="*60, "\n")
  cat(sprintf("Total scenarios: %d\n", nrow(all_scenarios)))
  cat(sprintf("Successfully processed: %d\n", successful_plots))
  cat(sprintf("Failed: %d\n", failed_plots))
  
  if (length(failed_sims) > 0) {
    cat(sprintf("Failed simulation IDs: %s\n", paste(failed_sims, collapse = ", ")))
  }
  
  cat("="*60, "\n")
  
  # Return summary
  return(list(
    total = nrow(all_scenarios),
    successful = successful_plots,
    failed = failed_plots,
    failed_sim_ids = failed_sims
  ))
}  # End of function process_stg_comp_mean_plots


plot_occupancy_dual <- function(pop_data, range_data, sim_id, plot_dir, start_year = 0,
                                formats = c("pdf", "png", "svg")) {
  
  suffix <- ifelse(start_year == 0, "full", "post100")
  cat(sprintf("  Creating occupancy plot (%s)...\n", suffix))
  
  occ_dir <- file.path(plot_dir, paste0("occupancy_sim", sim_id))
  dir.create(occ_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Define base_filename for output
  base_filename <- file.path(occ_dir, paste0("occupancy_sim", sim_id, "_", suffix))
  
  # Calculate occupancy summary
  if (!is.null(range_data) && "NOccupCells" %in% names(range_data)) {
    if (is.data.table(range_data)) {
      plot_range <- range_data[Year >= start_year]
      occ_summary <- plot_range[, .(mean_occ = mean(NOccupCells, na.rm = TRUE),
                                    sd_occ = sd(NOccupCells, na.rm = TRUE)), 
                                by = Year]
    } else {
      plot_range <- range_data %>% filter(Year >= start_year)
      occ_summary <- plot_range %>%
        group_by(Year) %>%
        summarise(
          mean_occ = mean(NOccupCells, na.rm = TRUE),
          sd_occ = sd(NOccupCells, na.rm = TRUE),
          .groups = 'drop'
        )
    }
  } else {
    # Calculate from pop data
    if (is.data.table(pop_data)) {
      plot_pop <- pop_data[Year >= start_year]
      occ_summary <- plot_pop[, .(occupied = sum(NInd > 0, na.rm = TRUE)), 
                              by = .(Year, Rep)][
                                , .(mean_occ = mean(occupied, na.rm = TRUE),
                                    sd_occ = sd(occupied, na.rm = TRUE)), 
                                by = Year]
    } else {
      plot_pop <- pop_data %>% filter(Year >= start_year)
      occ_summary <- plot_pop %>%
        group_by(Year, Rep) %>%
        summarise(occupied = sum(NInd > 0, na.rm = TRUE), .groups = 'drop') %>%
        group_by(Year) %>%
        summarise(
          mean_occ = mean(occupied, na.rm = TRUE),
          sd_occ = sd(occupied, na.rm = TRUE),
          .groups = 'drop'
        )
    }
  }
  
  # Calculate bounds
  occ_summary$upper_bound <- occ_summary$mean_occ + occ_summary$sd_occ
  occ_summary$lower_bound <- pmax(0, occ_summary$mean_occ - occ_summary$sd_occ)
  
  # Wrap plotting code ---
  plotting_code <- function() {
    par(mfrow = c(1, 1), mar = c(4, 5, 3, 2))  # more space on the left
    
    # Determine y-axis limits
    y_max <- max(occ_summary$upper_bound, na.rm = TRUE) * 1.1
    
    # Plot the mean line
    plot(occ_summary$Year, occ_summary$mean_occ,
         type = "l", lwd = 2, col = "darkgreen",
         xlab = "" , ylab = "", # "Year", ylab = "Mean landscape cell occupancy",
         main = paste("Mean number of occupied cells by year - Sim", sim_id, 
                      ifelse(start_year > 0, "(Post-100)", "")),
         xlim = c(start_year, 200),
         ylim = c(0, y_max),
         cex.lab = 2,  # Axis labels
         #cex.main = 1.5,  # Title
         cex.axis = 2)  #  Tick labels
    
    # Add the SD polygon with upper and lower bounds
    if (!all(is.na(occ_summary$sd_occ)) && any(occ_summary$sd_occ > 0)) {
      # Create the polygon coordinates
      x_coords <- c(occ_summary$Year, rev(occ_summary$Year))
      y_coords <- c(occ_summary$upper_bound, rev(occ_summary$lower_bound))
      
      # Remove any NA values that might cause issues
      valid_points <- !is.na(x_coords) & !is.na(y_coords)
      if (sum(valid_points) > 3) {  # Need at least 3 points for a polygon
        polygon(x_coords[valid_points], 
                y_coords[valid_points],
                col = rgb(0, 0.5, 0, alpha = 0.3), 
                border = NA)
      }
    }
    
    # Add vertical line at year 100 if showing full timeline
    if (start_year == 0) {
      abline(v = 100, col = "grey", lty = 3)
    }
    
    # Add horizontal line at y=0 for reference
    abline(h = 0, col = "grey", lty = 1, lwd = 0.5)
  } # End plotting_code
  # Save in multi-format ---
  save_plot_multi_format(plotting_code, base_filename,
                         width = 12, height = 8, dpi = 300, formats = formats)
  
  cat(sprintf("    ✓ Occupancy plot created (%s)\n", suffix))
}  # End of function plot_occupancy_dual


plot_spatial_distribution_by_stage <- function(pop_data, range_data, sim_id, plot_dir, 
                                               years_to_plot = NULL,
                                               rep_number = NULL,  
                                               seed_value = NULL,
                                               sample_size = NULL, # Keep for backward compatibility
                                               formats = c("pdf", "png", "svg")) {
  
  cat("  Creating spatial distribution plots by stage...\n")
  
  spatial_dir <- file.path(plot_dir, paste0("spatial_distribution_sim", sim_id))
  dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ============================================================================#
  # STEP 1: PREPARE DATA
  # ============================================================================#
  
  # Ensure all numeric columns are consistent type (double)
  if (is.data.table(pop_data)) {
    stage_cols <- c("NInd_stage1", "NInd_stage2", "NInd_stage3", "NInd_stage4")
    # Add NInd to columns if it exists
    if ("NInd" %in% names(pop_data)) {
      stage_cols <- c(stage_cols, "NInd")
    }
    for (col in stage_cols) {
      if (col %in% names(pop_data)) {
        set(pop_data, j = col, value = as.numeric(pop_data[[col]]))
      }
    }
    
    # Filter to years 100-200 only
    plot_data <- pop_data[Year >= 100 & Year <= 200]
  } else {
    plot_data <- pop_data %>% 
      mutate(across(c(starts_with("NInd_"), any_of("NInd")), as.numeric)) %>%
      filter(Year >= 100 & Year <= 200)
  }
  
  # Filter by replicate number if specified
  if (!is.null(rep_number)) {
    if (is.data.table(plot_data)) {
      # Check if rep_number exists in data
      if (!rep_number %in% unique(plot_data$Rep)) {
        stop(sprintf("Replicate %d not found in data. Available replicates: %s", 
                     rep_number, paste(sort(unique(plot_data$Rep)), collapse = ", ")))
      }
      plot_data <- plot_data[Rep == rep_number]
      cat(sprintf("    Filtering to replicate %d\n", rep_number))
    } else {
      if (!rep_number %in% unique(plot_data$Rep)) {
        stop(sprintf("Replicate %d not found in data. Available replicates: %s", 
                     rep_number, paste(sort(unique(plot_data$Rep)), collapse = ", ")))
      }
      plot_data <- plot_data %>% filter(Rep == rep_number)
      cat(sprintf("    Filtering to replicate %d\n", rep_number))
    }
  } else {
    # If no rep specified, randomly select one
    available_reps <- unique(plot_data$Rep)
    set.seed(seed_value)  # For reproducibility
    selected_rep <- sample(available_reps, 1)
    if (is.data.table(plot_data)) {
      plot_data <- plot_data[Rep == selected_rep]
    } else {
      plot_data <- plot_data %>% filter(Rep == selected_rep)
    }
    cat(sprintf("    No replicate specified. Randomly selected replicate %d\n", selected_rep))
  }
  
  # Define max year from data
  years_available <- unique(plot_data$Year)
  max_year <- max(years_available, na.rm = TRUE)
  
  # Select years to plot
  if (max_year >= 200) {
    years_to_plot <- c(100, 150, 200)
  } else {
    # calculate middle year between 100 and max
    mid_year <- round(mean(c(100, max_year)))
    years_to_plot <- c(100, mid_year, max_year)
    # ensure only available years
    years_to_plot <- years_to_plot[years_to_plot %in% years_available]
  }
  
  # Filter to selected years
  if (is.data.table(plot_data)) {
    plot_data <- plot_data[Year %in% years_to_plot]
  } else {
    plot_data <- plot_data %>% filter(Year %in% years_to_plot)
  }
  
  # Store original data before any sampling (for NInd sampling later)
  # Create a deep copy that works for both data.table and data.frame
  if (is.data.table(plot_data)) {
    # For data.table, create a new data.table from the existing one
    plot_data_original <- as.data.table(as.data.frame(plot_data))
  } else {
    # For data.frame, assignment creates a copy
    plot_data_original <- plot_data
  }
  
  cat(sprintf("    Using data from replicate %d\n", 
              ifelse(!is.null(rep_number), rep_number, selected_rep)))
  
  # ============================================================================#
  # STEP 2: DETERMINE STUDY AREA FROM ACTUAL DATA ONLY
  # ============================================================================#
  
  cat("    Determining study area boundary from actual data...\n")
  
  # Use ONLY the actual data extent, ignore range_data
  # Get the full extent from ALL years in pop_data for context
  if (is.data.table(pop_data)) {
    all_data_extent <- data.frame(
      xmin = min(pop_data$x, na.rm = TRUE),
      xmax = max(pop_data$x, na.rm = TRUE),
      ymin = min(pop_data$y, na.rm = TRUE),
      ymax = max(pop_data$y, na.rm = TRUE)
    )
  } else {
    all_data_extent <- data.frame(
      xmin = min(pop_data$x, na.rm = TRUE),
      xmax = max(pop_data$x, na.rm = TRUE),
      ymin = min(pop_data$y, na.rm = TRUE),
      ymax = max(pop_data$y, na.rm = TRUE)
    )
  }
  
  # Get the extent of the filtered plot data
  plot_data_extent <- data.frame(
    xmin = min(plot_data$x, na.rm = TRUE),
    xmax = max(plot_data$x, na.rm = TRUE),
    ymin = min(plot_data$y, na.rm = TRUE),
    ymax = max(plot_data$y, na.rm = TRUE)
  )
  
  cat(sprintf("    All data extent: X[%.1f, %.1f], Y[%.1f, %.1f]\n", 
              all_data_extent$xmin, all_data_extent$xmax,
              all_data_extent$ymin, all_data_extent$ymax))
  
  cat(sprintf("    Plot data (years %d-%d) extent: X[%.1f, %.1f], Y[%.1f, %.1f]\n", 
              min(years_to_plot), max(years_to_plot),
              plot_data_extent$xmin, plot_data_extent$xmax,
              plot_data_extent$ymin, plot_data_extent$ymax))
  
  # Use the broader extent to show context but ensure we capture the actual data
  study_area_boundary <- data.frame(
    xmin = min(all_data_extent$xmin, plot_data_extent$xmin),
    xmax = max(all_data_extent$xmax, plot_data_extent$xmax),
    ymin = min(all_data_extent$ymin, plot_data_extent$ymin),
    ymax = max(all_data_extent$ymax, plot_data_extent$ymax)
  )
  
  # Check if the range is reasonable
  x_range <- study_area_boundary$xmax - study_area_boundary$xmin
  y_range <- study_area_boundary$ymax - study_area_boundary$ymin
  
  cat(sprintf("    Data range: X=%.1f units, Y=%.1f units\n", x_range, y_range))
  
  # If range is very small (all points in same area), add minimum buffer
  if (x_range < 10) {
    cat("    Warning: Very narrow X range, adding buffer\n")
    study_area_boundary$xmin <- study_area_boundary$xmin - 5
    study_area_boundary$xmax <- study_area_boundary$xmax + 5
  }
  
  if (y_range < 10) {
    cat("    Warning: Very narrow Y range, adding buffer\n")
    study_area_boundary$ymin <- study_area_boundary$ymin - 5
    study_area_boundary$ymax <- study_area_boundary$ymax + 5
  }
  
  # Add 10% buffer for visibility
  x_buffer <- (study_area_boundary$xmax - study_area_boundary$xmin) * 0.1
  y_buffer <- (study_area_boundary$ymax - study_area_boundary$ymin) * 0.1
  
  study_area_boundary$xmin <- study_area_boundary$xmin - x_buffer
  study_area_boundary$xmax <- study_area_boundary$xmax + x_buffer
  study_area_boundary$ymin <- study_area_boundary$ymin - y_buffer
  study_area_boundary$ymax <- study_area_boundary$ymax + y_buffer
  
  cat(sprintf("    Final study area with buffer: X[%.1f, %.1f], Y[%.1f, %.1f]\n", 
              study_area_boundary$xmin, study_area_boundary$xmax,
              study_area_boundary$ymin, study_area_boundary$ymax))
  
  # Create factor for facet labels
  plot_data$Year_Label <- factor(plot_data$Year, 
                                 levels = sort(unique(plot_data$Year)))
  plot_data_original$Year_Label <- factor(plot_data_original$Year,
                                          levels = sort(unique(plot_data_original$Year)))
  
  # ============================================================================#
  # STEP 3: CREATE PLOT FUNCTION WITH SAMPLING OPTION
  # ============================================================================#
  
  stage_plot_with_hist <- function(data, stage_col, stage_name, 
                                   study_area, sim_id, file_suffix,
                                   apply_sampling = FALSE) {  # Add sampling parameter
    
    # Apply 10% sampling if requested (for NInd total)
    if (apply_sampling) {
      set.seed(seed_value)  # Use seed for reproducibility
      n_rows <- nrow(data)
      sample_size <- ceiling(n_rows * 0.1)  # 10% sample
      sample_idx <- sample(n_rows, sample_size)
      
      if (is.data.table(data)) {
        data <- data[sample_idx]
      } else {
        data <- data[sample_idx, ]
      }
      cat(sprintf("      Sampled %d cells (10%%) from %d total for %s\n", 
                  sample_size, n_rows, stage_name))
    }
    
    # Filter to occupied cells only
    if (is.data.table(data)) {
      occupied_data <- data[get(stage_col) > 0]
    } else {
      occupied_data <- data %>% filter(!!sym(stage_col) > 0)
    }
    
    if (nrow(occupied_data) == 0) {
      cat(sprintf("      No occupied cells for %s\n", stage_name))
      return(NULL)
    }
    
    # Histogram data for annotation
    if (is.data.table(occupied_data)) {
      hist_data <- occupied_data[, .(
        mean_ind   = as.numeric(mean(get(stage_col), na.rm = TRUE)),
        median_ind = as.numeric(median(get(stage_col), na.rm = TRUE)),
        max_ind    = as.numeric(max(get(stage_col), na.rm = TRUE)),
        q75_ind    = as.numeric(quantile(get(stage_col), 0.75, na.rm = TRUE)),
        n_cells    = .N
      ), by = Year_Label]
    } else {
      hist_data <- occupied_data %>%
        group_by(Year_Label) %>%
        summarise(
          mean_ind   = as.numeric(mean(!!sym(stage_col), na.rm = TRUE)),
          median_ind = as.numeric(median(!!sym(stage_col), na.rm = TRUE)),
          max_ind    = as.numeric(max(!!sym(stage_col), na.rm = TRUE)),
          q75_ind    = as.numeric(quantile(!!sym(stage_col), 0.75, na.rm = TRUE)),
          n_cells    = n(),
          .groups = "drop"
        )
    }
    
    # --- Build plots ---
    p_map <- ggplot() +
      geom_rect(data = study_area,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = "grey98", colour = "grey70", size = 0.5, alpha = 1) +
      geom_point(data = occupied_data,
                 aes(x = x, y = y),
                 colour = "black", size = 1, alpha = 0.9, shape = 1) +
      facet_wrap(~ Year_Label, ncol = 5) +
      scale_x_continuous(labels = scales::comma,
                         limits = c(study_area$xmin, study_area$xmax)) +
      scale_y_continuous(labels = scales::comma,
                         limits = c(study_area$ymin, study_area$ymax)) +
      theme_void() + 
      theme(
        strip.text = element_blank(),
        panel.border = element_rect(colour = "grey60", fill = NA, size = 0.5)
      ) +
      coord_fixed(ratio = 1)
    
    p_hist <- ggplot(hist_data, aes(x = Year_Label)) +
      geom_col(aes(y = mean_ind), fill = "black", alpha = 0.5, width = 0.6) +
      geom_errorbar(aes(ymin = median_ind, ymax = q75_ind), 
                    width = 0.2, colour = "darkgrey", size = 0.5) +
      geom_point(aes(y = max_ind), colour = "darkred", size = 2, shape = 17) +
      labs(x = NULL, y = "Individuals per cell",
           subtitle = sprintf("Mean (bars), median-Q75 (lines), max (triangles)\nTotal cells: %d%s", 
                              nrow(occupied_data),
                              ifelse(apply_sampling, " (10% sample)", ""))) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.subtitle = element_text(size = 10, colour = "grey40")
      )
    
    # --- Combine into one gtable object ---
    combined_plot <- gridExtra::arrangeGrob(
      p_map, p_hist, ncol = 1, heights = c(3, 1)
    )
    
    # --- Save using helper ---
    base_filename <- file.path(spatial_dir, 
                               paste0("spatial_", stage_col, "_sim", sim_id, "_", file_suffix))
    
    combined_grob <- gridExtra::arrangeGrob(p_map, p_hist, ncol = 1, heights = c(3, 1))
    save_grob_multi_format(combined_grob, base_filename, width = 16, height = 10, dpi = 300, formats = formats)
    return(TRUE)
  } # End stage_plot_with_hist
  
  # ============================================================================#
  # STEP 4: CREATE PLOTS FOR EACH STAGE (NO SAMPLING)
  # ============================================================================#
  
  stages_to_plot <- list(
    list(col = "NInd_stage4", name = "Stage 4 (Adults)", file = "stage4"),
    list(col = "NInd_stage3", name = "Stage 3 (Sub-adults)", file = "stage3"),
    list(col = "NInd_stage2", name = "Stage 2 (Juveniles)", file = "stage2"),
    list(col = "NInd_stage1", name = "Stage 1 (Seedlings)", file = "stage1")
  )
  
  for (stage_info in stages_to_plot) {
    cat(sprintf("    Creating %s plot...\n", stage_info$name))
    
    tryCatch({
      stage_plot_with_hist(
        data = plot_data_original,  # Use original unsampled data
        stage_col = stage_info$col,
        stage_name = stage_info$name,
        study_area = study_area_boundary,
        sim_id = sim_id,
        file_suffix = stage_info$file,
        apply_sampling = FALSE  # No sampling for individual stages
      )
    }, error = function(e) {
      cat(sprintf("      Error: %s\n", e$message))
    })
  }
  
  # ============================================================================#
  # STEP 5: CREATE COMBINED PLOTS
  # ============================================================================#
  
  # Add NInd plot if column exists
  if ("NInd" %in% names(plot_data_original)) {
    cat("    Creating NInd (total individuals per cell) plot with 10% sampling...\n")
    
    tryCatch({
      stage_plot_with_hist(
        data = plot_data_original,
        stage_col = "NInd",
        stage_name = "Total Individuals per Cell (NInd)",
        study_area = study_area_boundary,
        sim_id = sim_id,
        file_suffix = "NInd",
        apply_sampling = TRUE  # Apply 10% sampling for NInd
      )
    }, error = function(e) {
      cat(sprintf("      Error in NInd plot: %s\n", e$message))
    })
  }
  
  # Total individuals (sum of stages)
  cat("    Creating total individuals (sum of stages) plot...\n")
  
  if (is.data.table(plot_data_original)) {
    plot_data_original[, NInd_total := as.numeric(NInd_stage1) + as.numeric(NInd_stage2) + 
                         as.numeric(NInd_stage3) + as.numeric(NInd_stage4)]
  } else {
    plot_data_original <- plot_data_original %>%
      mutate(NInd_total = as.numeric(NInd_stage1) + as.numeric(NInd_stage2) + 
               as.numeric(NInd_stage3) + as.numeric(NInd_stage4))
  }
  
  tryCatch({
    stage_plot_with_hist(
      data = plot_data_original,
      stage_col = "NInd_total",
      stage_name = "Total individuals (all stages summed)",
      study_area = study_area_boundary,
      sim_id = sim_id,
      file_suffix = "total",
      apply_sampling = FALSE  # No sampling for stage sums
    )
  }, error = function(e) {
    cat(sprintf("      Error in total plot: %s\n", e$message))
  })
  
  # Reproductive stages
  cat("    Creating reproductive stages plot...\n")
  
  if (is.data.table(plot_data_original)) {
    plot_data_original[, NInd_repro := as.numeric(NInd_stage3) + as.numeric(NInd_stage4)]
  } else {
    plot_data_original <- plot_data_original %>%
      mutate(NInd_repro = as.numeric(NInd_stage3) + as.numeric(NInd_stage4))
  }
  
  tryCatch({
    stage_plot_with_hist(
      data = plot_data_original,
      stage_col = "NInd_repro",
      stage_name = "Reproductive stages (3+4)",
      study_area = study_area_boundary,
      sim_id = sim_id,
      file_suffix = "reproductive",
      apply_sampling = FALSE  # No sampling for stage combinations
    )
  }, error = function(e) {
    cat(sprintf("      Error in reproductive plot: %s\n", e$message))
  })
  
  cat("    ✓ Spatial distribution plots created\n")
}  # End of function plot_spatial_distribution_by_stage


#==============================================================================#
#   5. Per-simulation plotting drivers ----
#==============================================================================#

plot_single_simulation <- function(sim_num, all_scenarios, masterfolder,
                                   # plot specific
                                   include_temporal_plots = TRUE,
                                   include_spatial = FALSE, 
                                   rep_number = NULL,
                                   seed_value = NULL,
                                   include_occupancy_rasters = TRUE,
                                   template_dir = "./Calo_RangeShiftR/00__data_preparation/StudyArea",
                                   include_summary_report = TRUE,
                                   use_dynamic_years = TRUE,
                                   fix_coordinates = TRUE,
                                   include_full_timeline = TRUE,
                                   include_post100_timeline = TRUE,
                                   selected_years = NULL,
                                   occupancy_threshold = 0.6
) {
  
  # Extract simulation parameters
  sim_params <- all_scenarios[sim_num, ]
  sim_id <- sim_params$sim_id
  batch_num <- sim_params$batch_num
  ssp <- sim_params$ssp
  emig_prob <- sim_params$emig_prob
  
  cat("\n[Worker", Sys.getpid(), "] Starting plots for Sim", sim_id, 
      "- Batch:", batch_num, "SSP:", ssp, "Emig:", emig_prob, "\n")
  
  # Start timing
  start_time <- Sys.time()
  
  # Force garbage collection at start
  gc(verbose = FALSE)
  
  # Use user-provided selected_years (default = NULL: defer to create_all_plots_with_rasters logic)
  plot_result <- tryCatch({
    results <- create_all_plots_with_rasters(
      sim_id = sim_id,
      batch_num = batch_num,
      ssp = ssp,
      emig_prob = emig_prob,
      masterfolder = masterfolder,
      selected_years = selected_years,
      include_temporal_plots = include_temporal_plots,
      include_spatial = include_spatial, 
      rep_num = rep_number,
      seed_value = seed_value,
      include_occupancy_rasters = include_occupancy_rasters,
      template_dir = template_dir,
      include_summary_report = include_summary_report,
      use_dynamic_years = use_dynamic_years,
      fix_coordinates = fix_coordinates,
      include_full_timeline = include_full_timeline,
      include_post100_timeline = include_post100_timeline,
      occupancy_threshold = occupancy_threshold
    )
    
    plot_dir <- file.path(masterfolder, "Plots", paste0("sim_", sim_id))
    raster_dir <- file.path(masterfolder, "Output_Maps", paste0("sim_", sim_id))
    
    n_plot_files <- length(list.files(plot_dir, pattern = "\\.(pdf|png)$", recursive = TRUE))
    n_raster_files <- if (dir.exists(raster_dir)) {
      length(list.files(raster_dir, pattern = "\\.tif$", recursive = TRUE))
    } else {
      0
    }
    
    list(
      sim_id = sim_id,
      status = "success",
      plot_status = results$status,
      n_plots = n_plot_files,
      n_rasters = n_raster_files,
      occupancy_rasters_created = include_occupancy_rasters && n_raster_files > 0
    )
    
  }, error = function(e) {
    cat("  ✗ Plotting failed for sim", sim_id, ":", e$message, "\n")
    list(
      sim_id = sim_id,
      status = "failed",
      error = e$message
    )
  })
  
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "secs")
  plot_result$runtime <- as.numeric(runtime)
  
  gc(verbose = FALSE)
  
  cat("[Worker", Sys.getpid(), "] Completed Sim", sim_id, 
      "- Status:", plot_result$status, 
      "- Runtime:", round(plot_result$runtime, 1), "secs")
  
  if (plot_result$status == "success") {
    cat(" - Plots:", plot_result$n_plots)
    if (include_occupancy_rasters) {
      cat(" - Rasters:", plot_result$n_rasters)
    }
  }
  cat("\n")
  
  return(plot_result)
}  # End of function plot_single_simulation


run_parallel_plotting <- function(sim_indices, all_scenarios, masterfolder, 
                                  ncores = 4, sims_per_chunk = 2,
                                  include_temporal_plots = TRUE,
                                  include_spatial = TRUE, 
                                  rep_number = NULL,
                                  seed_value = NULL,
                                  include_occupancy_rasters = TRUE,
                                  template_dir = "./Calo_RangeShiftR/00__data_preparation/StudyArea",
                                  include_summary_report = TRUE,
                                  use_dynamic_years = TRUE,
                                  fix_coordinates = TRUE,
                                  include_full_timeline = TRUE,
                                  include_post100_timeline = TRUE,
                                  selected_years = NULL,
                                  occupancy_threshold = 0.6) {
  
  foreach::registerDoSEQ() # reset any leftover doSNOW/doParallel registration
  
  cat("\n" , strrep("=", 60), "\n")
  cat("STARTING PARALLEL PLOTTING\n")
  cat("Simulations to process:", length(sim_indices), "\n")
  cat("Cores:", ncores, "\n")
  cat("Sims per chunk:", sims_per_chunk, "\n")
  cat("Include occupancy rasters:", include_occupancy_rasters, "\n")
  cat("Template directory:", template_dir, "\n")
  if (!is.null(rep_number)) cat("Replicate number:", rep_number, "\n") 
  if (!is.null(seed_value)) cat("Seed value:", seed_value, "\n")  
  cat(strrep("=", 60), "\n")
  
  all_results <- list()
  chunks <- split(sim_indices, ceiling(seq_along(sim_indices) / sims_per_chunk))
  
  for (chunk_idx in seq_along(chunks)) {
    chunk_sims <- chunks[[chunk_idx]]
    cat("\n--- Processing chunk", chunk_idx, "of", length(chunks), "---\n")
    cat("  Simulations in chunk:", paste(chunk_sims, collapse = ", "), "\n")
    
    cl <- parallel::makeCluster(ncores, outfile = "")  # type = PSOCK to avoid snowglobal error, outfile = "" shows worker output
    
    parallel::clusterEvalQ(cl, {
      options(warn = 1)
      gc(reset = TRUE)
      NULL
    })
    
    doParallel::registerDoParallel(cl)
    
    chunk_results <- foreach::foreach( 
      sim_num = chunk_sims,
      .packages = c("data.table", "dplyr", "ggplot2", "gridExtra", 
                    "scales", "viridis", "svglite", "terra", "grid", "parallel", 
                    "doParallel", "foreach", "tidyr", "gridExtra", "patchwork"
      ),
      .export = c(
        # Main functions
        "plot_single_simulation",
        "create_all_plots_with_rasters", 
        "create_occupancy_rasters",
        "add_occupancy_rasters_to_workflow",
        "fix_sim_coords",
        "create_raster_from_template_interpolated",
        "format_axis_labels",
        "clean_label",
        "calculate_extinction_years",
        "plot_stage_dynamics_dual",
        "plot_stage_composition",
        "plot_occupancy_dual",
        "plot_spatial_distribution_by_stage",
        "save_plot_multi_format",
        "save_grob_multi_format",
        "all_scenarios",
        "masterfolder"
      ),
      .errorhandling = "pass",
      .verbose = FALSE
    ) %dopar% {
      result <- plot_single_simulation(
        sim_num, 
        all_scenarios, 
        masterfolder,
        include_temporal_plots = include_temporal_plots,
        include_spatial = include_spatial, 
        rep_number = rep_number, 
        seed_value = seed_value, 
        include_occupancy_rasters = include_occupancy_rasters,
        template_dir = template_dir,
        fix_coordinates = fix_coordinates,
        include_summary_report = include_summary_report,
        use_dynamic_years = use_dynamic_years,
        include_full_timeline = include_full_timeline,
        include_post100_timeline = include_post100_timeline,
        selected_years = selected_years,
        occupancy_threshold = occupancy_threshold
      )
      gc(verbose = FALSE, full = TRUE)
      return(result)
    }
    cat("\n--- DEBUG: Checking chunk_results structure ---\n")
    cat("Length of chunk_results:", length(chunk_results), "\n")
    cat("Class of chunk_results:", class(chunk_results), "\n")
    
    for (i in seq_along(chunk_results)) {
      cat("\nResult", i, ":\n")
      cat("  Class:", class(chunk_results[[i]]), "\n")
      if (inherits(chunk_results[[i]], "error")) {
        cat("  ERROR:", chunk_results[[i]]$message, "\n")
      } else if (is.list(chunk_results[[i]])) {
        cat("  Names:", paste(names(chunk_results[[i]]), collapse = ", "), "\n")
        if (!is.null(chunk_results[[i]]$status)) {
          cat("  Status:", chunk_results[[i]]$status, "\n")
        }
        if (!is.null(chunk_results[[i]]$error)) {
          cat("  Error message:", chunk_results[[i]]$error, "\n")
        }
      } else {
        cat("  Unexpected type:", utils::str(chunk_results[[i]]), "\n")
      }
    }
    cat("--- END DEBUG ---\n\n")
    
    on.exit({ # this forces sequential to close any half-dead workers
      try(parallel::stopCluster(cl), silent = TRUE)
      foreach::registerDoSEQ()
    }, add = TRUE)
    
    #parallel::stopCluster(cl)
    #doParallel::registerDoSEQ()
    all_results <- c(all_results, chunk_results)
    
    successful <- sum(sapply(chunk_results, function(x) {
      if (is.list(x) && !is.null(x$status)) {
        x$status == "success"
      } else {
        FALSE
      }
    }))
    failed <- sum(sapply(chunk_results, function(x) {
      if (is.list(x) && !is.null(x$status)) {
        x$status == "failed"
      } else if (inherits(x, "error")) {
        TRUE
      } else {
        FALSE
      }
    }))
    
    if (include_occupancy_rasters) {
      total_rasters <- sum(sapply(chunk_results, function(x) {
        if (x$status == "success" && !is.null(x$n_rasters)) x$n_rasters else 0
      }))
      cat("  Chunk", chunk_idx, "complete:", successful, "successful,", failed, "failed")
      cat(" - Total rasters created:", total_rasters, "\n")
    } else {
      cat("  Chunk", chunk_idx, "complete:", successful, "successful,", failed, "failed\n")
    }
    gc(verbose = FALSE, full = TRUE)
    Sys.sleep(2)
  }
  
  final_successful <- sum(sapply(all_results, function(x) x$status == "success"))
  final_failed <- sum(sapply(all_results, function(x) x$status == "failed"))
  
  cat("\n" , strrep("=", 60), "\n")
  cat("PARALLEL PLOTTING COMPLETE\n")
  cat("Total successful:", final_successful, "\n")
  cat("Total failed:", final_failed, "\n")
  
  if (include_occupancy_rasters) {
    total_plots <- sum(sapply(all_results, function(x) {
      if (x$status == "success" && !is.null(x$n_plots)) x$n_plots else 0
    }))
    total_rasters <- sum(sapply(all_results, function(x) {
      if (x$status == "success" && !is.null(x$n_rasters)) x$n_rasters else 0
    }))
    raster_success <- sum(sapply(all_results, function(x) {
      if (x$status == "success" && !is.null(x$occupancy_rasters_created)) {
        x$occupancy_rasters_created
      } else {
        FALSE
      }
    }))
    cat("Total plot files:", total_plots, "\n")
    cat("Total raster files:", total_rasters, "\n")
    cat("Simulations with rasters:", raster_success, "/", final_successful, "\n")
  }
  cat(strrep("=", 60), "\n")
  return(all_results)
}  # End of function run_parallel_plotting


create_all_plots_with_rasters <- function(sim_id, batch_num, ssp, emig_prob, masterfolder,
                                          # Section controls
                                          include_temporal_plots = TRUE,       # Controls section 4 (temporal plots)
                                          include_spatial = FALSE,             # Controls section 5 (spatial distribution)
                                          include_occupancy_rasters = FALSE,   # Superseded by stage 03; occupancy rasters now produced there
                                          include_summary_report = TRUE,       # Controls section 6 (summary report)
                                          # Temporal plot options
                                          include_full_timeline = TRUE,        # 0-200 year plots
                                          include_post100_timeline = TRUE,     # 100-200 year plots
                                          # Year selection
                                          selected_years = NULL,   # Beginning / middle / end of simulation time line
                                          use_dynamic_years = TRUE,            # Adapt years based on extinction (beginning, middle, end)
                                          # Other parameters
                                          occupancy_threshold = 0.6,
                                          fix_coordinates = TRUE,
                                          spatial_sample_size = NULL,
                                          rep_num = NULL,
                                          seed_value = NULL,
                                          template_dir = "./Calo_RangeShiftR/00__data_preparation/StudyArea") {
  
  cat("\n=== PLOTTING FOR SIM", sim_id, "===\n")
  cat("  Configuration:\n")
  cat("    Temporal plots:", include_temporal_plots, "\n")
  cat("    Spatial plots:", include_spatial, "\n")
  cat("    Occupancy rasters:", include_occupancy_rasters, "\n")
  cat("    Summary report:", include_summary_report, "\n")
  if (use_dynamic_years) {
    cat("    Using dynamic years based on extinction status\n")
  }
  
  # Setup main plot directory with sim_id
  plot_dir <- file.path(masterfolder, "Plots", paste0("sim_", sim_id))
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Results collector
  results <- list(
    status = list(),
    outputs = list(),
    extinction_info = list()
  )
  
  # ============================================================================#
  ## 1: Read in output .txt files ----
  # ============================================================================#
  
  cat("\n--- Reading data files ---\n")
  
  pop_file <- file.path(masterfolder, "Outputs", 
                        paste0("Batch", batch_num, "_Sim", sim_id, "_Land1_Pop.txt"))
  
  range_file <- file.path(masterfolder, "Outputs", 
                          paste0("Batch", batch_num, "_Sim", sim_id, "_Land1_Range.txt"))
  
  # Load population data (shared between plotting and raster creation)
  if (file.exists(pop_file)) {
    file_size_mb <- file.info(pop_file)$size / 1024^2
    cat(sprintf("Pop file size: %.1f MB\n", file_size_mb))
    
    pop_data <- tryCatch({
      if (file_size_mb > 1000) {
        essential_cols <- c("Year", "Rep", "x", "y", "NInd", "NJuvs", 
                            "NInd_stage1", "NInd_stage2", "NInd_stage3", "NInd_stage4")
        header <- fread(pop_file, nrows = 0)
        cols_to_read <- intersect(essential_cols, names(header))
        data <- fread(pop_file, select = cols_to_read)
      } else {
        data <- fread(pop_file)
      }
      if (!is.data.table(data)) setDT(data)
      cat("✓ Pop data loaded:", format(nrow(data), big.mark = ","), "rows\n")
      data
    }, error = function(e) {
      cat("✗ Failed to read pop data:", e$message, "\n")
      NULL
    })
  } else {
    cat("✗ Pop file not found:", pop_file, "\n")
    pop_data <- NULL
  }
  
  # Load range data
  range_data <- tryCatch({
    if (file.exists(range_file)) {
      data <- fread(range_file)
      cat("✓ Range data loaded:", format(nrow(data), big.mark = ","), "rows\n")
      data
    } else {
      NULL
    }
  }, error = function(e) {
    NULL
  })
  
  if (is.null(pop_data)) {
    cat("Cannot proceed without population data\n")
    return(NULL)
  }
  
  # ============================================================================#
  ## 2: Calculate extinction years ----
  # ============================================================================#
  
  cat("\n--- analysing population persistence/extinction ---\n")
  
  extinction_years <- calculate_extinction_years(range_data, final_sim_year = 200)
  results$extinction_info$years <- extinction_years
  
  # Report findings
  if (!is.null(extinction_years$n_persisted) && extinction_years$n_persisted > 0) {
    cat(sprintf("  %d replicates persisted to year 200\n", extinction_years$n_persisted))
  }
  if (!is.null(extinction_years$n_extinct) && extinction_years$n_extinct > 0) {
    cat(sprintf("  %d replicates went extinct\n", extinction_years$n_extinct))
    if (!is.na(extinction_years$mean)) {
      cat(sprintf("  Mean extinction year (for extinct reps): %.1f ± %.1f (SE)\n", 
                  extinction_years$mean, extinction_years$se))
    }
  }
  
  # ============================================================================#
  ## 2b: Dynamic year selection (using extinction timing) ----
  # ============================================================================#
  # Compute dynamic selected_years if requested, using extinction_years$extinction_timing
  if (use_dynamic_years) {
    # Extinction timing is a vector: c(beginning, middle, end)
    # Fallback logic: always have c(0, middle, 200), with middle being extinction_years$extinction_timing[2] if available
    dyn_years <- c(0, 200)
    if (!is.null(extinction_years$extinction_timing) && length(extinction_years$extinction_timing) >= 3) {
      # Use provided dynamic timing (should be c(0, middle, 200))
      dyn_years <- extinction_years$extinction_timing
    } else if (!is.null(extinction_years$extinction_timing) && length(extinction_years$extinction_timing) == 1) {
      # Only middle available (rare), so [0, middle, 200]
      dyn_years <- c(0, extinction_years$extinction_timing[1], 200)
    }
    # Uniqueness and ordering
    selected_years <- sort(unique(dyn_years))
    cat("  Dynamic year selection for raster/spatial plotting:", paste(selected_years, collapse = ", "), "\n")
  }
  
  # ============================================================================#
  ## 3: Create occupancy rasters and fix coords ----
  # ============================================================================#
  
  if (include_occupancy_rasters) { # only run if defintion = TRUE
    cat("\n--- Creating occupancy rasters ---\n")
    
    # Load template for coordinate fixing if needed
    pop_data_for_rasters <- pop_data
    
    if (fix_coordinates) { # RS outputs a false origin - needs transforming to ESRI:102022 coords first
      cat("  Applying coordinate alignment fix...\n")
      
      # Load template
      template_100m <- terra::rast(file.path(template_dir, "CAZ_100m.tif")) #"calo_current_CAZ_100m.tif")) # shares study area shape, extent, CRS
      
      # Fix coordinates
      pop_data_for_rasters <- fix_sim_coords(
        pop_data, 
        template_100m, 
        sim_resolution = 100
      )
    } # End fix_coordinates
    
    # Create the output occupancy rasters
    raster_result <- add_occupancy_rasters_to_workflow(
      sim_id = sim_id,
      batch_num = batch_num,
      ssp = ssp,
      emig_prob = emig_prob,
      masterfolder = masterfolder,
      target_years = selected_years,
      occupancy_threshold = occupancy_threshold,
      template_dir = template_dir,
      pop_data = pop_data_for_rasters
    )
    
    # Add raster results to main results
    results$status$occupancy_rasters <- raster_result$status
    if (raster_result$status == "success") {
      results$outputs$occupancy_metadata <- raster_result$metadata
    }
  } # End include_occupancy_rasters
  
  # ============================================================================#
  ## 4: Create temporal plots ----
  # ============================================================================#
  
  if (include_temporal_plots) {
    cat("\n--- Creating temporal plots (dual time windows) ---\n")
    
    # Full timeline (0-200)
    cat("\n  Full timeline (0-200):\n")
    
    tryCatch({
      extinction_info_full <- plot_stage_dynamics_dual(pop_data, sim_id, plot_dir, start_year = 0)
      results$extinction_info$stages_full <- extinction_info_full
      results$status$stage_dynamics_full <- "success"
    }, error = function(e) {
      results$status$stage_dynamics_full <- paste("failed:", e$message)
    })
    
    tryCatch({
      plot_stage_composition(pop_data, sim_id, plot_dir, NULL, start_year = 0)
      results$status$stage_composition_full <- "success"
    }, error = function(e) {
      results$status$stage_composition_full <- paste("failed:", e$message)
    })
    
    tryCatch({
      plot_occupancy_dual(pop_data, range_data, sim_id, plot_dir, start_year = 0)
      results$status$occupancy_full <- "success"
    }, error = function(e) {
      results$status$occupancy_full <- paste("failed:", e$message)
    })
    
    # Post-100 timeline (100-200)
    if (include_post100_timeline) {
      cat("\n  Post-100 timeline (100-200):\n")
      
      tryCatch({
        extinction_info_post100 <- plot_stage_dynamics_dual(pop_data, sim_id, plot_dir, start_year = 100)
        results$extinction_info$stages_post100 <- extinction_info_post100
        results$status$stage_dynamics_post100 <- "success"
      }, error = function(e) {
        results$status$stage_dynamics_post100 <- paste("failed:", e$message)
      })
      
      tryCatch({
        plot_stage_composition(pop_data, sim_id, plot_dir, NULL, start_year = 100)
        results$status$stage_composition_post100 <- "success"
      }, error = function(e) {
        results$status$stage_composition_post100 <- paste("failed:", e$message)
      })
      
      tryCatch({
        plot_occupancy_dual(pop_data, range_data, sim_id, plot_dir, start_year = 100)
        results$status$occupancy_post100 <- "success"
      }, error = function(e) {
        results$status$occupancy_post100 <- paste("failed:", e$message)
      })
    }
  } else {
    cat("\n--- Skipping temporal plots (disabled) ---\n")
    results$status$temporal_plots <- "skipped"
  }
  
  # ============================================================================#
  ## 5: Spatial plots: occupancy per stage/year ----
  # ============================================================================#
  
  if (include_spatial && !is.null(pop_data)) {
    cat("\n--- Creating spatial distribution plots ---\n")
    
    tryCatch({
      plot_spatial_distribution_by_stage(
        pop_data = pop_data,
        range_data = range_data,
        sim_id = sim_id,
        plot_dir = plot_dir,
        years_to_plot = selected_years,
        rep_number = rep_num,
        seed_value = seed_value
      )
      results$status$spatial_distribution <- "success"
    }, error = function(e) {
      cat("    ✗ Spatial distribution plots failed:", e$message, "\n")
      results$status$spatial_distribution <- paste("failed:", e$message)
    })
  } else if (!include_spatial) {
    cat("\n--- Skipping spatial plots (disabled) ---\n")
    results$status$spatial_plots <- "skipped"
  }
  
  # ============================================================================#
  ## 6: Summary report for writing up -----
  # ============================================================================#
  
  if (include_summary_report) {
    cat("\n--- Writing summary report ---\n")
    
    report_file <- file.path(plot_dir, paste0("summary_report_sim", sim_id, ".txt"))
    
    # Create appropriate extinction text based on persistence/extinction
    extinction_section <- if (!is.null(extinction_years$by_rep)) {
      lines <- c()
      
      if (extinction_years$n_persisted == nrow(extinction_years$by_rep)) {
        lines <- c(lines, "  All populations persisted to year 200 (no extinction)")
      } else if (extinction_years$n_extinct == nrow(extinction_years$by_rep)) {
        lines <- c(lines, sprintf("  All populations went extinct (mean: %.1f ± %.1f years SE)", 
                                  extinction_years$mean, extinction_years$se))
      } else {
        lines <- c(lines, sprintf("  Mixed outcomes: %d persisted, %d extinct", 
                                  extinction_years$n_persisted, extinction_years$n_extinct))
        if (extinction_years$n_extinct > 0) {
          lines <- c(lines, sprintf("  Mean extinction (extinct reps only): %.1f ± %.1f years (SE)", 
                                    extinction_years$mean, extinction_years$se))
        }
      }
      
      # Report dynamic target years if available
      if (!is.null(extinction_years$extinction_timing)) {
        lines <- c(lines, sprintf("  Dynamic target years for raster/spatial plots: %s", 
                                  paste(extinction_years$extinction_timing, collapse = ", ")))
      }
      
      # Details by replicate
      lines <- c(lines, "  By replicate:")
      for (i in 1:nrow(extinction_years$by_rep)) {
        rep_data <- extinction_years$by_rep[i, ]
        if (rep_data$status == "persisted") {
          lines <- c(lines, sprintf("    Rep %d: Persisted (survived to year 200)", rep_data$Rep))
        } else {
          lines <- c(lines, sprintf("    Rep %d: Extinct at year %d", 
                                    rep_data$Rep, rep_data$extinction_year))
        }
      }
      lines
    } else {
      "  No extinction data available"
    }
    
    # Count outputs
    all_plots <- list.files(plot_dir, pattern = "\\.(pdf|png)$", recursive = TRUE)
    plot_counts <- list(
      total = length(all_plots),
      pdfs = sum(grepl("\\.pdf$", all_plots)),
      pngs = sum(grepl("\\.png$", all_plots))
      # Optionally add counts for full_timeline/post100/spatial by filename if needed
    )
    
    # Count raster files if created
    raster_counts <- list(total = 0, files_100m = 0, files_1km = 0)
    if (include_occupancy_rasters) {
      raster_dir <- file.path(masterfolder, "Output_Maps", paste0("sim_", sim_id))
      if (dir.exists(raster_dir)) {
        all_rasters <- list.files(raster_dir, pattern = "\\.tif$", recursive = TRUE)
        raster_counts <- list(
          total = length(all_rasters),
          files_100m = sum(grepl("_100m_", all_rasters)),
          files_1km = sum(grepl("_1km_", all_rasters))
        )
      }
    }
    
    # Build report
    report <- c(
      paste("=" , strrep("=", 60)),
      paste("Simulation", sim_id, "Plotting summary"),
      paste("Generated:", Sys.time()),
      paste("=" , strrep("=", 60)),
      "",
      "Simulation parameters:",
      paste("  Batch:", batch_num),
      paste("  Climate Scenario:", ssp),
      paste("  Emigration Probability:", emig_prob),
      "",
      "Population persistence/extinction:",
      extinction_section,
      "",
      "Data summary:",
      paste("  Population data rows:", format(nrow(pop_data), big.mark = ",")),
      paste("  Range data rows:", ifelse(!is.null(range_data), 
                                         format(nrow(range_data), big.mark = ","), 
                                         "Not available")),
      "",
      "Plot generation status:",
      sapply(names(results$status), function(x) {
        status <- results$status[[x]]
        symbol <- ifelse(grepl("success", status), "✓", "✗")
        paste("  ", symbol, clean_label(x), ":", status)
      }),
      "",
      "Output summary:",
      paste("  Total plots created:", plot_counts$total),
      paste("    PDF files:", plot_counts$pdfs),
      paste("    PNG files:", plot_counts$pngs),
      if (include_occupancy_rasters) c(
        paste("  Total raster files:", raster_counts$total),
        paste("    100m resolution:", raster_counts$files_100m),
        paste("    1km resolution:", raster_counts$files_1km)
      ) else NULL,
      "",
      "File locations:",
      paste("  Main directory:", plot_dir),
      if (include_occupancy_rasters) paste("  Raster directory:", file.path(masterfolder, "Output_Maps", paste0("sim_", sim_id))) else NULL
    )
    
    writeLines(report, report_file)
    cat("  ✓ summary report written\n")
    
  } else {
    cat("\n--- Skipping summary report (disabled) ---\n")
    results$status$summary_report <- "skipped"
  }
  
  # Clean up memory
  if (nrow(pop_data) > 1e6) {
    rm(pop_data)
    gc(verbose = FALSE)
  }
  
  cat("\n=== PLOTTING COMPLETE FOR SIM", sim_id, "===\n")
  
  return(results)
}  # End of function create_all_plots_with_rasters


#==============================================================================#
#   6. Legacy per-sim occupancy rasters (optional; see note) ----
#==============================================================================#

# LEGACY / OPTIONAL: per-simulation occupancy rasters produced during plotting.
# Superseded for publication by 03__occupancy_probability.R (functions_occupancy.R).
# Retained only because create_all_plots_with_rasters() can call them when
# include_occupancy_rasters = TRUE (off by default). Safe to delete if unused.
create_raster_from_template_interpolated <- function(template, values, coords, name,
                                                     interpolate_gaps = TRUE) {
  .require_ns("terra")
  stopifnot(inherits(template, "SpatRaster"))
  stopifnot(length(values) == nrow(coords))
  stopifnot(all(c("x","y") %in% names(coords) | ncol(coords) >= 2))
  
  r <- template
  terra::values(r) <- NA
  if (!is.null(name)) names(r) <- name
  
  print("DEBUG: class of values argument to rasterizer:")
  print(class(values))
  print("DEBUG: str of values argument to rasterizer:")
  print(str(values))
  print("DEBUG: head of values argument to rasterizer:")
  print(head(values))
  
  xy <- if (all(c("x","y") %in% names(coords))) as.matrix(coords[, c("x","y")]) else as.matrix(coords[,1:2])
  cells <- terra::cellFromXY(r, xy)
  valid <- !is.na(cells)
  
  if (any(valid)) {
    df <- data.frame(cell = cells[valid], value = values[valid])
    print("DEBUG: class of df$value before aggregate:")
    print(class(df$value))
    agg <- stats::aggregate(value ~ cell, data = df, FUN = mean, na.rm = TRUE)
    print("DEBUG: class of agg$value after aggregate:")
    print(class(agg$value))
    print("DEBUG: str of agg$value after aggregate:")
    print(str(agg$value))
    # Force to numeric vector if not already
    agg$value <- as.numeric(agg$value)
    r[agg$cell] <- agg$value
  }
  
  before_na <- sum(is.na(terra::values(r)))
  if (isTRUE(interpolate_gaps)) {
    na_mask <- is.na(r)
    r_focal <- terra::focal(r, w = 3, fun = function(x, ...) mean(x, na.rm = TRUE))
    r[na_mask] <- r_focal[na_mask]
  }
  after_na <- sum(is.na(terra::values(r)))
  attr(r, "filled_cells_by_focal") <- before_na - after_na
  r
}  # End of function create_raster_from_template_interpolated


fix_sim_coords <- function(pop_data, template_raster, sim_resolution = 100) {
  
  cat("Fixing simulation coordinates...\n")
  
  # Get template raster properties
  template_ext <- terra::ext(template_raster)
  #template_res <- terra::res(template_raster)[1]  
  
  # Get the origin of the template raster (lower-left corner)
  template_origin_x <- template_ext$xmin
  template_origin_y <- template_ext$ymin
  
  template_origin_x <- terra::xmin(template_raster)
  template_origin_y <- terra::ymin(template_raster)
  template_res <- terra::res(template_raster)[1] # Assumes square cells
  
  # Important change: Snaps coordinates to cell centre instead of edges
  # This ensures every coordinate maps to exactly one cell
  
  # Calculate cell centres for the simulation coordinates
  # First, find which cell each coordinate belongs to
  cell_col <- floor((pop_data$x - template_origin_x) / template_res)
  cell_row <- floor((pop_data$y - template_origin_y) / template_res)
  
  # Now calculate the centre of each cell (not the edge)
  # This is the key fix - we ensure coordinates are at cell centres
  pop_data$x_aligned <- template_origin_x + (cell_col * template_res) + (template_res / 2)
  pop_data$y_aligned <- template_origin_y + (cell_row * template_res) + (template_res / 2)
  
  # Optional: Check for any coordinates outside the template extent
  # and handle them appropriately
  outside_extent <- pop_data$x_aligned < template_ext$xmin | 
    pop_data$x_aligned > template_ext$xmax |
    pop_data$y_aligned < template_ext$ymin | 
    pop_data$y_aligned > template_ext$ymax
  
  if (sum(outside_extent) > 0) {
    cat("  Warning:", sum(outside_extent), "coordinates fall outside template extent\n")
    # You can choose to remove these or handle differently
    # pop_data <- pop_data[!outside_extent, ]
  }
  
  # Replace original coordinates with aligned ones
  pop_data$x <- pop_data$x_aligned
  pop_data$y <- pop_data$y_aligned
  
  # Clean up temporary columns
  pop_data$x_aligned <- NULL
  pop_data$y_aligned <- NULL
  
  cat("  ✓ Aligned", nrow(pop_data), "coordinates to cell centres\n")
  
  return(pop_data)
}  # End of function fix_sim_coords


create_occupancy_rasters <- function(sim_id, batch_num, masterfolder, 
                                     target_years = c(100, 150, 200),
                                     occupancy_threshold = 0.7,  # 7 out of 10 replicates
                                     output_dir = NULL,
                                     template_paths = NULL,
                                     template_dir = "./Calo_RangeShiftR/00__data_preparation/StudyArea",
                                     pop_data = NULL, # Accept pre-loaded data
                                     fix_grid_lines = TRUE,  # as there were systematic grid lines being mapped as 'no data'
                                     extinction_info = NULL, # passes extinction info to map start-middle-end dynamically
                                     refugia_year = NULL, # highlights last remaining populations
                                     should_highlight = FALSE # only if checked
) { 
  
  cat("\n=== CREATING OCCUPANCY RASTERS FOR SIM", sim_id, "===\n")
  
  # Report extinction context
  if (!is.null(extinction_info)) {
    if (extinction_info$n_extinct > 0 && extinction_info$n_persisted == 0) {
      cat("  Context: All populations went extinct (mean year:", 
          round(extinction_info$mean), ")\n")
      cat("  Creating rasters to capture extinction dynamics\n")
    } else if (extinction_info$n_extinct > 0) {
      cat("  Context: Mixed outcomes -", extinction_info$n_extinct, "extinct,", 
          extinction_info$n_persisted, "persisted\n")
    }
  }
  
  # Setup output directory
  output_dir <- file.path(masterfolder, "Output_Maps", paste0("sim_", sim_id))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load templates
  cat("Loading template rasters...\n")
  template_100m_file <- file.path(template_dir, "CAZ_100m.tif")
  template_1km_file <- file.path(template_dir, "CAZ_1km_aggregated.tif")
  
  # Load the actual rasters
  template_100m <- terra::rast(template_100m_file)
  template_1km <- terra::rast(template_1km_file)
  
  cat("✓ Template rasters loaded:\n")
  cat("  100m template:", dim(template_100m)[1], "x", dim(template_100m)[2], "cells\n")
  cat("  1km template:", dim(template_1km)[1], "x", dim(template_1km)[2], "cells\n")
  
  
  # Load population data only if not provided
  if (is.null(pop_data)) {
    cat("Reading population data from file...\n")
    pop_file <- file.path(masterfolder, "Outputs", 
                          paste0("Batch", batch_num, "_Sim", sim_id, "_Land1_Pop.txt"))
    
    essential_cols <- c("Year", "Rep", "x", "y", "NInd")
    pop_data <- fread(pop_file, select = essential_cols)
    pop_data <- pop_data[Year %in% target_years]
  } else {
    if (!is.data.table(pop_data)) setDT(pop_data)
    pop_data <- pop_data[Year %in% target_years]
  }
  
  # Process each year
  results <- list()
  
  for (year in target_years) {
    cat("\nProcessing Year", year, "...\n")
    
    # Check if this is the refugia year
    is_refugia_year <- should_highlight && !is.null(refugia_year) && year == refugia_year
    if (is_refugia_year) {
      cat("  >>> Refugia year - creating extra visuals <<<\n")
    }
    
    # Extract data for current year
    year_data <- pop_data[Year == year]
    
    if (nrow(year_data) == 0) {
      cat("  ✗ No data for year", year, "\n")
      next
    }
    
    # Create occupancy indicator (1 if occupied, 0 if not)
    year_data[, occupied := as.integer(NInd > 0)]
    
    # Round coordinates to ensure exact cell alignment
    # This is a secondary safety measure to avoid grid lines appearing in output rasters
    template_res <- terra::res(template_100m)[1]
    year_data[, `:=`(
      x = round(x / template_res) * template_res,
      y = round(y / template_res) * template_res
    )]
    
    # Calculate occupancy statistics per cell across replicates
    occupancy_stats <- year_data[, .(
      n_reps = .N,                                                     # Number of replicates with data
      n_occupied = sum(occupied),                                      # Number of times occupied
      prob_occupied = mean(occupied),                                  # Probability of occupancy
      se_occupied = sqrt(mean(occupied) * (1 - mean(occupied)) / .N),  # Standard error for proportion
      mean_inds = mean(NInd),                                          # Mean individuals when present
      sd_inds = sd(NInd),                                              # SD of individuals
      max_inds = max(NInd)                                             # For refugia identification
    ), by = .(x, y)]
    
    # Create binary occupancy based on threshold
    occupancy_stats[, binary_occupied := as.integer(prob_occupied >= occupancy_threshold)]
    
    # If refugia year, identify core refugia areas
    if (is_refugia_year) {
      # Identify refugia: cells with high occupancy probability and population
      refugia_threshold <- quantile(occupancy_stats[prob_occupied > 0]$prob_occupied, 0.75, na.rm = TRUE)
      occupancy_stats[, is_refugia := as.integer(prob_occupied >= refugia_threshold)]
      
      cat("    Refugia cells identified:", sum(occupancy_stats$is_refugia), "\n")
      cat("    Refugia threshold:", round(refugia_threshold, 3), "\n")
    }
    
    # Create rasters with interpolation function at both resolutions
    tryCatch({
      # Outputs at 100m spatial resolution
      cat("  Creating 100m resolution rasters with gap filling...\n")
      
      # Get all study area cells from template
      all_xy <- as.data.table(terra::xyFromCell(template_100m, 1:terra::ncell(template_100m)))
      setnames(all_xy, c("x", "y"))
      
      # Occupancy stats for all study area cells via merge
      grid_stats <- merge(all_xy, occupancy_stats, by = c("x", "y"), all.x = TRUE)
      
      # NA's should become absences (never occupied)
      grid_stats[is.na(prob_occupied), prob_occupied := 0]
      grid_stats[is.na(binary_occupied), binary_occupied := 0]
      
      prob_100m <- create_raster_from_template_interpolated(
        template_100m,
        as.numeric(grid_stats$prob_occupied),
        grid_stats[, .(x, y)],
        paste0("occupancy_prob_100m_year", year),
        interpolate_gaps = FALSE
      )
      
      # prob_100m <- create_raster_from_template_interpolated_debug( # debug version
      #   template_100m, 
      #   occupancy_stats$prob_occupied,
      #   occupancy_stats[, .(x, y)],
      #   paste0("occupancy_prob_100m_year", year),
      #   interpolate_gaps = TRUE  # Try FALSE if it fails
      # )
      
      binary_100m <- create_raster_from_template_interpolated(
        template_100m,
        occupancy_stats$binary_occupied,
        occupancy_stats[, .(x, y)],
        paste0("occupancy_binary_100m_year", year),
        interpolate_gaps = FALSE  # Don't interpolate binary data
      )
      gc()
      
      se_100m <- create_raster_from_template_interpolated(
        template_100m,
        occupancy_stats$se_occupied,
        occupancy_stats[, .(x, y)],
        paste0("occupancy_se_100m_year", year),
        interpolate_gaps = TRUE
      )
      gc()
      # Create refugia raster if TRUE
      if (is_refugia_year) {
        refugia_100m <- create_raster_from_template_interpolated(
          template_100m,
          occupancy_stats$is_refugia,
          occupancy_stats[, .(x, y)],
          paste0("refugia_100m_year", year),
          interpolate_gaps = FALSE
        )
        
        # Save refugia raster
        refugia_file <- file.path(output_dir, 
                                  paste0("refugia_100m_sim", sim_id, "_year", year, ".tif"))
        terra::writeRaster(refugia_100m, refugia_file, overwrite = TRUE, datatype = "INT1U")
        cat("    ✓ Saved refugia raster:", basename(refugia_file), "\n")
      }
      
      # 1km version outputs
      cat("  Creating 1km resolution rasters...\n")
      
      # For 1km rasters, we need to aggregate the data spatially
      # Convert coordinates to 1km grid cells and aggregate
      coords_1km <- terra::cellFromXY(template_1km, as.matrix(occupancy_stats[, .(x, y)]))
      gc()
      # Add 1km cell IDs to occupancy stats
      occupancy_stats[, cell_1km := coords_1km]
      gc()
      # Aggregate to 1km resolution - calculate mean within each 1km cell
      occupancy_1km <- occupancy_stats[!is.na(cell_1km), .(
        prob_occupied_1km = mean(prob_occupied, na.rm = TRUE),
        se_occupied_1km = sqrt(mean(prob_occupied * (1 - prob_occupied) / n_reps, na.rm = TRUE)),
        binary_occupied_1km = as.integer(mean(prob_occupied, na.rm = TRUE) >= occupancy_threshold),
        n_cells_100m = .N  # Number of 100m cells contributing to each 1km cell
      ), by = cell_1km]
      
      # Create 1km rasters
      prob_1km <- template_1km
      terra::values(prob_1km) <- NA
      names(prob_1km) <- paste0("occupancy_prob_1km_year", year)
      prob_1km[occupancy_1km$cell_1km] <- occupancy_1km$prob_occupied_1km
      
      gc()
      binary_1km <- template_1km
      terra::values(binary_1km) <- NA
      names(binary_1km) <- paste0("occupancy_binary_1km_year", year)
      binary_1km[occupancy_1km$cell_1km] <- occupancy_1km$binary_occupied_1km
      gc()
      se_1km <- template_1km
      terra::values(se_1km) <- NA
      names(se_1km) <- paste0("occupancy_se_1km_year", year)
      se_1km[occupancy_1km$cell_1km] <- occupancy_1km$se_occupied_1km
      
      # Quality check for the 'grid lines' 'no data' issue
      na_count <- sum(is.na(terra::values(prob_100m)))
      total_cells <- terra::ncell(prob_100m)
      na_percent <- (na_count / total_cells) * 100
      
      if (na_percent > 5) {
        cat("  Warning: ", round(na_percent, 2), "% of cells are NA\n")
        cat("    Consider checking coordinate alignment\n")
      } else {
        cat("  ✓ Raster quality check passed (", round(na_percent, 2), "% NA cells)\n")
      }
      
      # Save rasters with resolution indicators
      # 100m resolution files
      prob_100m_file <- file.path(output_dir, paste0("occupancy_probability_100m_sim", sim_id, "_year", year, ".tif"))
      binary_100m_file <- file.path(output_dir, paste0("occupancy_binary_100m_sim", sim_id, "_year", year, ".tif"))
      se_100m_file <- file.path(output_dir, paste0("occupancy_se_100m_sim", sim_id, "_year", year, ".tif"))
      
      # 1km resolution files
      prob_1km_file <- file.path(output_dir, paste0("occupancy_probability_1km_sim", sim_id, "_year", year, ".tif"))
      binary_1km_file <- file.path(output_dir, paste0("occupancy_binary_1km_sim", sim_id, "_year", year, ".tif"))
      se_1km_file <- file.path(output_dir, paste0("occupancy_se_1km_sim", sim_id, "_year", year, ".tif"))
      
      # Write 100m rasters
      terra::writeRaster(prob_100m, prob_100m_file, overwrite = TRUE)
      terra::writeRaster(binary_100m, binary_100m_file, overwrite = TRUE, datatype = "INT1U")
      terra::writeRaster(se_100m, se_100m_file, overwrite = TRUE)
      
      # Write 1km rasters
      terra::writeRaster(prob_1km, prob_1km_file, overwrite = TRUE)
      terra::writeRaster(binary_1km, binary_1km_file, overwrite = TRUE, datatype = "INT1U")
      terra::writeRaster(se_1km, se_1km_file, overwrite = TRUE)
      
      gc()
      cat("  ✓ Saved rasters:\n")
      cat("    100m Resolution:\n")
      cat("      Probability:", basename(prob_100m_file), "\n")
      cat("      Binary:", basename(binary_100m_file), "\n") 
      cat("      Standard Error:", basename(se_100m_file), "\n")
      cat("    1km Resolution:\n")
      cat("      Probability:", basename(prob_1km_file), "\n")
      cat("      Binary:", basename(binary_1km_file), "\n")
      cat("      Standard Error:", basename(se_1km_file), "\n")
      
      # Store results
      results[[as.character(year)]] <- list(
        year = year,
        n_cells_100m = nrow(occupancy_stats),
        n_occupied_100m = sum(occupancy_stats$binary_occupied),
        mean_prob_100m = mean(occupancy_stats$prob_occupied),
        n_cells_1km = nrow(occupancy_1km),
        n_occupied_1km = sum(occupancy_1km$binary_occupied_1km),
        mean_prob_1km = mean(occupancy_1km$prob_occupied_1km),
        files = list(
          probability_100m = prob_100m_file,
          binary_100m = binary_100m_file,
          se_100m = se_100m_file,
          probability_1km = prob_1km_file,
          binary_1km = binary_1km_file,
          se_1km = se_1km_file,
          is_refugia_year = is_refugia_year,
          n_refugia_cells = if (is_refugia_year) sum(occupancy_stats$is_refugia) else NA
        )
      )
      
    }, error = function(e) {
      cat("  ✗ Failed to create rasters for year", year, ":", e$message, "\n")
    })
    
    # Clean up year data
    rm(year_data, occupancy_stats)
    if (exists("occupancy_1km")) rm(occupancy_1km)
    gc(verbose = FALSE)
  }
  
  # Save summary metadata
  metadata <- list(
    sim_id = sim_id,
    batch_num = batch_num,
    target_years = target_years,
    occupancy_threshold = occupancy_threshold,
    extinction_info = extinction_info,
    refugia_highlighted = should_highlight,
    refugia_year = refugia_year,
    n_replicates = length(unique(pop_data$Rep)),
    template_files = list(
      template_100m = template_100m_file,
      template_1km = template_1km_file
    ),
    extent_100m = as.vector(terra::ext(template_100m)),
    extent_1km = as.vector(terra::ext(template_1km)),
    resolution_100m = terra::res(template_100m),
    resolution_1km = terra::res(template_1km),
    processing_date = Sys.time(),
    results = results
  )
  
  metadata_file <- file.path(output_dir, paste0("occupancy_metadata_sim", sim_id, ".rds"))
  saveRDS(metadata, metadata_file)
  
  cat("\n✓ Occupancy raster creation complete for sim", sim_id, "\n")
  if (should_highlight) {
    cat("✓ Refugia analysis included for extinction dynamics\n")
  }
  
  return(metadata)
}  # End of function create_occupancy_rasters


add_occupancy_rasters_to_workflow <- function(sim_id, batch_num, ssp, emig_prob, masterfolder,
                                              target_years = c(100, 150, 200),
                                              occupancy_threshold = 0.6,
                                              template_dir = "./Calo_RangeShiftR/00__data_preparation/StudyArea",
                                              pop_data = NULL, # Accept pre-loaded data
                                              extinction_info = NULL,  # takes extinction information
                                              highlight_refugia = TRUE) {  # emphasise final locations prior to extinction
  
  cat("\n--- Adding occupancy rasters to workflow ---\n")
  should_highlight <- FALSE
  refugia_year <- NULL
  
  if (highlight_refugia && !is.null(extinction_info)) {
    # Check if any populations went extinct
    if (extinction_info$n_extinct > 0) {
      should_highlight <- TRUE
      # Identify the refugia year (latest year in target_years)
      refugia_year <- max(target_years)
      cat("  Will highlight refugial zones at year", refugia_year, "\n")
    }
  }
  
  result <- tryCatch({
    metadata <- create_occupancy_rasters(
      sim_id = sim_id,
      batch_num = batch_num, 
      masterfolder = masterfolder,
      target_years = target_years,
      occupancy_threshold = occupancy_threshold,
      template_dir = template_dir,
      pop_data = pop_data,  # Pass pre-loaded data
      extinction_info = extinction_info,
      refugia_year = refugia_year,
      should_highlight = should_highlight
    )
    
    if (!is.null(metadata)) {
      list(status = "success", metadata = metadata)
    } else {
      list(status = "failed", error = "Failed to create rasters")
    }
    
  }, error = function(e) {
    cat("✗ Occupancy raster creation failed:", e$message, "\n")
    list(status = "failed", error = e$message)
  })
  
  return(result)
}  # End of function add_occupancy_rasters_to_workflow


#==============================================================================#
#                        ----  End of functions_results.R ----
#==============================================================================#
