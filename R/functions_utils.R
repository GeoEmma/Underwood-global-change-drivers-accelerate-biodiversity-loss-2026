# # ---
# title: "functions_utils.R"
# manuscript: "Title: Climate change, defaunation and disease interact to accelerate biodiversity loss in a tropical hotspot"
# corresponding_author: "Emma L Underwood"
# coauthors: "Kerry A Brown, Rebekka Allgayer, Mark Mulligan, Nigel Walford, Jette Wolff"
# last update: "2026-08-19"
# output: Generic helpers used across the workflow (package loading, memory
#         management, plot saving and legends).
# requires: base R; ggplot2/grid for the legend helpers  (loaded centrally in config.R)
# ---

# Functions are grouped by pipeline stage. Bodies are unchanged from the
# original working code; only their location, headers and section markers
# have been standardised for the public repository.


#==============================================================================#
#   1. Package loading ----
#==============================================================================#

# Kept for reference; config.R is the primary package loader.
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos = 'http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}  # End of function install.load.package


#==============================================================================#
#   2. Fonts (optional) ----
#==============================================================================#

# Not called on load. Call explicitly from a driver if the Callisto MT
# manuscript font is installed; otherwise the fallback (serif) is used.
setup_plot_fonts <- function(preferred_font = "Callisto MT", fallback_font = "serif") {
  
  # Try to load extrafont
  font_to_use <- fallback_font
  
  tryCatch({
    # Load extrafont if available
    if (!require("extrafont", quietly = TRUE)) {
      install.packages("extrafont")
      library(extrafont)
    }
    
    # Import fonts (only needs to be done once per system)
    if (length(fonts()) == 0) {
      cat("Importing system fonts (this may take a minute)...\n")
      font_import(prompt = FALSE)
      loadfonts(quiet = TRUE)
    }
    
    # Check if preferred font is available
    available_fonts <- fonts()
    if (preferred_font %in% available_fonts) {
      font_to_use <- preferred_font
      cat("Using font:", preferred_font, "\n")
    } else {
      cat("Font", preferred_font, "not found. Available fonts similar to Callisto:\n")
      
      # Look for similar fonts
      callisto_like <- c("Callisto MT", "Calisto MT", "Georgia", "Cambria", 
                         "Book Antiqua", "Palatino", "Palatino Linotype")
      available_similar <- callisto_like[callisto_like %in% available_fonts]
      
      if (length(available_similar) > 0) {
        font_to_use <- available_similar[1]
        cat("Using similar font:", font_to_use, "\n")
      } else {
        cat("Using fallback font:", fallback_font, "\n")
      }
    }
  }, error = function(e) {
    cat("Could not load custom fonts, using:", fallback_font, "\n")
  })
  
  return(font_to_use)
}  # End of function setup_plot_fonts


#==============================================================================#
#   3. Memory and execution helpers ----
#==============================================================================#

check_memory <- function(label, force_gc = TRUE) {
  if (force_gc) {
    gc(verbose = FALSE)
  }
  mem_info <- gc()
  used_mb <- sum(mem_info[,2])
  max_mb <- sum(mem_info[,6])
  cat(sprintf("Memory [%s]: Used %.1f MB, Max %.1f MB\n", label, used_mb, max_mb))
  invisible(list(used = used_mb, max = max_mb))
}  # End of function check_memory


cleanup_objects <- function(..., force_gc = TRUE, aggressive = TRUE) {
  objects_to_remove <- as.character(substitute(list(...)))[-1]
  existing_objects <- objects_to_remove[sapply(objects_to_remove, exists, 
                                               envir = parent.frame())]
  
  if (length(existing_objects) > 0) {
    # Remove objects
    rm(list = existing_objects, envir = parent.frame())
    
    if (force_gc) {
      if (aggressive) {
        # Multiple gc passes for thorough cleanup
        for(i in 1:3) {
          invisible(gc(verbose = FALSE, full = TRUE))
        }
      } else {
        invisible(gc(verbose = FALSE))
      }
    }
    
    # Optional: Clear tempdir periodically
    if (aggressive) {
      temp_files <- list.files(tempdir(), full.names = TRUE)
      if (length(temp_files) > 100) {  # Threshold
        unlink(temp_files[1:50])  # Remove oldest files
      }
    }
  }
}  # End of function cleanup_objects


safe_execute <- function(code, description) {
  tryCatch({
    result <- code
    return(result)
  }, error = function(e) {
    cat("Error in", description, ":", e$message, "\n")
    return(NULL)
  })
}  # End of function safe_execute


.require_ns <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' must be installed to use this function.")
  }
}  # End of function .require_ns


#==============================================================================#
#   4. Labels and plot saving ----
#==============================================================================#

format_axis_labels <- function(x) {
  format(x, scientific = FALSE, big.mark = ",", trim = TRUE)
}  # End of function format_axis_labels


clean_label <- function(text) {
  gsub("_", " ", text)
}  # End of function clean_label


save_plot_multi_format <- function(plot_function, base_filename, 
                                   width = 8, height = 6, dpi = 300,
                                   formats = c("pdf", "png", "svg")) {
  
  # Save in each requested format
  for (format in formats) {
    if (format == "pdf") {
      pdf(paste0(base_filename, ".pdf"), width = width, height = height)
      plot_function()
      dev.off()
    } else if (format == "png") {
      png(paste0(base_filename, ".png"), 
          width = width * dpi, height = height * dpi, type = "cairo", res = dpi) # add cairo as X11 graphics are no longer working
      plot_function()
      dev.off()
    } else if (format == "svg") {
      svglite(paste0(base_filename, ".svg"), width = width, height = height)
      plot_function()
      dev.off()
    }
  } # End save_plot_multi_format
  
  # Report what was created
  cat("  Saved:", paste(paste0(base_filename, ".", formats), collapse = ", "), "\n")
}  # End of function save_plot_multi_format


save_grob_multi_format <- function(grob_object, base_filename,
                                   width = 8, height = 6, dpi = 300,
                                   formats = c("pdf", "png", "svg")) {
  
  for (format in formats) {
    if (format == "pdf") {
      pdf(paste0(base_filename, ".pdf"), width = width, height = height)
      grid::grid.draw(grob_object)
      dev.off()
    } else if (format == "png") {
      png(paste0(base_filename, ".png"), 
          width = width * dpi, height = height * dpi, res = dpi,
          type = "cairo")
      grid::grid.draw(grob_object)
      dev.off()
    } else if (format == "svg") {
      if (requireNamespace("svglite", quietly = TRUE)) {
        svglite::svglite(paste0(base_filename, ".svg"), width = width, height = height)
        grid::grid.draw(grob_object)
        dev.off()
      } else {
        cat("    svglite package not available, skipping SVG format\n")
      }
    }
  }
  
  cat(sprintf("    Files saved: %s\n", 
              paste(paste0(basename(base_filename), ".", formats), collapse = ", ")))
}  # End of function save_grob_multi_format


extract_legend <- function(plot) {
  require(ggplot2)
  # Get the grob
  g <- ggplotGrob(plot)
  # Find the legend
  legend_idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  if (length(legend_idx) > 0) {
    legend <- g$grobs[[legend_idx]]
    return(legend)
  }
  return(NULL)
}  # End of function extract_legend


save_all_legends <- function(plot_dir) {
  
  cat("  Saving plot legends separately...\n")
  
  # Create legend subdirectory
  legend_dir <- file.path(plot_dir, "Legends")
  dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ========================================#
  # 1. STAGE COMPOSITION LEGEND
  # ========================================#
  stage_colors <- c("lightblue", "skyblue", "steelblue", "darkblue", "purple")
  stage_names <- c("Seeds", "Seedlings", "Juveniles", "Sub-adults", "Adults")
  
  # Create dummy plot for legend
  p_stage <- ggplot(data.frame(
    x = 1:5, 
    y = 1, 
    stage = factor(stage_names, levels = stage_names)
  )) +
    geom_bar(aes(x = x, y = y, fill = stage), stat = "identity") +
    scale_fill_manual(
      values = setNames(stage_colors, stage_names),
      name = "Life stage (0-4)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.key.size = unit(1.2, "cm"),
      #legend.background = element_rect(fill = "white", color = "black", size = 0.5)
    )
  
  # Extract and save legend
  legend_grob <- extract_legend(p_stage)
  
  if (!is.null(legend_grob)) {
    # PDF
    pdf(file.path(legend_dir, paste0("stage_composition_legend.pdf")), 
        width = 4, height = 3)
    grid.draw(legend_grob)
    dev.off()
    
    # PNG
    png(file.path(legend_dir, paste0("stage_composition_legend.png")), 
        width = 1200, height = 900, res = 300, type = "cairo")
    grid.draw(legend_grob)
    dev.off()
    cat("    ✓ Stage composition legend saved\n")
  } else {
    cat("    No legend grob was extracted, file not saved.\n")
  }
  
  # ========================================#
  # 2. OCCUPANCY PLOT LEGEND
  # ========================================#
  #  green line = mean, shaded = SD
  
  # Create info card instead of extracting legend
  pdf(file.path(legend_dir, paste0("occupancy_legend.pdf")), 
      width = 8, height = 3)
  
  grid.newpage()
  grid.text("Occupancy plot", 0.5, 0.7, gp=gpar(fontsize=14, fontface="bold"))
  grid.text("Dark green line = Mean occupied cells", 0.5, 0.5, gp=gpar(fontsize=12))
  grid.text("Light green area = ± Standard deviation", 0.5, 0.3, gp=gpar(fontsize=12))
  grid.text("Grey dashed line = Year 100 (End spin-up / start of simulation - 2010 data)", 0.5, 0.2, gp=gpar(fontsize=12))
  
  # # superseding and using grid package instead
  # # Create empty plot area with defined x/y coords
  # plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  # text(0.5, 0.7, "Occupancy plot", cex = 1.5, font = 2)
  # text(0.5, 0.5, "Dark green line = Mean occupied cells", cex = 1.2)
  # text(0.5, 0.3, "Light green area = ± Standard deviation", cex = 1.2)
  # text(0.5, 0.1, "Grey dashed line = Year 100 (End spin-up / start of simulation - 2010 data)", cex = 1.2)
  dev.off()
  
  cat("    ✓ Occupancy legend info saved\n")
  
  # ========================================#
  # 3. SPATIAL DISTRIBUTION LEGEND
  # ========================================#
  # Create legend for spatial plots
  
  pdf(file.path(legend_dir, paste0("spatial_legend.pdf")), 
      width = 8, height = 3)
  
  grid.newpage()
  grid.text("Spatial distribution", 0.5, 0.9, gp=gpar(fontsize=14, fontface="bold"))
  grid.text("Black circles (outlined) = Occupied 100 x 100 cells", 0.5, 0.7, gp=gpar(fontsize=12))
  grid.text("Grey background = Study area boundary", 0.5, 0.5, gp=gpar(fontsize=12))
  grid.text("Histogram: Mean population per cell (bars), Median-Q75 (lines), Max (red triangles)", 0.5, 0.3, gp=gpar(fontsize=12))
  
  # superseding and using grid package instead
  # plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  # text(0.5, 0.8, "Spatial distribution", cex = 1.5, font = 2)
  # text(0.5, 0.6, "Black dots = Occupied cells", cex = 1.2)
  # text(0.5, 0.4, "Grey background = Study area boundary", cex = 1.2)
  # text(0.5, 0.2, "Histogram: Mean population per cell (bars), Median-Q75 (lines), Max (red triangles)", cex = 1)
  
  dev.off()
  
  cat("    ✓ Spatial distribution legend saved\n")
  
  # ========================================#
  # 4. STAGE DYNAMICS LEGEND
  # ========================================#
  # Create legend for spatial dynamics plots
  # blue line is mean total individuals per stage (across replicates),
  # shaded blue is interquartile range 25-75,
  # orange is median 
  
  pdf(file.path(legend_dir, paste0("stage_dynamics_legend.pdf")), 
      width = 8, height = 3)
  
  grid.newpage()
  grid.text("Stage dynamics (total individuals by year and stage)", 0.5, 0.9, gp=gpar(fontsize=14, fontface="bold"))
  grid.text("Solid blue line = Mean total population", 0.5, 0.7, gp=gpar(fontsize=12))
  grid.text("Light blue = ± Interquartile range", 0.5, 0.5, gp=gpar(fontsize=12))
  grid.text("Dashed orange line: Median total population", 0.5, 0.3, gp=gpar(fontsize=12))
  
  dev.off()
  
  cat("    ✓ Stage dynamics legend saved\n")
  
  return(legend_dir)
}  # End of function save_all_legends


#==============================================================================#
#   5. Extinction timing (shared: simulation + results) ----
#==============================================================================#

calculate_extinction_years <- function(range_data, final_sim_year = 200) {
  if (is.null(range_data) || nrow(range_data) == 0) {
    return(list(mean = NA, se = NA, by_rep = NULL, persisted = NA))
  }
  
  is_dt <- FALSE
  if (requireNamespace("data.table", quietly = TRUE)) {
    is_dt <- data.table::is.data.table(range_data)
  }
  
  if (is_dt) {
    max_years_by_rep <- range_data[, .(max_year = max(Year)), by = Rep]
  } else {
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      stop("Either 'data.table' or 'dplyr' must be available for aggregation.")
    }
    max_years_by_rep <- dplyr::summarise(
      dplyr::group_by(range_data, Rep),
      max_year = max(.data$Year),
      .groups = "drop"
    )
  }
  
  max_years_by_rep$status <- ifelse(max_years_by_rep$max_year >= final_sim_year, "persisted", "extinct")
  max_years_by_rep$extinction_year <- ifelse(max_years_by_rep$status == "extinct", max_years_by_rep$max_year, NA)
  
  extinct_reps   <- max_years_by_rep[max_years_by_rep$status == "extinct", ]
  persisted_reps <- max_years_by_rep[max_years_by_rep$status == "persisted", ]
  
  if (nrow(extinct_reps) > 0) {
    mean_extinction <- mean(extinct_reps$extinction_year, na.rm = TRUE)
    se_extinction   <- stats::sd(extinct_reps$extinction_year, na.rm = TRUE) / sqrt(nrow(extinct_reps))
  } else {
    mean_extinction <- NA_real_
    se_extinction   <- NA_real_
  }
  
  list(
    mean = mean_extinction,
    se = se_extinction,
    n_extinct = nrow(extinct_reps),
    n_persisted = nrow(persisted_reps),
    by_rep = max_years_by_rep
  )
}  # End of function calculate_extinction_years

#==============================================================================#
#                        ----  End of functions_utils.R ----
#==============================================================================#
