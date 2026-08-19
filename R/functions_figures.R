# # ---
# title: "functions_figures.R"
# manuscript: "Title: Climate change, defaunation and disease interact to accelerate biodiversity loss in a tropical hotspot"
# corresponding_author: "Emma L Underwood"
# coauthors: "Kerry A Brown, Rebekka Allgayer, Mark Mulligan, Nigel Walford, Jette Wolff"
# last update: "2026-08-19"
# output: Helpers for the publication figures built in 04__publication_figures.R.
# requires: data.table, dplyr, tidyr, ggplot2, patchwork, cowplot, scales, terra  (loaded centrally in config.R)
# ---


#==============================================================================#
#   1. Scenario recoding ----
#==============================================================================#

recode_climate     <- function(x) dplyr::recode(x,
                                                "current" = "Current", "ssp1" = "Low", "ssp5" = "High", .default = x)  # End of function recode_climate


recode_defaunation <- function(ep) dplyr::case_when(
  ep == 0.20 ~ "Low", ep == 0.10 ~ "Medium", ep == 0.05 ~ "High",
  TRUE       ~ as.character(ep))  # End of function recode_defaunation


#==============================================================================#
#   2. Data loading ----
#==============================================================================#

read_all_range_files <- function(all_scenarios, outputs_dir) {
  
  cat("\nReading", nrow(all_scenarios), "Range files...\n")
  
  rbindlist(lapply(seq_len(nrow(all_scenarios)), function(i) {
    
    row   <- all_scenarios[i, ]
    fpath <- file.path(outputs_dir,
                       sprintf("Batch%d_Sim%d_Land1_Range.txt",
                               row$batch_num, row$sim_id))
    
    if (!file.exists(fpath)) {
      cat(sprintf("  MISSING: Batch%d Sim%d\n", row$batch_num, row$sim_id))
      return(NULL)
    }
    
    dt <- fread(fpath, showProgress = FALSE)
    dt[, `:=`(batch_num = row$batch_num, sim_id = row$sim_id)]
    dt
    
  }), fill = TRUE)
}  # End of function read_all_range_files


#==============================================================================#
#   3. Theme and axis helpers ----
#==============================================================================#

theme_journal <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border       = element_rect(colour = "grey80", fill = NA, size = 0.4),
      plot.title         = element_text(face = "bold", size = base + 1),
      plot.subtitle      = element_text(size = base - 1, colour = "grey40"),
      axis.title         = element_text(size = base - 1),
      axis.text          = element_text(size = base - 2),
      legend.title       = element_text(size = base - 1, face = "bold"),
      legend.text        = element_text(size = base - 2)
    )
}  # End of function theme_journal


dual_x <- function(x) paste0(x, "\n(", x + 1910, ")")  # End of function dual_x


add_defaunation_strip <- function(df) {
  df |> mutate(
    defaunation_label = factor(
      DEFAUNATION_STRIP[as.character(defaunation)],
      levels = unname(DEFAUNATION_STRIP)
    )
  )
}  # End of function add_defaunation_strip


#==============================================================================#
#   4. Spatial boundary extraction ----
#==============================================================================#

get_boundary_df <- function(sim_id_val, target_year, inv, threshold) {
  
  row <- inv[inv$sim_id == sim_id_val & inv$year == target_year, ]
  if (nrow(row) == 0) return(NULL)
  
  r      <- terra::rast(row$path[1])
  binary <- terra::ifel(r >= threshold, 1L, NA)
  if (all(is.na(terra::values(binary)))) return(NULL)
  
  poly <- tryCatch(terra::as.polygons(binary, dissolve = TRUE), error = function(e) NULL)
  if (is.null(poly) || nrow(poly) == 0) return(NULL)
  
  df           <- as.data.frame(terra::geom(poly, df = TRUE))
  df$sim_id    <- sim_id_val
  df$used_year <- target_year
  df
}  # End of function get_boundary_df


#==============================================================================#
#                        ----  End of functions_figures.R ----
#==============================================================================#
