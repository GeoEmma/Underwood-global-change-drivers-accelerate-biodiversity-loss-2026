# # ---
# title: "functions_occupancy.R"
# manuscript: "Title: Climate change, defaunation and disease interact to accelerate biodiversity loss in a tropical hotspot"
# corresponding_author: "Emma L Underwood"
# coauthors: "Kerry A Brown, Rebekka Allgayer, Mark Mulligan, Nigel Walford, Jette Wolff"
# last update: "2026-08-19"
# output: Functions that turn raw RangeShiftR population output into occupancy
#         probability rasters (GeoTiff) and matching point vectors (GeoPackage).
#         Sourced by 03__occupancy_probability.R.
# requires: data.table, terra  (loaded centrally in config.R)
# ---

# These functions summarise, across replicates, the proportion of runs in which
# each 100 m cell was occupied at a given simulation year, then write one raster
# and one point layer per year. The point layers exist because occupied cells
# are sparse (~200 of ~15 million); they are invisible at map scale as a raster
# but symbolise more cleanly as points in QGIS.

#==============================================================================#
#                 1. Coordinate transform ----
#==============================================================================#

# Convert RangeShiftR grid coordinates to real-world cell centres.
#
# RangeShiftR resets the landscape to a false origin of (0, 0). The x and y
# columns in the Pop file are 0-based cell indices (column and row number), not
# metres. This was confirmed against the study landscape: max(x) = 2389 matches
# ncol(template) - 1 = 2389, and max(y) = 6379 matches nrow(template) - 1.
# The projected cell centre is therefore:
#     centre = origin + index * resolution + resolution / 2
# Using metres directly (an earlier approach) collapsed all records to within a
# few kilometres of the origin, so the index-based formula below is the correct
# one for these outputs.
rs_coords_to_projected <- function(dt, template_r) {

  stopifnot(all(c("x", "y") %in% names(dt)),
            inherits(template_r, "SpatRaster"))

  res_val <- terra::res(template_r)[1]   # 100 m
  ox      <- terra::xmin(template_r)
  oy      <- terra::ymin(template_r)
  ext_val <- terra::ext(template_r)

  # Index -> projected cell centre
  dt[, x := ox + x * res_val + res_val / 2]
  dt[, y := oy + y * res_val + res_val / 2]

  # Guard against floating-point drift at the extent edges
  dt[, x := pmin(pmax(x, ext_val$xmin), ext_val$xmax)]
  dt[, y := pmin(pmax(y, ext_val$ymin), ext_val$ymax)]

  invisible(dt)
} # End of function rs_coords_to_projected


#==============================================================================#
#                 2. Determine analysis years ----
#==============================================================================#

# Choose the years to map for one simulation, with an exit for populations that crash (early-extinction fallback).
#
# Reads the small Range file to test whether adults (NInd_stage4) are still
# present at each target year. If a target year is not viable, the last viable
# year before it is used instead; targets never reached return NA.
get_key_years <- function(sim_id,
                          all_scenarios,
                          masterfolder,
                          target_years = c(110L, 140L, 170L),
                          adult_col    = "NInd_stage4") {

  # Local alias so data.table does not resolve sim_id as a column name
  .sid <- sim_id
  scen <- all_scenarios[which(all_scenarios$sim_id == .sid), ]

  if (nrow(scen) == 0L)
    stop("sim_id ", .sid, " not found in all_scenarios")

  rpath <- file.path(masterfolder, "Outputs",
                     sprintf("Batch%d_Sim%d_Land1_Range.txt",
                             scen$batch_num[1L], .sid))

  if (!file.exists(rpath))
    stop("Range file not found:\n  ", rpath)

  # Range file is small; read only the columns needed
  rng <- fread(rpath, select = c("Rep", "Year", adult_col),
               showProgress = FALSE)
  setnames(rng, adult_col, "NAdults")

  # Mean adults per year across replicates; viable = mean adults > 1
  yr_adults <- rng[, .(mean_adults = mean(NAdults, na.rm = TRUE)), by = Year]
  viable    <- sort(yr_adults[mean_adults > 1, Year])

  actual_years <- vapply(target_years, function(ty) {
    if (ty %in% viable) return(ty)
    candidates <- viable[viable < ty]
    if (length(candidates) == 0L) return(NA_integer_)
    max(candidates)
  }, integer(1L))

  names(actual_years) <- paste0("t", seq_along(target_years))

  cat(sprintf("  Sim %02d key years: %s\n", .sid,
              paste(mapply(function(nm, yr)
                sprintf("%s=%s", nm, ifelse(is.na(yr), "NA", yr)),
                names(actual_years), actual_years),
                collapse = "  ")))

  actual_years
} # End of function get_key_years


#==============================================================================#
#                 3. Compute and write occupancy rasters ----
#==============================================================================#

# Build occupancy probability rasters (and optional point layers) for one simulation.
#
# Memory approach: read only five columns, filter to the required years before
# any grouping, collapse to one flag per replicate/cell, then average across
# replicates. Only cells with occupancy probability > 0 are written.
#
# terra::writeRaster() does not recompute GDAL statistics by default, which can
# leave QGIS or gdalinfo reporting a mean of -9999. terra::global() is called
# before writing to confirm real values are present.
compute_occ_prob_rasters <- function(sim_id,
                                      all_scenarios,
                                      masterfolder,
                                      template_r,
                                      key_years,
                                      fix_coords    = TRUE,
                                      occ_col       = "NInd",
                                      n_reps        = 10L,
                                      raster_dir    = NULL,
                                      points_dir    = NULL,
                                      points_format = "gpkg") {

  .sid <- sim_id
  scen <- all_scenarios[which(all_scenarios$sim_id == .sid), ]

  if (nrow(scen) == 0L)
    stop("sim_id ", .sid, " not found in all_scenarios")

  ppath <- file.path(masterfolder, "Outputs",
                     sprintf("Batch%d_Sim%d_Land1_Pop.txt",
                             scen$batch_num[1L], .sid))

  if (!file.exists(ppath))
    stop("Pop file not found:\n  ", ppath)

  # Unique viable years (drop NA fallbacks and any duplicates they create)
  years_needed <- unique(stats::na.omit(as.integer(key_years)))

  if (length(years_needed) == 0L) {
    message("  Sim ", .sid, ": extinct before all target years; no rasters produced.")
    return(stats::setNames(vector("list", length(key_years)), names(key_years)))
  }

  cat(sprintf("\nSim %02d: reading Pop file for years {%s}\n",
              .sid, paste(years_needed, collapse = ", ")))

  ## 3.1 Read minimal columns and filter to target years ----
  pop <- fread(ppath,
               select       = c("Rep", "Year", "x", "y", occ_col),
               showProgress = FALSE)
  pop <- pop[Year %in% years_needed]
  gc()

  cat(sprintf("  %s rows retained after year filter\n",
              format(nrow(pop), big.mark = ",")))

  if (nrow(pop) == 0L) {
    message("  No rows found for years: ", paste(years_needed, collapse = ", "))
    return(stats::setNames(vector("list", length(key_years)), names(key_years)))
  }

  ## 3.2 Convert coordinates ----
  if (fix_coords) rs_coords_to_projected(pop, template_r)

  setnames(pop, occ_col, "occ_val")

  ## 3.3 Occupancy probability per cell and year ----
  # Per replicate/year/cell: occupied if any individual present.
  # Per year/cell: occupied replicates / total replicates.
  occ <- pop[,
    .(occupied = as.integer(max(occ_val, na.rm = TRUE) > 0L)),
    by = .(Rep, Year, x, y)
  ][,
    .(n_occ = sum(occupied)),
    by = .(Year, x, y)
  ][,
    occ_prob := n_occ / n_reps
  ][occ_prob > 0]

  rm(pop); gc()

  cat(sprintf("  %s occupied cell-year records\n",
              format(nrow(occ), big.mark = ",")))

  ## 3.4 One raster (and point layer) per time point ----
  raster_list <- vector("list", length(key_years))
  names(raster_list) <- names(key_years)

  for (nm in names(key_years)) {

    yr <- key_years[[nm]]

    if (is.na(yr)) {
      cat(sprintf("  %s: extinct before this time point; skipped\n", nm))
      raster_list[[nm]] <- NULL
      next
    }

    yr_data <- occ[Year == yr]
    cat(sprintf("  %s (year %d): %d occupied cells\n", nm, yr, nrow(yr_data)))

    # Empty raster matching the template, NA background
    r <- terra::rast(template_r)
    terra::values(r) <- NA_real_
    names(r) <- sprintf("OccProb_Sim%02d_Year%03d", .sid, yr)

    if (nrow(yr_data) > 0L) {
      cell_ids <- terra::cellFromXY(r, as.matrix(yr_data[, .(x, y)]))
      valid    <- !is.na(cell_ids)

      if (sum(valid) < nrow(yr_data))
        cat(sprintf("    Note: %d cells fell outside the raster extent and were dropped\n",
                    sum(!valid)))

      if (any(valid))
        r[cell_ids[valid]] <- yr_data$occ_prob[valid]
    }

    # Confirm real values exist before writing
    r_stats <- terra::global(r, c("min", "max", "mean"), na.rm = TRUE)
    cat(sprintf("    Raster stats - min: %.3f  max: %.3f  mean: %.4f\n",
                r_stats$min, r_stats$max, r_stats$mean))

    # Provenance stored as raster metadata tags
    terra::metags(r) <- c(
      sim_id    = as.character(.sid),
      batch     = as.character(scen$batch_num[1L]),
      ssp       = as.character(scen$ssp),
      emig_prob = as.character(scen$emig_prob),
      time_pt   = nm,
      year      = as.character(yr)
    )

    raster_list[[nm]] <- r

    ### Write raster GeoTiff ----
    # >>> QGIS INPUT: these GeoTiffs are loaded directly in QGIS for the
    #     spatial occupancy figures (manuscript Fig 4, Fig S4) and are also
    #     read by 04__publication_figures.R. Do not rename them; the
    #     OccProb_Sim<XX>_Year<YYY>.tif pattern is parsed downstream.
    if (!is.null(raster_dir)) {
      dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)

      out_tif <- file.path(raster_dir,
                           sprintf("OccProb_Sim%02d_Year%03d.tif", .sid, yr))

      terra::writeRaster(
        r, out_tif,
        overwrite = TRUE,
        datatype  = "FLT4S",
        gdal      = c("COMPRESS=LZW", "PREDICTOR=3")  # PREDICTOR=3 for floats
      )

      # Read back to confirm the written statistics
      r_written     <- terra::rast(out_tif)
      written_stats <- terra::global(r_written, c("min", "max", "mean"), na.rm = TRUE)
      cat(sprintf("    Written - min: %.3f  max: %.3f  mean: %.4f  [%s]\n",
                  written_stats$min, written_stats$max, written_stats$mean,
                  basename(out_tif)))
      rm(r_written)
    }

    ### Write point vector ----
    # >>> QGIS INPUT: point layers for symbolising sparse occupied cells at
    #     map scale (manuscript Fig 4, Fig S4). One GeoPackage per simulation,
    #     one layer per year.
    if (!is.null(points_dir)) {
      fmt <- tolower(points_format)

      if (fmt == "gpkg") {
        out_pts    <- file.path(points_dir,
                                sprintf("OccPoints_Sim%02d.gpkg", .sid))
        layer_name <- sprintf("Sim%02d_Year%03d", .sid, yr)
      } else {
        out_pts    <- file.path(points_dir,
                                sprintf("OccPoints_Sim%02d_Year%03d.shp", .sid, yr))
        layer_name <- NULL
      }

      occ_raster_to_points(r, out_pts, layer_name)
    }
  }

  invisible(raster_list)
} # End of function compute_occ_prob_rasters


#==============================================================================#
#                 4. Single-simulation wrapper ----
#==============================================================================#

# Run occupancy raster and point production for one simulation end to end.
run_occ_rasters_single <- function(sim_id,
                                   all_scenarios,
                                   masterfolder,
                                   template_r,
                                   raster_dir    = NULL,
                                   points_dir    = NULL,
                                   points_format = "gpkg",
                                   target_years  = c(110L, 140L, 170L),
                                   fix_coords    = TRUE,
                                   occ_col       = "NInd",
                                   n_reps        = 10L) {

  # Warn if the template is not at the 100 m Pop-file resolution
  actual_res <- terra::res(template_r)[1]
  if (actual_res != 100)
    warning(sprintf(
      "Template resolution is %d m but Pop files are 100 m. Check the template.",
      actual_res))

  cat(sprintf("\n%s\nProcessing Sim %02d\n%s\n",
              strrep("-", 45), sim_id, strrep("-", 45)))

  t0 <- proc.time()

  ky <- get_key_years(sim_id, all_scenarios, masterfolder, target_years)

  rs <- compute_occ_prob_rasters(
    sim_id        = sim_id,
    all_scenarios = all_scenarios,
    masterfolder  = masterfolder,
    template_r    = template_r,
    key_years     = ky,
    fix_coords    = fix_coords,
    occ_col       = occ_col,
    n_reps        = n_reps,
    raster_dir    = raster_dir,
    points_dir    = points_dir,
    points_format = points_format
  )

  elapsed <- round((proc.time() - t0)[["elapsed"]] / 60, 2)
  cat(sprintf("  Sim %02d complete in %.2f min\n", sim_id, elapsed))

  invisible(list(key_years = ky, rasters = rs))
} # End of function run_occ_rasters_single


#==============================================================================#
#                 5. Batch wrapper ----
#==============================================================================#

# Process several simulations in sequence, freeing memory between each.
# For parallel processing, wrap run_occ_rasters_single() in foreach/doParallel.
run_occ_rasters_batch <- function(sim_ids,
                                  all_scenarios,
                                  masterfolder,
                                  template_r,
                                  raster_dir    = NULL,
                                  points_dir    = NULL,
                                  points_format = "gpkg",
                                  target_years  = c(110L, 140L, 170L),
                                  fix_coords    = TRUE,
                                  occ_col       = "NInd",
                                  n_reps        = 10L) {

  t_total <- proc.time()
  cat(sprintf("\n%s\nBatch occupancy rasters: %d simulations\n%s\n",
              strrep("=", 55), length(sim_ids), strrep("=", 55)))

  results <- vector("list", length(sim_ids))
  names(results) <- as.character(sim_ids)

  for (i in seq_along(sim_ids)) {
    sid <- sim_ids[[i]]
    cat(sprintf("\n[%d / %d]", i, length(sim_ids)))

    results[[as.character(sid)]] <- tryCatch(
      run_occ_rasters_single(
        sim_id        = sid,
        all_scenarios = all_scenarios,
        masterfolder  = masterfolder,
        template_r    = template_r,
        raster_dir    = raster_dir,
        points_dir    = points_dir,
        points_format = points_format,
        target_years  = target_years,
        fix_coords    = fix_coords,
        occ_col       = occ_col,
        n_reps        = n_reps
      ),
      error = function(e) {
        message(sprintf("  ERROR in Sim %02d: %s", sid, conditionMessage(e)))
        list(key_years = NULL, rasters = NULL, error = conditionMessage(e))
      }
    )

    gc()  # release memory between simulations
  }

  elapsed_total <- round((proc.time() - t_total)[["elapsed"]] / 60, 1)
  cat(sprintf("\n%s\nBatch complete: %d sims in %.1f min\n%s\n",
              strrep("=", 55), length(sim_ids), elapsed_total, strrep("=", 55)))

  invisible(results)
} # End of function run_occ_rasters_batch


#==============================================================================#
#                 6. Occupancy points export ----
#==============================================================================#

# Convert the non-NA cells of an occupancy raster to a point vector.
#
# Attributes written per point:
#   occ_prob   occupancy probability (0-1)
#   occ_class  High >= 0.8 | Med >= 0.5 | Low >= 0.2 | Rare < 0.2
#   sim_id, ssp, emig_prob, year   from the raster metadata tags, if present
occ_raster_to_points <- function(r,
                                 out_path,
                                 layer_name = NULL) {

  stopifnot(inherits(r, "SpatRaster"))

  occ_df <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)

  if (nrow(occ_df) == 0L) {
    message("  No occupied cells; point file not written.")
    return(invisible(NULL))
  }

  # terra names the value column after the raster layer; standardise it
  names(occ_df)[3] <- "occ_prob"

  # Occupancy class, for categorised styling in QGIS
  occ_df$occ_class <- as.character(cut(
    occ_df$occ_prob,
    breaks         = c(0, 0.2, 0.5, 0.8, 1.0),
    labels         = c("Rare", "Low", "Med", "High"),
    include.lowest = TRUE,
    right          = TRUE
  ))

  # Carry provenance from the raster metadata tags
  tags <- tryCatch(terra::metags(r), error = function(e) character(0))
  for (field in c("sim_id", "ssp", "emig_prob", "year")) {
    occ_df[[field]] <- if (field %in% names(tags)) tags[[field]] else NA_character_
  }

  pts <- terra::vect(occ_df, geom = c("x", "y"), crs = terra::crs(r))

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  ext <- tolower(tools::file_ext(out_path))

  if (ext == "gpkg") {
    lyr <- if (!is.null(layer_name)) layer_name else
      tools::file_path_sans_ext(basename(out_path))
    terra::writeVector(pts, out_path,
                       filetype  = "GPKG",
                       layer     = lyr,
                       overwrite = TRUE)
  } else {
    # Shapefile field names are truncated to 10 characters by GDAL
    terra::writeVector(pts, out_path, overwrite = TRUE)
  }

  cat(sprintf("  %d points -> %s\n", nrow(pts), basename(out_path)))
  invisible(pts)
} # End of function occ_raster_to_points


# Export point files for every time point in a raster list (one call per sim).
occ_rasters_to_points_all <- function(raster_list,
                                      out_dir,
                                      sim_id,
                                      format = "gpkg") {

  stopifnot(is.list(raster_list))
  format <- tolower(format)

  out_list <- vector("list", length(raster_list))
  names(out_list) <- names(raster_list)

  for (nm in names(raster_list)) {

    r <- raster_list[[nm]]
    if (is.null(r)) {
      cat(sprintf("  %s: NULL (extinct); skipped\n", nm))
      next
    }

    # Year from metadata, or parsed from the layer name as a fallback
    tags <- tryCatch(terra::metags(r), error = function(e) character(0))
    yr   <- if ("year" %in% names(tags)) tags[["year"]] else
              sub(".*Year([0-9]+).*", "\\1", names(r)[1])

    if (format == "gpkg") {
      out_path   <- file.path(out_dir, sprintf("OccPoints_Sim%02d.gpkg", sim_id))
      layer_name <- sprintf("Sim%02d_Year%s", sim_id, yr)
    } else {
      out_path   <- file.path(out_dir,
                              sprintf("OccPoints_Sim%02d_Year%s.shp", sim_id, yr))
      layer_name <- NULL
    }

    out_list[[nm]] <- occ_raster_to_points(r, out_path, layer_name)
  }

  invisible(out_list)
} # End of function occ_rasters_to_points_all

#==============================================================================#
#                        ----  End of functions_occupancy.R ----
#==============================================================================#
