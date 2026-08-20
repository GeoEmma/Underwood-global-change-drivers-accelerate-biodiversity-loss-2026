# # ---
# title: "04__publication_figures"
# manuscript: "Title: Climate change, defaunation and disease interact to accelerate biodiversity loss in a tropical hotspot"
# corresponding_author: "Emma L Underwood"
# coauthors: "Kerry A Brown, Rebekka Allgayer, Mark Mulligan, Nigel Walford, Jette Wolff"
# last update: "2026-08-19"
# output: Main manuscript figures produced in R
#         - Figure 3 (stressor hierarchy: occupancy trajectories + decline rate)
#         - Figure 4 (spatial refugia; final version composed in QGIS)
#         - Figure 5 (pathogen tipping point: years-to-extinction heatmap)
# note: These are the analytical originals. Final manuscript figures were
#       finalised in QGIS/Canva and differ visually (see README figure table).
# prerequisite: 01__simulations.R, and 03__occupancy_probability.R for the
#               spatial figures (reads OccRasters/OccPoints).
# ---

#==============================================================================#
#                           0. Workspace set up ----
#==============================================================================#

source("config.R")
source(file.path("R", "functions_utils.R"))
source(file.path("R", "functions_figures.R"))

# =============================================================================#
#                       1.  Load and harmonise all data ----
# =============================================================================#

# 1a.  Scenario mapping ----
all_scenarios <- fread(
  file.path(basefolder, "01__simulations/Master_simulation_mapping.csv")
)

# Recode SSP and emigration probability to display labels
# NOTE: recode_climate() moved to R/functions_figures.R

# NOTE: recode_defaunation() moved to R/functions_figures.R


all_scenarios <- all_scenarios |>
  filter(ssp != "ssp3") |>          # explicit SSP3 removal (not included in manuscript)
  mutate(
    climate     = factor(recode_climate(ssp),           levels = CLIMATE_LEVELS),
    defaunation = factor(recode_defaunation(emig_prob),  levels = DEFAUNATION_LEVELS),
    
    # Pathogen introduction year per simulation.
    # Batch 1: no pathogen.
    # Batch 2: current-climate sims (ssp == "current") are introduced at year 100
    #          (0 years of climate pre-conditioning); all SSP sims in Batch 2
    #          are introduced at year 101 (1 year of exposure).
    # Batch 3: year 131 (31 years pre-conditioning).
    # Batch 4: year 161 (61 years pre-conditioning).
    pathogen_year = dplyr::case_when(
      batch_num == 1L                            ~ NA_integer_,
      batch_num == 2L & ssp == "current"         ~ 100L,
      batch_num == 2L                            ~ 101L,
      batch_num == 3L                            ~ 131L,
      batch_num == 4L                            ~ 161L,
      TRUE                                       ~ NA_integer_
    )
  )

# Validation print — check this looks correct before running figures
cat("\n--- Scenario mapping summary (SSP3 excluded) ---\n")
print(
  all_scenarios |>
    count(batch_num, climate, defaunation) |>
    tidyr::pivot_wider(names_from = climate, values_from = n),
  n = 20
)


# =============================================================================#
#           1b.  Read ALL Range files (all batches) ----
# =============================================================================#

# NOTE: read_all_range_files() is found in R/functions_figures.R


range_raw <- read_all_range_files(all_scenarios, outputs_dir)
cat(sprintf("Read %s rows across %d simulations.\n",
            format(nrow(range_raw), big.mark = ","),
            length(unique(paste(range_raw$batch_num, range_raw$sim_id)))))


# 1c.  Mean across replicates, years >= 100 only ----
range_mean <- range_raw[Year >= 100,
                        .(mean_occ    = mean(NOccupCells,   na.rm = TRUE),
                          sd_occ      = sd(NOccupCells,     na.rm = TRUE),
                          mean_adults = mean(NInd_stage4,   na.rm = TRUE),
                          n_reps      = .N),
                        by = .(batch_num, sim_id, Year)
][, `:=`(
  se_occ = sd_occ / sqrt(n_reps),
  ci_lo  = pmax(0, mean_occ - 1.96 * sd_occ / sqrt(n_reps)),
  ci_hi  =         mean_occ + 1.96 * sd_occ / sqrt(n_reps)
)]

# Join scenario labels
range_mean <- range_mean |>
  left_join(
    all_scenarios |> select(batch_num, sim_id, climate, defaunation, pathogen_year),
    by = c("batch_num", "sim_id")
  )


# 1d.  Relative change from baseline (inner_join = same temporal resolution) --
baseline_sim_id <- all_scenarios |>
  filter(batch_num == 1, climate == "Current", defaunation == "Low") |>
  pull(sim_id)

stopifnot(length(baseline_sim_id) == 1)

baseline_occ <- range_mean |>
  filter(sim_id == baseline_sim_id) |>
  select(Year, baseline_occ = mean_occ)

range_rel <- range_mean |>
  inner_join(baseline_occ, by = "Year") |>
  mutate(
    rel_change = case_when(
      mean_occ == 0 | baseline_occ == 0 ~ NA_real_,
      TRUE ~ (mean_occ - baseline_occ) / baseline_occ * 100
    )
  )


# 1e.  Per-replicate extinction timing (needed for Figure 3 SE) ----
# For each rep in each pathogen sim, find the last year NInd_stage4 > 0.
# Done on range_raw so there is rep-level variance.
ext_per_rep <- range_raw[NInd_stage4 > 0 & batch_num > 1,
                         .(last_adult_year = max(Year)),
                         by = .(batch_num, sim_id, Rep)
]

# Sims where a rep survived to year 200: last_adult_year is a lower bound.
# Flag as censored to note this in the figure.
ext_per_rep[, censored := (last_adult_year >= 200L)]

# Add scenario labels and pathogen year
ext_per_rep <- ext_per_rep |>
  left_join(
    all_scenarios |> select(batch_num, sim_id, climate, defaunation, pathogen_year),
    by = c("batch_num", "sim_id")
  ) |>
  mutate(years_to_ext = last_adult_year - pathogen_year)

# Aggregate to mean ± SE per scenario
ext_summary <- ext_per_rep |>
  group_by(batch_num, sim_id, climate, defaunation, pathogen_year) |>
  summarise(
    mean_yte   = mean(years_to_ext, na.rm = TRUE),
    se_yte     = sd(years_to_ext,   na.rm = TRUE) / sqrt(n()),
    n_reps     = n(),
    n_censored = sum(censored),
    .groups    = "drop"
  ) |>
  mutate(
    precond_label = dplyr::case_when(
      batch_num == 2 & pathogen_year == 100L ~ "Minimal\n(0 yr — current)",
      batch_num == 2                         ~ "Minimal\n(1 yr — SSP)",
      batch_num == 3                         ~ "Moderate\n(31 yr)",
      batch_num == 4                         ~ "Severe\n(61 yr)"
    ),
    precond_label = factor(precond_label,
                           levels = c("Minimal\n(0 yr — current)", "Minimal\n(1 yr — SSP)",
                                      "Moderate\n(31 yr)", "Severe\n(61 yr)"))
  )


# 1f.  Refugial areas — last viable year per Batch 1 sim ----------------------#
refugia_b1 <- range_mean |>
  filter(batch_num == 1, mean_adults > 0) |>
  group_by(batch_num, sim_id, climate, defaunation) |>
  summarise(last_adult_year = max(Year), .groups = "drop")

cat("\nBatch 1 refugial summary:\n")
print(refugia_b1 |> select(sim_id, climate, defaunation, last_adult_year))


# =============================================================================#
# 2.  Shared theme and helpers ----
# =============================================================================#

# NOTE: theme_journal() moved to R/functions_figures.R


# Dual x-axis: simulation year + calendar year in brackets
# NOTE: dual_x() moved to R/functions_figures.R


# Ordered facet labels
DEFAUNATION_STRIP <- c(
  Low    = "Low defaunation",
  Medium = "Medium defaunation",
  High   = "High defaunation"
)

# NOTE: add_defaunation_strip() moved to R/functions_figures.R


# =============================================================================#
# FIGURE 3: Occupancy trends ----
# =============================================================================#
# (a)  Absolute occupancy — shows the three climate lines per defaunation panel
# (b)  Relative change    — makes the defaunation-dominance argument explicit
# (c)  Decline-rate grid — shows Climate x Defaunation interactions in one cell
# =============================================================================#


b1 <- range_rel |>
  filter(batch_num == 1) |>
  add_defaunation_strip()

# --- Panel (a): absolute occupancy -------------------------------------------#
fig3a <- ggplot(b1, aes(x = Year, y = mean_occ, colour = climate, fill = climate)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.22, colour = NA) +
  geom_line(size = 0.9, na.rm = TRUE) +
  facet_wrap(~ defaunation_label, ncol = 1, scales = "free_y",
             strip.position = "left") +
  scale_colour_manual(values = CLIMATE_COLOURS, name = "Climate",
                      labels = c("Current", "Low (SSP1-2.6)", "High (SSP5-8.5)")) +
  scale_fill_manual(  values = CLIMATE_COLOURS, name = "Climate",
                      labels = c("Current", "Low (SSP1-2.6)", "High (SSP5-8.5)")) +
  scale_x_continuous(breaks = seq(100, 200, 20), labels = dual_x,
                     expand = c(0.01, 0)) +
  scale_y_continuous(labels = comma) +
  labs(x = NULL, y = "Mean occupied cells",
       title = "(a)  Mean range occupancy",
       subtitle = "Shaded bands = 95% CI across 10 replicates") +
  theme_journal() +
  theme(
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 9),
    strip.placement   = "outside",
    panel.spacing     = unit(0.7, "lines"),
    legend.position   = "none"
  )


# --- Panel (b): relative change from baseline --------------------------------
fig3b <- ggplot(
  b1 |> filter(!is.na(rel_change)),
  aes(x = Year, y = rel_change, colour = climate)
) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", size = 0.4) +
  geom_line(size = 0.9, na.rm = TRUE) +
  facet_wrap(~ defaunation_label, ncol = 1, scales = "free_y",
             strip.position = "left") +
  scale_colour_manual(values = CLIMATE_COLOURS, name = "Climate",
                      labels = c("Current", "Low (SSP1-2.6)", "High (SSP5-8.5)")) +
  scale_x_continuous(breaks = seq(100, 200, 20), labels = dual_x,
                     expand = c(0.01, 0)) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  labs(x = NULL, y = "Change from baseline (%)",
       title = "(b)  Relative change from baseline",
       subtitle = "Baseline: Current climate x Low defaunation (Sim 1)") +
  theme_journal() +
  theme(
    strip.text.y.left = element_blank(),   # labels already in panel (a)
    strip.placement   = "outside",
    panel.spacing     = unit(0.7, "lines"),
    legend.position   = "none"
  )


# --- Panel (c): decline-rate interaction heatmap -----------------------------
# Metric: % change in mean occupied cells from year 100 to last available year,
# expressed per decade so short and long sims are comparable.
decline_rate <- b1 |>
  group_by(sim_id, climate, defaunation) |>
  summarise(
    occ_yr100  = mean_occ[Year == min(Year)],
    occ_final  = mean_occ[Year == max(Year)],
    n_years    = max(Year) - min(Year),
    .groups    = "drop"
  ) |>
  mutate(
    pct_per_decade = (occ_final - occ_yr100) / occ_yr100 / (n_years / 10) * 100
  )

# Text colour: white on dark fills (strong decline), dark on light
text_colour <- ifelse(decline_rate$pct_per_decade < -4, "white", "grey15")

fig3c <- ggplot(
  decline_rate,
  aes(x = climate, y = defaunation, fill = pct_per_decade)
) +
  # Slightly inset tiles (width/height < 1) to give clean separation between cells
  geom_tile(colour = "white", size = 2.5, width = 0.92, height = 0.92) +
  geom_text(
    # "per decade" removed — value is self-explanatory with axis context
    aes(label = sprintf("%.1f%%", pct_per_decade)),
    colour = text_colour, size = 4, fontface = "bold"
  ) +
  # Ecological colour scale:
  #   Green  = range expanding (positive value, only Current x Low)
  #   White  = no net change
  #   Red    = rapid range decline
  scale_fill_gradient2(
    low      = "#d7191c",   # red   = rapid decline (most negative value)
    mid      = "#f7f7f7",   # white  = no change
    high     = "#1a9641",   # green = expansion / slowest decline
    midpoint = 0,
    name     = "Decline rate\n(% per decade)",
    labels   = function(x) paste0(x, "%")
  ) +
  # x-axis at top to match the column reading direction of panels (a) and (b)
  scale_x_discrete(
    position = "top",
    limits   = CLIMATE_LEVELS,
    labels   = c("Current", "Low\n(SSP1)", "High\n(SSP5)")
  ) +
  # Defaunation: Low at top, High at bottom — matches row order of panels (a)/(b)
  scale_y_discrete(limits = rev(DEFAUNATION_LEVELS)) +
  labs(
    x     = NULL,
    y     = "Defaunation level",
    title = "(c)  Decline rate\nClimate x Defaunation"
  ) +
  # Square aspect ratio so tiles are square cubes
  coord_fixed(ratio = 1) +
  theme_journal() +
  theme(
    panel.grid      = element_blank(),
    panel.border    = element_rect(fill = NA, colour = "grey80", size = 0.4),
    axis.text.x     = element_text(lineheight = 0.85, size = 9),
    axis.text.x.top = element_text(margin = margin(b = 4)),
    legend.position = "none"   # legend saved separately below
  )

# --- Save Figure 3c legend as a standalone file ------------------------------
fig3c_legend_plot <- ggplot(
  decline_rate,
  aes(x = climate, y = defaunation, fill = pct_per_decade)
) +
  geom_tile() +
  scale_fill_gradient2(
    low      = "#d7191c",
    mid      = "#f7f7f7",
    high     = "#1a9641", 
    midpoint = 0,
    name     = "Decline rate\n(% per decade)",
    labels   = function(x) paste0(x, "%"),
    breaks   = c(2, 0, -2.5, -5, -7.5, -10)
  ) +
  theme_void() +
  theme(
    legend.position   = "right",
    legend.title      = element_text(size = 10, face = "bold", lineheight = 1.2),
    legend.text       = element_text(size = 9),
    legend.key.height = unit(1.4, "cm"),
    legend.key.width  = unit(0.45, "cm")
  )

ggsave(file.path(plot_dir, "Figure1c_legend.pdf"),
       cowplot::get_legend(fig3c_legend_plot),
       width = 2.5, height = 4, dpi = 300)
ggsave(file.path(plot_dir, "Figure1c_legend.png"),
       cowplot::get_legend(fig3c_legend_plot),
       width = 2.5, height = 4, dpi = 300, bg = "white")
cat("Figure 3c legend saved separately.\n")


# --- Shared x-axis label strip -----------------------------------------------
x_label <- ggplot() +
  labs(x = "Simulation year (calendar year)") +
  theme_void() +
  theme(axis.title.x = element_text(size = 10, face = "plain"))


# --- Assemble Figure 3 -------------------------------------------------------
# Panels (a) and (b) are tall (3 facets each); (c) is a small square.
# Layout: [a | b | c] with a shared legend row at bottom.

legend_plot <- ggplot(b1, aes(x = Year, y = mean_occ, colour = climate)) +
  geom_line() +
  scale_colour_manual(values = CLIMATE_COLOURS, name = "Climate scenario",
                      labels = c("Current", "Low (SSP1-2.6)", "High (SSP5-8.5)")) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title    = element_text(size = 10, face = "bold"),
        legend.text     = element_text(size = 9))
shared_legend <- cowplot::get_legend(legend_plot)

fig3_top <- cowplot::plot_grid(
  fig3a, fig3b, fig3c,
  nrow        = 1,
  rel_widths  = c(1.15, 1.15, 0.7),
  align       = "h",
  axis        = "tb"
)

fig3 <- cowplot::plot_grid(
  fig3_top, shared_legend,
  ncol        = 1,
  rel_heights = c(1, 0.07)
)

ggsave(file.path(plot_dir, "Figure1_stressor_hierarchy.pdf"),
       fig3, width = 15, height = 10, dpi = 300)
ggsave(file.path(plot_dir, "Figure1_stressor_hierarchy.png"),
       fig3, width = 15, height = 10, dpi = 300, bg = "white")
cat("Figure 3 saved.\n")


# =============================================================================#
# Figure 4: Spatial refugia ----
# =============================================================================#
# All 9 Batch 1 scenarios shown as overlaid range boundary polygons on one map.
# Colour = climate scenario; linetype = defaunation level.
# Background = ecoregion polygons for landscape context.
# Boundary drawn at occupancy probability >= 0.10.
# Year used per sim = last year with adults > 0 (from refugia_b1).
# =============================================================================#

BOUNDARY_THRESHOLD <- 0.10

# Load spatial context
template_r <- terra::rast(file.path(template_dir, "CAZ_100m.tif"))
veg        <- if (file.exists(veg_path)) terra::vect(veg_path) else NULL
pa         <- if (file.exists(pa_path))  terra::vect(pa_path)  else NULL

# Raster inventory
inv_all <- {
  tifs   <- list.files(occ_raster_dir, pattern = "\\.tif$", full.names = TRUE)
  fnames <- basename(tifs)
  data.frame(
    sim_id = suppressWarnings(as.integer(sub(".*Sim([0-9]+).*",  "\\1", fnames))),
    year   = suppressWarnings(as.integer(sub(".*Year([0-9]+).*", "\\1", fnames))),
    path   = tifs,
    stringsAsFactors = FALSE
  ) |> filter(!is.na(sim_id), !is.na(year))
}

# Extract boundary polygon from one raster at given threshold
# NOTE: get_boundary_df() moved to R/functions_figures.R


# For each Batch 1 sim, pick raster year closest to (and <= ) last_adult_year
BATCH1_SIMS <- all_scenarios |> filter(batch_num == 1) |> pull(sim_id)

boundary_list <- lapply(BATCH1_SIMS, function(sid) {
  
  last_yr    <- refugia_b1$last_adult_year[refugia_b1$sim_id == sid]
  if (length(last_yr) == 0) return(NULL)
  
  sim_inv    <- inv_all[inv_all$sim_id == sid, ]
  candidates <- sim_inv$year[sim_inv$year <= last_yr]
  if (length(candidates) == 0) candidates <- sim_inv$year
  best_yr    <- max(candidates)
  
  get_boundary_df(sid, best_yr, inv_all, BOUNDARY_THRESHOLD)
})

boundary_df <- rbindlist(Filter(Negate(is.null), boundary_list), fill = TRUE) |>
  left_join(
    all_scenarios |> select(sim_id, climate, defaunation),
    by = "sim_id"
  )

ext_r <- terra::ext(template_r)
xlims <- c(ext_r$xmin, ext_r$xmax)
ylims <- c(ext_r$ymin, ext_r$ymax)

# Ecoregion background
# FIX: fill = veg_coords$fill_col outside aes() is a vector (one value per
# polygon coordinate row). ggplot2 >= 3.5.0 requires constants outside aes()
# to be strictly length 1 — a vector of 42,000+ values causes check_aesthetics()
# to error. Fix: move fill_col inside aes() and add scale_fill_identity(guide =
# "none") so the hex strings are used directly without generating a legend.
veg_layer <- if (!is.null(veg)) {
  veg_coords <- as.data.frame(terra::geom(veg, df = TRUE))
  geom_ids   <- sort(unique(veg_coords$geom))
  fills      <- grey(seq(0.80, 0.93, length.out = max(length(geom_ids), 2)))
  veg_coords$fill_col <- fills[match(veg_coords$geom, geom_ids)]
  list(
    geom_polygon(data = veg_coords,
                 aes(x = x, y = y,
                     group    = interaction(geom, part),
                     fill     = fill_col),   # inside aes() — fixes the error
                 colour      = grey(0.68),
                 size        = 0.15,
                 inherit.aes = FALSE),
    scale_fill_identity(guide = "none")  # use hex values as-is; no legend
  )
} else {
  list(annotate("rect", xmin = xlims[1], xmax = xlims[2],
                ymin = ylims[1], ymax = ylims[2],
                fill = "#EFEFEF", colour = NA))
}

# Protected areas layer (placeholder if missing)
pa_layer <- if (!is.null(pa)) {
  pa_coords <- as.data.frame(terra::geom(pa, df = TRUE))
  list(
    geom_polygon(data = pa_coords,
                 aes(x = x, y = y, group = interaction(geom, part)),
                 fill = NA, colour = "#2d6a4f",
                 linewidth = 0.45, linetype = "dashed", inherit.aes = FALSE)
  )
} else {
  list(
    annotate("rect",
             xmin = xlims[1] + diff(xlims) * 0.12,
             xmax = xlims[1] + diff(xlims) * 0.55,
             ymin = ylims[1] + diff(ylims) * 0.55,
             ymax = ylims[1] + diff(ylims) * 0.82,
             fill = NA, colour = "#2d6a4f", linewidth = 0.45, linetype = "dashed"),
    annotate("text",
             x = xlims[1] + diff(xlims) * 0.34,
             y = ylims[1] + diff(ylims) * 0.68,
             label = "[Protected area\nplaceholder]",
             colour = "#2d6a4f", size = 2.2, fontface = "italic")
  )
}

# --- Load Current climate rasters for the filled baseline layer --------------
# For each defaunation level, the "Current climate" sim gives us the baseline
# occupancy probability surface to display as a grey-scale filled raster.
# SSP1 and SSP5 boundaries are then drawn as coloured outlines on top.

current_sims <- all_scenarios |>
  filter(batch_num == 1, climate == "Current") |>
  select(sim_id, defaunation)

current_raster_list <- lapply(seq_len(nrow(current_sims)), function(i) {
  
  sid <- current_sims$sim_id[i]
  def <- as.character(current_sims$defaunation[i])
  
  last_yr    <- refugia_b1$last_adult_year[refugia_b1$sim_id == sid]
  if (length(last_yr) == 0) return(NULL)
  
  sim_inv    <- inv_all[inv_all$sim_id == sid, ]
  candidates <- sim_inv$year[sim_inv$year <= last_yr]
  if (length(candidates) == 0) candidates <- sim_inv$year
  best_yr    <- max(candidates)
  
  row <- sim_inv[sim_inv$year == best_yr, ]
  if (nrow(row) == 0) return(NULL)
  
  r  <- terra::rast(row$path[1])
  df <- as.data.frame(r, xy = TRUE)
  names(df)[3] <- "occ_prob"
  df <- df[!is.na(df$occ_prob) & df$occ_prob > 0, ]
  df$defaunation <- factor(def, levels = DEFAUNATION_LEVELS)
  
  # Convert occ_prob to grey hex string so we can piggyback on the
  # scale_fill_identity() already present in veg_layer — no dual fill scale.
  # grey(0.75) = light grey (low prob), grey(0.10) = near-black (high prob).
  df$fill_colour <- grey(0.75 - df$occ_prob * 0.65)
  df
})

current_raster_df <- do.call(rbind, Filter(Negate(is.null), current_raster_list))

# SSP boundary lines only — Current is shown as filled raster above
ssp_boundaries <- boundary_df |> filter(climate != "Current")

# Build Figure 4 --------------------------------------------------------------#
fig4 <- ggplot() +
  veg_layer +    # ecoregion polygons + scale_fill_identity(guide="none")
  pa_layer +
  # Current climate: semi-transparent grey-scale filled raster.
  # Darker = higher occupancy probability.
  # Uses fill_colour (grey hex) — identity scale already in veg_layer handles it.
  geom_raster(
    data  = current_raster_df,
    aes(x = x, y = y, fill = fill_colour),
    alpha = 0.60
  ) +
  # SSP1 and SSP5: coloured boundary outlines overlaid on the baseline raster.
  geom_path(
    data  = ssp_boundaries,
    aes(x      = x,
        y      = y,
        group  = interaction(sim_id, geom, part),
        colour = climate),
    size  = 0.75,
    alpha = 0.92
  ) +
  facet_wrap(
    ~ defaunation,
    ncol           = 3,
    strip.position = "top"
  ) +
  annotate("text",
           x = xlims[2] - diff(xlims) * 0.06,
           y = ylims[2] - diff(ylims) * 0.03,
           label = "N ^", size = 3.0, fontface = "bold", colour = "grey25") +
  scale_colour_manual(
    values = c(
      "Low"  = unname(CLIMATE_COLOURS["Low"]),
      "High" = unname(CLIMATE_COLOURS["High"])
    ),
    name   = "SSP scenario boundary",
    labels = c("Low (SSP1-2.6)", "High (SSP5-8.5)")
  ) +
  coord_equal(xlim = xlims, ylim = ylims, expand = FALSE) +
  labs(
    title    = "Final viable range extent — all Batch 1 scenarios",
    subtitle = paste0(
      "Panels = defaunation level  |  ",
      "Grey fill = Current climate occupancy probability (darker = higher certainty)  |  ",
      "Coloured outlines = SSP scenario range boundaries (occ. prob. >= ",
      BOUNDARY_THRESHOLD, ")  |  ",
      "Each scenario shown at last year with adults > 0"
    ),
    x = "Easting (m)",
    y = "Northing (m)"
  ) +
  theme_journal() +
  theme(
    panel.grid       = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    strip.text       = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey95", colour = NA),
    legend.position  = "bottom",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9)
  )

ggsave(file.path(plot_dir, "Figure2_spatial_refugia.pdf"),
       fig4, width = 14, height = 16, dpi = 300)
ggsave(file.path(plot_dir, "Figure2_spatial_refugia.png"),
       fig4, width = 14, height = 16, dpi = 300, bg = "white")
cat("Figure 4 saved.\n")

# =============================================================================#
# FIGURE 5: Pathogen tipping points ----
# =============================================================================#
# Heatmap of mean years to extinction after pathogen introduction.
# Rows = defaunation level; columns = climate scenario.
# Facets = pathogen pre-conditioning period (Batches 2, 3, 4).
# Cell colour = years to extinction (red = fast, blue = slow/survived).
# Cell text = mean ± SE; asterisk if any replicates censored (survived to 200).
# =============================================================================#

# Cell text: "mean\n+/- SE"  with * if censored replicates present
ext_summary <- ext_summary |>
  mutate(
    cell_label = sprintf("%.1f\n+/- %.1f%s",
                         mean_yte,
                         se_yte,
                         ifelse(n_censored > 0, "*", ""))
  )

# Text colour: white on dark tiles (< 10 yr = red end of scale), dark on light
ext_summary <- ext_summary |>
  mutate(
    text_col = ifelse(mean_yte < 10, "white", "grey15")
  )

# Climate labels for x-axis (Batch 2 includes Current; Batches 3/4 don't)
# We keep the full grid and leave missing combinations blank naturally.
fig5 <- ggplot(
  ext_summary,
  aes(x = climate, y = defaunation, fill = mean_yte)
) +
  # Slightly inset tiles give a separated "square blocks" appearance
  # matching Figure 3c style — no external packages needed
  geom_tile(colour = "white", size = 2.0, width = 0.90, height = 0.90) +
  geom_text(
    aes(label = cell_label, colour = I(text_col)),
    size = 2.8, lineheight = 1, fontface = "bold"
  ) +
  facet_wrap(
    ~ precond_label,
    nrow   = 1,
    ncol   = 4,
    scales = "fixed"
  ) +
  # Ecological colour scale for survival time:
  #   Blue  = long survival = relatively safe
  #   Cream = intermediate
  #   Red   = rapid extinction = urgent
  # Midpoint at ~10 yr (roughly the middle of the observed 3–18 yr range)
  scale_fill_gradient2(
    low      = "#d73027",   # red   = fast extinction (< 5 yr)
    mid      = "#ffffbf",   # cream = ~10 yr
    high     = "#2166ac",   # blue  = long survival (> 15 yr)
    midpoint = 10,
    name     = "Years to\nextinction",
    na.value = "grey92",
    limits   = c(0, max(ext_summary$mean_yte, na.rm = TRUE) * 1.05)
  ) +
  # Defaunation: Low at bottom, High at top (bottom = less stress, ascending)
  scale_y_discrete(limits = rev(DEFAUNATION_LEVELS)) +
  scale_x_discrete(
    limits = CLIMATE_LEVELS,
    labels = c("Current", "Low\n(SSP1-2.6)", "High\n(SSP5-8.5)")
  ) +
  # Square aspect ratio so tiles are square — matches Figure 3c
  coord_fixed(ratio = 1) +
  labs(
    title    = "Years to extinction after pathogen introduction",
    subtitle = paste0(
      "Facets = pathogen arrival timing (climate pre-conditioning duration)  |  ",
      "Batch 2 current-climate sims: year 100 (0 yr); Batch 2 SSP sims: year 101 (1 yr)  |  ",
      "Cell: mean +/- SE across 10 replicates  |  ",
      "* replicate(s) survived year 200  |  grey = scenario absent from batch"
    ),
    x = "Climate scenario",
    y = "Defaunation level"
  ) +
  theme_journal(base = 10) +
  theme(
    panel.grid        = element_blank(),
    strip.text        = element_text(face = "bold", size = 9),
    strip.background  = element_rect(fill = "grey95", colour = NA),
    axis.text.x       = element_text(lineheight = 0.85),
    legend.position   = "none"    # legend saved separately
  )

# --- Save Figure 5 legend separately -----------------------------------------
fig5_legend_plot <- ggplot(ext_summary, aes(x = climate, y = defaunation,
                                            fill = mean_yte)) +
  geom_tile() +
  scale_fill_gradient2(
    low      = "#d73027",
    mid      = "#ffffbf",
    high     = "#2166ac",
    midpoint = 10,
    name     = "Years to\nextinction",
    na.value = "grey92",
    breaks   = c(0, 5, 10, 15, 18),
    limits   = c(0, max(ext_summary$mean_yte, na.rm = TRUE) * 1.05)
  ) +
  theme_void() +
  theme(
    legend.position   = "right",
    legend.title      = element_text(size = 10, face = "bold", lineheight = 1.2),
    legend.text       = element_text(size = 9),
    legend.key.height = unit(1.4, "cm"),
    legend.key.width  = unit(0.45, "cm")
  )

ggsave(file.path(plot_dir, "Figure3_legend.pdf"),
       cowplot::get_legend(fig5_legend_plot),
       width = 2.5, height = 4, dpi = 300)
ggsave(file.path(plot_dir, "Figure3_legend.png"),
       cowplot::get_legend(fig5_legend_plot),
       width = 2.5, height = 4, dpi = 300, bg = "white")
cat("Figure 5 legend saved separately.\n")

ggsave(file.path(plot_dir, "Figure3_pathogen_tipping_point.pdf"),
       fig5, width = 14, height = 5, dpi = 300)
ggsave(file.path(plot_dir, "Figure3_pathogen_tipping_point.png"),
       fig5, width = 14, height = 5, dpi = 300, bg = "white")
cat("Figure 5 saved.\n")



# =============================================================================#
# Session summary
# =============================================================================#

cat("\n", strrep("=", 55), "\n")
cat("All figures saved to:", plot_dir, "\n")
cat(strrep("=", 55), "\n")
cat("Figure 3: Figure1_stressor_hierarchy.pdf/.png\n")
cat("Figure 4: Figure2_spatial_refugia.pdf/.png\n")
cat("Figure 5: Figure3_pathogen_tipping_point.pdf/.png\n")
cat("\nKey values to record for Figure 4 annotation:\n")

# Print key quantitative findings for Figure 4 annotation
cat("\n--- Decline rates (Figure 3c) ---\n")
print(decline_rate |>
        select(climate, defaunation, pct_per_decade) |>
        arrange(defaunation, climate))

cat("\n--- Relative change at year 100 (defaunation damage before climate) ---\n")
b1_yr100 <- range_rel |>
  filter(batch_num == 1, Year == 100) |>
  select(sim_id, climate, defaunation, rel_change) |>
  arrange(defaunation, climate)
print(b1_yr100)

cat("\n--- Extinction timing range across Batches 2-4 ---\n")
cat(sprintf("  Min mean years to extinction: %.1f\n", min(ext_summary$mean_yte, na.rm=TRUE)))
cat(sprintf("  Max mean years to extinction: %.1f\n", max(ext_summary$mean_yte, na.rm=TRUE)))
cat(sprintf("  Range compression (Severe vs Minimal, High x High):\n"))
high_high <- ext_summary |>
  filter(climate == "High", defaunation == "High") |>
  select(precond_label, mean_yte, se_yte)
print(high_high)


#==============================================================================#
#                        ----  End of workflow ----
#==============================================================================#