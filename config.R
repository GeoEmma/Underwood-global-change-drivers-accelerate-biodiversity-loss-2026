# # ---
# title: "config.R"
# manuscript: "Title: Climate change, defaunation and disease interact to accelerate biodiversity loss in a tropical hotspot"
# corresponding_author: "Emma L Underwood"
# coauthors: "Kerry A Brown, Rebekka Allgayer, Mark Mulligan, Nigel Walford, Jette Wolff"
# last update: "2026-08-19"
# output: Central configuration sourced by every numbered script.
#         Sets graphics options, loads packages, and resolves namespace
#         conflicts. Working directory is provided by the RStudio project
#         (*.Rproj) via here::here(); no setwd() is used.
# ---

#==============================================================================#
#                           Graphics / headless option ----
#==============================================================================#

# Use Cairo for bitmap devices so plotting works on headless HPC nodes
# (avoids the "unable to load R_X11.so / libXt.so.6" error).
options(bitmapType = "cairo")

#==============================================================================#
#                           Packages ----
#==============================================================================#

# Install (if needed) and load. For a fully pinned environment use renv:
#   renv::restore()   # installs the exact versions in renv.lock
.load_pkg <- function(x) {
  if (!requireNamespace(x, quietly = TRUE))
    install.packages(x, repos = "http://cran.us.r-project.org")
  suppressPackageStartupMessages(library(x, character.only = TRUE))
}  # End of function .load_pkg

.required_packages <- c(
  "here",        # project-root-relative paths
  "data.table",  # fast I/O and aggregation
  "terra",       # raster / vector handling
  "dplyr", "tidyr",
  "ggplot2", "patchwork", "cowplot", "scales", "viridis",
  "readr",
  "foreach", "doParallel"   # parallel simulation / plotting
)
invisible(lapply(.required_packages, .load_pkg))

# RangeShiftR is not on CRAN. Install once, pinned to the released tag:
#   pak::pak("RangeShifter/RangeShiftR-pkg/RangeShiftR@v3.0.0")
# It is loaded by the scripts that need it (01, and 03 for the template CRS).

#==============================================================================#
#                           Namespace conflicts ----
#==============================================================================#

# Force the tidyverse verbs to win against terra/stats equivalents.
select <- dplyr::select
filter <- dplyr::filter

#==============================================================================#
#                           Optional: manuscript font ----
#==============================================================================#

# Our manuscript uses "Callisto MT". It is applied in Canva during figure
# finalisation, so R does not need it. To use it in R outputs, install the font
# and uncomment (requires functions_utils.R to be sourced first):
#   plot_font <- setup_plot_fonts("Callisto MT", "serif")

#==============================================================================#
#                        ----  End of config.R ----
#==============================================================================#
