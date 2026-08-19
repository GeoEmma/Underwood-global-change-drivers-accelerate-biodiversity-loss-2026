# # ---
# title: "01b__simulation_validation"
# manuscript: "Title: Climate change, defaunation and disease interact to accelerate biodiversity loss in a tropical hotspot"
# corresponding_author: "Emma L Underwood"
# coauthors: "Kerry A Brown, Rebekka Allgayer, Mark Mulligan, Nigel Walford, Jette Wolff"
# last update: "2026-08-19"
# output: A per-simulation run summary (simulation_validation.csv): replicate
#         count, whether adults persisted to year 200, and mean extinction year
#         for replicates that went extinct. Read for reporting / QA.
# prerequisite: 01__simulations.R
# ---

#==============================================================================#
#                           0. Workspace set up ----
#==============================================================================#

source("config.R")
source(file.path("R", "functions_utils.R"))
source(file.path("R", "functions_simulation.R"))   # calculate_extinction_years()

## Paths ----
masterfolder <- here::here("outputs", "simulations", "Master")
range_dir    <- file.path(masterfolder, "Outputs")

#==============================================================================#
#                           1. Load scenario mapping ----
#==============================================================================#

all_scenarios <- fread(
  here::here("outputs", "simulations", "Master_simulation_mapping.csv")
)

#==============================================================================#
#                           2. Summarise each simulation ----
#==============================================================================#

# Null-coalescing helper for optional list elements
`%||%` <- function(a, b) if (is.null(a)) b else a

summarise_one <- function(i) {
  row   <- all_scenarios[i, ]
  rfile <- file.path(range_dir,
                     sprintf("Batch%d_Sim%d_Land1_Range.txt", row$batch_num, row$sim_id))
  if (!file.exists(rfile)) {
    return(data.frame(sim_id = row$sim_id, batch_num = row$batch_num,
                      ssp = row$ssp, emig_prob = row$emig_prob,
                      status = "range_file_missing",
                      n_reps = NA_integer_, n_persisted = NA_integer_,
                      n_extinct = NA_integer_, mean_extinction_year = NA_real_,
                      se_extinction_year = NA_real_))
  }

  rng <- fread(rfile, showProgress = FALSE)
  ext <- calculate_extinction_years(rng, final_sim_year = 200)

  data.frame(
    sim_id               = row$sim_id,
    batch_num            = row$batch_num,
    ssp                  = row$ssp,
    emig_prob            = row$emig_prob,
    status               = "ok",
    n_reps               = if (!is.null(ext$by_rep)) nrow(ext$by_rep) else NA_integer_,
    n_persisted          = ext$n_persisted %||% NA_integer_,
    n_extinct            = ext$n_extinct   %||% NA_integer_,
    mean_extinction_year = ext$mean        %||% NA_real_,
    se_extinction_year   = ext$se          %||% NA_real_
  )
}  # End of function summarise_one


validation <- do.call(rbind, lapply(seq_len(nrow(all_scenarios)), summarise_one))

#==============================================================================#
#                           3. Write output ----
#==============================================================================#

out_csv <- file.path(masterfolder, "simulation_validation.csv")
write.csv(validation, out_csv, row.names = FALSE)
cat("Wrote", out_csv, "with", nrow(validation), "rows\n")
print(validation)

#==============================================================================#
#                        ----  End of script ----
#==============================================================================#
