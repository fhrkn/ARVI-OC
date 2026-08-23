# Aggregate MBO outputs across scenarios and runs.
#
# This script is intentionally read-only with respect to optimization. It only
# loads saved bundles and builds summary tables for reporting/export.

library(here)
library(gt)

results_root <- here::here("scenario-prespecification")

bundle_files <- list.files(
  results_root,
  pattern = "^parameters_chosen_models.*\\.rds$",
  full.names = TRUE
)
names(bundle_files) <- gsub(
  "\\.rds$",
  "",
  basename(bundle_files)
)

targets_files <- list.files(
  results_root,
  pattern = "^parameters_targets.*\\.rds$",
  full.names = TRUE
)
names(targets_files) <- gsub(
  "\\.rds$",
  "",
  basename(targets_files)
)

if (length(bundle_files) == 0L) {
  stop("No mbo_result_bundle_*.rds files found under MBO_scenarios/results/.")
}

bundle_data <- purrr::map(bundle_files, readRDS)
targets <- purrr::map(targets_files, readRDS)


# 1. Summary of parameters chosen for the control group and the final treatment effect for each target treatment effect size.
# Find unique parameter sets (ctrl group parameters are the same for difffernt mu_treatment_effect values)
DGMs <- c("LBM", "OMP", "PEF")
scenarios_sets <- setNames(c("ACTT1", "ACTT2"), c("ACTT1", "ACTT2"))

parameters_DGM_df <- list()
for (dgm in DGMs) {
  DGM_parameters <- bundle_data[
    grepl(dgm, names(bundle_data)) & !grepl("taooh", names(bundle_data))
  ]
  if (length(DGM_parameters) == 0L) {
    warning(paste0("No parameter sets found for DGM ", dgm, "."))
    parameters_DGM_df[[dgm]] <- rep(NULL, times = 3)
    next
  }
  ctrl_parameter_sets <- purrr::map(
    DGM_parameters,
    function(x) {
      params <- x[[1]]
      params$drift_start <- paste0(
        "rpois(",
        round(environment(params$drift_start)$lambda, 4),
        ")"
      )
      unlist(lapply(params, function(x) {
        paste0(if (is.numeric(x)) round(x, 4) else x, collapse = ", ")
      }))
    }
  )

  parameters_DGM_df[[dgm]] <- as.data.frame(do.call(
    cbind,
    ctrl_parameter_sets
  ))
}

parameter_table <- do.call(rbind, parameters_DGM_df)
parameter_table <- cbind(
  Parameter = rownames(parameter_table),
  parameter_table
)
names(parameter_table) <- c("Parameter", scenarios_sets)
markov.misc_commit <- setNames(
  as.data.frame(as.list(c(
    "markov.misc_commit:f9bfe57",
    rep(NA, length(scenarios_sets))
  ))),
  names(parameter_table)
)
parameter_table <- rbind(
  parameter_table,
  markov.misc_commit
)
parameter_table |> View()

parameter_table |>
  gt() |>
  gtsave(
    filename = here(
      "scenario-prespecification",
      "parameters_chosen_models.html"
    )
  )

treatment_effect_params <- list()
for (dgm in DGMs) {
  treatment_effect_params[[dgm]] <- purrr::map_dfr(
    scenarios_sets,
    function(scenario) {
      treatment_effect_param_data <- bundle_data[
        grepl(dgm, names(bundle_data)) &
          grepl("taooh", names(bundle_data)) &
          grepl(scenario, names(bundle_data))
      ]
      if (length(treatment_effect_param_data) == 0L) {
        warning(paste0(
          "No treatment effect parameter sets found for DGM ",
          dgm,
          " and scenario ",
          scenario,
          "."
        ))
        return(data.frame(
          DGM = dgm,
          scenario = NA,
          target_taooh = NA,
          mu_treatment_effect = NA
        ))
      }

      scenario_treatment_effects <- purrr::map(
        treatment_effect_param_data,
        function(x) {
          if (dgm == "LBM") {
            round(x$parameters$mu_treatment_effect, 4)
          }
        }
      )
      names <- strsplit(names(treatment_effect_param_data), "_")
      trt_effect_df <- data.frame(
        DGM = dgm,
        scenario = unlist(lapply(names, function(x) x[4])),
        target_taooh = unlist(lapply(names, function(x) x[7])),
        mu_treatment_effect = unlist(scenario_treatment_effects)
      )
      trt_effect_df
    }
  )

  # if (dgm == "LBM") {
  #   for (i in scenarios_sets) {
  #     ctrl_parameter_sets[grepl(
  #       i,
  #       names(ctrl_parameter_sets)
  #     )][[1]]["mu_treatment_effect"] <- paste0(
  #       c(
  #         ctrl_parameter_sets[grepl(
  #           i,
  #           names(ctrl_parameter_sets)
  #         )][[1]]["mu_treatment_effect"],
  #         treatment_effect_params[[i]]
  #       ),
  #       collapse = ", "
  #     )
  #   }
  # }
}

gt(do.call(rbind, treatment_effect_params)) |>
  gtsave(
    filename = here(
      "scenario-prespecification",
      "parameters_targets.html"
    )
  )
