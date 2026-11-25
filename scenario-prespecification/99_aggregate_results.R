# Aggregate MBO outputs across scenarios and runs.
#
# This script is intentionally read-only with respect to optimization. It only
# loads saved bundles and builds summary tables for reporting/export.

library(here)
library(gt)

results_root <- here::here("scenario-prespecification")

# Record the markov.misc version used (current commit of the local repo)
if (FALSE) {
  markov_misc_commit <- system2(
    "git",
    c(
      "-C",
      shQuote(path.expand("~/github/markov.misc")),
      "log",
      "-1",
      "--format=%h"
    ),
    stdout = TRUE
  )
  if (length(markov_misc_commit) == 0L) markov_misc_commit <- "unknown"
}
# Hardcoded stable commit at time of publication: 8b64521
markov_misc_commit <- "8b64521"


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
  stop(
    "No parameters_chosen_models*.rds files found under scenario-prespecification/."
  )
}

bundle_data <- purrr::map(bundle_files, readRDS)
targets <- purrr::map(targets_files, readRDS)


# Parse file names: parameters_chosen_models_<scenario>_<DGM|control>_...
parse_name <- function(nm) {
  parts <- strsplit(sub("^parameters_(chosen_models|targets)_", "", nm), "_")
  data.frame(
    name = nm,
    scenario = vapply(parts, `[`, "", 1),
    dgm = vapply(
      parts,
      function(x) setdiff(x[-1], c("control", "taooh", "final"))[1],
      ""
    ),
    is_control = vapply(
      parts,
      function(x) x[2] == "control" || x[3] == "control",
      NA
    ),
    target_taooh = vapply(
      parts,
      function(x) {
        if ("taooh" %in% x) x[which(x == "taooh") + 1] else NA_character_
      },
      ""
    )
  )
}

# Collapse any parameter value to a single display string
format_value <- function(x) {
  if (is.null(x)) {
    return("NULL")
  }
  if (is.function(x)) {
    lambda <- environment(x)$lambda
    return(
      if (is.null(lambda)) {
        "function"
      } else {
        paste0("rpois(", signif(lambda, 6), ")")
      }
    )
  }
  if (is.list(x)) {
    return(paste(vapply(x, format_value, ""), collapse = "; "))
  }
  paste0(if (is.numeric(x)) signif(x, 6) else format(x), collapse = ", ")
}

# Fields that are not parameters of the DGM (metadata), or are free to choose (e.g. n_patients)
drop_fields <- c(
  "merged_fit",
  "created",
  "r_version",
  "rms_version",
  "n_patients",
  "state2_day14",
  "actt2_thresholds",
  "mortality_day28",
  "n"
)

flatten_params <- function(params) {
  params <- params[setdiff(names(params), drop_fields)]
  vapply(params, format_value, "")
}

# Treatment-effect field(s) per DGM, and their value in the null scenario
treatment_fields <- list(
  LBM = "mu_treatment_effect",
  PEF = "hazard_ratio",
  OMP = c("tx_coef", "trt_coef")
)
null_values <- c(LBM = 0, PEF = 1, OMP = 0)
effect_field <- c(
  LBM = "mu_treatment_effect",
  PEF = "hazard_ratio",
  OMP = "trt_coef"
)

file_info <- parse_name(names(bundle_data))

# Every bundle should have a "parameters" element. Control bundles still
# holding candidates are flagged; their first candidate is used for the checks.
get_parameters <- function(nm) {
  x <- bundle_data[[nm]]
  if (!is.null(x$parameters)) {
    return(x$parameters)
  }
  warning(
    nm,
    ": no 'parameters' element; found candidate(s) ",
    paste(names(x), collapse = ", "),
    if (length(x) > 1L) ". Several candidates: using the first.",
    call. = FALSE
  )
  x[[1]]
}
all_parameters <- setNames(
  lapply(file_info$name, get_parameters),
  file_info$name
)

# Parameters without treatment-effect fields and metadata, as comparable strings
comparable <- function(params, dgm) {
  params <- params[setdiff(
    names(params),
    c(drop_fields, treatment_fields[[dgm]])
  )]
  vapply(params, format_value, "")
}

differences <- function(a, b) {
  fields <- union(names(a), names(b))
  fields[
    vapply(fields, function(f) !identical(a[f], b[f]), NA) |
      !(fields %in% names(a)) |
      !(fields %in% names(b))
  ]
}

# Checks: (1) scenarios differing only by taooh share all non-treatment
# parameters, (2) control bundles match the corresponding final bundles.
for (scen in unique(file_info$scenario)) {
  for (dgm in unique(file_info$dgm)) {
    grp <- file_info[file_info$scenario == scen & file_info$dgm == dgm, ]
    if (nrow(grp) < 2L) {
      next
    }
    ref <- comparable(all_parameters[[grp$name[1]]], dgm)
    for (nm in grp$name[-1]) {
      diff_fields <- differences(ref, comparable(all_parameters[[nm]], dgm))
      if (length(diff_fields) > 0L) {
        warning(
          scen,
          " ",
          dgm,
          ": ",
          nm,
          " differs from ",
          grp$name[1],
          " in: ",
          paste(diff_fields, collapse = ", "),
          call. = FALSE
        )
      }
    }
  }
}

# 1. One parameter set per trial and DGM, taken from the final (taooh) bundles
# with the treatment-effect parameters removed.
final_info <- file_info[!file_info$is_control, ]
final_info <- final_info[!duplicated(final_info[c("scenario", "dgm")]), ]

parameter_long <- purrr::pmap_dfr(
  final_info,
  function(name, scenario, dgm, ...) {
    p <- comparable(all_parameters[[name]], dgm)
    data.frame(
      DGM = dgm,
      Parameter = names(p),
      scenario = scenario,
      Value = unname(p)
    )
  }
)

parameter_table <- parameter_long |>
  tidyr::pivot_wider(names_from = scenario, values_from = Value) |>
  dplyr::arrange(factor(DGM, levels = c("LBM", "PEF", "OMP")))

parameter_table |>
  gt(groupname_col = "DGM", rowname_col = "Parameter") |>
  sub_missing(missing_text = "") |>
  tab_source_note(paste0("markov.misc commit: ", markov_misc_commit)) |>
  gtsave(
    filename = here(
      "scenario-prespecification",
      "parameters_chosen_models.html"
    )
  )

# 2. Treatment effect parameter: one row per scenario, incl. the null scenario
treatment_effect_params <- purrr::pmap_dfr(
  final_info[c("scenario", "dgm")],
  function(scenario, dgm) {
    grp <- file_info[
      file_info$scenario == scenario &
        file_info$dgm == dgm &
        !file_info$is_control,
    ]
    data.frame(
      DGM = dgm,
      Trial = scenario,
      Scenario = c("null", paste0("taooh_", grp$target_taooh)),
      Parameter = effect_field[[dgm]],
      Value = c(
        null_values[[dgm]],
        vapply(
          grp$name,
          USE.NAMES = FALSE,
          function(nm) {
            round(all_parameters[[nm]][[effect_field[[dgm]]]], 4)
          },
          0
        )
      )
    )
  }
) |>
  dplyr::arrange(factor(DGM, levels = c("LBM", "PEF", "OMP")), Trial)

treatment_effect_params |>
  gt(groupname_col = "DGM") |>
  gtsave(
    filename = here("scenario-prespecification", "parameters_targets.html")
  )
