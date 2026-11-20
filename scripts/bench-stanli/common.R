# Shared setup for the stanli-vs-CmdStan benchmarks: the bundled panels'
# Stan data, built exactly as countsynth_fit() builds it, written once to
# JSON so both engines read byte-identical data.

pkgload::load_all("/Users/afranks/Dropbox/course/countsynth", quiet = TRUE)

TEST_DIR <- "/Users/afranks/Dropbox/course/test_countsynth"
OUT_DIR <- "/Users/afranks/Dropbox/course/countsynth/scripts/bench-stanli/out"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CASES <- list(
  # NB, fixed dispersion, censoring, te regression, rank shrinkage: the
  # heaviest code path through joint.stan.
  fertility_te_race = list(config = "configs/fertility_te_race.yml",
                           type = "race", rank = 4L),
  # Poisson, heavy censoring (the 47%-censored "other" group), rank
  # shrinkage, legacy treatment hierarchy.
  mortality_full = list(config = "configs/mortality_full.yml",
                        type = "race", rank = 5L),
  # NB, flat Dirichlet weights, legacy hierarchy: the small, plain model.
  fertility_full = list(config = "configs/fertility_full.yml",
                        type = "total", rank = 3L),
  # Poisson, cut mode: what stage 2's many short fits look like.
  mortality_con = list(config = "configs/mortality_con.yml",
                       type = "con_noncon", rank = 3L)
)

# The data list countsynth_fit() would hand to $sample().
case_stan_data <- function(name, model_treated = TRUE) {
  cs <- CASES[[name]]
  withr::local_dir(TEST_DIR)
  config <- read_countsynth_config(cs$config)
  data <- countsynth_data(config, type = cs$type)
  te_design <- if (model_treated && !is.null(config$model$treatment_effects)) {
    build_te_design(config$model$treatment_effects, data)
  } else {
    NULL
  }
  sd <- stan_data_joint(
    data,
    rank = cs$rank,
    model_treated = model_treated,
    outcome_distribution = config$model$outcome_distribution,
    nb_disp = config$model$nb_disp,
    sample_disp = config$model$sample_disp,
    adjust_for_missingness = config$model$adjust_for_missingness,
    gen_ypred = TRUE,
    te_design = te_design,
    time_fac_shape = gamma_shape_from_pct(
      config$model$factor_variation_pct, DEFAULT_TIME_FAC_SHAPE
    ),
    time_fe_shape = gamma_shape_from_pct(
      config$model$time_level_variation_pct, DEFAULT_TIME_FE_SHAPE
    ),
    rank_shrinkage = config$model$rank_shrinkage,
    time_level = config$model[["time_level"]] %||% "centered"
  )
  list(sd = sd, config = config, data = data)
}

case_json <- function(name, model_treated = TRUE) {
  path <- file.path(OUT_DIR, sprintf("%s_%s.json", name,
                                     if (model_treated) "joint" else "base"))
  if (!file.exists(path)) {
    cmdstanr::write_stan_json(case_stan_data(name, model_treated)$sd, path)
  }
  path
}

stan_file <- function(name = "joint") {
  file.path("/Users/afranks/Dropbox/course/countsynth/inst/stan",
            paste0(name, ".stan"))
}
INCLUDE_DIR <- "/Users/afranks/Dropbox/course/countsynth/inst/stan/include"
