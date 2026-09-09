test_that("schema requires exactly one of outcomes / prefixes", {
  expect_error(
    bpnmf_schema("u", "t", "tr"),
    "exactly one"
  )
  expect_error(
    bpnmf_schema(
      "u", "t", "tr",
      outcomes = list(bpnmf_outcome("x", "lab")),
      outcomes_from_prefixes = bpnmf_prefixes("p_")
    ),
    "exactly one"
  )
})

test_that("cut mode requires model_treated", {
  expect_error(
    bpnmf_model_opts(model_treated = FALSE, inference_mode = "cut"),
    "model_treated"
  )
})

test_that("sample_disp requires NB", {
  expect_error(
    bpnmf_model_opts(outcome_distribution = "Poisson", sample_disp = TRUE),
    "NB"
  )
})

test_that("stage2_mcmc seed is rejected", {
  expect_error(
    bpnmf_cut_opts(stage2_mcmc = list(random_seed = 1)),
    "seed authority"
  )
  expect_error(
    bpnmf_cut_opts(stage2_mcmc = list(seed = 1)),
    "seed authority"
  )
})

test_that("stage1_method defaults to sampling and validates its ADVI knobs", {
  expect_equal(bpnmf_cut_opts()$stage1_method, "sample")
  expect_null(bpnmf_cut_opts()$stage1_variational)

  opts <- bpnmf_cut_opts(
    stage1_method = "variational",
    stage1_variational = list(algorithm = "meanfield", draws = 500L)
  )
  expect_equal(opts$stage1_method, "variational")
  expect_equal(opts$stage1_variational$draws, 500L)

  expect_error(bpnmf_cut_opts(stage1_method = "advi"), "stage1_method")
  # The seed is mcmc.seed's job, so it is not an ADVI knob.
  expect_error(
    bpnmf_cut_opts(
      stage1_method = "variational", stage1_variational = list(seed = 1)
    ),
    "Unknown key"
  )
  # Tuning a method you are not using is a mistake worth flagging.
  expect_warning(
    bpnmf_cut_opts(stage1_variational = list(iter = 100)),
    "ignored"
  )
})

test_that("aggregate unit needs exactly one selector", {
  expect_error(bpnmf_aggregate_unit("All"), "exactly one")
  expect_error(
    bpnmf_aggregate_unit(
      "All",
      include_all_units = TRUE, include_treated_units = TRUE
    ),
    "exactly one"
  )
  spec <- bpnmf_aggregate_unit("All", include_all_units = TRUE)
  expect_s3_class(spec, "bpnmf_aggregate_unit")
})

test_that("figures normalization matches the Python validator", {
  expect_identical(bpnmf_output_opts(figures = TRUE)$figures, sort(FIGURE_NAMES))
  expect_identical(bpnmf_output_opts(figures = FALSE)$figures, character())
  expect_identical(bpnmf_output_opts(figures = "all")$figures, sort(FIGURE_NAMES))
  expect_identical(bpnmf_output_opts(figures = "none")$figures, character())
  expect_identical(
    bpnmf_output_opts(figures = c("unit_fit", "ppc"))$figures,
    c("unit_fit", "ppc")
  )
  expect_error(bpnmf_output_opts(figures = "bogus"), "unknown figure")
})

test_that("adapt_delta bounds are enforced", {
  expect_error(bpnmf_mcmc_opts(adapt_delta = 1), "between 0 and 1")
  expect_error(bpnmf_mcmc_opts(adapt_delta = 0), "between 0 and 1")
})

write_yaml_config <- function(text) {
  path <- withr::local_tempfile(fileext = ".yaml", .local_envir = parent.frame())
  writeLines(text, path)
  path
}

BASE_YAML <- '
data:
  input_file: "data.csv"
  output_dir: "results"
  schema:
    unit_col: state
    time_col: time
    treatment_col: exposed
    outcomes_from_prefixes:
      outcome_prefix: "births_"
      denominator_prefix: "pop_"
model:
  outcome_distribution: "NB"
  types:
    total:
      groups: ["total"]
      ranks_to_test: [5]
      total_all: true
mcmc:
  num_warmup: 500
  num_samples: 1000
  thinning: 5
  target_accept: 0.85
  random_seed: 123
output:
  figures: "all"
  draws_format: parquet
cut:
  num_stage1_draws: 10
  stage2_mcmc:
    num_warmup: 200
'

test_that("YAML loader accepts the Python schema and translates names", {
  cfg <- read_bpnmf_config(write_yaml_config(BASE_YAML))
  expect_s3_class(cfg, "bpnmf_config")
  expect_equal(cfg$mcmc$iter_warmup, 500L)
  expect_equal(cfg$mcmc$iter_sampling, 1000L)
  expect_equal(cfg$mcmc$thin, 5L)
  expect_equal(cfg$mcmc$adapt_delta, 0.85)
  expect_equal(cfg$mcmc$seed, 123L)
  expect_identical(cfg$output$figures, sort(FIGURE_NAMES))
  expect_equal(cfg$output$draws_format, "parquet")
  expect_equal(cfg$cut$num_stage1_draws, 10L)
  expect_equal(cfg$cut$stage2_mcmc$num_warmup, 200)
  expect_equal(cfg$model$types$total$total_all, TRUE)
})

test_that("YAML loader reads the cut stage-1 method", {
  yaml <- sub(
    "  num_stage1_draws: 10",
    paste(
      "  num_stage1_draws: 10",
      "  stage1_method: variational",
      "  stage1_variational:",
      "    algorithm: meanfield",
      "    draws: 400",
      sep = "\n"
    ),
    BASE_YAML,
    fixed = TRUE
  )
  cfg <- read_bpnmf_config(write_yaml_config(yaml))
  expect_equal(cfg$cut$stage1_method, "variational")
  expect_equal(cfg$cut$stage1_variational$algorithm, "meanfield")
  expect_equal(cfg$cut$stage1_variational$draws, 400)

  expect_error(
    read_bpnmf_config(write_yaml_config(sub(
      "  stage1_variational:", "  stage1_variational:\n    bogus: 1",
      yaml, fixed = TRUE
    ))),
    "Unknown config key"
  )
})

test_that("unknown YAML keys are rejected at every level", {
  expect_error(
    read_bpnmf_config(write_yaml_config(sub(
      "mcmc:", "parallel:\n  analysis_workers: -1\nmcmc:", BASE_YAML
    ))),
    "Unknown config key"
  )
  expect_error(
    read_bpnmf_config(write_yaml_config(sub(
      "  figures: \"all\"", "  figures: \"all\"\n  filename_pattern: \"x\"",
      BASE_YAML
    ))),
    "Unknown config key"
  )
})

test_that("quoted string booleans are rejected", {
  expect_error(
    read_bpnmf_config(write_yaml_config(sub(
      "model:", "model:\n  sample_disp: \"false\"", BASE_YAML
    ))),
    "boolean"
  )
})

test_that("the Python repo's shipped configs load when present", {
  py_configs <- Sys.glob(
    file.path(
      "~", "hierarchical-bayesian-NMF-refactor", "configs",
      c("base_config.yaml", "fertility_config.yaml", "fertility_smoke_test.yaml",
        "fertility_cut_smoke_test.yaml", "education_joint.yaml",
        "education_cut.yaml", "test_config.yaml")
    )
  )
  skip_if(length(py_configs) == 0, "Python repo configs not present")
  for (path in py_configs) {
    cfg <- read_bpnmf_config(path)
    expect_s3_class(cfg, "bpnmf_config")
  }
})

test_that("time_aggregation parses, and the old `aggregation` key still loads", {
  with_data <- function(block) {
    sub("  schema:", paste0(block, "\n  schema:"), BASE_YAML, fixed = TRUE)
  }

  n <- read_bpnmf_config(write_yaml_config(with_data(
    "  time_aggregation:\n    enabled: true\n    n_periods: 7"
  )))
  expect_equal(n$time_aggregation$n_periods, 7L)
  expect_null(n$time_aggregation$period)

  p <- read_bpnmf_config(write_yaml_config(with_data(
    "  time_aggregation:\n    enabled: true\n    period: quarterly"
  )))
  expect_equal(p$time_aggregation$period, "quarterly")
  expect_null(p$time_aggregation$n_periods)

  expect_warning(
    old <- read_bpnmf_config(write_yaml_config(with_data(
      "  aggregation:\n    enabled: true\n    period: yearly"
    ))),
    "deprecated"
  )
  expect_equal(old$time_aggregation$period, "yearly")

  expect_error(
    read_bpnmf_config(write_yaml_config(with_data(
      paste0("  aggregation:\n    enabled: true\n",
             "  time_aggregation:\n    enabled: true")
    ))),
    "both"
  )
  expect_error(
    read_bpnmf_config(write_yaml_config(with_data(
      "  time_aggregation:\n    enabled: true\n    periods: 7"
    ))),
    "Unknown"
  )
})

test_that("plain named lists work anywhere a bpnmf_* object is expected", {
  cfg <- bpnmf_config(
    input_file = "data.csv", output_dir = "results",
    schema = list(
      unit_col = "state", time_col = "time", treatment_col = "exposed",
      outcomes_from_prefixes = list(
        outcome_prefix = "births_", denominator_prefix = "pop_"
      )
    ),
    model = list(types = list(total = list(groups = "total", ranks_to_test = 3))),
    mcmc = list(iter_warmup = 500, convergence = list(ess_min = 200)),
    output = list(
      figures = TRUE,
      aggregate_units = list(
        list(unit = "All treated", include_treated_units = TRUE)
      )
    ),
    time_aggregation = list(enabled = TRUE, n_periods = 3)
  )
  # Every level is coerced through its real constructor, so the result is
  # indistinguishable from the object-built form.
  expect_s3_class(cfg$schema, "bpnmf_schema")
  expect_s3_class(cfg$schema$outcomes_from_prefixes, "bpnmf_prefixes")
  expect_s3_class(cfg$model$types$total, "bpnmf_type")
  expect_s3_class(cfg$mcmc$convergence, "bpnmf_convergence")
  expect_s3_class(cfg$output$aggregate_units[[1]], "bpnmf_aggregate_unit")
  expect_s3_class(cfg$time_aggregation, "bpnmf_time_aggregation")
  expect_equal(cfg$mcmc$convergence$ess_min, 200)
  expect_equal(cfg$time_aggregation$n_periods, 3)
  # Defaults the list did not mention are still filled in.
  expect_identical(cfg$mcmc$gate_params, c("mu_ctrl", "te"))
})

test_that("list form and constructor form produce identical configs", {
  args <- list(input_file = "d.csv", output_dir = "o")
  built <- do.call(bpnmf_config, c(args, list(
    schema = bpnmf_schema(
      "state", "time", "exposed",
      outcomes_from_prefixes = bpnmf_prefixes("births_", "pop_")
    ),
    model = bpnmf_model_opts(types = list(total = bpnmf_type("total", 3)))
  )))
  listed <- do.call(bpnmf_config, c(args, list(
    schema = list(
      unit_col = "state", time_col = "time", treatment_col = "exposed",
      outcomes_from_prefixes = list(
        outcome_prefix = "births_", denominator_prefix = "pop_"
      )
    ),
    model = list(types = list(total = list(groups = "total", ranks_to_test = 3)))
  )))
  expect_equal(built, listed)
})

test_that("coercion catches typos and still runs the real validators", {
  schema <- bpnmf_schema("u", "t", "x", outcomes_from_prefixes = bpnmf_prefixes("b_"))
  expect_error(
    bpnmf_config("d", "o", schema = list(
      unit_col = "u", time_col = "t", treatment_col = "x", unit_column = "oops"
    )),
    'Unknown name "unit_column" in schema'
  )
  expect_error(
    bpnmf_model_opts(types = list(a = list(groups = "g", ranks = 3))),
    'Unknown name "ranks" in types\\$a'
  )
  # A named list must not become a way to skip validation.
  expect_error(
    bpnmf_config("d", "o", schema = schema, mcmc = list(adapt_delta = 1.5)),
    "strictly between 0 and 1"
  )
  expect_error(
    bpnmf_mcmc_opts(convergence = list(ess_min = -5)), "not >= 0"
  )
  expect_error(
    bpnmf_config("d", "o", schema = list("u", "t", "x")),
    "must be a .*bpnmf_schema.* object or a named list"
  )
})

test_that("factor_variation_pct maps to the Gamma shape and defaults to 20", {
  # time_fac ~ Gamma(a, a) enters the log-rate as log(time_fac), so
  # sd(log) ~ 1/sqrt(a); asking for a p% swing means sd(log) = log1p(p/100).
  expect_equal(gamma_shape_from_pct(NULL), 20)
  expect_equal(gamma_shape_from_pct(NULL, DEFAULT_TIME_FE_SHAPE), 1)
  # Exact inversion: Var[log x] = trigamma(shape) = log1p(pct/100)^2.
  for (p in c(5, 25, 100, 300)) {
    expect_equal(trigamma(gamma_shape_from_pct(p)), log1p(p / 100)^2,
                 tolerance = 1e-8)
  }
  # ~25% is the historical Gamma(20, 20), which is why that is the default.
  expect_equal(gamma_shape_from_pct(25), 20, tolerance = 0.05)
  # The historical time_fe prior, Gamma(1, 1), is a ~260% swing -- the
  # delta-method 1/sd^2 approximation would have called it 100%.
  expect_equal(gamma_shape_from_pct(100 * (exp(sqrt(trigamma(1))) - 1)), 1,
               tolerance = 1e-6)
  # Monotone: a wider expected swing is a looser (smaller-shape) prior.
  shapes <- vapply(c(5, 10, 25, 50, 100), gamma_shape_from_pct, numeric(1))
  expect_true(all(diff(shapes) < 0))
  expect_error(gamma_shape_from_pct(0), "not >= ")
  expect_error(gamma_shape_from_pct(-5), "not >= ")

  expect_null(bpnmf_model_opts()$time_level_variation_pct)
  expect_equal(
    bpnmf_model_opts(time_level_variation_pct = 80)$time_level_variation_pct, 80
  )

  expect_null(bpnmf_model_opts()$factor_variation_pct)
  expect_equal(bpnmf_model_opts(factor_variation_pct = 50)$factor_variation_pct, 50)
  expect_error(bpnmf_model_opts(factor_variation_pct = -1), "not >= ")
})

test_that("factor_variation_pct round-trips through YAML", {
  with_model <- function(block) {
    sub("model:", paste0("model:\n", block), BASE_YAML, fixed = TRUE)
  }
  absent <- read_bpnmf_config(write_yaml_config(BASE_YAML))
  expect_null(absent$model$factor_variation_pct)

  set <- read_bpnmf_config(write_yaml_config(
    with_model("  factor_variation_pct: 60")
  ))
  expect_equal(set$model$factor_variation_pct, 60)
  lvl <- read_bpnmf_config(write_yaml_config(
    with_model("  time_level_variation_pct: 80")
  ))
  expect_equal(lvl$model$time_level_variation_pct, 80)
  expect_equal(
    trigamma(gamma_shape_from_pct(set$model$factor_variation_pct)),
    log1p(0.6)^2
  )
  expect_error(
    read_bpnmf_config(write_yaml_config(with_model("  factor_variation: 60"))),
    "Unknown"
  )
})

test_that("rank_shrinkage opts validate and accept the TRUE shorthand", {
  o <- bpnmf_rank_shrinkage_opts()
  expect_s3_class(o, "bpnmf_rank_shrinkage_opts")
  expect_equal(o$group_mass_prior, c(2, 1))
  expect_equal(o$unit_sd_prior, 1)

  # The mass prior is a (shape, rate) pair, the spread prior one scale.
  expect_error(bpnmf_rank_shrinkage_opts(group_mass_prior = 2), "length 2")
  expect_error(
    bpnmf_rank_shrinkage_opts(group_mass_prior = c(2, 0)), "not >= "
  )
  expect_error(bpnmf_rank_shrinkage_opts(unit_sd_prior = c(1, 2)), "length 1")
  expect_error(bpnmf_rank_shrinkage_opts(unit_sd_prior = -1), "not >= ")

  expect_null(bpnmf_model_opts()$rank_shrinkage)
  expect_equal(bpnmf_model_opts(rank_shrinkage = TRUE)$rank_shrinkage, o)
  expect_null(bpnmf_model_opts(rank_shrinkage = FALSE)$rank_shrinkage)
  # A named list of the constructor's arguments, as elsewhere in the config.
  expect_equal(
    bpnmf_model_opts(rank_shrinkage = list(unit_sd_prior = 3))$rank_shrinkage$unit_sd_prior,
    3
  )
  expect_error(
    bpnmf_model_opts(rank_shrinkage = list(unit_sd = 3)), "Unknown"
  )
})

test_that("rank_shrinkage round-trips through YAML", {
  with_model <- function(block) {
    sub("model:", paste0("model:\n", block), BASE_YAML, fixed = TRUE)
  }
  absent <- read_bpnmf_config(write_yaml_config(BASE_YAML))
  expect_null(absent$model$rank_shrinkage)

  flag <- read_bpnmf_config(write_yaml_config(
    with_model("  rank_shrinkage: true")
  ))
  expect_equal(flag$model$rank_shrinkage, bpnmf_rank_shrinkage_opts())
  expect_null(
    read_bpnmf_config(write_yaml_config(with_model("  rank_shrinkage: false")))$model$rank_shrinkage
  )

  set <- read_bpnmf_config(write_yaml_config(with_model(
    "  rank_shrinkage:\n    group_mass_prior: [3, 2]\n    unit_sd_prior: 0.5"
  )))
  expect_equal(set$model$rank_shrinkage$group_mass_prior, c(3, 2))
  expect_equal(set$model$rank_shrinkage$unit_sd_prior, 0.5)
  # Unspecified keys keep the constructor defaults.
  partial <- read_bpnmf_config(write_yaml_config(
    with_model("  rank_shrinkage:\n    unit_sd_prior: 0.5")
  ))
  expect_equal(partial$model$rank_shrinkage$group_mass_prior, c(2, 1))

  expect_error(
    read_bpnmf_config(write_yaml_config(
      with_model("  rank_shrinkage:\n    group_mass: [3, 2]")
    )),
    "Unknown"
  )
})
