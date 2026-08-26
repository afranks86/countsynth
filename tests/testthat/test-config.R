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
