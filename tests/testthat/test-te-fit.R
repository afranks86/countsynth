# Sampler-level checks for the treatment-effect regression. The load-bearing
# one is the reconstruction test: Stan's `te` must equal the R-side design
# times the sampled coefficients, which is the only thing tying the formula
# the user wrote to the linear predictor the likelihood actually used.

te_fit_data <- function() {
  groups <- c("g1", "g2")
  units <- c("A", "B", "C", "D")
  times <- seq(as.Date("2020-01-01"), by = "month", length.out = 6)
  withr::with_seed(7, {
    grid <- expand.grid(
      unit = units, time = times, group = groups, stringsAsFactors = FALSE
    )
    grid$treatment <- as.integer(
      (grid$unit == "C" & grid$time >= times[4]) |
        (grid$unit == "D" & grid$time >= times[5])
    )
    grid$denominator <- stats::runif(nrow(grid), 5000, 30000)
    grid$outcome <- stats::rpois(nrow(grid), 50 * grid$denominator / 1e4)
  })
  arrays <- build_model_arrays(grid, groups)
  arrays$df <- grid
  arrays$type <- "both"
  structure(arrays, class = c("countsynth_data", "list"))
}

te_fit_config <- function(data, te_opts, inference_mode = NULL) {
  countsynth_config(
    input_file = "unused.csv", output_dir = tempfile(),
    schema = countsynth_schema(
      unit_col = "unit", time_col = "time", treatment_col = "treatment",
      outcomes = list(countsynth_outcome("outcome_g1", "g1"))
    ),
    model = countsynth_model_opts(
      outcome_distribution = "NB",
      types = list(both = countsynth_type(groups = data$groups, ranks_to_test = 2)),
      inference_mode = inference_mode,
      treatment_effects = te_opts
    ),
    mcmc = countsynth_mcmc_opts(
      auto_parallelism = FALSE, chains = 1, parallel_chains = 1,
      iter_warmup = 100, iter_sampling = 100, thin = 1, seed = 99,
      progress = FALSE
    )
  )
}

# te = X %*% beta + sum_j Z_j[, q] * u_{j,q}[level] + treatment_kt_z * scale
reconstruct_te <- function(design, te_draws, kt, it_scale) {
  W <- te_wide_matrix(te_draws)
  n_draws <- nrow(W)
  out <- matrix(0, n_draws, design$n_exposed)
  if (length(design$x_names) > 0) {
    beta <- W[, te_coef_key("(fixed)", design$x_names, NA), drop = FALSE]
    out <- out + beta %*% t(design$X)
  }
  for (t in design$terms) {
    for (q in seq_along(t$predictors)) {
      u <- W[, te_coef_key(t$label, t$predictors[[q]], t$levels), drop = FALSE]
      out <- out + u[, t$index, drop = FALSE] *
        matrix(t$Z[, q], n_draws, design$n_exposed, byrow = TRUE)
    }
  }
  out + kt * it_scale
}

test_that("joint fit's te equals the R-side design reconstruction", {
  skip_on_cran()
  skip_if_no_cmdstan()
  data <- te_fit_data()
  opts <- countsynth_te_opts(
    formula = ~ 1 + event_time + (1 + event_time | group) + (1 | unit)
  )
  fit <- countsynth_fit(data, rank = 2, config = te_fit_config(data, opts))

  expect_equal(fit$stan_data$te_reg, 1L)
  expect_false(is.null(fit$te_design))

  te_stan <- extract_var_matrix(fit$fit, "te")
  recon <- reconstruct_te(
    fit$te_design, countsynth_te_draws(fit),
    extract_var_matrix(fit$fit, "treatment_kt_z"),
    as.numeric(extract_var_matrix(fit$fit, "treatment_it_scale"))
  )
  # CmdStan writes 6 significant digits to its CSV, so agreement is checked
  # relative to the scale of te rather than to machine epsilon.
  expect_lt(max(abs(te_stan - recon)) / max(abs(te_stan)), 1e-4)
})

test_that("the regression replaces the legacy hierarchy in the sampler", {
  skip_on_cran()
  skip_if_no_cmdstan()
  data <- te_fit_data()
  opts <- countsynth_te_opts(formula = ~ 1 + event_time)
  fit <- countsynth_fit(data, rank = 2, config = te_fit_config(data, opts))
  vars <- fit$fit$metadata()$stan_variables
  expect_true(all(c("te_beta", "te") %in% vars))
  # Legacy blocks are declared zero-size, and CmdStan omits zero-size
  # parameters from its output entirely.
  expect_false(any(
    c("category_treatment_effect", "unit_treatment_effect_z",
      "unit_category_te_z", "treatment_category_scale") %in% vars
  ))
})

test_that("without a formula the legacy treatment block is used", {
  skip_on_cran()
  skip_if_no_cmdstan()
  data <- te_fit_data()
  fit <- countsynth_fit(data, rank = 2, config = te_fit_config(data, NULL))
  expect_equal(fit$stan_data$te_reg, 0L)
  expect_null(fit$te_design)
  legacy <- posterior::as_draws_matrix(
    fit$fit$draws(variables = "category_treatment_effect")
  )
  expect_equal(ncol(legacy), length(data$groups))
  expect_error(countsynth_te_draws(fit), "no treatment-effect regression")
})

test_that("te draws and summaries have the documented shape", {
  skip_on_cran()
  skip_if_no_cmdstan()
  data <- te_fit_data()
  opts <- countsynth_te_opts(formula = ~ 1 + event_time + (1 + event_time | group))
  fit <- countsynth_fit(data, rank = 2, config = te_fit_config(data, opts))

  te <- countsynth_te_draws(fit)
  expect_true(all(
    c(".draw", ".chain", ".iteration", "term", "predictor", "level", "value")
      %in% names(te)
  ))
  # 2 fixed + 2 predictors x 2 group levels.
  n_coef <- length(unique(paste(te$term, te$predictor, te$level)))
  expect_equal(n_coef, 2L + 4L)
  expect_true(all(is.na(te$level[te$term == "(fixed)"])))
  expect_setequal(unique(te$term), c("(fixed)", "group"))

  tbl <- countsynth_te_coef_table(fit)
  expect_equal(nrow(tbl), n_coef)
  expect_true(all(tbl$p_positive >= 0 & tbl$p_positive <= 1))
  expect_true(all(tbl$lower_95 <= tbl$median & tbl$median <= tbl$upper_95))
})

test_that("regression summaries and plots cover population and by-level", {
  skip_on_cran()
  skip_if_no_cmdstan()
  data <- te_fit_data()
  opts <- countsynth_te_opts(formula = ~ 1 + event_time + (1 + event_time | group))
  fit <- countsynth_fit(data, rank = 2, config = te_fit_config(data, opts))

  pop <- countsynth_te_regression_summary(fit, predictor = "event_time")
  expect_equal(sort(unique(pop$x)), c(0, 1, 2))
  expect_equal(unique(pop$level), "population")
  expect_true(all(pop$lower_95 <= pop$lower_67))
  expect_true(all(pop$upper_67 <= pop$upper_95))

  by_group <- countsynth_te_regression_summary(
    fit,
    predictor = "event_time", by = "group"
  )
  expect_setequal(as.character(by_group$level), data$groups)
  # The percent scale transforms each draw before summarizing, so it agrees
  # with transforming the log-scale summary only up to the interpolation
  # between the two central draws (an even draw count here).
  pct <- countsynth_te_regression_summary(
    fit,
    predictor = "event_time", scale = "percent"
  )
  expect_equal(pct$median, 100 * (exp(pop$median) - 1), tolerance = 1e-3)
  expect_equal(order(pct$median), order(pop$median))

  expect_s3_class(
    countsynth_te_regression_plot(fit, predictor = "event_time"), "ggplot"
  )
  expect_s3_class(
    countsynth_te_regression_plot(fit, predictor = "event_time", by = "group"),
    "ggplot"
  )
  expect_s3_class(countsynth_te_coef_plot(fit, terms = "all"), "ggplot")
  expect_s3_class(plot(fit, which = "te_regression"), "ggplot")
  expect_error(
    countsynth_te_regression_plot(fit, predictor = "event_time", by = "nosuch"),
    "not a random-effect grouping term"
  )
  expect_error(
    countsynth_te_regression_plot(fit, predictor = "group"),
    "categorical"
  )
})

test_that("cut mode pools coefficient draws with matching provenance", {
  skip_on_cran()
  skip_if_no_cmdstan()
  data <- te_fit_data()
  opts <- countsynth_te_opts(formula = ~ 1 + event_time + (1 + event_time | group))
  config <- te_fit_config(data, opts, inference_mode = "cut")
  config$cut <- countsynth_cut_opts(
    num_stage1_draws = 2, stage2_draws_per_component = 5,
    stage2_mcmc = list(num_warmup = 100, num_samples = 100)
  )
  fit <- suppressWarnings(countsynth_cut_fit(data, rank = 2, config = config))

  te <- countsynth_te_draws(fit)
  expect_true(all(
    c("cut_component", "stage1_draw", "stage1_chain", "stage1_iteration")
      %in% names(te)
  ))
  # Coefficient draws are thinned on the same indices as the te draws, so the
  # two frames index the same pooled stage-2 draws.
  expect_setequal(unique(te$.draw), unique(fit$draws$.draw))
  expect_setequal(unique(te$cut_component), unique(fit$draws$cut_component))
  expect_equal(length(unique(te$.draw)), 2L * 5L)
  expect_s3_class(
    countsynth_te_regression_plot(fit, predictor = "event_time", by = "group"),
    "ggplot"
  )
})
