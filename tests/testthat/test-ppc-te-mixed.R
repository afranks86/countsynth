# Mixed predictive check of the treatment-effect shrinkage. The load-bearing
# tests are the base reconstructions: what the check keeps fixed must be
# exactly mu_ctrl plus the parts of te that do NOT belong to one unit, or the
# redrawn pieces are being added on top of their own posterior values.

mixed_data <- function() {
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

mixed_config <- function(data, te_opts, inference_mode = NULL) {
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

# Exposed-cell indices the way the Stan data define them.
mixed_exposed <- function(data) {
  D <- length(data$units)
  N <- length(data$times)
  exp_cell <- which(!flatten_kdn(data$control_idx_array))
  sub <- kdn_from_flat(exp_cell, D, N)
  list(
    exp_cell = exp_cell, exp_k = sub$k, exp_d = sub$d,
    exp_kd = (sub$d - 1L) * length(data$groups) + sub$k,
    n_exposed = length(exp_cell), K = length(data$groups), D = D
  )
}

mixed_formula <- ~ 1 + event_time + (1 | group) + (1 | unit) +
  (1 + event_time | group:unit)

test_that("only random effects whose levels sit inside one unit are redrawn", {
  data <- mixed_data()
  design <- build_te_design(countsynth_te_opts(formula = mixed_formula), data)
  comps <- te_unit_components(design, mixed_exposed(data))
  labels <- vapply(comps, `[[`, character(1), "label")
  expect_equal(labels, c(
    "(1 | unit)", "(1 | group:unit)", "(event_time | group:unit)",
    "exposed cell"
  ))
  # (1 | group) holds te_re_z[1:2] and te_re_scale[1]; skipping it must not
  # shift the offsets of the terms after it.
  expect_equal(comps[[1]]$z_cols, sprintf("te_re_z[%d]", 3:4))
  expect_equal(comps[[1]]$scale_col, "te_re_scale[2]")
  expect_equal(comps[[2]]$z_cols, sprintf("te_re_z[%d]", 5:8))
  expect_equal(comps[[3]]$z_cols, sprintf("te_re_z[%d]", 9:12))
  expect_equal(comps[[3]]$scale_col, "te_re_scale[4]")
})

test_that("a covariate shared across units stays at its posterior draw", {
  data <- mixed_data()
  # `own` is one-to-one with the units, so its levels are per-unit
  # coefficients under another name. `parity` alternates by month, so each of
  # its levels covers both treated units (C from month 4, D from month 5).
  covs <- expand.grid(
    unit = data$units, time = data$times, stringsAsFactors = FALSE
  )
  covs$own <- tolower(covs$unit)
  covs$parity <- ifelse(as.integer(format(covs$time, "%m")) %% 2 == 1, "odd", "even")
  covs$time <- as.character(covs$time)
  design <- build_te_design(
    countsynth_te_opts(
      formula = ~ 1 + (1 | own) + (1 | parity), covariates = covs
    ),
    data
  )
  comps <- te_unit_components(design, mixed_exposed(data))
  expect_equal(
    vapply(comps, `[[`, character(1), "label"),
    c("(1 | own)", "exposed cell")
  )
})

# A draws matrix built by hand from known coefficients: te is the full linear
# predictor, so base must come back as mu + the non-unit pieces alone.
test_that("base keeps the fixed surface and shared terms, nothing else", {
  data <- mixed_data()
  design <- build_te_design(countsynth_te_opts(formula = mixed_formula), data)
  ex <- mixed_exposed(data)
  comps <- te_unit_components(design, ex)
  n_draws <- 7
  n_exp <- ex$n_exposed
  utot <- sum(vapply(design$terms, function(t) {
    length(t$levels) * length(t$predictors)
  }, integer(1)))
  qtot <- sum(vapply(design$terms, function(t) length(t$predictors), integer(1)))

  withr::with_seed(11, {
    beta <- matrix(stats::rnorm(n_draws * ncol(design$X)), n_draws)
    re_scale <- matrix(stats::rexp(n_draws * qtot), n_draws)
    re_z <- matrix(stats::rnorm(n_draws * utot), n_draws)
    it_scale <- stats::rexp(n_draws, 10)
    kt_z <- matrix(stats::rnorm(n_draws * n_exp), n_draws)
    mu <- matrix(stats::rnorm(n_draws * n_exp, 5), n_draws)
  })
  fixed_and_shared <- beta %*% t(design$X)
  te <- fixed_and_shared + kt_z * it_scale
  uo <- 0L
  qo <- 0L
  for (t in design$terms) {
    for (q in seq_along(t$predictors)) {
      u <- re_z[, uo + seq_along(t$levels), drop = FALSE] * re_scale[, qo + q]
      piece <- u[, t$index, drop = FALSE] * rep(t$Z[, q], each = n_draws)
      te <- te + piece
      if (t$label == "group") fixed_and_shared <- fixed_and_shared + piece
      uo <- uo + length(t$levels)
    }
    qo <- qo + length(t$predictors)
  }
  m <- cbind(te, re_scale, re_z, it_scale, kt_z)
  colnames(m) <- c(
    sprintf("te[%d]", seq_len(n_exp)), sprintf("te_re_scale[%d]", seq_len(qtot)),
    sprintf("te_re_z[%d]", seq_len(utot)), "treatment_it_scale[1]",
    sprintf("treatment_kt_z[%d]", seq_len(n_exp))
  )

  parts <- te_mixed_parts(m, comps, mu)
  expect_equal(parts$base, mu + fixed_and_shared, ignore_attr = TRUE)
  expect_equal(parts$scale[, 1], re_scale[, 2])
  expect_equal(parts$scale[, 4], it_scale)
})

test_that("legacy hierarchy keeps only the group effect, in cut's naming", {
  data <- mixed_data()
  ex <- mixed_exposed(data)
  comps <- te_unit_components(NULL, ex)
  expect_equal(
    vapply(comps, `[[`, character(1), "label"),
    c("(1 | unit)", "(1 | group:unit)", "exposed cell")
  )
  n_draws <- 5
  K <- ex$K
  D <- ex$D
  withr::with_seed(12, {
    cat_eff <- matrix(stats::rnorm(n_draws * K), n_draws)
    unit_scale <- stats::rexp(n_draws)
    unit_z <- matrix(stats::rnorm(n_draws * D), n_draws)
    uc_scale <- stats::rexp(n_draws)
    uc_z <- matrix(stats::rnorm(n_draws * K * D), n_draws) # column (d-1)*K+k
    it_scale <- stats::rexp(n_draws, 10)
    kt_z <- matrix(stats::rnorm(n_draws * ex$n_exposed), n_draws)
  })
  te <- cat_eff[, ex$exp_k, drop = FALSE] +
    unit_z[, ex$exp_d, drop = FALSE] * unit_scale +
    uc_z[, ex$exp_kd, drop = FALSE] * uc_scale + kt_z * it_scale
  kd <- seq_len(K * D)
  m <- cbind(te, unit_scale, unit_z, uc_scale, uc_z, it_scale, kt_z)
  # Cut stage 2 declares treatment_it_scale as a scalar, the joint model as
  # an array of one: both spellings must resolve.
  colnames(m) <- c(
    sprintf("te[%d]", seq_len(ex$n_exposed)), "treatment_unit_scale[1]",
    sprintf("unit_treatment_effect_z[%d]", seq_len(D)), "unit_category_scale[1]",
    sprintf("unit_category_te_z[%d,%d]", (kd - 1) %% K + 1, (kd - 1) %/% K + 1),
    "treatment_it_scale", sprintf("treatment_kt_z[%d]", seq_len(ex$n_exposed))
  )
  parts <- te_mixed_parts(m, comps, rep(2, ex$n_exposed))
  expect_equal(parts$base, 2 + cat_eff[, ex$exp_k], ignore_attr = TRUE)
})

# A source with one group, `n_units` treated units of six exposed cells each,
# and a model whose unit-level scale is `tau_model`; observations come from
# unit effects with sd `tau_true`.
synthetic_source <- function(n_units, tau_model, tau_true, n_draws = 400,
                             seed = 21) {
  n_exp <- 6L * n_units
  unit_idx <- rep(seq_len(n_units), each = 6L)
  withr::with_seed(seed, {
    b <- stats::rnorm(n_units, 0, tau_true)
    outcome <- stats::rpois(n_exp, exp(log(1000) + b[unit_idx]))
  })
  list(
    base = matrix(log(1000), n_draws, n_exp),
    scale = cbind(rep(tau_model, n_draws), rep(1e-4, n_draws)),
    phi = NULL,
    components = list(
      list(label = "(1 | unit)", level = unit_idx, n_levels = n_units, mult = NULL),
      list(label = "exposed cell", level = seq_len(n_exp), n_levels = n_exp, mult = NULL)
    ),
    unit = sprintf("u%02d", unit_idx), group = rep("g", n_exp),
    outcome = outcome, exp_d = unit_idx
  )
}

test_that("p-values are spread out when the unit scale is right", {
  src <- synthetic_source(40, tau_model = 0.1, tau_true = 0.1)
  res <- withr::with_seed(1, te_mixed_check(src, seq_along(src$outcome)))
  expect_equal(nrow(res$pvals), 40)
  extreme <- mean(res$pvals$pval < 0.05 | res$pvals$pval > 0.95)
  expect_lt(extreme, 0.25)
  expect_gt(stats::sd(res$pvals$pval), 0.2)
})

test_that("over-shrinkage piles the p-values at 0 and 1", {
  # Unit effects of 10% against a model that allows 0.5%: the plain
  # posterior predictive would still look fine here, because each unit's
  # posterior effect would sit on its data.
  src <- synthetic_source(40, tau_model = 0.005, tau_true = 0.1)
  res <- withr::with_seed(1, te_mixed_check(src, seq_along(src$outcome)))
  extreme <- mean(res$pvals$pval < 0.05 | res$pvals$pval > 0.95)
  expect_gt(extreme, 0.6)
})

test_that("suppressed exposed cells are left out of both totals", {
  src <- synthetic_source(3, tau_model = 0.1, tau_true = 0.1, n_draws = 50)
  src$outcome[1:6] <- NA # all of u01's exposed cells
  src$outcome[7] <- NA # one of u02's
  fake <- structure(list(), class = "fake_mixed_fit")
  local_mocked_bindings(te_mixed_source = function(x) src)
  res <- countsynth_ppc_te_mixed(fake, seed = 3)
  expect_equal(res$pvals$unit, c("u02", "u03"))
  expect_s3_class(res$plot, "ggplot")
  expect_s3_class(res$calibration, "ggplot")
  expect_equal(res$redrawn, c("(1 | unit)", "exposed cell"))
  # Same seed, same p-values: the report must not change between renders.
  again <- countsynth_ppc_te_mixed(fake, seed = 3)
  expect_equal(again$pvals, res$pvals)
})

# The reconstruction the R-side tests above assume, checked against what
# Stan actually sampled.
expected_base <- function(design, te_draws, mu_exp) {
  W <- te_wide_matrix(te_draws)
  n_draws <- nrow(W)
  out <- mu_exp
  if (length(design$x_names) > 0) {
    beta <- W[, te_coef_key("(fixed)", design$x_names, NA), drop = FALSE]
    out <- out + beta %*% t(design$X)
  }
  for (t in design$terms) {
    if (t$label != "group") next
    for (q in seq_along(t$predictors)) {
      u <- W[, te_coef_key(t$label, t$predictors[[q]], t$levels), drop = FALSE]
      out <- out + u[, t$index, drop = FALSE] *
        matrix(t$Z[, q], n_draws, design$n_exposed, byrow = TRUE)
    }
  }
  out
}

test_that("joint fit: base matches the sampled fixed and group effects", {
  skip_on_cran()
  skip_if_no_cmdstan()
  data <- mixed_data()
  opts <- countsynth_te_opts(
    formula = ~ 1 + event_time + (1 | group) + (1 | unit)
  )
  fit <- countsynth_fit(data, rank = 2, config = mixed_config(data, opts))
  src <- te_mixed_source(fit)
  mu_exp <- extract_var_matrix(fit$fit, "mu_ctrl")[, fit$stan_data$exp_cell]
  want <- expected_base(fit$te_design, countsynth_te_draws(fit), mu_exp)
  # CmdStan writes 6 significant digits, so compare relative to the scale.
  expect_lt(max(abs(src$base - want)) / max(abs(want)), 1e-4)
  expect_equal(
    vapply(src$components, `[[`, character(1), "label"),
    c("(1 | unit)", "exposed cell")
  )

  res <- countsynth_ppc_te_mixed(fit)
  expect_setequal(res$pvals$unit, c("C", "D"))
  expect_true(all(res$pvals$pval >= 0 & res$pvals$pval <= 1))

  out <- withr::local_tempdir()
  suppressWarnings(countsynth_report(
    countsynth_draws(fit), out,
    figures = "ppc", fit = fit, print_tables = FALSE, html_tables = FALSE,
    progress = FALSE
  ))
  diag <- file.path(out, "model_diagnostics")
  expect_true(file.exists(file.path(diag, "ppc_te_mixed.png")))
  expect_true(file.exists(file.path(diag, "ppc_te_mixed_calibration.png")))
  pv <- utils::read.csv(file.path(diag, "ppc_pvalues.csv"))
  expect_equal(sum(pv$check_type == "te_mixed"), nrow(res$pvals))
})

test_that("joint fit, legacy hierarchy: base is mu_ctrl plus the group effect", {
  skip_on_cran()
  skip_if_no_cmdstan()
  data <- mixed_data()
  fit <- countsynth_fit(data, rank = 2, config = mixed_config(data, NULL))
  src <- te_mixed_source(fit)
  sd <- fit$stan_data
  want <- extract_var_matrix(fit$fit, "mu_ctrl")[, sd$exp_cell] +
    extract_var_matrix(fit$fit, "category_treatment_effect")[, sd$exp_k]
  expect_lt(max(abs(src$base - want)) / max(abs(want)), 1e-4)
})

test_that("cut fit saves the check's inputs on the pooled draws", {
  skip_on_cran()
  skip_if_no_cmdstan()
  data <- mixed_data()
  opts <- countsynth_te_opts(
    formula = ~ 1 + event_time + (1 | group) + (1 | unit)
  )
  config <- mixed_config(data, opts, inference_mode = "cut")
  config$cut <- countsynth_cut_opts(
    num_stage1_draws = 2, stage2_draws_per_component = 5,
    stage2_mcmc = list(num_warmup = 100, num_samples = 100)
  )
  fit <- suppressWarnings(countsynth_cut_fit(data, rank = 2, config = config))
  expect_equal(nrow(fit$te_mixed$base), length(unique(fit$draws$.draw)))

  # mu_ctrl on the exposed cells, per pooled draw, from the draws frame.
  ex <- mixed_exposed(data)
  cells <- build_cell_table(data)[ex$exp_cell, ]
  d <- fit$draws
  key <- paste(d$unit, d$group, d$time)
  mu_exp <- t(vapply(sort(unique(d$.draw)), function(i) {
    rows <- d$.draw == i
    d$mu[rows][match(paste(cells$unit, cells$group, cells$time), key[rows])]
  }, numeric(ex$n_exposed)))
  want <- expected_base(fit$te_design, countsynth_te_draws(fit), mu_exp)
  expect_lt(max(abs(fit$te_mixed$base - want)) / max(abs(want)), 1e-4)

  res <- countsynth_ppc_te_mixed(fit)
  expect_setequal(res$pvals$unit, c("C", "D"))
})
