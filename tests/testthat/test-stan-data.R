fixture_data <- function() {
  fx <- jsonlite::fromJSON(fixture_path("logdensity_parity.json"))
  d <- fx$data
  K <- d$K
  D <- d$D
  N <- d$N
  list(
    fx = fx,
    data = structure(
      list(
        Y = array(d$Y, dim = c(K, D, N)),
        denominators = array(d$denominators, dim = c(K, D, N)),
        control_idx_array = array(d$control == 1, dim = c(K, D, N)),
        missing_idx_array = array(d$missing == 1, dim = c(K, D, N)),
        groups = paste0("g", seq_len(K)),
        units = LETTERS[seq_len(D)],
        times = seq(as.Date("2020-01-01"), by = "month", length.out = N),
        df = NULL,
        type = "test"
      ),
      class = c("countsynth_data", "list")
    )
  )
}

test_that("stan_data_joint index sets follow the Python mask logic", {
  fd <- fixture_data()
  data <- fd$data
  KDN <- prod(dim(data$Y))
  missing_flat <- flatten_kdn(data$missing_idx_array)
  control_flat <- flatten_kdn(data$control_idx_array)

  # model_treated = TRUE: likelihood on all non-missing cells; censoring
  # adjustment set = all cells. Observed cells carry no censoring term of
  # their own (the Python model's "not in 1..9" term counted it twice).
  sd <- stan_data_joint(data, rank = 2, model_treated = TRUE)
  expect_equal(as.integer(sd$obs_cell), which(!missing_flat))
  expect_equal(as.integer(sd$cens_cell), which(missing_flat))
  expect_null(sd$notcens_cell)
  expect_equal(as.integer(sd$exp_cell), which(!control_flat))
  expect_equal(sd$n_exposed, sum(!control_flat))

  # Exposed-cell order must be ascending flat order (numpy reshape(-1)).
  fx_exposed <- jsonlite::fromJSON(fixture_path("flatten_parity.json"))$exposed_cells_1based
  expect_equal(as.integer(sd$exp_cell), as.integer(fx_exposed))

  # model_treated = FALSE: likelihood and adjustment on control cells only.
  sd0 <- stan_data_joint(data, rank = 2, model_treated = FALSE)
  expect_equal(as.integer(sd0$obs_cell), which(!missing_flat & control_flat))
  expect_equal(as.integer(sd0$cens_cell), which(missing_flat & control_flat))

  # exp_k/exp_d agree with kdn_from_flat.
  sub <- kdn_from_flat(as.integer(sd$exp_cell), sd$D, sd$N)
  expect_identical(as.integer(sd$exp_k), sub$k)
  expect_identical(as.integer(sd$exp_d), sub$d)

  # Dispersion levels cover every cell: groups by default, units on request.
  sub_all <- kdn_from_flat(seq_len(KDN), sd$D, sd$N)
  expect_equal(sd$n_disp, sd$K)
  expect_identical(as.integer(sd$cell_disp), sub_all$k)
  expect_identical(as.integer(sd$obs_disp), sub_all$k[sd$obs_cell])
  expect_identical(as.integer(sd$cens_disp), sub_all$k[sd$cens_cell])
  sdu <- stan_data_joint(data, rank = 2, model_treated = TRUE, disp_by = "unit")
  expect_equal(sdu$n_disp, sdu$D)
  expect_identical(as.integer(sdu$cell_disp), sub_all$d)
  expect_identical(as.integer(sdu$obs_disp), sub_all$d[sdu$obs_cell])
  expect_identical(as.integer(sdu$cens_disp), sub_all$d[sdu$cens_cell])
})

test_that("disp_prior_rate enters the log density as -rate * sum(sqrt(disp))", {
  skip_on_cran()
  skip_if_not_installed("stanli")
  skip_if_not(stanli::stanli_available(), "stanli runtime missing")
  data <- fixture_data()$data
  model_at <- function(rate) {
    test_env <- parent.frame()
    sd <- stan_data_joint(data, rank = 2, sample_disp = TRUE, disp_prior_rate = rate)
    json <- withr::local_tempfile(fileext = ".json", .local_envir = test_env)
    cmdstanr::write_stan_json(sd, json)
    stanli::stanli_model(
      file = stan_file_path("joint"), data = json,
      include_paths = system.file("stan", "include", package = "countsynth")
    )
  }
  expect_equal(stan_data_joint(data, rank = 2)$disp_prior_rate, 100)
  m100 <- model_at(100)
  m10 <- model_at(10)
  # One point, dispersion pinned, everything else random but shared.
  disp <- c(0.01, 0.04)
  q <- stanli::unconstrain(m100, list(disp = disp), partial = TRUE, seed = 1)
  shift <- stanli::log_prob_grad(m10, q)$lp - stanli::log_prob_grad(m100, q)$lp
  expect_equal(shift, (100 - 10) * sum(sqrt(disp)), tolerance = 1e-10)
})

# The joint program under stanli on the fixture panel with a sampled
# dispersion, at one random point where only the dispersion coordinates move:
# eval_at(u) sets every one of them to u on the unconstrained scale, u_at(d)
# is the u that puts the dispersion at d.
stanli_disp_point <- function(env = parent.frame()) {
  sd <- stan_data_joint(fixture_data()$data, rank = 2, sample_disp = TRUE)
  json <- withr::local_tempfile(fileext = ".json", .local_envir = env)
  cmdstanr::write_stan_json(sd, json)
  m <- stanli::stanli_model(
    file = stan_file_path("joint"), data = json,
    include_paths = system.file("stan", "include", package = "countsynth")
  )
  q_at <- function(d) {
    stanli::unconstrain(m, list(disp = rep(d, sd$n_disp)), partial = TRUE, seed = 1)
  }
  q <- q_at(1e-6)
  idx <- which(q != q_at(1e-7))
  list(
    sd = sd,
    idx = idx,
    u_at = function(d) q_at(d)[idx[1]],
    eval_at = function(u) {
      q[idx] <- u
      stanli::log_prob_grad(m, q)
    }
  )
}

test_that("the censored NB mass stays exact as a sampled dispersion nears its floor", {
  skip_on_cran()
  skip_if_not_installed("stanli")
  skip_if_not(stanli::stanli_available(), "stanli runtime missing")
  p <- stanli_disp_point()
  sd <- p$sd
  expect_gt(sd$n_cens, 0)
  expect_length(p$idx, sd$n_disp)
  floor <- sd$disp_floor
  disp <- c(1e-8, 1e-9, 5e-10, 2e-10, 1.1e-10, 1.01e-10)
  lp <- vapply(disp, function(d) p$eval_at(p$u_at(d))$lp, numeric(1))
  # Less the dispersion's own prior and the Jacobian of its logit onto
  # (floor, 1), the log density is the Poisson's plus the NB's first-order
  # departure from it, a straight line in disp this close to zero. Before the
  # stable forms, phi past ~1e10 scattered it over several nats per censored
  # cell; the floor stops phi at 1e10.
  rest <- lp - sd$n_disp * (-0.5 * log(disp) - sd$disp_prior_rate * sqrt(disp) +
    log(disp - floor) + log1p(-disp) - log1p(-floor))
  expect_lt(max(abs(stats::resid(stats::lm(rest ~ disp)))), 1e-6)
})

test_that("the gradient on a sampled dispersion stays exact down to its floor", {
  skip_on_cran()
  skip_if_not_installed("stanli")
  skip_if_not(stanli::stanli_available(), "stanli runtime missing")
  p <- stanli_disp_point()
  # Logit-scale points that, unfloored, put disp at 1e-6 ... 1e-20. Without
  # the floor the NB2 phi-partials cancel to rounding past phi ~1e12 and the
  # chain rule multiplies them by phi: autodiff here read 72 at 1e-16 and
  # 2,347 from 1e-17 on, against a finite difference of 1.
  u <- stats::qlogis(10^-(6:20))
  h <- 1e-4
  err <- vapply(u, function(v) {
    autodiff <- sum(p$eval_at(v)$grad[p$idx])
    finite <- (p$eval_at(v + h)$lp - p$eval_at(v - h)$lp) / (2 * h)
    autodiff - finite
  }, numeric(1))
  expect_lt(max(abs(err)), 1e-3)
})

test_that("a sampled dispersion starts near the Poisson, apart per chain", {
  data <- fixture_data()$data
  expect_null(disp_init(stan_data_joint(data, rank = 2), seed = 1))
  expect_null(disp_init(
    stan_data_joint(data, rank = 2, outcome_distribution = "Poisson"), seed = 1
  ))

  sd <- stan_data_joint(data, rank = 2, sample_disp = TRUE, disp_prior_rate = 20)
  init <- disp_init(sd, seed = 7)
  starts <- lapply(1:4, init)
  for (s in starts) {
    expect_named(s, "disp")
    expect_length(s$disp, sd$n_disp)
    expect_true(all(s$disp >= exp(-1) / 400 & s$disp <= exp(1) / 400))
  }
  expect_false(identical(starts[[1]]$disp, starts[[2]]$disp))
  expect_identical(init(3), disp_init(sd, seed = 7)(3))

  # One level stays a length-1 array, which write_stan_json keeps a vector.
  sd1 <- sd
  sd1$n_disp <- 1L
  expect_true(is.array(disp_init(sd1, seed = 1)(1)$disp))

  # A rate strong enough to start below the floor starts 10x above it.
  strong <- stan_data_joint(data, rank = 2, sample_disp = TRUE, disp_prior_rate = 1e6)
  expect_true(all(disp_init(strong, seed = 1)(1)$disp == 10 * strong$disp_floor))
})

test_that("countsynth_fit starts each chain's dispersion from disp_init", {
  skip_on_cran()
  skip_if_no_cmdstan()
  tmp <- withr::local_tempdir()
  config <- make_test_config(write_test_csv(make_test_data(), file.path(tmp, "panel.csv")), tmp)
  config$model$sample_disp <- TRUE
  data <- countsynth_data(config, type = "both")
  fit <- suppressWarnings(countsynth_fit(data, rank = 2, config = config))
  used <- fit$fit$init()
  expect_length(used, 2)
  want <- disp_init(fit$stan_data, config$mcmc$seed)
  for (c in 1:2) {
    expect_equal(as.numeric(used[[c]]$disp), as.numeric(want(c)$disp), tolerance = 1e-8)
  }
})

test_that("stan_data_stage2 subsets are within-exposed positions", {
  fd <- fixture_data()
  data <- fd$data
  mu_flat <- flatten_kdn(array(fd$fx$data$mu_ctrl_fixed, dim = dim(data$Y)))
  # A distinct concentration per cell, so the exposed subset is checkable.
  phi_flat <- seq_along(mu_flat) + 0.5
  sd2 <- stan_data_stage2(
    data, mu_flat,
    phi_flat = phi_flat, outcome_distribution = "NB"
  )
  control_flat <- flatten_kdn(data$control_idx_array)
  missing_flat <- flatten_kdn(data$missing_idx_array)
  exp_cell <- which(!control_flat)
  expect_identical(as.integer(sd2$exp_cell), exp_cell)
  expect_identical(
    as.integer(sd2$obs_e), which(!missing_flat[exp_cell])
  )
  expect_identical(
    as.integer(sd2$cens_e), which(missing_flat[exp_cell])
  )
  expect_identical(as.numeric(sd2$phi_exposed), phi_flat[exp_cell])
  expect_error(
    stan_data_stage2(data, mu_flat, phi_flat = NULL, outcome_distribution = "NB"),
    "phi_flat"
  )
})

test_that("validate_cut_data catches empty stage likelihoods", {
  fd <- fixture_data()
  data <- fd$data
  expect_silent(validate_cut_data(data))
  all_control <- data
  all_control$control_idx_array[] <- TRUE
  expect_error(validate_cut_data(all_control), "stage 2")
})

test_that("time_fac_shape reaches the Stan data with the historical default", {
  data <- fixture_data()$data
  expect_equal(stan_data_joint(data, rank = 2)$time_fac_shape, 20)
  expect_equal(
    stan_data_joint(data, rank = 2, time_fac_shape = 6.08)$time_fac_shape, 6.08
  )
  # The config route: model opts -> shape, default and explicit.
  expect_equal(stan_data_joint(data, rank = 2)$time_fe_shape, 1)
  expect_equal(
    stan_data_joint(data, rank = 2, time_fe_shape = 10)$time_fe_shape, 10
  )
  expect_equal(
    gamma_shape_from_pct(countsynth_model_opts()$factor_variation_pct), 20
  )
  expect_equal(
    gamma_shape_from_pct(
      countsynth_model_opts()$time_level_variation_pct, DEFAULT_TIME_FE_SHAPE
    ),
    1
  )
})

test_that("rank_shrinkage reaches the Stan data and is off by default", {
  data <- fixture_data()$data
  off <- stan_data_joint(data, rank = 4)
  expect_equal(off$rank_shrink, 0L)
  # Placeholders, not meaningful values: with rank_shrink = 0 every parameter
  # the block declares is zero-size, so Stan never reads them. They must
  # still be present and positive or the data list would fail validation.
  expect_true(all(vapply(
    off[c("group_mass_shape", "group_mass_rate", "unit_sd_scale")],
    function(x) is.numeric(x) && x > 0, logical(1)
  )))

  on <- stan_data_joint(
    data, rank = 4,
    rank_shrinkage = countsynth_rank_shrinkage_opts(
      group_mass_prior = c(3, 1.5), unit_sd_prior = 0.5
    )
  )
  expect_equal(on$rank_shrink, 1L)
  expect_equal(on$group_mass_shape, 3)
  expect_equal(on$group_mass_rate, 1.5)
  expect_equal(on$unit_sd_scale, 0.5)

  # Adding the fields must not disturb the rest of the data list.
  expect_equal(off[names(off) != "rank_shrink"][["y"]], on[["y"]])
})

test_that("rank_shrinkage is disabled at rank 1 rather than aborting", {
  data <- fixture_data()$data
  # A ranks_to_test sweep including 1 must still run: with one component
  # there is nothing to shrink and both hyperparameters would be sampled
  # straight from their priors.
  expect_warning(
    sd <- stan_data_joint(
      data, rank = 1, rank_shrinkage = countsynth_rank_shrinkage_opts()
    ),
    "ignored at rank 1"
  )
  expect_equal(sd$rank_shrink, 0L)
})
