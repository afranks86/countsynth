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
      class = c("bpnmf_data", "list")
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
  # adjustment set = all cells.
  sd <- stan_data_joint(data, rank = 2, model_treated = TRUE)
  expect_equal(as.integer(sd$obs_cell), which(!missing_flat))
  expect_equal(as.integer(sd$cens_cell), which(missing_flat))
  expect_equal(as.integer(sd$notcens_cell), which(!missing_flat))
  expect_equal(as.integer(sd$exp_cell), which(!control_flat))
  expect_equal(sd$n_exposed, sum(!control_flat))

  # Exposed-cell order must be ascending flat order (numpy reshape(-1)).
  fx_exposed <- jsonlite::fromJSON(fixture_path("flatten_parity.json"))$exposed_cells_1based
  expect_equal(as.integer(sd$exp_cell), as.integer(fx_exposed))

  # model_treated = FALSE: likelihood and adjustment on control cells only.
  sd0 <- stan_data_joint(data, rank = 2, model_treated = FALSE)
  expect_equal(as.integer(sd0$obs_cell), which(!missing_flat & control_flat))
  expect_equal(as.integer(sd0$cens_cell), which(missing_flat & control_flat))
  expect_equal(as.integer(sd0$notcens_cell), which(!missing_flat & control_flat))

  # exp_k/exp_d agree with kdn_from_flat.
  sub <- kdn_from_flat(as.integer(sd$exp_cell), sd$D, sd$N)
  expect_identical(as.integer(sd$exp_k), sub$k)
  expect_identical(as.integer(sd$exp_d), sub$d)

  # cell_unit covers every cell
  expect_equal(length(sd$cell_unit), KDN)
  expect_identical(as.integer(sd$cell_unit), kdn_from_flat(seq_len(KDN), sd$D, sd$N)$d)
})

test_that("stan_data_stage2 subsets are within-exposed positions", {
  fd <- fixture_data()
  data <- fd$data
  mu_flat <- flatten_kdn(array(fd$fx$data$mu_ctrl_fixed, dim = dim(data$Y)))
  sd2 <- stan_data_stage2(
    data, mu_flat,
    phi_unit = fd$fx$data$phi_fixed, outcome_distribution = "NB"
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
  expect_error(
    stan_data_stage2(data, mu_flat, phi_unit = NULL, outcome_distribution = "NB"),
    "phi_unit"
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
    gamma_shape_from_pct(bpnmf_model_opts()$factor_variation_pct), 20
  )
  expect_equal(
    gamma_shape_from_pct(
      bpnmf_model_opts()$time_level_variation_pct, DEFAULT_TIME_FE_SHAPE
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
    rank_shrinkage = bpnmf_rank_shrinkage_opts(
      group_mass_prior = c(3, 1.5), unit_sd_prior = 0.5
    )
  )
  expect_equal(on$rank_shrink, 1L)
  expect_equal(on$group_mass_shape, 3)
  expect_equal(on$group_mass_rate, 1.5)
  expect_equal(on$unit_sd_scale, 0.5)

  # Adding the fields must not disturb the rest of the data list, which the
  # log-density parity fixture depends on.
  expect_equal(off[names(off) != "rank_shrink"][["y"]], on[["y"]])
})

test_that("rank_shrinkage is disabled at rank 1 rather than aborting", {
  data <- fixture_data()$data
  # A ranks_to_test sweep including 1 must still run: with one component
  # there is nothing to shrink and both hyperparameters would be sampled
  # straight from their priors.
  expect_warning(
    sd <- stan_data_joint(
      data, rank = 1, rank_shrinkage = bpnmf_rank_shrinkage_opts()
    ),
    "ignored at rank 1"
  )
  expect_equal(sd$rank_shrink, 0L)
})
