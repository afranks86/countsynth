# The common time level is centered by default: log(time_fe) sums to zero over
# periods within each group, so unit_fe_mu owns the level. The checks that
# matter are that the constraint actually holds in the draws, and that the
# centered prior has the marginal spread time_level_variation_pct promises --
# the sum-to-zero constraint shrinks each element's variance by (N - 1) / N,
# which the model has to undo.

tl_data <- function() {
  df <- make_test_data()
  arrays <- build_model_arrays(df, c("g1", "g2"))
  arrays$df <- df
  arrays$type <- "both"
  structure(arrays, class = c("countsynth_data", "list"))
}

test_that("time_level reaches the Stan data, centered by default", {
  data <- tl_data()
  sd <- stan_data_joint(data, rank = 2)
  expect_equal(sd$center_time, 1L)
  # One knob for both parameterizations: the centered sd is the Gamma's own
  # sd(log), so the historical Gamma(1, 1) maps to 1.28.
  expect_equal(sd$time_fe_sd, sqrt(trigamma(1)))
  expect_equal(
    stan_data_joint(data, rank = 2, time_level = "uncentered")$center_time, 0L
  )
  # And a requested swing means the same spread either way.
  shape <- gamma_shape_from_pct(40, DEFAULT_TIME_FE_SHAPE)
  expect_equal(
    stan_data_joint(data, rank = 2, time_fe_shape = shape)$time_fe_sd,
    log1p(0.4)
  )
  expect_error(stan_data_joint(data, rank = 2, time_level = "pinned"))
})

test_that("time_level validates and defaults to centered", {
  expect_equal(countsynth_model_opts()$time_level, "centered")
  expect_equal(
    countsynth_model_opts(time_level = "uncentered")$time_level, "uncentered"
  )
  expect_error(countsynth_model_opts(time_level = "pinned"), "should be one of")
})

test_that("the centered time level sums to zero and has the promised spread", {
  skip_on_cran()
  skip_if_no_cmdstan()
  # Prior only: every likelihood term removed, so the draws are the prior.
  sd <- stan_data_joint(
    tl_data(), rank = 2, model_treated = FALSE,
    adjust_for_missingness = FALSE, gen_ypred = FALSE
  )
  empty <- as.array(integer(0))
  for (f in c("obs_cell", "y", "obs_unit", "cens_cell", "cens_unit",
              "notcens_cell", "notcens_unit")) {
    sd[[f]] <- empty
  }
  sd$n_obs <- 0L
  sd$n_cens <- 0L
  sd$n_notcens <- 0L
  fit <- countsynth_stan_model("joint")$sample(
    data = sd, chains = 1, iter_warmup = 300, iter_sampling = 1500, seed = 5,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE
  )
  lt <- log(posterior::draws_of(
    posterior::as_draws_rvars(fit$draws(variables = "time_fe"))$time_fe
  ))                                             # [draw, n, k]

  # The constraint, to CmdStan's CSV precision.
  expect_lt(max(abs(apply(lt, c(1, 3), sum))), 1e-5)
  # The marginal sd of each period's log level. Without the sqrt(N / (N - 1))
  # correction this would come out at 1.28 * sqrt(4 / 5) = 1.15 for N = 5,
  # so the tolerance separates the two comfortably.
  expect_equal(mean(apply(lt, c(2, 3), stats::sd)), sd$time_fe_sd, tolerance = 0.05)
})
