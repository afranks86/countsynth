# The stanli backend has to be a drop-in: the same Stan programs, the same
# data, and fits that every downstream reader (draws frame, gate, cut stage
# 2) accepts unchanged.

skip_if_no_stanli <- function() {
  testthat::skip_if_not_installed("stanli")
  testthat::skip_if_not(stanli::stanli_available(), "stanli runtime missing")
}

stanli_test_config <- function(tmp, ...) {
  df <- make_test_data()
  path <- write_test_csv(df, file.path(tmp, "panel.csv"))
  config <- make_test_config(path, tmp, ...)
  config$mcmc$backend <- "stanli"
  config
}

test_that("both backends compute the same log density and gradient", {
  skip_on_cran()
  skip_if_no_cmdstan()
  skip_if_no_stanli()
  tmp <- withr::local_tempdir()
  config <- stanli_test_config(tmp)
  data <- countsynth_data(config, type = "both")
  sd <- stan_data_joint(data, rank = 2)
  json <- withr::local_tempfile(fileext = ".json")
  cmdstanr::write_stan_json(sd, json)

  sl <- stanli::stanli_model(
    file = stan_file_path("joint"), data = json,
    include_paths = system.file("stan", "include", package = "countsynth")
  )
  q <- withr::with_seed(1, stats::runif(sl$n_unconstrained, -1, 1))
  got <- stanli::log_prob_grad(sl, q)

  # CmdStan's standalone log_prob method at the same unconstrained point.
  exe <- countsynth_stan_model("joint")$exe_file()
  qfile <- withr::local_tempfile(fileext = ".json")
  out <- withr::local_tempfile(fileext = ".csv")
  jsonlite::write_json(list(params_r = q), qfile, digits = NA)
  status <- system2(exe, c(
    "log_prob", paste0("unconstrained_params=", qfile), "jacobian=1",
    "data", paste0("file=", json),
    "output", paste0("file=", out), "sig_figs=18"
  ), stdout = FALSE, stderr = FALSE)
  expect_equal(status, 0)
  want <- unlist(utils::read.csv(out, comment.char = "#")[1, ])

  expect_equal(got$lp, unname(want[1]), tolerance = 1e-10)
  expect_equal(got$grad, unname(want[-1]), tolerance = 1e-10)
})

test_that("a stanli joint fit flows through draws and the gate", {
  skip_on_cran()
  skip_if_no_stanli()
  tmp <- withr::local_tempdir()
  config <- stanli_test_config(tmp)
  data <- countsynth_data(config, type = "both")
  fit <- suppressWarnings(countsynth_fit(data, rank = 2, config = config))

  expect_s3_class(fit$fit, "stanli_cstanfit")
  draws <- countsynth_draws(fit)
  expect_equal(length(unique(draws$.draw)), 2L * 150L)
  gate <- convergence_gate(fit)
  expect_true(is.logical(gate$converged))
  expect_identical(fit$fit$metadata()$max_treedepth, 10L)
})

test_that("cut mode runs both stages on stanli", {
  skip_on_cran()
  skip_if_no_stanli()
  tmp <- withr::local_tempdir()
  config <- stanli_test_config(tmp, inference_mode = "cut")
  config$cut <- countsynth_cut_opts(
    num_stage1_draws = 2, stage2_draws_per_component = 5,
    stage2_mcmc = list(num_warmup = 100, num_samples = 100)
  )
  data <- countsynth_data(config, type = "both")
  fit <- suppressWarnings(countsynth_cut_fit(data, rank = 2, config = config))

  expect_s3_class(fit$stage1$fit, "stanli_cstanfit")
  expect_equal(length(unique(fit$draws$cut_component)), 2L)
  # The per-component quota is per chain: 2 components x 2 chains x 5.
  expect_equal(length(unique(fit$draws$.draw)), 2L * 2L * 5L)
})

test_that("unknown sampler options are refused, not ignored", {
  skip_if_no_stanli()
  model <- countsynth_stan_model("joint", backend = "stanli")
  expect_error(model$sample(data = list(), sig_figs = 18), "sig_figs")
})

test_that("mcmc.backend is a validated config field, CmdStan by default", {
  expect_identical(countsynth_mcmc_opts()$backend, "cmdstan")
  expect_identical(parse_yaml_mcmc(list(backend = "stanli"))$backend, "stanli")
  expect_identical(parse_yaml_mcmc(list())$backend, "cmdstan")
  expect_error(countsynth_mcmc_opts(backend = "numpyro"), "should be one of")
})
