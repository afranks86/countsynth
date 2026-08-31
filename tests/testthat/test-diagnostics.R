test_that("convergence status bands match the Python gate", {
  th <- bpnmf_convergence() # 1.01 / 1.05 / 400 / 0.25
  expect_equal(convergence_status(1.0, 1000, th), "PASS")
  expect_equal(convergence_status(1.009, 400, th), "PASS")
  expect_equal(convergence_status(1.01, 1000, th), "WARN")
  expect_equal(convergence_status(1.0, 399, th), "WARN")
  expect_equal(convergence_status(1.049, 100, th), "WARN")
  expect_equal(convergence_status(1.05, 1000, th), "FAIL")
  expect_equal(convergence_status(1.0, 99, th), "FAIL")
})

test_that("custom thresholds are honored", {
  th <- bpnmf_convergence(
    rhat_warn = 1.1, rhat_fail = 1.2, ess_min = 100, ess_fail_fraction = 0.5
  )
  expect_equal(convergence_status(1.15, 1000, th), "WARN")
  expect_equal(convergence_status(1.05, 60, th), "WARN")
  expect_equal(convergence_status(1.05, 49, th), "FAIL")
})

test_that("divergences are summarized as a rate over retained transitions", {
  fake_fit <- function(divergent) {
    sd <- posterior::as_draws_array(array(
      divergent,
      dim = c(length(divergent) / 2L, 2L, 1L),
      dimnames = list(NULL, NULL, "divergent__")
    ))
    list(sampler_diagnostics = function() sd)
  }
  d <- divergence_summary(fake_fit(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0)))
  expect_equal(d$count, 1)
  expect_equal(d$transitions, 10)
  expect_equal(d$fraction, 0.1)

  clean <- divergence_summary(fake_fit(rep(0, 10)))
  expect_equal(clean$count, 0)
  expect_equal(clean$fraction, 0)
})

test_that("divergence_fail_fraction sets the gate's tolerance", {
  expect_equal(bpnmf_convergence()$divergence_fail_fraction, 0.01)
  # 0 restores the old rule: any divergence at all fails.
  strict <- bpnmf_convergence(divergence_fail_fraction = 0)
  expect_true(0 <= strict$divergence_fail_fraction)
  expect_false(1e-9 <= strict$divergence_fail_fraction)
  expect_error(bpnmf_convergence(divergence_fail_fraction = 1.5))
})

test_that("gate_params prefix matching errors on zero matches", {
  expect_error(
    match_gate_params(c("mu_ctrl", "te"), "bogus"),
    "matched no posterior variables"
  )
  gated <- match_gate_params(c("mu_ctrl", "te", "state_fe_z"), c("mu", "te"))
  expect_identical(gated, c(TRUE, TRUE, FALSE))
})

test_that("a variational fit is reported as ungated, not as failed", {
  # `converged = NA` is the whole point: ADVI produces one stream of draws
  # from an approximating family, so R-hat / ESS / divergences do not exist
  # and must not be silently reported as either a pass or a failure.
  g <- variational_gate()
  expect_identical(g$method, "variational")
  expect_true(is.na(g$converged))
  expect_false(isTRUE(g$converged))
  expect_false(isFALSE(g$converged))
  expect_true(all(vapply(
    g[c("rhat_max", "ess_bulk_min", "ess_tail_min", "divergences",
        "divergence_fraction")],
    is.na, logical(1)
  )))
})

test_that("NA gate fields round-trip through the convergence JSON as null", {
  path <- withr::local_tempfile(fileext = ".json")
  write_convergence_json(variational_gate(), path)
  back <- jsonlite::read_json(path)
  expect_null(back$converged)
  expect_null(back$rhat_max)
  expect_identical(back$method, "variational")
})

test_that("is_variational_fit sees through the bpnmf_fit wrapper", {
  vb <- structure(list(), class = c("CmdStanVB", "CmdStanFit", "R6"))
  mcmc <- structure(list(), class = c("CmdStanMCMC", "CmdStanFit", "R6"))
  expect_true(is_variational_fit(vb))
  expect_false(is_variational_fit(mcmc))
  expect_true(is_variational_fit(structure(list(fit = vb), class = "bpnmf_fit")))
  expect_false(is_variational_fit(structure(list(fit = mcmc), class = "bpnmf_fit")))
})

test_that("chain-based diagnostics refuse a variational fit", {
  vb <- structure(
    list(config = bpnmf_config(
      input_file = "x.csv", output_dir = tempdir(),
      schema = bpnmf_schema("u", "t", "tr", outcomes = list(bpnmf_outcome("o", "g"))),
      model = bpnmf_model_opts(types = list(g = bpnmf_type("g", 2)))
    )),
    class = "bpnmf_fit"
  )
  vb$fit <- structure(list(), class = c("CmdStanVB", "CmdStanFit", "R6"))
  expect_error(parameter_diagnostics(vb), "need multiple MCMC chains")
  expect_error(bpnmf_trace_plot(vb), "need MCMC chains")
  expect_true(is.na(convergence_gate(vb)$converged))
})

test_that("ADVI tuning arguments are checked against the settable knobs", {
  expect_identical(check_variational_args(NULL), list())
  expect_identical(
    check_variational_args(list(eta = 0.1)), list(eta = 0.1)
  )
  expect_error(check_variational_args(list(adapt_delta = 0.9)), "Unknown ADVI")
  expect_error(check_variational_args(list(seed = 1)), "mcmc.seed")
})

test_that("only run-specific ADVI complaints warn, not CmdStan's standard banner", {
  banner <- c(
    "This is Automatic Differentiation Variational Inference.",
    "(EXPERIMENTAL ALGORITHM: expect frequent deprecations).",
    "This procedure has not been thoroughly tested and may be unstable."
  )
  quiet <- list(output = function() cat(banner, sep = "\n"))
  expect_silent(warn_advi_quality(quiet))

  noisy <- list(output = function() {
    cat(c(banner, "Informational Message: The ELBO at a previous iteration is",
          "larger than the ELBO upon convergence!",
          "The variational approximation may be poor."), sep = "\n")
  })
  expect_warning(warn_advi_quality(noisy), "may be poor")
})

test_that("gate failure bullets name the criterion that tripped", {
  th <- bpnmf_convergence(
    rhat_warn = 1.1, rhat_fail = 1.5, ess_min = 100,
    ess_fail_fraction = 0.25, divergence_fail_fraction = 0.01
  )
  # The fertility_yearly failure: R-hat and ESS both blown, no divergences.
  gate <- list(
    rhat_max = 1.8258, ess_bulk_min = 5.865, ess_tail_min = 15.886,
    divergences = 0, divergence_fraction = 0, converged = FALSE
  )
  b <- gate_failure_bullets(gate, th)
  expect_length(b, 2)
  expect_match(b[1], "max R-hat 1.83 \\(fails at 1.5\\)")
  expect_match(b[2], "min ESS 5.87 \\(fails below 25")

  # Divergences alone.
  div <- list(
    rhat_max = 1.01, ess_bulk_min = 900, ess_tail_min = 900,
    divergences = 40L, divergence_fraction = 0.02, converged = FALSE
  )
  db <- gate_failure_bullets(div, th)
  expect_length(db, 1)
  expect_match(db, "40 divergences = 2.00% \\(fails above 1.00%\\)")

  # A WARN-level R-hat with ESS above the fail floor trips no branch; the
  # numbers still have to be reported rather than an empty bullet list.
  edge <- list(
    rhat_max = 1.2, ess_bulk_min = 60, ess_tail_min = 60,
    divergences = 0, divergence_fraction = 0, converged = FALSE
  )
  eb <- gate_failure_bullets(edge, th)
  expect_length(eb, 1)
  expect_match(eb, "max R-hat 1.2, min ESS 60, 0 divergences")

  # ESS uses the smaller of bulk and tail.
  tail_bad <- list(
    rhat_max = 1.01, ess_bulk_min = 900, ess_tail_min = 10,
    divergences = 0, divergence_fraction = 0, converged = FALSE
  )
  expect_match(gate_failure_bullets(tail_bad, th), "min ESS 10")
})

test_that("gate_worst_parameters is silent without a fit", {
  expect_equal(gate_worst_parameters(NULL), character())
})
