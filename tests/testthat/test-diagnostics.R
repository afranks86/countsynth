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

test_that("the gate defaults to mu_ctrl and te, with an 'all' opt-out", {
  # Every real Stan variable, so an over-matching prefix would show up here.
  bases <- c(
    "mu_ctrl", "te", "time_fac", "time_fe", "unit_weight", "disp", "phi_unit",
    "state_fe_mu", "state_fe_sigma", "state_fe_z", "state_category_scale",
    "treatment_it_scale", "treatment_state_scale", "treatment_category_scale",
    "treatment_kt_z", "state_treatment_effect_z", "state_category_te_z",
    "category_treatment_effect"
  )
  gated <- match_gate_params(bases, DEFAULT_GATE_PARAMS)
  # "te" must not reach treatment_* (which begins "tr"), or the default would
  # silently pull the factor-adjacent scales back into the gate.
  expect_identical(bases[gated], c("mu_ctrl", "te"))

  expect_true(all(match_gate_params(bases, "all")))
  expect_true(all(match_gate_params(bases, NULL)))
})

test_that("every config path lands on the default gate_params", {
  expect_identical(bpnmf_mcmc_opts()$gate_params, c("mu_ctrl", "te"))
  # The YAML loader passes NULL explicitly when the key is absent, which would
  # otherwise skip the signature default.
  expect_identical(bpnmf_mcmc_opts(gate_params = NULL)$gate_params, c("mu_ctrl", "te"))
  expect_identical(bpnmf_mcmc_opts(gate_params = "all")$gate_params, "all")
  expect_identical(bpnmf_mcmc_opts(gate_params = "te")$gate_params, "te")
  expect_identical(parse_yaml_mcmc(NULL)$gate_params, c("mu_ctrl", "te"))
  expect_identical(parse_yaml_mcmc(list())$gate_params, c("mu_ctrl", "te"))
  expect_identical(
    parse_yaml_mcmc(list(gate_params = list("te")))$gate_params, "te"
  )
  expect_identical(parse_yaml_mcmc(list(gate_params = "all"))$gate_params, "all")
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

test_that("gate failure bullets name the criterion and the level it reached", {
  th <- bpnmf_convergence(
    rhat_warn = 1.1, rhat_fail = 1.5, ess_min = 100,
    ess_fail_fraction = 0.25, divergence_fail_fraction = 0.01
  )
  gate <- function(rhat, bulk, tail, div, frac) {
    list(rhat_max = rhat, ess_bulk_min = bulk, ess_tail_min = tail,
         divergences = div, divergence_fraction = frac, converged = FALSE)
  }

  # fertility_yearly: both criteria blown past their fail thresholds.
  hard <- gate_failure_bullets(gate(1.8258, 5.865, 15.886, 0L, 0), th)
  expect_match(hard[1], "max R-hat 1.83, at or above the fail threshold 1.5")
  expect_match(hard[2], "min ESS 5.87, below the fail floor 25")
  expect_false(any(grepl("clean PASS", hard)))

  # infant_mortality/total: nothing reaches a fail threshold, but ESS is under
  # ess_min, and the gate demands a clean PASS -- so it fails. Saying only
  # "FAILED" next to an R-hat of 1.05 explains nothing.
  soft <- gate_failure_bullets(gate(1.0531, 71.05, 506.8, 2L, 0.0005), th)
  expect_false(any(grepl("R-hat", soft))) # 1.05 is under the warn threshold
  expect_match(soft[1], "min ESS 71, below ess_min 100")
  expect_match(soft[2], "2 divergences = 0.05%, within the 1.00% allowance")
  expect_match(soft[3], "requires a clean PASS")

  # infant_mortality/race: R-hat exactly at the warn threshold.
  warn <- gate_failure_bullets(gate(1.1029, 27.21, 118.5, 0L, 0), th)
  expect_match(warn[1], "max R-hat 1.1, at or above the warn threshold 1.1")
  expect_match(warn[2], "min ESS 27.2, below ess_min 100")
  expect_match(warn[3], "requires a clean PASS")

  # Divergences past their own threshold are a fail, not a warn.
  div <- gate_failure_bullets(gate(1.01, 900, 900, 40L, 0.02), th)
  expect_match(div[1], "40 divergences = 2.00%, above divergence_fail_fraction 1.00%")
  expect_false(any(grepl("clean PASS", div)))

  # ESS uses the smaller of bulk and tail.
  expect_match(
    gate_failure_bullets(gate(1.01, 900, 10, 0L, 0), th)[1], "min ESS 10"
  )

  # A cut manifest carries no top-level R-hat/ESS; say so rather than nothing.
  expect_match(
    gate_failure_bullets(list(converged = FALSE), th),
    "no passing status"
  )
})

test_that("gate_worst_parameters is silent without a fit", {
  expect_equal(gate_worst_parameters(NULL), character())
})
