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
