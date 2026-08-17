# Statistical parity vs the Python implementation on the bundled fertility
# smoke configuration. Slow (~20 min: full Stan fit); run explicitly with
# BPNMF_RUN_STATISTICAL_PARITY=true. Compares post-treatment estimands
# against frozen Python golden summaries (inst/python-golden).

test_that("R posterior estimands match the Python golden summaries", {
  skip_on_cran()
  skip_if_no_cmdstan()
  skip_if_not(
    identical(Sys.getenv("BPNMF_RUN_STATISTICAL_PARITY"), "true"),
    "Set BPNMF_RUN_STATISTICAL_PARITY=true to run the slow parity fit"
  )

  golden_dir <- system.file("python-golden", package = "bpnmf")
  golden <- utils::read.csv(
    file.path(golden_dir, "post_treatment_summary.csv")
  )

  cfg <- bpnmf_example_config(
    output_dir = file.path(withr::local_tempdir(), "parity"),
    output = bpnmf_output_opts(figures = FALSE, print_tables = FALSE)
  )
  dat <- bpnmf_data(cfg)
  fit <- bpnmf_fit(dat, config = cfg)
  pt <- bpnmf_post_treatment_summary(bpnmf_draws(fit))

  m <- merge(pt, golden, by = c("unit", "group"), suffixes = c("_r", "_py"))
  expect_equal(nrow(m), nrow(golden))

  # Same estimand from two samplers: means must agree well within the
  # posterior interval half-width (observed ~0.08 at fixture creation).
  ci_halfwidth <- (m$excess_pct_upper_95_py - m$excess_pct_lower_95_py) / 2
  rel <- abs(m$excess_pct_mean_r - m$excess_pct_mean_py) / ci_halfwidth
  expect_lt(max(rel), 0.5)

  # Expected (counterfactual) totals within 1%.
  expect_lt(
    max(abs(m$expected_mean_r / m$expected_mean_py - 1)),
    0.01
  )
})
