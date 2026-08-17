# Python golden summaries

Reference outputs from the Python `bayesian_panel_nmf` implementation on the
shipped fertility smoke-test configuration (`fertility_smoke_test.yaml`:
total births, rank 3, bimonthly, NB, 4 chains x 200+200, seed 8675309),
generated 2026-08-17 with numpyro 0.21 / jax 0.10.2.

Used by `tests/testthat/test-statistical-parity.R` to check that the R/Stan
implementation reproduces the same post-treatment estimands within Monte
Carlo error. Smoke settings do not converge (by design); the estimand
comparison is still informative because both posteriors target the same
model. Observed R-vs-Python agreement at creation time: max
|excess_pct_mean| difference 0.098pp, max difference / CI-halfwidth 0.079.
