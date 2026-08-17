test_that("flat index ordering matches numpy reshape(-1)", {
  fx <- jsonlite::fromJSON(fixture_path("flatten_parity.json"))
  K <- fx$K
  D <- fx$D
  N <- fx$N
  Y <- unflatten_kdn(fx$Y_flat, K, D, N)
  expect_equal(flatten_kdn(Y), fx$Y_flat)

  control <- unflatten_kdn(as.logical(fx$control_flat), K, D, N)
  expect_identical(which(!flatten_kdn(control)), fx$exposed_cells_1based)

  sub <- kdn_from_flat(seq_len(K * D * N), D, N)
  expect_equal(flat_idx(sub$k, sub$d, sub$n, D, N), seq_len(K * D * N))
})

test_that("chain quotas match Python _chain_quotas", {
  fx <- jsonlite::fromJSON(fixture_path("cut_mechanics.json"), simplifyVector = FALSE)
  for (case in fx$quotas) {
    expect_identical(
      chain_quotas(case$m, case$chains),
      as.integer(unlist(case$result)),
      info = sprintf("m=%d chains=%d", case$m, case$chains)
    )
  }
})

test_that("strided subsample indices match numpy linspace().round()", {
  fx <- jsonlite::fromJSON(fixture_path("cut_mechanics.json"), simplifyVector = FALSE)
  for (case in fx$strided) {
    expect_identical(
      strided_indices(case$n, case$quota),
      as.integer(unlist(case$result_1based)),
      info = sprintf("n=%d quota=%d", case$n, case$quota)
    )
  }
  expect_error(strided_indices(5, 6), "only 5")
})

test_that("logsumexp is stable and correct", {
  x <- c(1000, 1000.5, 999)
  expect_equal(logsumexp(x), 1000.5 + log(sum(exp(x - 1000.5))))
  expect_equal(logsumexp(c(-Inf, -Inf)), -Inf)
})

test_that("unit_slug matches Python _slug", {
  expect_equal(unit_slug("New York"), "new_york")
  expect_equal(unit_slug("Texas"), "texas")
})
