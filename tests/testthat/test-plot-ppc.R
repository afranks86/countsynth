# PPC handling of cells with no observed control-period data.

make_ppc_draws <- function(suppress = TRUE) {
  units <- c("Kentucky", "Texas")
  groups <- c("total", "otherraceeth")
  times <- seq(as.Date("2016-01-01"), by = "year", length.out = 8)
  grid <- expand.grid(
    unit = units, time = times, group = groups, stringsAsFactors = FALSE
  )
  grid$treatment <- as.integer(grid$time >= as.Date("2022-01-01"))
  grid$denominator <- 50000
  grid$outcome <- 100
  if (suppress) {
    grid$outcome[grid$unit == "Kentucky" & grid$group == "otherraceeth"] <- NA
  }
  draws <- withr::with_seed(3, dplyr::bind_rows(lapply(1:20, function(d) {
    g <- grid
    g$.draw <- d
    g$.chain <- 1L
    g$.iteration <- d
    g$mu <- log(100) + stats::rnorm(nrow(g), 0, 0.05)
    g$mu_treated <- g$mu
    g$ypred <- stats::rpois(nrow(g), 100)
    g
  })))
  class(draws) <- c("bpnmf_draws", class(draws))
  draws
}

test_that("a cell with no observed control data is dropped, not scored p = 1", {
  draws <- make_ppc_draws()
  expect_warning(res <- bpnmf_ppc_abs(draws), "no observed control-period data")
  # The old path let max(all-NA, na.rm = TRUE) return -Inf, which beat every
  # predicted statistic and reported a perfect p = 1 for an empty cell.
  expect_false(any(
    res$pvals$unit == "Kentucky" & res$pvals$group == "otherraceeth"
  ))
  expect_equal(nrow(res$pvals), 3)
  expect_true(all(is.finite(res$plot$data$diff_in_diff)))
})

test_that("the warning names the dropped cells once, not once per draw", {
  draws <- make_ppc_draws()
  warns <- character()
  withCallingHandlers(
    bpnmf_ppc_abs(draws),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warns, 1)
  expect_match(warns, "Kentucky/otherraceeth")
  expect_false(any(grepl("no non-missing arguments", warns)))
})

test_that("rmse and acf checks drop the same cell", {
  draws <- make_ppc_draws()
  rmse <- suppressWarnings(bpnmf_ppc_rmse(draws))
  expect_equal(nrow(rmse$pvals), 3)
  expect_true(all(is.finite(rmse$plot$data$diff_in_diff)))

  acf <- suppressWarnings(bpnmf_ppc_acf(draws, lag = 2))
  expect_false(any(
    acf$pvals$unit == "Kentucky" & acf$pvals$group == "otherraceeth"
  ))
})

test_that("fully observed panels are untouched and warn about nothing", {
  draws <- make_ppc_draws(suppress = FALSE)
  expect_no_warning(res <- bpnmf_ppc_abs(draws))
  expect_equal(nrow(res$pvals), 4)
})
