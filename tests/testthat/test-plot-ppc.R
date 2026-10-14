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
  class(draws) <- c("countsynth_draws", class(draws))
  draws
}

test_that("a cell with no observed control data is dropped, not scored p = 1", {
  draws <- make_ppc_draws()
  expect_warning(res <- countsynth_ppc_abs(draws), "no observed control-period data")
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
    countsynth_ppc_abs(draws),
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
  rmse <- suppressWarnings(countsynth_ppc_rmse(draws))
  expect_equal(nrow(rmse$pvals), 3)
  expect_true(all(is.finite(rmse$plot$data$diff_in_diff)))

  acf <- suppressWarnings(countsynth_ppc_acf(draws, lag = 2))
  expect_false(any(
    acf$pvals$unit == "Kentucky" & acf$pvals$group == "otherraceeth"
  ))
})

test_that("fully observed panels are untouched and warn about nothing", {
  draws <- make_ppc_draws(suppress = FALSE)
  expect_no_warning(res <- countsynth_ppc_abs(draws))
  expect_equal(nrow(res$pvals), 4)
})

# Degenerate cells in the unit-correlation check: a unit whose observed series
# never moves has residual -exp(mu), a smooth curve that correlates
# near-perfectly with every other such unit.
make_flat_unit_draws <- function() {
  units <- c(sprintf("busy%d", 1:6), sprintf("flat%d", 1:5), "sparse1")
  times <- seq(as.Date("2010-01-01"), by = "year", length.out = 12)
  grid <- expand.grid(
    unit = units, time = times, group = "total", stringsAsFactors = FALSE
  )
  grid$treatment <- 0L
  grid$denominator <- 50000
  withr::with_seed(5, {
    grid$outcome <- stats::rpois(nrow(grid), 100)
    grid$outcome[startsWith(grid$unit, "flat")] <- 0
    # zeros in 9 of 12 periods, but not flat
    sparse <- grid$unit == "sparse1"
    grid$outcome[sparse] <- rep(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3), each = 1)
    draws <- dplyr::bind_rows(lapply(1:40, function(d) {
      g <- grid
      g$.draw <- d
      g$.chain <- 1L
      g$.iteration <- d
      g$mu <- ifelse(startsWith(g$unit, "flat") | g$unit == "sparse1",
                     log(1) + 0.1 * as.integer(format(g$time, "%Y")) %% 3,
                     log(100)) + stats::rnorm(nrow(g), 0, 0.05)
      g$mu_treated <- g$mu
      g$ypred <- stats::rpois(nrow(g), exp(g$mu))
      g
    }))
  })
  class(draws) <- c("countsynth_draws", class(draws))
  draws
}

test_that("the unit-correlation check drops flat and mostly-zero cells", {
  draws <- make_flat_unit_draws()
  expect_warning(res <- countsynth_ppc_unit_corr(draws), "no variation")
  expect_equal(nrow(res$pvals), 1)
  # All five all-zero units and the mostly-zero one are gone, so the observed
  # spectral norm is no longer dominated by -exp(mu) curves.
  expect_gt(res$pvals$pval, 0.05)
})

test_that("flat residual columns are what inflate the observed spectral norm", {
  trend <- -exp(seq(0, 1, length.out = 12))
  noise <- withr::with_seed(7, matrix(stats::rnorm(12 * 6), 12, 6))
  # Six noisy units alone, versus the same six plus five units whose residual
  # is the same deterministic curve -- the all-zero-outcome case.
  expect_lt(
    countsynth:::spectral_norm(noise),
    countsynth:::spectral_norm(cbind(noise, matrix(trend, 12, 5)))
  )
})

test_that("a mostly-zero cell with some variation is kept at a loose threshold", {
  draws <- make_flat_unit_draws()
  draws <- draws[!startsWith(draws$unit, "flat"), ]
  class(draws) <- c("countsynth_draws", class(draws))
  expect_warning(
    tight <- countsynth_ppc_unit_corr(draws, max_zero_frac = 0.5),
    "mostly zeros"
  )
  expect_silent(loose <- countsynth_ppc_unit_corr(draws, max_zero_frac = 1))
  expect_true(is.data.frame(tight$pvals) && is.data.frame(loose$pvals))
})
