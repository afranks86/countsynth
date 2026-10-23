# Pooled units under staggered adoption: a member enters the pooled effect on
# its own start date, not the earliest member's.

# A control, B exposed from month 2, C only from month 3, 1 group, 3 times,
# 4 draws. Every exposed cell carries exactly +20%. Pooling B and C from B's
# start would put C's unexposed month 2 into the counterfactual and dilute
# that to (3 * 120 + 100) / 400 - 1 = +15%.
make_staggered_draws <- function() {
  times <- seq(as.Date("2021-01-01"), by = "month", length.out = 3)
  grid <- expand.grid(
    unit = c("A", "B", "C"), time = times, stringsAsFactors = FALSE
  )
  grid$group <- "total"
  start <- c(A = NA, B = 2, C = 3)
  grid$treatment <- as.integer(
    !is.na(start[grid$unit]) & match(grid$time, times) >= start[grid$unit]
  )
  grid$outcome <- ifelse(grid$treatment == 1, 120, 100)
  grid$denominator <- 1000
  draws <- dplyr::bind_rows(lapply(1:4, function(d) {
    g <- grid
    g$.draw <- d
    g$.chain <- 1L
    g$.iteration <- d
    g$mu <- log(100)
    g$mu_treated <- g$mu + ifelse(g$treatment == 1, log(1.2), 0)
    g$ypred <- 100
    g
  }))
  class(draws) <- c("countsynth_draws", class(draws))
  draws
}

pool_treated <- function(draws) {
  add_aggregate_units(
    draws,
    list(countsynth_aggregate_unit("All treated", include_treated_units = TRUE))
  )
}

test_that("the pooled unit's own rows stay the calendar-time total", {
  pooled <- pool_treated(make_staggered_draws())
  expect_equal(attr(pooled, "aggregate_unit_sources"), list(`All treated` = c("B", "C")))
  own <- pooled[pooled$unit == "All treated" & pooled$.draw == 1, ]
  own <- own[order(own$time), ]
  # Exposed from B's start, and totals both members whatever their exposure.
  expect_equal(own$treatment, c(0, 1, 1))
  expect_equal(own$outcome, c(200, 220, 240))
})

test_that("the summary table counts each member only once it is exposed", {
  pooled <- pool_treated(make_staggered_draws())
  st <- countsynth_summary_table(pooled, "All treated")
  # B's months 2-3 and C's month 3: three exposed cells, not four.
  expect_equal(st$`Denominator-Years`, 3000L)
  expect_equal(st$Observed, 360L)
  expect_equal(st$`Obs Rate`, 120)
  expect_equal(st$`Exp Rate`, 100)
  expect_match(st$`Pct Change CI`, "^20.0% \\(20.0%, 20.0%\\)")
})

test_that("post_treatment_summary and the interval plot agree", {
  pooled <- pool_treated(make_staggered_draws())
  pt <- countsynth_post_treatment_summary(pooled)
  row <- pt[pt$unit == "All treated", ]
  expect_equal(row$observed, 360)
  expect_equal(row$expected_mean, 300)
  expect_equal(row$excess_pct_mean, 20, tolerance = 1e-10)
  # Periods in which any member was exposed.
  expect_equal(row$n_periods, 2)
  # The real units are untouched by the rebuild.
  expect_equal(pt$observed[pt$unit == "B"], 240)
  expect_equal(pt$observed[pt$unit == "C"], 120)

  p <- countsynth_interval_plot(pooled, estimand = "ratio", method = "mu")
  pooled_row <- p$data[p$data$unit == "All treated", ]
  expect_equal(pooled_row$median, 20, tolerance = 1e-10)
})

test_that("a member's suppressed exposed cell is imputed and flagged", {
  draws <- make_staggered_draws()
  hole <- draws$unit == "C" & draws$time == as.Date("2021-03-01")
  draws$outcome[hole] <- NA_real_
  pooled <- pool_treated(draws)
  st <- countsynth_summary_table(pooled, "All treated")
  expect_true(st$Imputed)
  # Filled with exp(mu_treated) = 120, so the total is unchanged.
  expect_equal(st$Observed, 360L)
  pt <- countsynth_post_treatment_summary(pooled)
  expect_true(pt$observed_imputed[pt$unit == "All treated"])
})

test_that("pooled groups of the pooled unit respect adoption too", {
  draws <- make_staggered_draws()
  small <- draws
  small$group <- "small"
  both_groups <- dplyr::bind_rows(draws, small)
  class(both_groups) <- class(draws)
  reporting <- add_aggregate_groups(
    pool_treated(both_groups),
    list(countsynth_aggregate_group("All", include_all_groups = TRUE))
  )
  expect_equal(names(attr(reporting, "aggregate_unit_sources")), "All treated")
  st <- countsynth_summary_table(reporting, "All treated")
  cross <- st[st$Group == "All*", ]
  expect_equal(cross$Observed, 720L)
  expect_match(cross$`Pct Change CI`, "^20.0%")
})

test_that("missing members fall back to the pooled rows, with a warning", {
  pooled <- pool_treated(make_staggered_draws())
  alone <- pooled[pooled$unit == "All treated", ]
  expect_warning(
    st <- countsynth_summary_table(alone, "All treated"),
    "missing member"
  )
  expect_match(st$`Pct Change CI`, "^15.0%")
})
