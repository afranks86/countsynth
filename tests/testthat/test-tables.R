# Table math on a hand-built draws frame with known answers
# (make_draws_frame(), in setup.R).

test_that("auto_detect_target picks the unit with most treated periods", {
  expect_equal(auto_detect_target(make_draws_frame()), "B")
})

test_that("compute_quantiles gives exact quantiles of ypred", {
  q <- compute_quantiles(make_draws_frame())
  expect_equal(nrow(q), 6) # 2 units x 3 times x 1 group
  row <- q[q$unit == "A" & q$time == as.Date("2021-01-01"), ]
  expect_equal(row$ypred_mean, mean(101:104))
  expect_equal(row$ypred_median, median(101:104))
  expect_equal(row$ypred_lower, quantile(101:104, 0.025, names = FALSE))
})

test_that("post-treatment summary matches hand-computed estimands", {
  pt <- countsynth_post_treatment_summary(make_draws_frame())
  expect_equal(nrow(pt), 1)
  expect_equal(pt$unit, "B")
  expect_equal(pt$n_periods, 2)
  expect_equal(pt$observed, 120 + 130)
  # expected = sum over 2 periods of exp(mu) per draw, then mean over draws
  per_draw_expected <- 2 * (101:104)
  expect_equal(pt$expected_mean, mean(per_draw_expected))
  # excess = treated - expected = 0.2 * expected per draw
  expect_equal(pt$excess_mean, mean(0.2 * per_draw_expected))
  expect_equal(pt$excess_pct_mean, 20, tolerance = 1e-10)
})

test_that("summary table computes person-year-weighted rates", {
  st <- countsynth_summary_table(make_draws_frame(), "B", rate_normalizer = 1000)
  expect_equal(nrow(st), 1)
  # 2 post periods x denominator 1000 (years = 1 without period bounds)
  expect_equal(st$`Person-Years`, 2000L)
  expect_equal(st$Observed, 250L)
  # pct change is exactly 20% in every draw
  expect_match(st$`Pct Change CI`, "^20.0% \\(20.0%, 20.0%\\)")
  # every draw has treated > untreated -> p = 0 -> significant
  expect_match(st$Group, "\\*$")
})

test_that("expected_vs_observed detail computes gaps", {
  detail <- countsynth_expected_vs_observed(make_draws_frame(), "B")
  expect_equal(nrow(detail), 6)
  row <- detail[detail$unit == "A" & detail$time == as.Date("2021-01-01"), ]
  expect_equal(row$gap, 100 - mean(101:104))
  expect_true(all(detail$treated_unit == (detail$unit == "B")))
})

test_that("aggregate units pool with sums and logsumexp", {
  draws <- make_draws_frame()
  agg <- add_aggregate_units(
    draws,
    list(countsynth_aggregate_unit("Both", include_all_units = TRUE))
  )
  both <- agg[agg$unit == "Both", ]
  expect_equal(nrow(both), 12) # 4 draws x 3 times x 1 group
  row <- both[both$.draw == 1 & both$time == as.Date("2021-01-01"), ]
  expect_equal(row$outcome, 210) # A 100 + B 110
  expect_equal(row$mu, logsumexp(c(log(101), log(101))))
  expect_equal(row$treatment, 0)
  row2 <- both[both$.draw == 1 & both$time == as.Date("2021-02-01"), ]
  expect_equal(row2$treatment, 1) # max over units
})

test_that("aggregate unit collision requires overwrite", {
  draws <- make_draws_frame()
  expect_error(
    add_aggregate_units(
      draws, list(countsynth_aggregate_unit("A", include_all_units = TRUE))
    ),
    "collides"
  )
  agg <- add_aggregate_units(
    draws,
    list(countsynth_aggregate_unit("A", include_all_units = TRUE, overwrite = TRUE))
  )
  # A replaced by the aggregate: 12 rows for A (4 draws x 3 times)
  expect_equal(sum(agg$unit == "A"), 12)
  expect_equal(
    agg$outcome[agg$unit == "A" & agg$.draw == 1 & agg$time == as.Date("2021-01-01")],
    210
  )
})

test_that("strict aggregate spec errors on missing units", {
  draws <- make_draws_frame()
  expect_error(
    add_aggregate_units(
      draws,
      list(countsynth_aggregate_unit(
        "X",
        include_units = c("A", "ZZ"), strict = TRUE
      ))
    ),
    "missing unit"
  )
})

test_that("PPC statistics are computed on control cells", {
  draws <- make_draws_frame()
  r <- countsynth_ppc_rmse(draws)
  # Unit B is treated -> PPC restricted to treated_units per Python logic;
  # B's control period rows only.
  expect_true(all(r$pvals$pval >= 0 & r$pvals$pval <= 1))
  expect_s3_class(r$plot, "ggplot")
})

test_that("interval plot effects are exact for the synthetic frame", {
  draws <- make_draws_frame()
  df <- draws[!is.na(draws$treatment) & draws$treatment == 1, ]
  df$years <- years_per_row(df)
  eff <- compute_draw_effects(
    df,
    estimand = "ratio", method = "mu", rate_normalizer = 1000,
    agg_cols = c(".draw", "unit")
  )
  expect_equal(unique(round(eff$causal_effect, 10)), 20)
  p <- countsynth_interval_plot(draws, estimand = "ratio", method = "mu")
  expect_s3_class(p, "ggplot")
})

test_that("dodged interval segments stay horizontal and track their points", {
  # geom_segment + position_dodge dodges y but not yend, slanting every
  # interval; the linerange layers must sit at the same y as the median point.
  draws <- make_draws_frame()
  g2 <- draws
  g2$group <- "other"
  g2$mu_treated <- g2$mu + ifelse(g2$treatment == 1, log(1.5), 0)
  draws <- dplyr::bind_rows(draws, g2)
  class(draws) <- c("countsynth_draws", class(draws))

  layers <- ggplot2::ggplot_build(
    countsynth_interval_plot(draws, estimand = "ratio", method = "mu")
  )$data
  ci_95 <- layers[[2]]
  ci_67 <- layers[[3]]
  points <- layers[[4]]

  expect_gt(length(unique(points$y)), 1) # groups really are dodged apart
  expect_equal(ci_95$y, points$y)
  expect_equal(ci_67$y, points$y)
  expect_false("yend" %in% names(ci_95))
})
