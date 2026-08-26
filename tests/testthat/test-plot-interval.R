# Interval-plot banding: aggregate units are split off rather than ranked in.

make_interval_draws <- function() {
  # 1 control + 3 treated units, 1 group, 4 times, 20 draws, known effects.
  units <- c("A", "B", "C", "D")
  times <- seq(as.Date("2021-01-01"), by = "month", length.out = 4)
  grid <- expand.grid(unit = units, time = times, stringsAsFactors = FALSE)
  grid$group <- "total"
  grid$treatment <- as.integer(
    grid$unit %in% c("B", "C", "D") & grid$time >= times[3]
  )
  grid$outcome <- 100
  grid$denominator <- 1000
  bump <- c(A = 0, B = log(1.1), C = log(1.3), D = log(0.9))
  draws <- dplyr::bind_rows(lapply(1:20, function(d) {
    g <- grid
    g$.draw <- d
    g$.chain <- 1L
    g$.iteration <- d
    g$mu <- log(100)
    g$mu_treated <- g$mu + ifelse(g$treatment == 1, bump[g$unit], 0)
    g$ypred <- 100
    g
  }))
  class(draws) <- c("bpnmf_draws", class(draws))
  draws
}

test_that("separate_units puts named units in their own band", {
  draws <- make_interval_draws()
  agg <- add_aggregate_units(
    draws,
    list(bpnmf_aggregate_unit(unit = "All treated", include_treated_units = TRUE))
  )
  p <- bpnmf_interval_plot(agg, separate_units = "All treated")

  expect_s3_class(p$facet, "FacetGrid")
  bands <- stats::setNames(as.character(p$data$.band), as.character(p$data$unit))
  expect_equal(bands[["All treated"]], "separate")
  expect_equal(unname(bands[c("B", "C", "D")]), rep("units", 3))
  # Both bands still sorted by median effect (D < B < C).
  expect_equal(
    levels(p$data$unit)[levels(p$data$unit) != "All treated"],
    c("D", "B", "C")
  )
})

test_that("separate_units is inert when it names no present unit", {
  draws <- make_interval_draws()
  expect_s3_class(bpnmf_interval_plot(draws, separate_units = "Nope")$facet, "FacetNull")
  expect_s3_class(bpnmf_interval_plot(draws)$facet, "FacetNull")
  expect_false(".band" %in% names(bpnmf_interval_plot(draws)$data))
})

test_that("interval_aggregates gates whether the aggregate reaches the plot", {
  draws <- make_interval_draws()
  spec <- list(
    bpnmf_aggregate_unit(unit = "All treated", include_treated_units = TRUE)
  )
  for (on in c(FALSE, TRUE)) {
    dir <- withr::local_tempdir()
    bpnmf_report(
      draws, dir,
      target_unit = "C", figures = "interval",
      aggregate_units = spec, interval_aggregates = on, print_tables = FALSE
    )
    expect_true(file.exists(file.path(dir, "figs", "interval.png")))
  }
  # The gate is on the plot data, which the PNG hides -- check it directly.
  agg <- add_aggregate_units(draws, spec)
  expect_false("All treated" %in% bpnmf_interval_plot(draws)$data$unit)
  expect_true(
    "All treated" %in%
      bpnmf_interval_plot(agg, separate_units = "All treated")$data$unit
  )
})

test_that("output opts and YAML carry interval_aggregates", {
  expect_false(bpnmf_output_opts()$interval_aggregates)
  expect_true(bpnmf_output_opts(interval_aggregates = TRUE)$interval_aggregates)
  parsed <- parse_yaml_output(list(interval_aggregates = TRUE))
  expect_true(parsed$interval_aggregates)
})
