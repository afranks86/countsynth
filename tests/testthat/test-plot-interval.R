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
  class(draws) <- c("countsynth_draws", class(draws))
  draws
}

test_that("separate_units puts named units in their own band", {
  draws <- make_interval_draws()
  agg <- add_aggregate_units(
    draws,
    list(countsynth_aggregate_unit(unit = "All treated", include_treated_units = TRUE))
  )
  p <- countsynth_interval_plot(agg, separate_units = "All treated")

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

test_that("aggregate units are split off by default, and can be opted out of", {
  agg <- add_aggregate_units(
    make_interval_draws(),
    list(countsynth_aggregate_unit(unit = "All treated", include_treated_units = TRUE))
  )
  # No separate_units argument at all: the frame names its own aggregates.
  auto <- countsynth_interval_plot(agg)
  expect_s3_class(auto$facet, "FacetGrid")
  expect_equal(
    as.character(auto$data$.band[auto$data$unit == "All treated"]), "separate"
  )
  # character() ranks everything in one band.
  flat <- countsynth_interval_plot(agg, separate_units = character())
  expect_s3_class(flat$facet, "FacetNull")
  expect_true("All treated" %in% flat$data$unit)
})

test_that("separate_units is inert when it names no present unit", {
  draws <- make_interval_draws()
  expect_s3_class(countsynth_interval_plot(draws, separate_units = "Nope")$facet, "FacetNull")
  expect_s3_class(countsynth_interval_plot(draws)$facet, "FacetNull")
  expect_false(".band" %in% names(countsynth_interval_plot(draws)$data))
})

test_that("interval_aggregates gates whether the aggregate reaches the plot", {
  draws <- make_interval_draws()
  spec <- list(
    countsynth_aggregate_unit(unit = "All treated", include_treated_units = TRUE)
  )
  for (on in c(FALSE, TRUE)) {
    dir <- withr::local_tempdir()
    countsynth_report(
      draws, dir,
      target_unit = "C", figures = "interval",
      aggregate_units = spec, interval_aggregates = on, print_tables = FALSE
    )
    expect_true(file.exists(file.path(dir, "figs", "interval.png")))
  }
  # The gate is on the plot data, which the PNG hides -- check it directly.
  agg <- add_aggregate_units(draws, spec)
  expect_false("All treated" %in% countsynth_interval_plot(draws)$data$unit)
  expect_true("All treated" %in% countsynth_interval_plot(agg)$data$unit)
})

test_that("a configured aggregate unit becomes the default headline unit", {
  # Reproduces the yearly-fertility setup: aggregate_units configured,
  # target_unit left unset. C has the largest single-unit effect and D the
  # longest name; neither should win over the pooled unit.
  draws <- make_interval_draws()
  spec <- list(
    countsynth_aggregate_unit(unit = "All treated", include_treated_units = TRUE)
  )
  dir <- withr::local_tempdir()
  res <- countsynth_report(
    draws, dir,
    figures = character(), aggregate_units = spec, print_tables = FALSE
  )
  expect_equal(res$target_unit, "All treated")
  expect_true("All treated" %in% res$per_unit$unit)
  # ... and it reaches the by-unit table, not just post_treatment_summary.
  by_unit <- utils::read.csv(file.path(dir, "figs", "summary_table_by_unit.csv"))
  expect_true("All treated" %in% by_unit$Unit)

  # Without any aggregate the old rule still applies: most treated periods.
  plain <- countsynth_report(
    draws, withr::local_tempdir(),
    figures = character(), print_tables = FALSE
  )
  expect_true(plain$target_unit %in% c("B", "C", "D"))
})

test_that("aggregate_unit_names round-trips and survives row subsetting", {
  draws <- make_interval_draws()
  expect_equal(aggregate_unit_names(draws), character())
  agg <- add_aggregate_units(
    draws,
    list(countsynth_aggregate_unit(unit = "All treated", include_treated_units = TRUE))
  )
  expect_equal(aggregate_unit_names(agg), "All treated")
  # The report subsets to post-treatment rows before plotting, so the marker
  # has to survive that.
  expect_equal(aggregate_unit_names(agg[agg$treatment == 1, ]), "All treated")
  expect_equal(aggregate_unit_names(dplyr::filter(agg, .data$unit != "A")), "All treated")
})

test_that("output opts and YAML carry interval_aggregates", {
  expect_true(countsynth_output_opts()$interval_aggregates)
  expect_false(countsynth_output_opts(interval_aggregates = FALSE)$interval_aggregates)
  expect_true(parse_yaml_output(list())$interval_aggregates)
  expect_false(parse_yaml_output(list(interval_aggregates = FALSE))$interval_aggregates)
})
