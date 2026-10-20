# Group-axis pooling: the implied total behind a by-race / by-age split.

# 2 units (A control, B treated in periods 2-3), 2 groups of deliberately
# different size, 3 times, 4 draws. "big" is 10x "small" in both count and
# exposure, and carries a +10% effect against the small group's +50%, so an
# exposure-weighted total (+13.6%) is distinguishable from the unweighted
# mean of the two percentages (+30%).
make_group_draws <- function() {
  units <- c("A", "B")
  times <- seq(as.Date("2021-01-01"), by = "month", length.out = 3)
  grid <- expand.grid(
    unit = units, time = times, group = c("big", "small"),
    stringsAsFactors = FALSE
  )
  grid$treatment <- as.integer(grid$unit == "B" & grid$time >= times[2])
  big <- grid$group == "big"
  grid$outcome <- ifelse(big, 1000, 100)
  grid$denominator <- ifelse(big, 10000, 1000)
  dplyr::bind_rows(lapply(1:4, function(d) {
    g <- grid
    g$.draw <- d
    g$.chain <- 1L
    g$.iteration <- d
    is_big <- g$group == "big"
    g$mu <- ifelse(is_big, log(1000 + 10 * d), log(100 + d))
    g$mu_treated <- g$mu + ifelse(
      g$treatment == 1, ifelse(is_big, log(1.1), log(1.5)), 0
    )
    g$ypred <- exp(g$mu)
    structure(g, class = c("countsynth_draws", class(g)))
  }))
}

pooled_spec <- function(...) {
  list(countsynth_aggregate_group("All", include_all_groups = TRUE, ...))
}

test_that("aggregate groups pool with sums and logsumexp", {
  agg <- add_aggregate_groups(make_group_draws(), pooled_spec())
  all_rows <- agg[agg$group == "All", ]
  expect_equal(nrow(all_rows), 24) # 4 draws x 3 times x 2 units

  row <- all_rows[all_rows$.draw == 1 & all_rows$unit == "A" &
    all_rows$time == as.Date("2021-01-01"), ]
  expect_equal(row$outcome, 1100) # 1000 + 100
  expect_equal(row$denominator, 11000) # 10000 + 1000
  expect_equal(row$mu, logsumexp(c(log(1010), log(101))))
  expect_equal(row$treatment, 0)

  # treatment is max'd over the pooled axis, same as on the unit axis
  treated <- all_rows[all_rows$.draw == 1 & all_rows$unit == "B" &
    all_rows$time == as.Date("2021-02-01"), ]
  expect_equal(treated$treatment, 1)

  # The real groups are untouched and still there.
  expect_setequal(unique(agg$group), c("big", "small", "All"))
  expect_equal(sum(agg$group == "big"), 24)
})

test_that("the pooled effect is exposure-weighted, not an average of groups", {
  agg <- add_aggregate_groups(make_group_draws(), pooled_spec())
  treated <- agg[agg$group == "All" & agg$unit == "B" &
    agg$time == as.Date("2021-02-01"), ]
  ratio <- exp(treated$mu_treated - treated$mu)
  # (1.1 * 1000 + 1.5 * 100) / 1100, constant across draws since both groups
  # scale identically with the draw index.
  expect_equal(ratio, rep(1250 / 1100, 4))
  # Emphatically not the mean of the two group percentages.
  expect_false(isTRUE(all.equal(ratio[[1]], mean(c(1.1, 1.5)))))
})

test_that("pooling one group is refused rather than silently renaming it", {
  single <- make_draws_frame() # one group, "total"
  expect_warning(
    out <- add_aggregate_groups(single, pooled_spec()),
    "at least two groups"
  )
  expect_equal(unique(out$group), "total")
  expect_equal(aggregate_group_names(out), character())

  # Same refusal when an exclusion whittles a real multi-group fit down to one.
  expect_warning(
    add_aggregate_groups(
      make_group_draws(),
      list(countsynth_aggregate_group(
        "All", include_all_groups = TRUE, exclude_groups = "small"
      ))
    ),
    "at least two groups"
  )
})

test_that("aggregate groups compose with aggregate units", {
  draws <- make_group_draws()
  with_units <- add_aggregate_units(
    draws, list(countsynth_aggregate_unit("Both", include_all_units = TRUE))
  )
  both <- add_aggregate_groups(with_units, pooled_spec())

  # Both markers survive the second pass -- bind_rows drops attributes.
  expect_equal(aggregate_unit_names(both), "Both")
  expect_equal(aggregate_group_names(both), "All")

  # The cross cell exists and totals all four source cells.
  cell <- both[both$unit == "Both" & both$group == "All" &
    both$.draw == 1 & both$time == as.Date("2021-01-01"), ]
  expect_equal(nrow(cell), 1)
  expect_equal(cell$outcome, 2 * 1100)
  expect_equal(cell$mu, logsumexp(log(c(1010, 1010, 101, 101))))

  # Pooling in the other order agrees (logsumexp is associative).
  other <- add_aggregate_units(
    add_aggregate_groups(draws, pooled_spec()),
    list(countsynth_aggregate_unit("Both", include_all_units = TRUE))
  )
  swapped <- other[other$unit == "Both" & other$group == "All" &
    other$.draw == 1 & other$time == as.Date("2021-01-01"), ]
  expect_equal(swapped$mu, cell$mu)
  expect_equal(swapped$outcome, cell$outcome)
})

test_that("a suppressed cell is imputed once and stays flagged through both axes", {
  draws <- make_group_draws()
  hole <- draws$unit == "A" & draws$group == "small" &
    draws$time == as.Date("2021-01-01")
  draws$outcome[hole] <- NA_real_

  units_first <- add_aggregate_units(
    draws, list(countsynth_aggregate_unit("Both", include_all_units = TRUE))
  )
  # The unit pass fills A's missing cell from exp(mu_treated) and flags it.
  pooled_unit <- units_first[units_first$unit == "Both" &
    units_first$group == "small" & units_first$.draw == 1 &
    units_first$time == as.Date("2021-01-01"), ]
  expect_true(pooled_unit$outcome_imputed)
  expect_false(is.na(pooled_unit$outcome))

  # The group pass must not reset the flag just because the hole is now full.
  both <- add_aggregate_groups(units_first, pooled_spec())
  cross <- both[both$unit == "Both" & both$group == "All" &
    both$.draw == 1 & both$time == as.Date("2021-01-01"), ]
  expect_true(cross$outcome_imputed)
  expect_false(is.na(cross$outcome))

  # An untouched cell stays unflagged.
  clean <- both[both$unit == "Both" & both$group == "All" &
    both$.draw == 1 & both$time == as.Date("2021-03-01"), ]
  expect_false(clean$outcome_imputed)
})

test_that("group selection honors include_groups, exclude_groups and strict", {
  draws <- make_group_draws()
  draws <- dplyr::bind_rows(draws, transform(draws, group = "third"))
  class(draws) <- c("countsynth_draws", class(draws))

  picked <- add_aggregate_groups(
    draws,
    list(countsynth_aggregate_group("Two", include_groups = c("big", "small")))
  )
  cell <- picked[picked$group == "Two" & picked$unit == "A" &
    picked$.draw == 1 & picked$time == as.Date("2021-01-01"), ]
  expect_equal(cell$outcome, 1100) # "third" left out

  expect_warning(
    add_aggregate_groups(
      draws,
      list(countsynth_aggregate_group(
        "Two",
        include_groups = c("big", "small", "nope")
      ))
    ),
    "missing group"
  )
  expect_error(
    add_aggregate_groups(
      draws,
      list(countsynth_aggregate_group(
        "Two",
        include_groups = c("big", "small", "nope"), strict = TRUE
      ))
    ),
    "missing group"
  )
})

test_that("aggregate group collision requires overwrite", {
  draws <- make_group_draws()
  expect_error(
    add_aggregate_groups(
      draws, list(countsynth_aggregate_group("big", include_all_groups = TRUE))
    ),
    "collides"
  )
  over <- add_aggregate_groups(
    draws,
    list(countsynth_aggregate_group(
      "big",
      include_all_groups = TRUE, overwrite = TRUE
    ))
  )
  expect_equal(sum(over$group == "big"), 24) # 4 draws x 3 times x 2 units
  cell <- over[over$group == "big" & over$unit == "A" & over$.draw == 1 &
    over$time == as.Date("2021-01-01"), ]
  expect_equal(cell$outcome, 1100) # the pool, not the original 1000
})

test_that("aggregate_group_names round-trips and survives row subsetting", {
  draws <- make_group_draws()
  expect_equal(aggregate_group_names(draws), character())
  agg <- add_aggregate_groups(draws, pooled_spec())
  expect_equal(aggregate_group_names(agg), "All")
  expect_equal(aggregate_group_names(agg[agg$treatment == 1, ]), "All")
  expect_equal(
    aggregate_group_names(dplyr::filter(agg, .data$unit != "A")), "All"
  )
})

test_that("the spec needs exactly one include selector", {
  expect_error(
    countsynth_aggregate_group("All"), "exactly one include selector"
  )
  expect_error(
    countsynth_aggregate_group(
      "All",
      include_all_groups = TRUE, include_groups = "big"
    ),
    "exactly one include selector"
  )
  expect_s3_class(
    countsynth_aggregate_group("All", include_all_groups = TRUE),
    "countsynth_aggregate_group"
  )
})

test_that("output opts and YAML carry aggregate_groups", {
  expect_null(countsynth_output_opts()$aggregate_groups)
  opts <- countsynth_output_opts(aggregate_groups = pooled_spec())
  expect_s3_class(opts$aggregate_groups[[1]], "countsynth_aggregate_group")

  parsed <- parse_yaml_output(list(
    aggregate_groups = list(
      list(group = "All races", include_all_groups = TRUE),
      list(group = "Non-Hispanic", include_groups = list("nhwhite", "nhblack"))
    )
  ))
  expect_length(parsed$aggregate_groups, 2)
  expect_equal(parsed$aggregate_groups[[1]]$group, "All races")
  expect_true(parsed$aggregate_groups[[1]]$include_all_groups)
  expect_equal(
    parsed$aggregate_groups[[2]]$include_groups, c("nhwhite", "nhblack")
  )
  expect_error(
    parse_yaml_output(list(aggregate_groups = list(list(group = "All", bogus = 1)))),
    "bogus"
  )
})

test_that("the report carries the pooled group into tables and figures", {
  draws <- make_group_draws()
  dir <- withr::local_tempdir()
  res <- countsynth_report(
    draws, dir,
    target_unit = "B", figures = c("unit_fit", "raw_rate", "interval"),
    aggregate_groups = pooled_spec(), print_tables = FALSE, html_tables = FALSE
  )
  # Group labels carry a significance star, so match the stem.
  expect_true(any(startsWith(res$summary$Group, "All")))
  expect_true("All" %in% res$per_unit$group)

  # Its own per-group figure directory, raw_rate included -- that plot is
  # built from a frame without the aggregate units, so it needs the pooled
  # group added separately or it would be handed no rows.
  expect_true(file.exists(file.path(dir, "figures", "All", "fit_B.png")))
  expect_true(file.exists(file.path(dir, "figures", "All", "raw_rate.png")))

  # interval.png is gated by interval_aggregates; check the plot data, which
  # the PNG hides.
  agg <- add_aggregate_groups(draws, pooled_spec())
  expect_false("All" %in% countsynth_interval_plot(draws)$data$group)
  expect_true("All" %in% countsynth_interval_plot(agg)$data$group)

  # group_comparison stays on the fitted groups: a total faceted beside its
  # own components reads as a peer comparison when it is not.
  gc <- countsynth_group_comparison_plot(draws)
  expect_setequal(unique(gc$data$group), c("big", "small"))
})
