# The pooled-row reductions in aggregate_over() are hand-vectorized (see the
# comment above group_max()), so the semantics a grouped summarise gave for
# free -- max()'s NA rule, logsumexp()'s handling of a non-finite maximum --
# are now this package's to get right, and are pinned here.

test_that("group_sum / group_max / group_logsumexp reduce per group", {
  gid <- c(1L, 1L, 2L, 2L, 3L, 1L)
  x <- c(1, 2, 10, 20, 5, 3)
  expect_equal(group_sum(x, gid, 3L), c(6, 30, 5))
  expect_equal(group_max(x, gid, 3L), c(3, 20, 5))
  expect_equal(
    group_logsumexp(x, gid, 3L),
    c(log(sum(exp(c(1, 2, 3)))), log(sum(exp(c(10, 20)))), 5)
  )
})

test_that("group_max carries max()'s NA rule, group_logsumexp drops NAs", {
  gid <- c(1L, 1L, 2L, 2L, 3L, 3L)
  x <- c(1, NA, 7, 2, NA, NA)
  # max(): a group holding any NA is NA. A group holding only NAs is too.
  expect_equal(group_max(x, gid, 3L), c(NA, 7, NA))
  # logsumexp over the non-missing values; an all-NA group has nothing to
  # pool and stays NA rather than collapsing to -Inf.
  expect_equal(group_logsumexp(x, gid, 3L), c(1, log(exp(7) + exp(2)), NA))
})

test_that("group_logsumexp passes a non-finite maximum straight through", {
  gid <- c(1L, 1L, 2L, 2L, 3L, 3L)
  # All -Inf pools to -Inf (a cell contributing no count at all); any +Inf
  # dominates the sum. Both match logsumexp()'s own early return.
  x <- c(-Inf, -Inf, Inf, 4, -Inf, 2)
  expect_equal(group_logsumexp(x, gid, 3L), c(-Inf, Inf, 2))
})

test_that("pooling matches a grouped summarise, NA imputation included", {
  draws <- make_draws_frame()
  # Suppress one unit's count in one period, so the pooled cell has to be
  # filled from mean(exp(mu_treated)) rather than summed to NA.
  draws$outcome[draws$unit == "A" & draws$time == unique(draws$time)[2]] <- NA

  agg <- aggregate_over(draws, "unit", "both", c("A", "B"))
  reference <- draws |>
    dplyr::group_by(.data$.draw, .data$time, .data$group) |>
    dplyr::summarise(
      ypred = sum(.data$ypred),
      denominator = sum(.data$denominator),
      treatment = max(.data$treatment),
      mu = logsumexp(.data$mu),
      mu_treated = logsumexp(.data$mu_treated),
      .groups = "drop"
    )
  joined <- dplyr::left_join(
    agg, reference,
    by = c(".draw", "time", "group"), suffix = c("", "_ref")
  )
  for (col in c("ypred", "denominator", "treatment", "mu", "mu_treated")) {
    expect_equal(joined[[col]], joined[[paste0(col, "_ref")]], info = col)
  }
  # The suppressed cell is filled, flagged, and filled with the same number
  # for every draw -- compute_quantiles() groups by outcome and would split
  # the draws apart if the fill varied across them.
  filled <- agg[agg$time == unique(draws$time)[2], ]
  expect_true(all(filled$outcome_imputed))
  expect_equal(length(unique(filled$outcome)), 1L)
  expect_false(anyNA(agg$outcome))
})

test_that("observed_cell_frame collapses draws without moving the figures", {
  draws <- make_draws_frame()
  draws$outcome[draws$unit == "A" & draws$time == unique(draws$time)[3]] <- NA
  cells <- observed_cell_frame(draws)
  expect_equal(nrow(cells), nrow(dplyr::distinct(draws, unit, time, group)))
  # mu/mu_treated become the log of the cell's posterior mean rate, not one
  # draw's value: aggregate_over() fills suppressed cells from
  # mean(exp(mu_treated)), and that has to survive the collapse.
  one_cell <- draws[draws$unit == "A" & draws$time == unique(draws$time)[1], ]
  expect_equal(
    cells$mu_treated[cells$unit == "A" &
      cells$time == unique(draws$time)[1]],
    log(mean(exp(one_cell$mu_treated)))
  )

  # The payoff: pooling a reporting group over the cells gives the figures
  # the same numbers as pooling over the whole draws frame.
  spec <- list(countsynth_aggregate_group("all", include_all_groups = TRUE))
  draws$group <- rep(c("g1", "g2"), length.out = nrow(draws))
  cells <- observed_cell_frame(draws)
  keys <- c("unit", "time", "group", "outcome", "denominator", "treatment")
  from_draws <- observed_rows(add_aggregate_groups(draws, spec))[keys]
  from_cells <- observed_rows(add_aggregate_groups(cells, spec))[keys]
  ord <- function(d) as.data.frame(d[do.call(order, d[1:3]), ])
  expect_equal(ord(from_cells), ord(from_draws), ignore_attr = TRUE)
})

test_that("observed_cell_frame keeps the attributes reporting reads", {
  draws <- add_aggregate_units(
    make_draws_frame(),
    list(countsynth_aggregate_unit("both", include_all_units = TRUE))
  )
  cells <- observed_cell_frame(draws)
  # Load-bearing, not cosmetic: the collapsed frame is handed to the raw-rate
  # and group-comparison figures, which read has_denominator to decide
  # whether they are plotting a rate or a count at all.
  expect_s3_class(cells, "countsynth_draws")
  expect_identical(
    attr(cells, "has_denominator"), attr(draws, "has_denominator")
  )
  expect_identical(aggregate_unit_names(cells), aggregate_unit_names(draws))
  expect_identical(draws_has_denominator(cells), draws_has_denominator(draws))
})
