# The reporting frame's size is known from the cell counts before any of it
# is built. report_size_guard() is what turns the package's worst failure --
# a report that runs for hours over tens of millions of rows, is killed by
# the OS, and writes nothing -- into a message at second zero.

test_that("the guard stays quiet when the frame fits", {
  draws <- make_draws_frame()
  expect_silent(report_size_guard(draws, NULL, NULL, budget = 1e9))
})

test_that("the guard projects the aggregate rows it has not built yet", {
  draws <- make_draws_frame()
  base <- report_size_guard(draws, NULL, NULL, budget = Inf)
  expect_equal(base, nrow(draws))

  # One pooled unit covers every (draw, time, group); one pooled group covers
  # every (draw, time, unit), the pooled unit included. The group axis is the
  # expensive one, and the projection has to say so before the work starts.
  n_draws <- max(draws$.draw)
  n_times <- length(unique(draws$time))
  n_units <- length(unique(draws$unit))
  n_groups <- length(unique(draws$group))
  with_both <- report_size_guard(
    draws,
    list(countsynth_aggregate_unit("all", include_all_units = TRUE)),
    list(countsynth_aggregate_group("all", include_all_groups = TRUE)),
    budget = Inf
  )
  expect_equal(
    with_both,
    nrow(draws) + n_draws * n_times * n_groups +
      n_draws * n_times * (n_units + 1)
  )
})

test_that("the guard names a thinning factor when the frame will not fit", {
  draws <- make_draws_frame()
  expect_warning(
    report_size_guard(draws, NULL, NULL, budget = 1),
    "Projected reporting frame"
  )
  msg <- tryCatch(
    report_size_guard(draws, NULL, NULL, budget = 1),
    warning = function(w) conditionMessage(w)
  )
  # Actionable, not just alarming: the message has to carry the remedy.
  expect_match(msg, "Thin to 1 draw in")
  expect_match(msg, "mcmc.thinning")
})

test_that("the guard tolerates a frame it cannot project", {
  expect_null(report_size_guard(make_draws_frame()[0, ], NULL, NULL))
  no_draw <- make_draws_frame()
  no_draw$.draw <- NULL
  expect_null(report_size_guard(no_draw, NULL, NULL))
})
