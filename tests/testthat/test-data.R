make_wide_df <- function() {
  df <- make_test_data()
  df |>
    tidyr::pivot_wider(
      id_cols = c("unit", "time", "treatment"),
      names_from = "group",
      values_from = c("outcome", "denominator"),
      names_glue = "{ifelse(.value == 'outcome', 'births', 'pop')}_{group}"
    )
}

wide_schema <- function() {
  countsynth_schema(
    unit_col = "unit", time_col = "time", treatment_col = "treatment",
    outcomes_from_prefixes = countsynth_prefixes(
      outcome_prefix = "births_", denominator_prefix = "pop_"
    )
  )
}

wide_config <- function(types, ...) {
  countsynth_config(
    input_file = "unused.csv", output_dir = tempdir(),
    schema = wide_schema(),
    model = countsynth_model_opts(types = types),
    ...
  )
}

test_that("prefix resolution finds labels and validates denominators", {
  wide <- make_wide_df()
  cfg <- wide_config(
    list(both = countsynth_type(groups = c("g1", "g2"), ranks_to_test = 2))
  )
  dat <- countsynth_data(cfg, df = wide)
  expect_identical(dat$groups, c("g1", "g2"))
  expect_identical(dat$units, c("A", "B", "C"))
  expect_equal(dim(dat$Y), c(2, 3, 5))
  # denominators scaled per 10k
  raw <- make_test_data()
  a_g1_t1 <- raw$denominator[raw$unit == "A" & raw$group == "g1" &
    raw$time == min(raw$time)]
  expect_equal(dat$denominators[1, 1, 1], a_g1_t1 / 1e4)
  # treatment: unit C exposed from period 4
  expect_true(all(dat$control_idx_array[, 1:2, ]))
  expect_identical(unname(dat$control_idx_array[1, 3, ]), c(TRUE, TRUE, TRUE, FALSE, FALSE))
})

test_that("include filter and missing include labels error correctly", {
  wide <- make_wide_df()
  schema <- countsynth_schema(
    unit_col = "unit", time_col = "time", treatment_col = "treatment",
    outcomes_from_prefixes = countsynth_prefixes(
      "births_", "pop_",
      include = c("g1", "nope")
    )
  )
  cfg <- countsynth_config(
    input_file = "unused.csv", output_dir = tempdir(), schema = schema,
    model = countsynth_model_opts(
      types = list(t = countsynth_type(groups = "g1", ranks_to_test = 2))
    )
  )
  expect_error(countsynth_data(cfg, df = wide), "nope")
})

test_that("synthetic total group sums subgroups", {
  wide <- make_wide_df()
  cfg <- wide_config(
    list(total = countsynth_type(
      groups = "total", ranks_to_test = 2, total_all = TRUE
    ))
  )
  dat <- countsynth_data(cfg, df = wide)
  expect_identical(dat$groups, "total")
  raw <- make_test_data()
  expected_total <- sum(raw$outcome[raw$unit == "A" & raw$time == min(raw$time)])
  expect_equal(dat$Y[1, 1, 1], expected_total)
  # denominator also summed for the synthetic total
  expected_denom <- sum(raw$denominator[raw$unit == "A" & raw$time == min(raw$time)])
  expect_equal(dat$denominators[1, 1, 1], expected_denom / 1e4)
})

test_that("total requires total_from or total_all", {
  wide <- make_wide_df()
  cfg <- wide_config(
    list(total = countsynth_type(groups = "total", ranks_to_test = 2))
  )
  expect_error(countsynth_data(cfg, df = wide), "total_from")
})

test_that("total_from validates undefined and duplicate labels", {
  wide <- make_wide_df()
  cfg <- wide_config(
    list(total = countsynth_type(
      groups = "total", ranks_to_test = 2, total_from = c("g1", "zz")
    ))
  )
  expect_error(countsynth_data(cfg, df = wide), "undefined")
})

test_that("date filter is start-inclusive and end-exclusive", {
  wide <- make_wide_df()
  cfg <- wide_config(
    list(both = countsynth_type(groups = c("g1", "g2"), ranks_to_test = 2)),
    start_date = "2020-02-01", end_date = "2020-04-01"
  )
  dat <- countsynth_data(cfg, df = wide)
  expect_identical(
    dat$times, seq(as.Date("2020-02-01"), by = "month", length.out = 2)
  )
})

test_that("exclude_units drops units", {
  wide <- make_wide_df()
  cfg <- wide_config(
    list(both = countsynth_type(
      groups = c("g1", "g2"), ranks_to_test = 2, exclude_units = "B"
    ))
  )
  dat <- countsynth_data(cfg, df = wide)
  expect_identical(dat$units, c("A", "C"))
})

test_that("a zero denominator against a positive outcome is a hard error", {
  wide <- make_wide_df()
  wide$pop_g1[3] <- 0
  cfg <- wide_config(
    list(both = countsynth_type(groups = c("g1", "g2"), ranks_to_test = 2))
  )
  # births_g1[3] is a positive rpois() draw, so this is a real outcome >
  # denominator (0) mismatch, not the deterministic (0, 0) case below.
  expect_true(wide$births_g1[3] > 0)
  expect_error(countsynth_data(cfg, df = wide), "Invalid outcome/denominator")
})

test_that("negative outcomes and negative denominators are hard errors", {
  wide <- make_wide_df()
  cfg <- wide_config(
    list(both = countsynth_type(groups = c("g1", "g2"), ranks_to_test = 2))
  )
  neg_outcome <- wide
  neg_outcome$births_g1[3] <- -1
  expect_error(countsynth_data(cfg, df = neg_outcome), "Invalid outcome/denominator")

  neg_denom <- wide
  neg_denom$pop_g1[3] <- -1
  expect_error(countsynth_data(cfg, df = neg_denom), "Invalid outcome/denominator")
})

test_that("an outcome exceeding a positive denominator is not an error", {
  # The denominator is an exposure/offset, not a trial count, so it need not
  # bound the outcome -- e.g. a denominator reported in thousands while the
  # outcome is a raw count.
  wide <- make_wide_df()
  wide$pop_g1[3] <- wide$births_g1[3] - 1
  expect_true(wide$pop_g1[3] > 0) # a real, positive-but-smaller denominator
  cfg <- wide_config(
    list(both = countsynth_type(groups = c("g1", "g2"), ranks_to_test = 2))
  )
  dat <- expect_silent(countsynth_data(cfg, df = wide))
  expect_false(any(dat$excluded_idx_array))
})

test_that("a (0, 0) outcome/denominator cell is valid and excluded from the likelihood", {
  wide <- make_wide_df()
  wide$pop_g1[3] <- 0
  wide$births_g1[3] <- 0
  cfg <- wide_config(
    list(both = countsynth_type(groups = c("g1", "g2"), ranks_to_test = 2))
  )
  dat <- expect_silent(countsynth_data(cfg, df = wide))
  expect_equal(sum(dat$excluded_idx_array), 1)
  # unit C (row 3), group g1, first time period.
  expect_true(dat$excluded_idx_array[1, 3, 1])
  expect_false(any(dat$excluded_idx_array[2, , ])) # g2 untouched
  expect_equal(dat$Y[1, 3, 1], 0)
  expect_false(dat$missing_idx_array[1, 3, 1]) # excluded, not missing/suppressed

  sd <- stan_data_joint(dat, rank = 2)
  excluded_cell <- flat_idx(1, 3, 1, length(dat$units), length(dat$times))
  expect_false(excluded_cell %in% sd$obs_cell)
  expect_false(excluded_cell %in% sd$cens_cell)
  expect_false(excluded_cell %in% sd$notcens_cell)
})

test_that("a bad denominator outside the analysis window is not an error", {
  # Reported case: a 0 denominator in a period the config's end_date excludes.
  # Those rows never reach the model, so failing on them rejects a good run.
  wide <- make_wide_df()
  last <- max(as.Date(wide$time))
  wide$pop_g1[as.Date(wide$time) == last] <- 0
  types <- list(both = countsynth_type(groups = c("g1", "g2"), ranks_to_test = 2))

  expect_error(countsynth_data(wide_config(types), df = wide), "Invalid outcome/denominator")
  # end_date is exclusive, so this drops exactly the offending period.
  ok <- countsynth_data(
    wide_config(types, end_date = format(last)), df = wide
  )
  expect_false(last %in% as.Date(ok$df$time))
  expect_true(all(ok$df$denominator > 0))

  # ... and start_date on the other end.
  first <- min(as.Date(wide$time))
  wide2 <- make_wide_df()
  wide2$pop_g1[as.Date(wide2$time) == first] <- 0
  expect_error(countsynth_data(wide_config(types), df = wide2), "Invalid outcome/denominator")
  expect_silent(
    countsynth_data(wide_config(types, start_date = format(first + 1)), df = wide2)
  )
})

test_that("a bad denominator in an excluded unit or unmodeled group is ignored", {
  wide <- make_wide_df()
  wide$pop_g1[wide$unit == "C"] <- 0
  # Excluded units are dropped before the check.
  expect_silent(countsynth_data(
    wide_config(list(both = countsynth_type(
      groups = c("g1", "g2"), ranks_to_test = 2, exclude_units = "C"
    ))),
    df = wide
  ))
  # Same rows, unit not excluded: still a hard error.
  expect_error(
    countsynth_data(
      wide_config(list(both = countsynth_type(c("g1", "g2"), ranks_to_test = 2))),
      df = wide
    ),
    "Invalid outcome/denominator"
  )
  # A group the type does not model is never validated.
  expect_silent(countsynth_data(
    wide_config(list(one = countsynth_type(groups = "g2", ranks_to_test = 2))),
    df = wide
  ))
})

test_that("the denominator error names the periods and units involved", {
  wide <- make_wide_df()
  bad_time <- sort(unique(as.Date(wide$time)))[2]
  wide$pop_g1[as.Date(wide$time) == bad_time & wide$unit == "B"] <- 0
  expect_error(
    countsynth_data(
      wide_config(list(both = countsynth_type(c("g1", "g2"), ranks_to_test = 2))),
      df = wide
    ),
    format(bad_time),
    fixed = TRUE
  )
})

test_that("duplicate (group, unit, time) rows are a hard error", {
  wide <- make_wide_df()
  wide <- rbind(wide, wide[1, ])
  cfg <- wide_config(
    list(both = countsynth_type(groups = c("g1", "g2"), ranks_to_test = 2))
  )
  expect_error(countsynth_data(cfg, df = wide), "[Dd]uplicate")
})

test_that("unbalanced panel errors unless allowed, then marks missing", {
  wide <- make_wide_df()
  wide <- wide[-2, ]
  cfg <- wide_config(
    list(both = countsynth_type(groups = c("g1", "g2"), ranks_to_test = 2))
  )
  expect_error(countsynth_data(cfg, df = wide), "Unbalanced panel")
  cfg2 <- wide_config(
    list(both = countsynth_type(groups = c("g1", "g2"), ranks_to_test = 2)),
    allow_unbalanced_panel = TRUE
  )
  expect_warning(dat <- countsynth_data(cfg2, df = wide), "structurally")
  expect_equal(sum(dat$missing_idx_array), 2) # both groups of the dropped row
})

test_that("temporal aggregation sums outcomes, maxes treatment, means denominators", {
  wide <- make_wide_df()
  cfg <- wide_config(
    list(both = countsynth_type(groups = c("g1", "g2"), ranks_to_test = 2)),
    time_aggregation = countsynth_time_aggregation(enabled = TRUE, period = "quarterly")
  )
  dat <- countsynth_data(cfg, df = wide)
  raw <- make_test_data()
  # Q1 2020 = Jan+Feb+Mar for unit A group g1
  q1 <- raw[raw$unit == "A" & raw$group == "g1" &
    raw$time < as.Date("2020-04-01"), ]
  expect_equal(dat$Y[1, 1, 1], sum(q1$outcome))
  expect_equal(dat$denominators[1, 1, 1], mean(q1$denominator) / 1e4)
  # treated unit C: Q2 (Apr+May) has treatment (Apr 1 onward)
  expect_false(dat$control_idx_array[1, 3, 2])
  # period boundary columns exist for person-years
  expect_true(all(c("start_date", "end_date") %in% names(dat$df)))
  expect_equal(
    as.Date(dat$df$end_date[dat$df$start_date == as.Date("2020-01-01")][1]),
    as.Date("2020-03-31")
  )
})

test_that("no denominator means denominator 1 everywhere, aggregation included", {
  # Regression: temporal aggregation used to fabricate a literal denominator
  # of 1 for a denominator-less outcome, which then got divided by
  # denominator_scale downstream (1e-4) instead of staying 1 -- silently
  # rescaling every rate the moment aggregation was turned on.
  grid <- expand.grid(
    unit = c("A", "B"),
    time = seq(as.Date("2020-01-01"), by = "month", length.out = 6),
    stringsAsFactors = FALSE
  )
  grid$treatment <- 0L
  grid$outcome_x <- 3L
  cfg <- countsynth_config(
    input_file = "unused.csv", output_dir = tempdir(),
    schema = list(
      unit_col = "unit", time_col = "time", treatment_col = "treatment",
      outcomes = list(list(outcome_col = "outcome_x", label = "x"))
    ),
    model = countsynth_model_opts(
      types = list(t = countsynth_type(groups = "x", ranks_to_test = 1))
    )
  )
  flat <- countsynth_data(cfg, df = tibble::as_tibble(grid))
  expect_true(all(flat$denominators == 1))

  cfg$time_aggregation <- countsynth_time_aggregation(enabled = TRUE, period = "quarterly")
  agg <- countsynth_data(cfg, df = tibble::as_tibble(grid))
  expect_true(all(agg$denominators == 1))
})

test_that("n_periods aggregation combines N consecutive periods", {
  # 5 monthly points, blocks of 2 -> Jan+Feb, Mar+Apr, May (partial).
  wide <- make_wide_df()
  cfg <- wide_config(
    list(both = countsynth_type(groups = c("g1", "g2"), ranks_to_test = 2)),
    time_aggregation = countsynth_time_aggregation(enabled = TRUE, n_periods = 2)
  )
  expect_warning(dat <- countsynth_data(cfg, df = wide), "final block holds 1")
  expect_equal(dim(dat$Y)[3], 3)

  raw <- make_test_data()
  janfeb <- raw[raw$unit == "A" & raw$group == "g1" &
    raw$time < as.Date("2020-03-01"), ]
  expect_equal(dat$Y[1, 1, 1], sum(janfeb$outcome))
  expect_equal(dat$denominators[1, 1, 1], mean(janfeb$denominator) / 1e4)

  # Blocks are anchored to the panel start, not the calendar, and the short
  # trailing block records its real one-month exposure.
  bounds <- unique(dat$df[, c("start_date", "end_date")])
  bounds <- bounds[order(bounds$start_date), ]
  expect_equal(as.Date(bounds$start_date),
               as.Date(c("2020-01-01", "2020-03-01", "2020-05-01")))
  # Each block ends the day before the next input period starts; the trailing
  # block has no successor, so its width is extrapolated from the last gap.
  expect_equal(as.Date(bounds$end_date),
               as.Date(c("2020-02-29", "2020-04-30", "2020-05-30")))
})

test_that("n_periods works at daily resolution, where calendar bins cannot", {
  # 28 daily points in one month: `period` has nothing to bin on, n_periods
  # turns them into 4 weeks.
  days <- seq(as.Date("2021-03-01"), by = "day", length.out = 28)
  grid <- expand.grid(unit = c("A", "B"), time = days, stringsAsFactors = FALSE)
  grid$treatment <- as.integer(grid$unit == "B" & grid$time >= days[15])
  grid$births_total <- 10
  grid$pop_total <- 1000
  cfg <- countsynth_config(
    input_file = "unused.csv", output_dir = tempdir(),
    schema = wide_schema(),
    model = countsynth_model_opts(
      types = list(t = countsynth_type(groups = "total", ranks_to_test = 1))
    ),
    time_aggregation = countsynth_time_aggregation(enabled = TRUE, n_periods = 7)
  )
  dat <- countsynth_data(cfg, df = tibble::as_tibble(grid))
  expect_equal(dim(dat$Y)[3], 4)
  expect_equal(dat$Y[1, 1, 1], 70) # 7 days x 10
  wk <- unique(dat$df[, c("start_date", "end_date")])
  wk <- wk[order(wk$start_date), ]
  expect_equal(as.Date(wk$start_date), days[c(1, 8, 15, 22)])
  expect_equal(as.Date(wk$end_date), days[c(7, 14, 21, 28)])

  # Calendar monthly collapses the same panel to a single time point.
  cfg$time_aggregation <- countsynth_time_aggregation(enabled = TRUE, period = "monthly")
  expect_equal(dim(countsynth_data(cfg, df = tibble::as_tibble(grid))$Y)[3], 1)
})

test_that("period and n_periods are mutually exclusive", {
  expect_error(
    countsynth_time_aggregation(enabled = TRUE, period = "monthly", n_periods = 3),
    "only one of"
  )
  expect_null(countsynth_time_aggregation()$period)
  # Enabling without either keeps the historical bimonthly default.
  expect_equal(countsynth_time_aggregation(enabled = TRUE)$period, "bimonthly")
  expect_null(countsynth_time_aggregation(enabled = TRUE)$n_periods)
})

test_that("date auto-parsing handles multiple formats", {
  wide <- make_wide_df()
  wide$time <- format(as.Date(wide$time), "%m/%d/%Y")
  cfg <- wide_config(
    list(both = countsynth_type(groups = c("g1", "g2"), ranks_to_test = 2))
  )
  dat <- countsynth_data(cfg, df = wide)
  expect_s3_class(dat$times, "Date")
  expect_equal(min(dat$times), as.Date("2020-01-01"))
})
