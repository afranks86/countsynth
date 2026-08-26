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
  bpnmf_schema(
    unit_col = "unit", time_col = "time", treatment_col = "treatment",
    outcomes_from_prefixes = bpnmf_prefixes(
      outcome_prefix = "births_", denominator_prefix = "pop_"
    )
  )
}

wide_config <- function(types, ...) {
  bpnmf_config(
    input_file = "unused.csv", output_dir = tempdir(),
    schema = wide_schema(),
    model = bpnmf_model_opts(types = types),
    ...
  )
}

test_that("prefix resolution finds labels and validates denominators", {
  wide <- make_wide_df()
  cfg <- wide_config(
    list(both = bpnmf_type(groups = c("g1", "g2"), ranks_to_test = 2))
  )
  dat <- bpnmf_data(cfg, df = wide)
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
  schema <- bpnmf_schema(
    unit_col = "unit", time_col = "time", treatment_col = "treatment",
    outcomes_from_prefixes = bpnmf_prefixes(
      "births_", "pop_",
      include = c("g1", "nope")
    )
  )
  cfg <- bpnmf_config(
    input_file = "unused.csv", output_dir = tempdir(), schema = schema,
    model = bpnmf_model_opts(
      types = list(t = bpnmf_type(groups = "g1", ranks_to_test = 2))
    )
  )
  expect_error(bpnmf_data(cfg, df = wide), "nope")
})

test_that("synthetic total group sums subgroups", {
  wide <- make_wide_df()
  cfg <- wide_config(
    list(total = bpnmf_type(
      groups = "total", ranks_to_test = 2, total_all = TRUE
    ))
  )
  dat <- bpnmf_data(cfg, df = wide)
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
    list(total = bpnmf_type(groups = "total", ranks_to_test = 2))
  )
  expect_error(bpnmf_data(cfg, df = wide), "total_from")
})

test_that("total_from validates undefined and duplicate labels", {
  wide <- make_wide_df()
  cfg <- wide_config(
    list(total = bpnmf_type(
      groups = "total", ranks_to_test = 2, total_from = c("g1", "zz")
    ))
  )
  expect_error(bpnmf_data(cfg, df = wide), "undefined")
})

test_that("date filter is start-inclusive and end-exclusive", {
  wide <- make_wide_df()
  cfg <- wide_config(
    list(both = bpnmf_type(groups = c("g1", "g2"), ranks_to_test = 2)),
    start_date = "2020-02-01", end_date = "2020-04-01"
  )
  dat <- bpnmf_data(cfg, df = wide)
  expect_identical(
    dat$times, seq(as.Date("2020-02-01"), by = "month", length.out = 2)
  )
})

test_that("exclude_units drops units", {
  wide <- make_wide_df()
  cfg <- wide_config(
    list(both = bpnmf_type(
      groups = c("g1", "g2"), ranks_to_test = 2, exclude_units = "B"
    ))
  )
  dat <- bpnmf_data(cfg, df = wide)
  expect_identical(dat$units, c("A", "C"))
})

test_that("non-positive denominators are a hard error", {
  wide <- make_wide_df()
  wide$pop_g1[3] <- 0
  cfg <- wide_config(
    list(both = bpnmf_type(groups = c("g1", "g2"), ranks_to_test = 2))
  )
  expect_error(bpnmf_data(cfg, df = wide), "non-positive denominator")
})

test_that("duplicate (group, unit, time) rows are a hard error", {
  wide <- make_wide_df()
  wide <- rbind(wide, wide[1, ])
  cfg <- wide_config(
    list(both = bpnmf_type(groups = c("g1", "g2"), ranks_to_test = 2))
  )
  expect_error(bpnmf_data(cfg, df = wide), "[Dd]uplicate")
})

test_that("unbalanced panel errors unless allowed, then marks missing", {
  wide <- make_wide_df()
  wide <- wide[-2, ]
  cfg <- wide_config(
    list(both = bpnmf_type(groups = c("g1", "g2"), ranks_to_test = 2))
  )
  expect_error(bpnmf_data(cfg, df = wide), "Unbalanced panel")
  cfg2 <- wide_config(
    list(both = bpnmf_type(groups = c("g1", "g2"), ranks_to_test = 2)),
    allow_unbalanced_panel = TRUE
  )
  expect_warning(dat <- bpnmf_data(cfg2, df = wide), "structurally")
  expect_equal(sum(dat$missing_idx_array), 2) # both groups of the dropped row
})

test_that("temporal aggregation sums outcomes, maxes treatment, means denominators", {
  wide <- make_wide_df()
  cfg <- wide_config(
    list(both = bpnmf_type(groups = c("g1", "g2"), ranks_to_test = 2)),
    time_aggregation = bpnmf_time_aggregation(enabled = TRUE, period = "quarterly")
  )
  dat <- bpnmf_data(cfg, df = wide)
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

test_that("n_periods aggregation combines N consecutive periods", {
  # 5 monthly points, blocks of 2 -> Jan+Feb, Mar+Apr, May (partial).
  wide <- make_wide_df()
  cfg <- wide_config(
    list(both = bpnmf_type(groups = c("g1", "g2"), ranks_to_test = 2)),
    time_aggregation = bpnmf_time_aggregation(enabled = TRUE, n_periods = 2)
  )
  expect_warning(dat <- bpnmf_data(cfg, df = wide), "final block holds 1")
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
  cfg <- bpnmf_config(
    input_file = "unused.csv", output_dir = tempdir(),
    schema = wide_schema(),
    model = bpnmf_model_opts(
      types = list(t = bpnmf_type(groups = "total", ranks_to_test = 1))
    ),
    time_aggregation = bpnmf_time_aggregation(enabled = TRUE, n_periods = 7)
  )
  dat <- bpnmf_data(cfg, df = tibble::as_tibble(grid))
  expect_equal(dim(dat$Y)[3], 4)
  expect_equal(dat$Y[1, 1, 1], 70) # 7 days x 10
  wk <- unique(dat$df[, c("start_date", "end_date")])
  wk <- wk[order(wk$start_date), ]
  expect_equal(as.Date(wk$start_date), days[c(1, 8, 15, 22)])
  expect_equal(as.Date(wk$end_date), days[c(7, 14, 21, 28)])

  # Calendar monthly collapses the same panel to a single time point.
  cfg$time_aggregation <- bpnmf_time_aggregation(enabled = TRUE, period = "monthly")
  expect_equal(dim(bpnmf_data(cfg, df = tibble::as_tibble(grid))$Y)[3], 1)
})

test_that("period and n_periods are mutually exclusive", {
  expect_error(
    bpnmf_time_aggregation(enabled = TRUE, period = "monthly", n_periods = 3),
    "only one of"
  )
  expect_null(bpnmf_time_aggregation()$period)
  # Enabling without either keeps the historical bimonthly default.
  expect_equal(bpnmf_time_aggregation(enabled = TRUE)$period, "bimonthly")
  expect_null(bpnmf_time_aggregation(enabled = TRUE)$n_periods)
})

test_that("date auto-parsing handles multiple formats", {
  wide <- make_wide_df()
  wide$time <- format(as.Date(wide$time), "%m/%d/%Y")
  cfg <- wide_config(
    list(both = bpnmf_type(groups = c("g1", "g2"), ranks_to_test = 2))
  )
  dat <- bpnmf_data(cfg, df = wide)
  expect_s3_class(dat$times, "Date")
  expect_equal(min(dat$times), as.Date("2020-01-01"))
})
