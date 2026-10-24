# Suppressed and absent cells: which reach the likelihood, and on what range.

# Two units, six months of 50 deaths against 4000 births each.
monthly_wide <- function(months = 6, by = "month") {
  times <- seq(as.Date("2020-01-01"), by = by, length.out = months)
  w <- expand.grid(unit = c("A", "B"), time = times, stringsAsFactors = FALSE)
  w$treatment <- 0L
  w$deaths_x <- 50
  w$births_x <- 4000
  tibble::as_tibble(w)
}

mortality_config <- function(time_aggregation = countsynth_time_aggregation(),
                             model = list(), ...) {
  countsynth_config(
    input_file = "unused.csv", output_dir = tempdir(),
    schema = countsynth_schema(
      unit_col = "unit", time_col = "time", treatment_col = "treatment",
      outcomes = list(
        countsynth_outcome("deaths_x", "x", denominator_col = "births_x")
      )
    ),
    model = do.call(countsynth_model_opts, c(
      list(types = list(t = countsynth_type(groups = "x", ranks_to_test = 1))),
      model
    )),
    time_aggregation = time_aggregation,
    output = countsynth_output_opts(denominator_time_unit = "none"),
    ...
  )
}

quarterly <- countsynth_time_aggregation(enabled = TRUE, period = "quarterly")

suppress_cell <- function(w, unit, month) {
  w$deaths_x[w$unit == unit & w$time == as.Date(month)] <- NA
  w
}

test_that("a lone suppressed cell is censored on the suppression range", {
  w <- suppress_cell(monthly_wide(), "A", "2020-02-01")
  dat <- countsynth_data(mortality_config(), df = w)
  expect_true(dat$missing_idx_array[1, 1, 2])
  expect_equal(c(dat$cens_lower[1, 1, 2], dat$cens_upper[1, 1, 2]), c(1, 9))
  sd <- stan_data_joint(dat, rank = 1)
  expect_equal(sd$n_cens, 1)
  expect_equal(sd$n_icens, 0)
  expect_equal(c(sd$supp_lo, sd$supp_hi), c(1L, 9L))
})

test_that("a quarter pooling a suppressed month is censored on its summed range", {
  w <- suppress_cell(monthly_wide(), "A", "2020-02-01")
  expect_warning(
    dat <- countsynth_data(mortality_config(quarterly), df = w),
    "pools a suppressed period"
  )
  # Jan and Mar observed at 50 each, Feb somewhere in 1..9.
  expect_true(is.na(dat$Y[1, 1, 1]) || dat$missing_idx_array[1, 1, 1])
  expect_equal(c(dat$cens_lower[1, 1, 1], dat$cens_upper[1, 1, 1]), c(101, 109))
  sd <- stan_data_joint(dat, rank = 1)
  expect_equal(sd$n_cens, 0)
  expect_equal(sd$n_icens, 1)
  expect_equal(c(sd$icens_lo, sd$icens_hi), c(101L, 109L))
  # Births still sum over all three months.
  expect_equal(dat$denominators[1, 1, 1], 12000 / 1e4)
})

test_that("a fully suppressed quarter spans the range once per month", {
  w <- monthly_wide()
  for (m in c("2020-01-01", "2020-02-01", "2020-03-01")) {
    w <- suppress_cell(w, "A", m)
  }
  dat <- suppressWarnings(countsynth_data(mortality_config(quarterly), df = w))
  expect_equal(c(dat$cens_lower[1, 1, 1], dat$cens_upper[1, 1, 1]), c(3, 27))
})

test_that("suppression_range sets the range, alone and pooled", {
  w <- suppress_cell(monthly_wide(), "A", "2020-02-01")
  model <- list(suppression_range = c(1, 10))
  lone <- countsynth_data(mortality_config(model = model), df = w)
  sd <- stan_data_joint(lone, rank = 1)
  expect_equal(c(sd$supp_lo, sd$supp_hi), c(1L, 10L))
  expect_equal(sd$n_cens, 1)
  expect_equal(sd$n_icens, 0)

  pooled <- suppressWarnings(
    countsynth_data(mortality_config(quarterly, model = model), df = w)
  )
  expect_equal(
    c(pooled$cens_lower[1, 1, 1], pooled$cens_upper[1, 1, 1]), c(101, 110)
  )
})

test_that("suppression_range is validated and read from YAML", {
  expect_equal(countsynth_model_opts()$suppression_range, c(1L, 9L))
  expect_error(countsynth_model_opts(suppression_range = c(9, 1)), "lower <= upper")
  expect_error(countsynth_model_opts(suppression_range = c(-1, 9)))
  expect_error(countsynth_model_opts(suppression_range = 5))
  expect_equal(
    parse_yaml_model(list(suppression_range = list(0, 4)))$suppression_range,
    c(0L, 4L)
  )
})

test_that("a synthetic total sums its subgroups' ranges", {
  times <- seq(as.Date("2020-01-01"), by = "month", length.out = 4)
  w <- expand.grid(unit = c("A", "B"), time = times, stringsAsFactors = FALSE)
  w$treatment <- 0L
  w$deaths_a <- 20
  w$deaths_b <- 30
  w$births_a <- 1000
  w$births_b <- 2000
  w$deaths_a[w$unit == "A" & w$time == times[1]] <- NA
  cfg <- countsynth_config(
    input_file = "unused.csv", output_dir = tempdir(),
    schema = countsynth_schema(
      unit_col = "unit", time_col = "time", treatment_col = "treatment",
      outcomes_from_prefixes = countsynth_prefixes(
        outcome_prefix = "deaths_", denominator_prefix = "births_"
      )
    ),
    model = countsynth_model_opts(types = list(t = countsynth_type(
      groups = "total", total_from = c("a", "b"), ranks_to_test = 1
    )))
  )
  dat <- countsynth_data(cfg, df = tibble::as_tibble(w))
  # b's 30 is known; a's count is 1..9.
  expect_true(dat$missing_idx_array[1, 1, 1])
  expect_equal(c(dat$cens_lower[1, 1, 1], dat$cens_upper[1, 1, 1]), c(31, 39))
  expect_equal(stan_data_joint(dat, rank = 1)$n_icens, 1)
})

test_that("an absent month inside a block is an unbalanced panel", {
  w <- monthly_wide()
  w <- w[!(w$unit == "A" & w$time == as.Date("2020-02-01")), ]
  expect_error(
    suppressWarnings(countsynth_data(mortality_config(quarterly), df = w)),
    "missing input periods"
  )
  cfg <- mortality_config(quarterly, allow_unbalanced_panel = TRUE)
  warnings <- character()
  dat <- withCallingHandlers(
    countsynth_data(cfg, df = w),
    warning = function(cnd) {
      warnings <<- c(warnings, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_match(paste(warnings, collapse = "\n"), "left out of the likelihood")
  # A's Q1 is gone rather than summed over two months, and absent rather than
  # suppressed, so nothing censors it.
  expect_true(dat$excluded_idx_array[1, 1, 1])
  expect_false(dat$missing_idx_array[1, 1, 1])
  expect_equal(dat$Y[1, 2, 1], 150) # B's Q1 is untouched
  sd <- stan_data_joint(dat, rank = 1)
  expect_equal(sd$n_cens + sd$n_icens, 0)
})

test_that("a bin the date window cuts short carries its real coverage", {
  cfg <- mortality_config(quarterly, end_date = "2020-05-01")
  expect_warning(
    dat <- countsynth_data(cfg, df = monthly_wide()),
    "only partly covered"
  )
  q2 <- dat$df[dat$df$time == as.Date("2020-04-01"), ][1, ]
  expect_equal(as.Date(q2$start_date), as.Date("2020-04-01"))
  expect_equal(as.Date(q2$end_date), as.Date("2020-04-30"))
  # A full quarter keeps its calendar dates.
  q1 <- dat$df[dat$df$time == as.Date("2020-01-01"), ][1, ]
  expect_equal(as.Date(q1$end_date), as.Date("2020-03-31"))
})

test_that("input coarser than the bin keeps its own dates", {
  # Half-years in quarterly bins: nothing to combine, and labelling each as a
  # quarter would halve its person-time.
  w <- monthly_wide(months = 4, by = "6 months")
  expect_warning(
    dat <- countsynth_data(mortality_config(quarterly), df = w),
    "longer than a quarterly bin"
  )
  jan <- dat$df[dat$df$time == as.Date("2020-01-01"), ][1, ]
  expect_equal(as.Date(jan$end_date), as.Date("2020-06-30"))
})

test_that("an un-aggregated panel carries each period's length", {
  # Without bounds every row counted as a whole year: a bimonthly panel's
  # person-years were 6x too large.
  w <- monthly_wide(months = 6, by = "2 months")
  dat <- countsynth_data(mortality_config(), df = w)
  expect_true(all(c("start_date", "end_date") %in% names(dat$df)))
  first <- dat$df[dat$df$time == as.Date("2020-01-01"), ][1, ]
  expect_equal(as.Date(first$end_date), as.Date("2020-02-29"))
  last <- dat$df[dat$df$time == max(dat$df$time), ][1, ]
  expect_equal(as.Date(last$end_date), as.Date("2020-12-31"))
  expect_equal(mean(years_per_row(dat$df)), 1 / 6, tolerance = 0.01)
})

test_that("a panel-wide gap does not stretch the period before it", {
  w <- monthly_wide(months = 6)
  w <- w[w$time != as.Date("2020-04-01"), ]
  dat <- countsynth_data(mortality_config(), df = w)
  mar <- dat$df[dat$df$time == as.Date("2020-03-01"), ][1, ]
  expect_equal(as.Date(mar$end_date), as.Date("2020-03-31"))
})

test_that("cut stage 2 splits exposed suppressed cells the same way", {
  w <- monthly_wide()
  w$treatment[w$unit == "B" & w$time >= as.Date("2020-04-01")] <- 1L
  w <- suppress_cell(w, "B", "2020-05-01")
  dat <- suppressWarnings(countsynth_data(mortality_config(quarterly), df = w))
  mu <- rep(log(150), length(dat$Y))
  sd2 <- stan_data_stage2(dat, mu, phi_unit = c(10, 10))
  expect_equal(sd2$n_cens, 0)
  expect_equal(sd2$n_icens, 1)
  expect_equal(c(sd2$icens_lo, sd2$icens_hi), c(101L, 109L))
  expect_equal(c(sd2$supp_lo, sd2$supp_hi), c(1L, 9L))
})

test_that("an absent member cell is not imputed into a pooled total", {
  # What countsynth_draws() emits for a cell with no row: every observed
  # column NA, and a prediction built on a placeholder exposure.
  times <- seq(as.Date("2021-01-01"), by = "month", length.out = 4)
  g <- expand.grid(unit = c("A", "B", "C"), time = times, stringsAsFactors = FALSE)
  g$treatment <- as.integer(g$unit != "A" & g$time >= times[3])
  g$outcome <- 100
  g$denominator <- 1000
  absent <- g$unit == "C" & g$time == times[2]
  g$outcome[absent] <- NA
  g$denominator[absent] <- NA
  g$treatment[absent] <- NA
  g$mu <- ifelse(absent, log(1000), log(100))
  g$mu_treated <- g$mu
  g$ypred <- exp(g$mu)
  frame <- dplyr::bind_rows(
    dplyr::mutate(g, group = "total"), dplyr::mutate(g, group = "other")
  )
  frame$.draw <- 1L
  class(frame) <- c("countsynth_draws", class(frame))
  attr(frame, "has_denominator") <- TRUE

  pooled <- add_aggregate_units(
    frame,
    list(countsynth_aggregate_unit("All treated", include_treated_units = TRUE))
  )
  cell <- pooled[pooled$unit == "All treated" & pooled$group == "total" &
    pooled$time == times[2], ]
  expect_true(is.na(cell$outcome))
  expect_false(cell$outcome_imputed)

  # The pooled group's raw rate no longer spikes where C is absent.
  cells <- add_aggregate_groups(
    observed_cell_frame(frame),
    list(countsynth_aggregate_group("All", include_all_groups = TRUE))
  )
  p <- countsynth_raw_rate_plot(cells, group = "All")
  treated <- p$data[p$data$treatment_group == "Treated", ]
  expect_equal(treated$y_value, rep(100, 4))

  # A suppressed cell -- outcome NA on a reported row -- is still filled.
  frame$outcome[frame$unit == "B" & frame$time == times[1]] <- NA
  pooled <- add_aggregate_units(
    frame,
    list(countsynth_aggregate_unit("All treated", include_treated_units = TRUE))
  )
  filled <- pooled[pooled$unit == "All treated" & pooled$group == "total" &
    pooled$time == times[1], ]
  expect_true(filled$outcome_imputed)
  expect_equal(filled$outcome, 200)
})
