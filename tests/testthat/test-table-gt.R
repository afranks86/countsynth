# gt rendering and the consolidated report table set.

test_that("bpnmf_gt_table renders one unit and every unit", {
  skip_if_not_installed("gt")
  draws <- make_draws_frame()

  one <- bpnmf_gt_table(draws)
  expect_s3_class(one, "gt_tbl")
  # Imputed becomes a source note, not a column.
  expect_false("Imputed" %in% names(one$`_data`))
  expect_true(grepl("^B", one$`_heading`$title))

  many <- bpnmf_gt_table(draws, by_unit = TRUE)
  expect_s3_class(many, "gt_tbl")
  expect_equal(many$`_heading`$title, "Post-treatment effect by unit")
  # Unit becomes the row-group column rather than a data column.
  expect_true("Unit" %in% many$`_boxhead`$var[many$`_boxhead`$type == "row_group"])

  # With nothing treated there is no unit to detect, so the abort comes from
  # target detection; by_unit gets past that and hits the empty-table guard.
  untreated <- draws[draws$treatment == 0, ]
  expect_error(bpnmf_gt_table(untreated), "No treated units")
  expect_error(bpnmf_gt_table(untreated, by_unit = TRUE), "No post-treatment rows")
})

test_that("bpnmf_summary_table_by_unit stacks per-unit tables", {
  draws <- make_draws_frame()
  by_unit <- bpnmf_summary_table_by_unit(draws)
  expect_true("Unit" %in% names(by_unit))
  expect_equal(unique(by_unit$Unit), "B") # A is never treated
  # Identical to the single-unit table apart from the Unit column.
  expect_equal(
    by_unit[setdiff(names(by_unit), "Unit")],
    bpnmf_summary_table(draws, "B")
  )
})

test_that("the report writes HTML tables and can be told not to", {
  skip_if_not_installed("gt")
  draws <- make_draws_frame()
  html <- c("summary_table.html", "summary_table_by_unit.html")

  on_dir <- withr::local_tempdir()
  bpnmf_report(draws, on_dir, figures = character(), print_tables = FALSE)
  expect_true(all(file.exists(file.path(on_dir, "figs", html))))

  off_dir <- withr::local_tempdir()
  bpnmf_report(
    draws, off_dir,
    figures = character(), print_tables = FALSE, html_tables = FALSE
  )
  expect_false(any(file.exists(file.path(off_dir, "figs", html))))
  # The CSVs are unaffected either way.
  expect_true(file.exists(file.path(off_dir, "figs", "summary_table_by_unit.csv")))
})

test_that("terminal output is one table by default, two on request", {
  draws <- make_draws_frame()
  # cli headers go to stderr and print() to stdout, so collect both.
  quiet <- function(...) {
    msgs <- character()
    out <- utils::capture.output(withCallingHandlers(
      bpnmf_report(draws, withr::local_tempdir(), figures = character(),
                   html_tables = FALSE, ...),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ))
    paste(c(out, msgs), collapse = "\n")
  }
  default_out <- quiet(print_tables = TRUE)
  expect_match(default_out, "Post-treatment effect by unit")
  expect_false(grepl("observed vs expected", default_out))

  both_out <- quiet(print_tables = TRUE, print_target_table = TRUE)
  expect_match(both_out, "observed vs expected")
  expect_match(both_out, "Post-treatment effect by unit")

  expect_equal(trimws(quiet(print_tables = FALSE)), "")
})

test_that("html_tables flows through the config layers", {
  expect_true(bpnmf_output_opts()$html_tables)
  expect_false(bpnmf_output_opts(html_tables = FALSE)$html_tables)
  expect_true(parse_yaml_output(list())$html_tables)
  expect_false(parse_yaml_output(list(html_tables = FALSE))$html_tables)
  # print_target_table now defaults off, since its rows are in the by-unit table.
  expect_false(bpnmf_output_opts()$print_target_table)
  expect_false(parse_yaml_output(list())$print_target_table)
})

test_that("denominator_may_be_affected flows through the config layers, on by default", {
  # Default TRUE: a user who doesn't already know about this concern has no
  # reason to go looking for a flag to turn it on, so the caution starts
  # showing and has to be deliberately silenced instead.
  expect_true(bpnmf_output_opts()$denominator_may_be_affected)
  expect_false(
    bpnmf_output_opts(denominator_may_be_affected = FALSE)$denominator_may_be_affected
  )
  expect_true(parse_yaml_output(list())$denominator_may_be_affected)
  expect_false(
    parse_yaml_output(
      list(denominator_may_be_affected = FALSE)
    )$denominator_may_be_affected
  )
})

test_that("the denominator caveat is a footnote on the expected columns", {
  skip_if_not_installed("gt")
  draws <- make_draws_frame()

  flagged <- bpnmf_gt_table(draws)
  # A real gt footnote anchored to the column it is about, not a loose source
  # note. Expected only, not Exp Rate: the model parameterizes the rate
  # directly, so only the count-scale counterfactual needs a denominator
  # multiplied back in, and that's the assumption being flagged.
  fn <- flagged[["_footnotes"]]
  expect_equal(fn$colname, "Expected")
  expect_true(all(fn$locname == "columns_columns"))
  expect_match(unlist(fn$footnotes), "conditional on the observed denominator",
               all = FALSE)

  quiet <- bpnmf_gt_table(draws, denominator_may_be_affected = FALSE)
  expect_equal(nrow(quiet[["_footnotes"]]), 0)

  # With a denominator, % Change belongs with the rates: it estimates a
  # change in the rate, holding the measured exposure fixed.
  spanners <- flagged[["_spanners"]]
  rate_vars <- unlist(
    spanners$vars[grepl("Rate per", unlist(spanners$spanner_label))]
  )
  expect_true("Pct Change CI" %in% rate_vars)

  # Flows through the full report, printed only alongside the table itself.
  out_dir <- withr::local_tempdir()
  msg <- utils::capture.output(
    bpnmf_report(
      draws, out_dir,
      figures = character(), html_tables = FALSE, print_tables = TRUE
    ),
    type = "message"
  )
  expect_match(paste(msg, collapse = "\n"), "Caution")
  expect_match(paste(msg, collapse = "\n"), "percent change in the rate")
})

test_that("with no denominator, rate columns are dropped and % change is a count", {
  skip_if_not_installed("gt")
  draws <- make_draws_frame()
  draws$denominator <- NULL

  # Regression: bpnmf_summary_table()/bpnmf_interval_plot() used to error
  # ("Column `denominator` not found") whenever no denominator was ever
  # configured -- they must instead treat that as denominator == 1.
  tbl <- expect_silent(bpnmf_summary_table(draws, "B"))
  expect_true(nrow(tbl) > 0)
  expect_s3_class(bpnmf_interval_plot(draws), "ggplot")

  # The exposure column would truncate to 0 and every rate would divide by
  # nothing, so the rate columns are omitted rather than reported as garbage.
  expect_false(any(
    c("Denominator-Years", "Obs Rate", "Exp Rate", "Rate Diff CI") %in%
      names(tbl)
  ))
  expect_true("Pct Change CI" %in% names(tbl))

  g <- bpnmf_gt_table(draws)
  spanners <- g[["_spanners"]]
  expect_equal(unlist(spanners$spanner_label), "Counts")
  expect_true("Pct Change CI" %in% unlist(spanners$vars))
  # Nothing is held fixed without a denominator, so there is no assumption
  # about it to caveat.
  expect_equal(nrow(g[["_footnotes"]]), 0)

  out_dir <- withr::local_tempdir()
  msg <- utils::capture.output(
    bpnmf_report(
      draws, out_dir,
      figures = character(), html_tables = FALSE, print_tables = TRUE
    ),
    type = "message"
  )
  expect_match(paste(msg, collapse = "\n"), "percent change in the raw count")
  expect_false(grepl("Caution", paste(msg, collapse = "\n")))
})

test_that("a denominator-free run renders every figure in the report", {
  # Reported case: mortality.yml with `denominator_prefix` commented out blew
  # up in bpnmf_group_comparison_plot() ("Column `denominator` not found").
  # Every rate-aware consumer has to tolerate a run with no denominator, so
  # exercise the whole figure set rather than the tables alone.
  draws <- make_draws_frame()
  draws$denominator <- NULL
  attr(draws, "has_denominator") <- FALSE

  out_dir <- withr::local_tempdir()
  expect_no_error(
    bpnmf_report(draws, out_dir, print_tables = FALSE, html_tables = FALSE)
  )
  figs <- list.files(file.path(out_dir, "figs"), recursive = TRUE)
  expect_true(all(
    c("raw_rate.png", "group_comparison.png", "interval.png") %in% figs
  ))

  # A rate is undefined without an exposure, so the rate plots fall back to
  # counts rather than plotting count * rate_multiplier under a rate axis.
  expect_equal(bpnmf_raw_rate_plot(draws)$labels$y, "Count")
  expect_equal(bpnmf_group_comparison_plot(draws)$labels$y, "Count")
})

test_that("bpnmf_draws stamps has_denominator and always emits the column", {
  # The column is always present (1 everywhere when unconfigured) so every
  # consumer sees one schema; the attribute is what carries the distinction.
  with_denom <- make_draws_frame()
  expect_true(draws_has_denominator(with_denom))

  # No flag stamped (hand-built frame): fall back to column presence.
  no_col <- with_denom
  no_col$denominator <- NULL
  attr(no_col, "has_denominator") <- NULL
  expect_false(draws_has_denominator(no_col))

  # Flag wins over column presence once stamped: a denominator-free run still
  # has a denominator column, filled with 1s.
  filled <- with_denom
  filled$denominator <- 1
  attr(filled, "has_denominator") <- FALSE
  expect_false(draws_has_denominator(filled))
})

test_that("rate scale and denominator noun are configurable, neutral by default", {
  skip_if_not_installed("gt")
  draws <- make_draws_frame()

  # The package can't know what a denominator counts, so no noun by default --
  # "person-years" was a demographic assumption baked into a general-purpose
  # package -- and the exposure column is named for what it is.
  # The noun is composed from the denominator's own label and the time unit
  # it is weighted by, so the time weighting stays visible instead of a bare
  # "Rate per 1,000" implying a rate per denominator.
  expect_equal(format_rate_label(1000), "Rate per 1,000 denominator-years")
  expect_equal(
    format_rate_label(1e5, "person"), "Rate per 100,000 person-years"
  )
  expect_equal(format_rate_label(1e5, "birth", "none"), "Rate per 100,000 births")
  expect_equal(exposure_column_name(), "Denominator-Years")
  expect_equal(exposure_column_name("person"), "Person-Years")
  expect_equal(exposure_column_name("birth", "none"), "Births")
  expect_true("Denominator-Years" %in% names(bpnmf_summary_table(draws, "B")))
  expect_true(
    "Person-Years" %in%
      names(bpnmf_summary_table(draws, "B", denominator_label = "person"))
  )

  spanner <- function(g) unlist(g[["_spanners"]]$spanner_label)
  expect_true("Rate per 1,000 denominator-years" %in% spanner(bpnmf_gt_table(draws)))
  expect_true(
    "Rate per 100,000 person-years" %in%
      spanner(bpnmf_gt_table(draws, rate_normalizer = 1e5,
                             denominator_label = "person"))
  )
  expect_equal(
    bpnmf_interval_plot(
      draws,
      estimand = "diff", rate_normalizer = 1e5, denominator_label = "person"
    )$labels$x,
    "Rate Difference (per 100,000 person-years)"
  )
  # The raw-rate plots don't time-weight at all (sum(outcome)/sum(denominator)),
  # so their label must not claim a per-time rate.
  expect_equal(
    bpnmf_raw_rate_plot(
      draws,
      rate_multiplier = 1e5, denominator_label = "birth"
    )$labels$y,
    "Rate per 100,000 births"
  )

  # The time unit is part of the quantity: weighting by months rather than
  # years scales the exposure, and hence the rates, by 12.
  yearly <- bpnmf_summary_table(draws, "B")
  monthly <- bpnmf_summary_table(draws, "B", denominator_time_unit = "month")
  expect_equal(monthly$`Denominator-Months`, yearly$`Denominator-Years` * 12L)
  # Obs Rate is rounded to 2dp in the table, so compare on that scale.
  expect_equal(monthly$`Obs Rate`, round(yearly$`Obs Rate` / 12, 2))
})

test_that("rate_normalizer and denominator_label reach the report from config", {
  # Regression: rate_normalizer was a hardcoded default on six function
  # signatures and was never threaded through bpnmf_report(), so a bpnmf_run()
  # pipeline was stuck at 1000 with no way to override it.
  expect_equal(bpnmf_output_opts()$rate_normalizer, 1000)
  expect_equal(bpnmf_output_opts()$denominator_label, "denominator")
  expect_equal(bpnmf_output_opts()$denominator_time_unit, "year")
  expect_error(bpnmf_output_opts(denominator_time_unit = "fortnight"))
  expect_equal(bpnmf_output_opts(rate_normalizer = 1e5)$rate_normalizer, 1e5)
  expect_equal(parse_yaml_output(list())$rate_normalizer, 1000)
  expect_equal(
    parse_yaml_output(list(rate_normalizer = 1e5))$rate_normalizer, 1e5
  )
  expect_equal(
    parse_yaml_output(list(denominator_label = "births"))$denominator_label,
    "births"
  )
  expect_error(bpnmf_output_opts(rate_normalizer = 0))
  expect_error(bpnmf_output_opts(denominator_label = 42))

  skip_if_not_installed("gt")
  out_dir <- withr::local_tempdir()
  bpnmf_report(
    make_draws_frame(), out_dir,
    figures = character(), print_tables = FALSE,
    rate_normalizer = 1e5, denominator_label = "person"
  )
  html <- paste(
    readLines(file.path(out_dir, "figs", "summary_table.html"), warn = FALSE),
    collapse = "\n"
  )
  expect_match(html, "Rate per 100,000 person-years")
})
