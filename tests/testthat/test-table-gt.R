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
