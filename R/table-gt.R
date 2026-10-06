# Publication-ready HTML rendering of the headline summary table. The tibble
# from countsynth_summary_table() is already display-shaped -- counts rounded, CIs
# pre-formatted as strings -- so this module is presentation only: column
# labels, spanners, grouping, and the footnotes that explain the `*` and the
# dagger. gt lives in Suggests, so every entry point checks for it first.

# Column labels carry their own units, so the headers can stay short.
GT_LABELS <- list(
  `Person-Years` = "Person-years",
  Observed = "Observed",
  Expected = "Expected",
  `Diff (95% CI)` = "Difference (95% CI)",
  `Obs Rate` = "Observed",
  `Exp Rate` = "Expected",
  `Rate Diff CI` = "Difference (95% CI)",
  `Pct Change CI` = "% change (95% CI)"
)

#' Render a summary table as a `gt` HTML table
#'
#' Formats [countsynth_summary_table()] (or [countsynth_summary_table_by_unit()], when
#' `by_unit` is `TRUE`) for display: counts and rates split under their own
#' spanners, one row group per unit, and footnotes for the significance star
#' and the imputation dagger. Requires the `gt` package.
#'
#' @param draws A `countsynth_draws` frame.
#' @param target_unit Unit to summarize. Ignored when `by_unit` is `TRUE`;
#'   auto-detected when `NULL` (an aggregate unit wins -- see
#'   [auto_detect_target()]).
#' @param by_unit Show every treated unit, grouped by unit, instead of one.
#' @param rate_normalizer Rates are per this many person-years.
#' @param title,subtitle Header text. `NULL` builds a default from the unit.
#' @return A `gt_tbl`. Print it to view, or save with [gt::gtsave()].
#' @export
#' @examples
#' \dontrun{
#' draws <- countsynth_draws(fit)
#' countsynth_gt_table(draws)                  # headline unit, one row per group
#' countsynth_gt_table(draws, by_unit = TRUE)  # every treated unit
#' gt::gtsave(countsynth_gt_table(draws), "summary.html")
#' }
countsynth_gt_table <- function(draws, target_unit = NULL, by_unit = FALSE,
                           rate_normalizer = 1000,
                           title = NULL, subtitle = NULL) {
  rlang::check_installed("gt", reason = "to render HTML summary tables")
  checkmate::assert_flag(by_unit)

  if (by_unit) {
    tbl <- countsynth_summary_table_by_unit(draws, rate_normalizer = rate_normalizer)
    title <- title %||% "Post-treatment effect by unit"
  } else {
    target_unit <- target_unit %||% auto_detect_target(draws)
    tbl <- countsynth_summary_table(draws, target_unit, rate_normalizer = rate_normalizer)
    title <- title %||% sprintf("%s — observed vs expected", target_unit)
  }
  if (nrow(tbl) == 0) {
    cli::cli_abort("No post-treatment rows to tabulate.")
  }

  # `Imputed` drives a footnote; it would be a redundant column next to the
  # dagger already in the Group label.
  any_imputed <- any(tbl$Imputed)
  tbl$Imputed <- NULL
  rate_label <- sprintf(
    "Rate per %s person-years", format(rate_normalizer, big.mark = ",")
  )

  g <- gt::gt(
    tbl,
    groupname_col = if (by_unit) "Unit" else NULL,
    rowname_col = "Group"
  )
  g <- gt::tab_header(g, title = title, subtitle = subtitle)
  g <- gt::cols_label(g, .list = GT_LABELS)
  g <- gt::tab_spanner(g, label = "Counts", columns = c(
    "Observed", "Expected", "Diff (95% CI)"
  ))
  g <- gt::tab_spanner(g, label = rate_label, columns = c(
    "Obs Rate", "Exp Rate", "Rate Diff CI", "Pct Change CI"
  ))
  g <- gt::fmt_number(g, columns = c("Person-Years", "Observed", "Expected"),
                      decimals = 0, use_seps = TRUE)
  g <- gt::cols_align(g, align = "right", columns = -1)
  # Source notes rather than footnotes: both marks live in the row label, so
  # anchoring a numbered footnote to some arbitrary column would misdirect.
  g <- gt::tab_source_note(
    g, "* two-sided posterior p < 0.05 for the treated/untreated contrast."
  )
  if (any_imputed) {
    g <- gt::tab_source_note(g, paste(
      "\u2020 post-treatment window includes a cell whose suppressed count",
      "was imputed when building the aggregate unit."
    ))
  }
  gt::tab_options(g, table.font.size = gt::px(13), data_row.padding = gt::px(4))
}

# Write the report's HTML tables. gt is optional, so a missing install is a
# warning that names the fix, not a failed run -- the CSVs already landed.
write_gt_tables <- function(draws, target_unit, figs_dir) {
  if (!requireNamespace("gt", quietly = TRUE)) {
    cli::cli_warn(c(
      "Skipping HTML tables: the {.pkg gt} package is not installed.",
      i = 'Install it with {.code install.packages("gt")}.'
    ))
    return(invisible(character()))
  }
  written <- character()
  specs <- list(
    list(file = "summary_table.html", by_unit = FALSE),
    list(file = "summary_table_by_unit.html", by_unit = TRUE)
  )
  for (spec in specs) {
    path <- file.path(figs_dir, spec$file)
    ok <- tryCatch(
      {
        gt::gtsave(
          countsynth_gt_table(draws, target_unit, by_unit = spec$by_unit),
          path
        )
        TRUE
      },
      error = function(e) {
        cli::cli_warn("Could not write {.file {spec$file}}: {conditionMessage(e)}")
        FALSE
      }
    )
    if (ok) written <- c(written, path)
  }
  invisible(written)
}
