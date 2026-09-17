# Publication-ready HTML rendering of the headline summary table. The tibble
# from bpnmf_summary_table() is already display-shaped -- counts rounded, CIs
# pre-formatted as strings -- so this module is presentation only: column
# labels, spanners, grouping, and the footnotes that explain the `*` and the
# dagger. gt lives in Suggests, so every entry point checks for it first.

# Column labels carry their own units, so the headers can stay short.
GT_LABELS <- list(
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
#' Formats [bpnmf_summary_table()] (or [bpnmf_summary_table_by_unit()], when
#' `by_unit` is `TRUE`) for display: counts and rates split under their own
#' spanners, one row group per unit, and footnotes for the significance star
#' and the imputation dagger. Requires the `gt` package.
#'
#' @param draws A `bpnmf_draws` frame.
#' @param target_unit Unit to summarize. Ignored when `by_unit` is `TRUE`;
#'   auto-detected when `NULL` (an aggregate unit wins -- see
#'   [auto_detect_target()]).
#' @param by_unit Show every treated unit, grouped by unit, instead of one.
#' @param title,subtitle Header text. `NULL` builds a default from the unit.
#' @inheritParams bpnmf_output_opts
#' @return A `gt_tbl`. Print it to view, or save with [gt::gtsave()].
#' @export
#' @examples
#' \dontrun{
#' draws <- bpnmf_draws(fit)
#' bpnmf_gt_table(draws)                  # headline unit, one row per group
#' bpnmf_gt_table(draws, by_unit = TRUE)  # every treated unit
#' gt::gtsave(bpnmf_gt_table(draws), "summary.html")
#' }
bpnmf_gt_table <- function(draws, target_unit = NULL, by_unit = FALSE,
                           rate_normalizer = 1000,
                           title = NULL, subtitle = NULL,
                           denominator_label = "denominator",
                           denominator_time_unit = "year",
                           denominator_may_be_affected = TRUE) {
  rlang::check_installed("gt", reason = "to render HTML summary tables")
  checkmate::assert_flag(by_unit)

  if (by_unit) {
    tbl <- bpnmf_summary_table_by_unit(
      draws,
      rate_normalizer = rate_normalizer,
      denominator_label = denominator_label,
      denominator_time_unit = denominator_time_unit
    )
    title <- title %||% "Post-treatment effect by unit"
  } else {
    target_unit <- target_unit %||% auto_detect_target(draws)
    tbl <- bpnmf_summary_table(
      draws, target_unit,
      rate_normalizer = rate_normalizer,
      denominator_label = denominator_label,
      denominator_time_unit = denominator_time_unit
    )
    title <- title %||% sprintf("%s — observed vs expected", target_unit)
  }
  if (nrow(tbl) == 0) {
    cli::cli_abort("No post-treatment rows to tabulate.")
  }

  # `Imputed` drives a footnote; it would be a redundant column next to the
  # dagger already in the Group label.
  any_imputed <- any(tbl$Imputed)
  tbl$Imputed <- NULL
  rate_label <- format_rate_label(
    rate_normalizer, denominator_label, denominator_time_unit
  )
  exposure_col <- exposure_column_name(
    denominator_label, denominator_time_unit
  )
  # A configured denominator makes mu's log(denominator) offset identical on
  # the treated and untreated side, so "Pct Change" estimates the percent
  # change in the RATE (holding that measured exposure fixed) and belongs
  # with the rate columns. Without one there is no measured exposure at all:
  # bpnmf_summary_table() drops the rate columns entirely, and the same
  # number is a change in the raw COUNT, so it belongs under Counts.
  has_denom <- draws_has_denominator(draws)
  count_cols <- c("Observed", "Expected", "Diff (95% CI)")
  rate_cols <- c("Obs Rate", "Exp Rate", "Rate Diff CI")

  g <- gt::gt(
    tbl,
    groupname_col = if (by_unit) "Unit" else NULL,
    rowname_col = "Group"
  )
  g <- gt::tab_header(g, title = title, subtitle = subtitle)
  g <- gt::cols_label(g, .list = GT_LABELS[names(GT_LABELS) %in% names(tbl)])
  if (has_denom) {
    g <- gt::tab_spanner(g, label = "Counts", columns = count_cols)
    g <- gt::tab_spanner(
      g,
      label = rate_label, columns = c(rate_cols, "Pct Change CI")
    )
  } else {
    g <- gt::tab_spanner(
      g,
      label = "Counts", columns = c(count_cols, "Pct Change CI")
    )
  }
  g <- gt::fmt_number(
    g,
    columns = intersect(c(exposure_col, "Observed", "Expected"), names(tbl)),
    decimals = 0, use_seps = TRUE
  )
  g <- gt::cols_align(g, align = "right", columns = -1)
  # Source notes rather than footnotes for these two: both marks live in the
  # row label, so anchoring a numbered footnote to some arbitrary column
  # would misdirect.
  g <- gt::tab_source_note(
    g, "* two-sided posterior p < 0.05 for the treated/untreated contrast."
  )
  if (any_imputed) {
    g <- gt::tab_source_note(g, paste(
      "\u2020 post-treatment window includes a cell whose suppressed count",
      "was imputed when building the aggregate unit."
    ))
  }
  # This one does have a column to point at, so it is a real footnote rather
  # than a source note. It belongs on Expected alone, not on Exp Rate: the
  # model parameterizes the rate directly, so the counterfactual RATE needs
  # no assumption about what the denominator would have been (it cancels --
  # Exp Rate is untreated_count / denom_val). Turning that rate into a
  # counterfactual COUNT is what requires multiplying by a denominator, and
  # the one used is the observed, possibly treatment-affected one. Only
  # meaningful with a denominator at all -- without one nothing is being held
  # fixed. gt's numbered marker avoids colliding with the literal * and
  # dagger already carried in the row labels.
  if (isTRUE(denominator_may_be_affected) && has_denom) {
    g <- gt::tab_footnote(
      g,
      footnote = paste(
        "Counterfactual count, conditional on the observed denominator: it",
        "multiplies the estimated counterfactual rate by the denominator as",
        "actually observed, which treatment may itself have changed (e.g.",
        "births as the denominator for an infant mortality rate, when the",
        "exposure could also change the number of births). The rate columns",
        "do not share this assumption -- the denominator cancels out of",
        "them. If the effect on the count is what you are after, fit the",
        "model with no denominator: the counterfactual count is then",
        "extrapolated in its own right rather than conditioned on an",
        "observed denominator treatment may have moved. Set",
        "output.denominator_may_be_affected: false to silence this."
      ),
      locations = gt::cells_column_labels(columns = "Expected")
    )
  }
  gt::tab_options(g, table.font.size = gt::px(13), data_row.padding = gt::px(4))
}

# Write the report's HTML tables. gt is optional, so a missing install is a
# warning that names the fix, not a failed run -- the CSVs already landed.
write_gt_tables <- function(draws, target_unit, figs_dir,
                            rate_normalizer = 1000,
                            denominator_label = "denominator",
                            denominator_time_unit = "year",
                            denominator_may_be_affected = TRUE) {
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
          bpnmf_gt_table(
            draws, target_unit, by_unit = spec$by_unit,
            rate_normalizer = rate_normalizer,
            denominator_label = denominator_label,
            denominator_time_unit = denominator_time_unit,
            denominator_may_be_affected = denominator_may_be_affected
          ),
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
