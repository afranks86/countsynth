# Figure/table orchestration with artifact-layout parity: figures and tables
# under <output_dir>/figs/ with the same filenames the Python package writes
# (fit_<unit>.png, gap_<unit>.png, raw_rate.png, interval.png,
# group_comparison.png, ppc/ppc_*.png + ppc_pvalues.csv,
# summary_table_by_unit.csv, expected_vs_observed.csv,
# post_treatment_summary.csv). Tables ALWAYS write; only figures are gated by
# the `figures` selection. The target unit's headline table is not written
# separately: it is the `Unit == target_unit` subset of
# summary_table_by_unit.csv.

# write.csv serializes doubles at full 15-17 significant digits, which makes
# the tables unreadable and implies precision the posterior does not have.
# Round to 6 significant digits on the way out; integer/character columns
# pass through untouched.
write_table_csv <- function(df, path, digits = 6) {
  num <- vapply(df, function(x) is.double(x) && !inherits(x, "Date"), logical(1))
  df[num] <- lapply(df[num], signif, digits = digits)
  utils::write.csv(df, path, row.names = FALSE)
}

# The `Imputed` flag only means something once a cell has actually been
# imputed; an all-FALSE column is a wasted column in an already-wide table.
# It stays in the CSV either way, where a stable schema matters more.
drop_unused_imputed <- function(tbl) {
  if ("Imputed" %in% names(tbl) && !any(tbl$Imputed)) tbl$Imputed <- NULL
  tbl
}

save_plot <- function(plot, path, width = 10, height = 6) {
  dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggplot2::ggsave(
    path, plot,
    width = width, height = height, dpi = 150, device = dev
  )
}

# PPC plots facet by (unit, group) or, for unit_corr, by group alone (see
# ppc_histogram()'s ncol logic in plot-ppc.R). A fixed canvas size squashes
# every row once unit/group counts grow past what fits at that size, so scale
# the canvas to the facet grid instead.
ppc_plot_dims <- function(n_facets, ncol,
                          per_facet_width = 3.6, per_facet_height = 1.9,
                          min_width = 9, min_height = 5) {
  n_facets <- max(n_facets, 1L)
  facet_nrow <- ceiling(n_facets / ncol)
  list(
    width = max(min_width, ncol * per_facet_width),
    height = max(min_height, facet_nrow * per_facet_height)
  )
}

#' Generate figures and tables from posterior draws
#'
#' Port of `reports.generate_reports`. Per-unit figures are nested under
#' `figs/<group>/` when more than one group is reported; cross-group figures
#' and all tables go directly under `figs/`.
#'
#' @param draws A `bpnmf_draws` frame.
#' @param output_dir Run directory; artifacts go to `<output_dir>/figs/`.
#' @param target_unit Headline unit (auto-detected when `NULL`).
#' @param groups Groups to render per-unit figures for (`NULL` = all).
#' @param figures Character vector of figure names to render (see
#'   [bpnmf_output_opts()]); `NULL` renders everything, `character()` renders
#'   tables only.
#' @param aggregate_units List of [bpnmf_aggregate_unit()] specs.
#' @param ppc_draws Alternative draws frame for the PPC suite only (cut mode:
#'   the stage-1 PPC frame).
#' @param ppc_units,ppc_exclude_units,ppc_acf_lags,ppc_unit_corr_max_time
#'   PPC options (see [bpnmf_output_opts()]).
#' @param fit_gap_per_unit Also render fit/gap for every treated unit
#'   (restricted to the `"total"` group).
#' @param interval_aggregates Include the `aggregate_units` in `interval.png`,
#'   split off into their own band above the individual units. Set `FALSE` to
#'   plot only the real units. No effect when `aggregate_units` is `NULL`.
#' @param print_tables Print the by-unit summary table to the terminal.
#' @param print_target_table Also print the target unit's own table. Its rows
#'   are the target unit's slice of the by-unit table, so this is off by
#'   default when `target_unit` is not set explicitly.
#' @param html_tables Also write `summary_table.html` and
#'   `summary_table_by_unit.html` via [bpnmf_gt_table()]. Needs the `gt`
#'   package; warns and skips when it is missing.
#' @param fit Optional `bpnmf_fit` / `bpnmf_cut_fit` the draws came from.
#'   Required for the `"te_regression"` figures, which read the
#'   treatment-effect design and coefficient draws rather than the draws
#'   frame; without it that figure is skipped.
#' @inheritParams bpnmf_output_opts
#' @return Invisible list with `summary`, `per_unit`, `detail`,
#'   `target_unit`, `figs_dir`, `treated_units`.
#' @export
bpnmf_report <- function(draws, output_dir, target_unit = NULL, groups = NULL,
                         figures = NULL, aggregate_units = NULL,
                         ppc_draws = NULL, ppc_units = NULL,
                         ppc_exclude_units = NULL, ppc_acf_lags = NULL,
                         ppc_unit_corr_max_time = NULL,
                         fit_gap_per_unit = FALSE,
                         interval_aggregates = TRUE, print_tables = TRUE,
                         print_target_table = FALSE, html_tables = TRUE,
                         fit = NULL, rate_normalizer = 1000,
                         denominator_label = "denominator",
                         denominator_time_unit = "year",
                         denominator_may_be_affected = TRUE) {
  selected <- figures %||% FIGURE_NAMES
  unknown <- setdiff(selected, FIGURE_NAMES)
  if (length(unknown) > 0) {
    cli::cli_abort(
      "figures contains unknown name{?s} {.val {unknown}};
       valid names are {.val {sort(FIGURE_NAMES)}}."
    )
  }
  figs_dir <- file.path(output_dir, "figs")
  dir.create(figs_dir, recursive = TRUE, showWarnings = FALSE)

  reporting <- if (!is.null(aggregate_units)) {
    add_aggregate_units(draws, aggregate_units)
  } else {
    draws
  }

  # Detect on `reporting`, not `draws`: an aggregate unit exists only in the
  # former, and when one is configured it is the headline unit by default.
  target_unit <- target_unit %||% auto_detect_target(reporting)
  if (is.null(target_unit)) {
    cli::cli_abort("No treated units in draws and no target_unit specified.")
  }

  all_groups <- unique(reporting$group)
  report_groups <- groups %||% all_groups
  bad <- setdiff(report_groups, all_groups)
  if (length(bad) > 0) {
    cli::cli_abort(
      "groups {.val {bad}} not present in draws (have {.val {all_groups}})."
    )
  }

  # From `reporting` so the aggregate unit appears in the by-unit table
  # alongside post_treatment_summary.csv, which is built from the same frame.
  treated_units <- identify_treated_units(reporting)
  fit_gap_units <- function(grp) {
    if (fit_gap_per_unit && grp == "total") {
      unique(c(target_unit, treated_units))
    } else {
      target_unit
    }
  }

  quantiles <- compute_quantiles(reporting)

  # Per-group unit-level figures.
  for (grp in report_groups) {
    grp_dir <- if (length(report_groups) > 1) {
      file.path(figs_dir, grp)
    } else {
      figs_dir
    }
    dir.create(grp_dir, recursive = TRUE, showWarnings = FALSE)
    for (unit in fit_gap_units(grp)) {
      slug <- unit_slug(unit)
      if ("unit_fit" %in% selected) {
        save_plot(
          bpnmf_unit_fit_plot(quantiles, unit, grp),
          file.path(grp_dir, sprintf("fit_%s.png", slug))
        )
      }
      if ("unit_gap" %in% selected) {
        save_plot(
          bpnmf_unit_gap_plot(quantiles, unit, grp),
          file.path(grp_dir, sprintf("gap_%s.png", slug))
        )
      }
    }
    if ("raw_rate" %in% selected) {
      save_plot(
        bpnmf_raw_rate_plot(
          draws,
          group = grp, separate_unit = target_unit,
          rate_multiplier = rate_normalizer,
          denominator_label = denominator_label
        ),
        file.path(grp_dir, "raw_rate.png")
      )
    }
  }

  # Cross-group figures.
  if ("interval" %in% selected) {
    # An aggregate unit pools the same draws as the units it covers, so it is
    # shown in its own band rather than ranked among them; `interval_aggregates
    # = FALSE` drops it from the figure entirely.
    save_plot(
      bpnmf_interval_plot(
        if (interval_aggregates) reporting else draws,
        estimand = "ratio", method = "mu",
        rate_normalizer = rate_normalizer,
        denominator_label = denominator_label,
        denominator_time_unit = denominator_time_unit,
        separate_units = if (interval_aggregates) NULL else character()
      ),
      file.path(figs_dir, "interval.png"),
      width = 10, height = 8
    )
  }
  if ("group_comparison" %in% selected) {
    save_plot(
      bpnmf_group_comparison_plot(
        draws,
        rate_multiplier = rate_normalizer,
        denominator_label = denominator_label
      ),
      file.path(figs_dir, "group_comparison.png"),
      width = 11, height = 7
    )
  }
  if ("te_regression" %in% selected && !is.null(fit) &&
    !is.null(fit$te_design)) {
    te_dir <- file.path(figs_dir, "te")
    dir.create(te_dir, recursive = TRUE, showWarnings = FALSE)
    te_figs <- te_report_figures(fit)
    for (nm in names(te_figs)) {
      save_plot(te_figs[[nm]], file.path(te_dir, paste0(nm, ".png")),
        width = 10, height = 6
      )
    }
    write_table_csv(
      bpnmf_te_coef_table(fit), file.path(te_dir, "te_coefficients.csv")
    )
  }
  if ("ppc" %in% selected) {
    ppc_dir <- file.path(figs_dir, "ppc")
    dir.create(ppc_dir, recursive = TRUE, showWarnings = FALSE)
    ppc_source <- ppc_draws %||% reporting
    if (!is.null(ppc_draws) && !is.null(aggregate_units)) {
      ppc_source <- add_aggregate_units(ppc_draws, aggregate_units)
    }
    ppc <- bpnmf_ppc_plots(
      ppc_source,
      acf_lags = ppc_acf_lags %||% 1,
      max_treat_date = ppc_unit_corr_max_time,
      ppc_units = ppc_units,
      ppc_exclude_units = ppc_exclude_units
    )
    for (nm in names(ppc$plots)) {
      is_unit_corr <- identical(nm, "ppc_unit_corr")
      check_type <- if (identical(nm, "ppc_abs_residual")) "abs" else sub("^ppc_", "", nm)
      n_facets <- sum(ppc$pvals$check_type == check_type)
      dims <- ppc_plot_dims(n_facets, ncol = if (is_unit_corr) 2L else 3L)
      save_plot(
        ppc$plots[[nm]], file.path(ppc_dir, paste0(nm, ".png")),
        width = dims$width, height = dims$height
      )
    }
    write_table_csv(
      ppc$pvals, file.path(ppc_dir, "ppc_pvalues.csv")
    )
  }

  # Tables: always written, never gated by `figures`.
  # Returned (and printed) but not written: its rows are the target unit's
  # slice of summary_table_by_unit.csv below.
  summary_tbl <- bpnmf_summary_table(
    reporting, target_unit,
    rate_normalizer = rate_normalizer,
    denominator_label = denominator_label,
    denominator_time_unit = denominator_time_unit
  )
  by_unit_tbl <- bpnmf_summary_table_by_unit(
    reporting, treated_units,
    rate_normalizer = rate_normalizer,
    denominator_label = denominator_label,
    denominator_time_unit = denominator_time_unit
  )
  write_table_csv(
    by_unit_tbl, file.path(figs_dir, "summary_table_by_unit.csv")
  )
  detail <- bpnmf_expected_vs_observed(reporting, target_unit)
  write_table_csv(
    detail, file.path(figs_dir, "expected_vs_observed.csv")
  )
  per_unit <- bpnmf_post_treatment_summary(reporting)
  write_table_csv(
    per_unit, file.path(figs_dir, "post_treatment_summary.csv")
  )
  if (html_tables) {
    write_gt_tables(
      reporting, target_unit, figs_dir,
      rate_normalizer = rate_normalizer,
      denominator_label = denominator_label,
      denominator_time_unit = denominator_time_unit,
      denominator_may_be_affected = denominator_may_be_affected
    )
  }

  if (print_tables && nrow(by_unit_tbl) > 0) {
    # One formatted table, not two: the by-unit table is the target table plus
    # every other unit, so printing both repeats the headline rows. The raw
    # post_treatment_summary frame is 14 numeric columns wide and unreadable
    # in a terminal -- it stays a CSV, for joining and plotting.
    if (print_target_table && nrow(summary_tbl) > 0) {
      cli::cli_h1("{target_unit} \u2014 observed vs expected")
      print(as.data.frame(drop_unused_imputed(summary_tbl)))
    }
    cli::cli_h1("Post-treatment effect by unit")
    print(as.data.frame(drop_unused_imputed(by_unit_tbl)), row.names = FALSE)
    # The terminal table can't carry the gt version's spanners or its
    # column-anchored footnote (see bpnmf_gt_table()), so say in words what
    # "Pct Change" estimates -- and, with a denominator, what Expected
    # assumes. Without a denominator nothing is being held fixed, so there is
    # no such assumption to caveat.
    if (draws_has_denominator(reporting)) {
      cli::cli_alert_info(paste(
        "Pct Change above is the estimated percent change in the rate",
        "(a configured denominator was held fixed on both sides)."
      ))
      if (isTRUE(denominator_may_be_affected)) {
        cli::cli_alert_warning(paste(
          "Caution: Expected is a counterfactual count, conditional on the",
          "observed denominator -- it multiplies the estimated counterfactual",
          "rate by the denominator as actually observed, which treatment may",
          "itself have changed. The rate columns do not share this",
          "assumption (the denominator cancels out of them). Set",
          "{.field output.denominator_may_be_affected} to FALSE to silence",
          "this."
        ))
      }
    } else {
      cli::cli_alert_info(paste(
        "Pct Change above is the estimated percent change in the raw count",
        "(no denominator was configured, so it can't be attributed to a",
        "rate change specifically); rate columns are omitted for the same",
        "reason."
      ))
    }
  }

  invisible(list(
    summary = summary_tbl, per_unit = per_unit, detail = detail,
    target_unit = target_unit, figs_dir = figs_dir,
    treated_units = treated_units
  ))
}
