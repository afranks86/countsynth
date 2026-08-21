# Figure/table orchestration with artifact-layout parity: figures and tables
# under <output_dir>/figs/ with the same filenames the Python package writes
# (fit_<unit>.png, gap_<unit>.png, raw_rate.png, interval.png,
# group_comparison.png, ppc/ppc_*.png + ppc_pvalues.csv, summary_table.csv,
# summary_table_by_unit.csv, expected_vs_observed.csv,
# post_treatment_summary.csv). Tables ALWAYS write; only figures are gated by
# the `figures` selection.

save_plot <- function(plot, path, width = 10, height = 6) {
  dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggplot2::ggsave(
    path, plot,
    width = width, height = height, dpi = 150, device = dev
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
#' @param print_tables Print summary tables to the terminal.
#' @return Invisible list with `summary`, `per_unit`, `detail`,
#'   `target_unit`, `figs_dir`, `treated_units`.
#' @export
bpnmf_report <- function(draws, output_dir, target_unit = NULL, groups = NULL,
                         figures = NULL, aggregate_units = NULL,
                         ppc_draws = NULL, ppc_units = NULL,
                         ppc_exclude_units = NULL, ppc_acf_lags = NULL,
                         ppc_unit_corr_max_time = NULL,
                         fit_gap_per_unit = FALSE, print_tables = TRUE) {
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

  target_unit <- target_unit %||% auto_detect_target(draws)
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

  treated_units <- identify_treated_units(draws)
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
        bpnmf_raw_rate_plot(draws, group = grp, separate_unit = target_unit),
        file.path(grp_dir, "raw_rate.png")
      )
    }
  }

  # Cross-group figures.
  if ("interval" %in% selected) {
    save_plot(
      bpnmf_interval_plot(draws, estimand = "ratio", method = "mu"),
      file.path(figs_dir, "interval.png"),
      width = 10, height = 8
    )
  }
  if ("group_comparison" %in% selected) {
    save_plot(
      bpnmf_group_comparison_plot(draws),
      file.path(figs_dir, "group_comparison.png"),
      width = 11, height = 7
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
      acf_lags = ppc_acf_lags %||% 6,
      max_treat_date = ppc_unit_corr_max_time,
      ppc_units = ppc_units,
      ppc_exclude_units = ppc_exclude_units
    )
    for (nm in names(ppc$plots)) {
      save_plot(ppc$plots[[nm]], file.path(ppc_dir, paste0(nm, ".png")),
        width = 11, height = 8
      )
    }
    utils::write.csv(
      ppc$pvals, file.path(ppc_dir, "ppc_pvalues.csv"),
      row.names = FALSE
    )
  }

  # Tables: always written, never gated by `figures`.
  summary_tbl <- bpnmf_summary_table(reporting, target_unit)
  utils::write.csv(
    summary_tbl, file.path(figs_dir, "summary_table.csv"),
    row.names = FALSE
  )
  by_unit_tbl <- dplyr::bind_rows(lapply(treated_units, function(tu) {
    tbl <- bpnmf_summary_table(reporting, tu)
    if (nrow(tbl) == 0) {
      return(tbl)
    }
    tibble::add_column(tbl, Unit = tu, .before = 1)
  }))
  utils::write.csv(
    by_unit_tbl, file.path(figs_dir, "summary_table_by_unit.csv"),
    row.names = FALSE
  )
  detail <- bpnmf_expected_vs_observed(reporting, target_unit)
  utils::write.csv(
    detail, file.path(figs_dir, "expected_vs_observed.csv"),
    row.names = FALSE
  )
  per_unit <- bpnmf_post_treatment_summary(reporting)
  utils::write.csv(
    per_unit, file.path(figs_dir, "post_treatment_summary.csv"),
    row.names = FALSE
  )

  if (print_tables && nrow(summary_tbl) > 0) {
    cli::cli_h1("{target_unit} \u2014 Observed vs Expected")
    print(as.data.frame(summary_tbl))
    if (nrow(per_unit) > 0) {
      cli::cli_h1("Post-treatment totals by unit (ranked by % excess)")
      print(as.data.frame(per_unit), digits = 4)
    }
  }

  invisible(list(
    summary = summary_tbl, per_unit = per_unit, detail = detail,
    target_unit = target_unit, figs_dir = figs_dir,
    treated_units = treated_units
  ))
}
