# Figure and table orchestration. Output is split three ways by what the
# reader wants from it, under <output_dir>/:
#
#   figures/            what happened -- fit_<unit>.png, gap_<unit>.png,
#                       raw_rate.png, interval.png, group_comparison.png,
#                       plus <group>/ per-group and te/ effect figures
#   tables/             the numbers behind it -- summary_table*.html,
#                       summary_table_by_unit.csv, expected_vs_observed.csv,
#                       post_treatment_summary.csv
#   model_diagnostics/  whether to believe the fit -- ppc_*.png,
#                       ppc_pvalues.csv, and rank_*.{png,csv} under shrinkage
#
# Tables ALWAYS write; only figures are gated by the `figures` selection. The
# target unit's headline table is not written separately: it is the
# `Unit == target_unit` subset of summary_table_by_unit.csv.

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

# Total physical memory, or NA when we cannot tell. It only scales a warning,
# so not knowing it just means reporting the projected size on its own.
system_memory_bytes <- function() {
  val <- tryCatch(
    {
      if (identical(Sys.info()[["sysname"]], "Darwin")) {
        as.numeric(system2("sysctl", c("-n", "hw.memsize"), stdout = TRUE))
      } else if (file.exists("/proc/meminfo")) {
        line <- grep(
          "^MemTotal:", readLines("/proc/meminfo", warn = FALSE),
          value = TRUE
        )
        as.numeric(sub("^MemTotal:[^0-9]*([0-9]+).*$", "\\1", line[1])) * 1024
      } else {
        NA_real_
      }
    },
    error = function(e) NA_real_,
    warning = function(w) NA_real_
  )
  if (length(val) != 1L || !is.finite(val)) NA_real_ else val
}

fmt_gb <- function(x) sprintf("%.1f GB", x / 2^30)
fmt_rows <- function(x) {
  if (x >= 1e6) sprintf("%.1fM", x / 1e6) else format(round(x), big.mark = ",")
}

# Everything the report will hold is known from the cell counts before any of
# it is built, so say so at second zero. The failure this exists to prevent is
# the worst one the package has: a report that reduces over tens of millions
# of rows, runs for hours, gets killed by the OS, and leaves nothing behind to
# say what it was doing.
report_size_guard <- function(draws, aggregate_units, aggregate_groups,
                              budget = NULL) {
  n_rows <- nrow(draws)
  if (n_rows == 0 || !".draw" %in% names(draws)) {
    return(invisible(NULL))
  }
  n_groups <- length(attr(draws, "groups") %||% unique(draws$group))
  n_units <- length(attr(draws, "units") %||% unique(draws$unit))
  n_times <- length(attr(draws, "times") %||% unique(draws$time))
  n_draws <- max(draws$.draw)
  n_unit_specs <- length(aggregate_units %||% list())
  n_group_specs <- length(aggregate_groups %||% list())
  # One synthetic unit per unit spec covers every (draw, time, group); one
  # synthetic group per group spec covers every (draw, time, unit), the
  # pooled units included, which is why the group axis is the expensive one.
  extra <- n_unit_specs * n_draws * n_times * n_groups +
    n_group_specs * n_draws * n_times * (n_units + n_unit_specs)
  projected <- n_rows + extra
  sample_n <- min(1000L, n_rows)
  per_row <- as.numeric(
    utils::object.size(utils::head(draws, sample_n))
  ) / sample_n
  projected_bytes <- projected * per_row
  # The draws frame stays live while the reporting frame is built beside it,
  # so the working set is roughly twice what the rows themselves weigh.
  peak <- 2 * projected_bytes
  ram <- system_memory_bytes()
  budget <- budget %||% if (is.na(ram)) Inf else 0.5 * ram
  if (peak <= budget) {
    return(invisible(projected))
  }
  thin <- ceiling(peak / budget)
  size_msg <- if (is.na(ram)) {
    sprintf(
      "Projected reporting frame: %s rows (~%s), peaking near %s.",
      fmt_rows(projected), fmt_gb(projected_bytes), fmt_gb(peak)
    )
  } else {
    sprintf(
      paste(
        "Projected reporting frame: %s rows (~%s), peaking near %s",
        "against %s of memory."
      ),
      fmt_rows(projected), fmt_gb(projected_bytes), fmt_gb(peak), fmt_gb(ram)
    )
  }
  made_of <- sprintf(
    "%s draw rows plus %s aggregate rows.", fmt_rows(n_rows), fmt_rows(extra)
  )
  thin_msg <- sprintf(
    "Thin to 1 draw in %d for ~%s rows, or drop the aggregate groups.",
    thin, fmt_rows(projected / thin)
  )
  cli::cli_warn(c(
    size_msg,
    i = made_of,
    i = thin_msg,
    i = "Thinning is set with {.field mcmc.thinning} and costs little when
         the chains are autocorrelated; aggregate groups with
         {.field output.aggregate_groups}."
  ))
  invisible(projected)
}

#' Generate figures and tables from posterior draws
#'
#' Port of `reports.generate_reports`. Per-unit figures are nested under
#' `figs/<group>/` when more than one group is reported; cross-group figures
#' and all tables go directly under `figs/`.
#'
#' @param draws A `countsynth_draws` frame.
#' @param output_dir Run directory; artifacts go to `<output_dir>/figs/`.
#' @param target_unit Headline unit (auto-detected when `NULL`).
#' @param groups Groups to render per-unit figures for (`NULL` = all).
#' @param figures Character vector of figure names to render (see
#'   [countsynth_output_opts()]); `NULL` renders everything, `character()` renders
#'   tables only.
#' @param aggregate_units List of [countsynth_aggregate_unit()] specs.
#' @param aggregate_groups List of [countsynth_aggregate_group()] specs.
#'   Applied after `aggregate_units`, so a pooled group also covers the
#'   pooled units. The pooled group flows into every table and gets its own
#'   `figures/<group>/` directory; it enters `interval.png` only when
#'   `interval_aggregates` is `TRUE`, and never enters
#'   `group_comparison.png`, which exists to contrast the fitted groups
#'   against each other.
#' @param ppc_draws Alternative draws frame for the PPC suite only (cut mode:
#'   the stage-1 PPC frame).
#' @param ppc_units,ppc_exclude_units,ppc_acf_lags,ppc_unit_corr_max_time
#'   PPC options (see [countsynth_output_opts()]).
#' @param fit_gap_per_unit Also render fit/gap for every treated unit
#'   (restricted to the `"total"` group).
#' @param interval_aggregates Include the `aggregate_units` in `interval.png`,
#'   split off into their own band above the individual units, and the
#'   `aggregate_groups` as an extra colored series. Set `FALSE` to plot only
#'   the real units and groups. No effect when both are `NULL`.
#' @param progress Report progress to the terminal as each stage of the
#'   report is built. On by default: a report over tens of millions of draw
#'   rows can run for a long time, and a silent run gives no way to tell a
#'   slow stage from a stuck one.
#' @param print_tables Print the by-unit summary table to the terminal.
#' @param print_target_table Also print the target unit's own table. Its rows
#'   are the target unit's slice of the by-unit table, so this is off by
#'   default when `target_unit` is not set explicitly.
#' @param html_tables Also write `summary_table.html` and
#'   `summary_table_by_unit.html` via [countsynth_gt_table()]. Needs the `gt`
#'   package; warns and skips when it is missing.
#' @param fit Optional `countsynth_fit` / `countsynth_cut_fit` the draws came from.
#'   Required for the `"te_regression"` and `"rank_shrinkage"` figures, which
#'   read the treatment-effect design and the factor block rather than the
#'   draws frame; without it those figures are skipped.
#' @inheritParams countsynth_output_opts
#' @return Invisible list with `summary`, `per_unit`, `detail`,
#'   `target_unit`, `figures_dir`, `tables_dir`, `diagnostics_dir`,
#'   `treated_units`.
#' @export
countsynth_report <- function(draws, output_dir, target_unit = NULL, groups = NULL,
                         figures = NULL, aggregate_units = NULL,
                         aggregate_groups = NULL,
                         ppc_draws = NULL, ppc_units = NULL,
                         ppc_exclude_units = NULL, ppc_acf_lags = NULL,
                         ppc_unit_corr_max_time = NULL,
                         fit_gap_per_unit = FALSE,
                         interval_aggregates = TRUE, progress = TRUE,
                         print_tables = TRUE,
                         print_target_table = FALSE, html_tables = TRUE,
                         fit = NULL, rate_normalizer = 1000,
                         denominator_label = "denominator",
                         rate_label = NULL,
                         denominator_time_unit = "year",
                         denominator_may_be_affected = TRUE) {
  checkmate::assert_flag(progress)
  selected <- figures %||% FIGURE_NAMES
  unknown <- setdiff(selected, FIGURE_NAMES)
  if (length(unknown) > 0) {
    cli::cli_abort(
      "figures contains unknown name{?s} {.val {unknown}};
       valid names are {.val {sort(FIGURE_NAMES)}}."
    )
  }
  # Three destinations rather than one figs/: what happened (figures), the
  # numbers behind it (tables), and whether to believe the fit at all
  # (model_diagnostics). The last is named for what it holds -- the rank
  # figures in it were never posterior predictive checks, so "ppc" was wrong.
  figures_dir <- file.path(output_dir, "figures")
  tables_dir <- file.path(output_dir, "tables")
  diagnostics_dir <- file.path(output_dir, "model_diagnostics")
  for (d in c(figures_dir, tables_dir, diagnostics_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  report_size_guard(draws, aggregate_units, aggregate_groups)

  if (progress) cli::cli_progress_step("Building reporting frame")
  reporting <- if (!is.null(aggregate_units)) {
    add_aggregate_units(draws, aggregate_units)
  } else {
    draws
  }
  # Groups second, over a frame that already holds the pooled units, so the
  # pooled group covers them too -- the "all treated units, all groups" cell.
  if (!is.null(aggregate_groups)) {
    reporting <- add_aggregate_groups(reporting, aggregate_groups)
  }
  # raw_rate.png pools treated vs control across units itself, so it must not
  # see the aggregate units -- each would double-count every unit it covers.
  # It does need the aggregate groups, or the pooled group's own figure
  # directory below would be handed a frame with no rows for it.
  #
  # Both it and group_comparison.png read observed columns only, so they work
  # from the collapsed cell frame: pooling groups on the draws-sized frame
  # repeated the most expensive reduction in the report to produce a figure
  # that never looks at a draw (see observed_cell_frame()).
  observed_cells <- observed_cell_frame(draws)
  raw_rate_source <- if (!is.null(aggregate_groups)) {
    add_aggregate_groups(observed_cells, aggregate_groups)
  } else {
    observed_cells
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

  if (progress) cli::cli_progress_step("Summarising posterior draws")
  quantiles <- compute_quantiles(reporting)

  # Per-group unit-level figures.
  for (grp in report_groups) {
    # Closed at the foot of the loop rather than left for the next step to
    # supersede: cli renders a step's completion line by re-evaluating its
    # template, and by then `grp` has already moved on to the next group.
    if (progress) cli::cli_progress_step("Rendering figures for group {grp}")
    grp_dir <- if (length(report_groups) > 1) {
      file.path(figures_dir, grp)
    } else {
      figures_dir
    }
    dir.create(grp_dir, recursive = TRUE, showWarnings = FALSE)
    for (unit in fit_gap_units(grp)) {
      slug <- unit_slug(unit)
      if ("unit_fit" %in% selected) {
        save_plot(
          countsynth_unit_fit_plot(quantiles, unit, grp),
          file.path(grp_dir, sprintf("fit_%s.png", slug))
        )
      }
      if ("unit_gap" %in% selected) {
        save_plot(
          countsynth_unit_gap_plot(quantiles, unit, grp),
          file.path(grp_dir, sprintf("gap_%s.png", slug))
        )
      }
    }
    if ("raw_rate" %in% selected) {
      save_plot(
        countsynth_raw_rate_plot(
          raw_rate_source,
          group = grp, separate_unit = target_unit,
          rate_multiplier = rate_normalizer,
          denominator_label = denominator_label,
          rate_label = rate_label
        ),
        file.path(grp_dir, "raw_rate.png")
      )
    }
    if (progress) cli::cli_progress_done()
  }

  # Cross-group figures.
  if (progress) cli::cli_progress_step("Rendering cross-group figures")
  if ("interval" %in% selected) {
    # An aggregate unit pools the same draws as the units it covers, so it is
    # shown in its own band rather than ranked among them; `interval_aggregates
    # = FALSE` drops it from the figure entirely.
    save_plot(
      countsynth_interval_plot(
        if (interval_aggregates) reporting else draws,
        estimand = "ratio", method = "mu",
        rate_normalizer = rate_normalizer,
        denominator_label = denominator_label,
        denominator_time_unit = denominator_time_unit,
        rate_label = rate_label,
        separate_units = if (interval_aggregates) NULL else character()
      ),
      file.path(figures_dir, "interval.png"),
      width = 10, height = 8
    )
  }
  if ("group_comparison" %in% selected) {
    save_plot(
      countsynth_group_comparison_plot(
        observed_cells,
        rate_multiplier = rate_normalizer,
        denominator_label = denominator_label,
        rate_label = rate_label
      ),
      file.path(figures_dir, "group_comparison.png"),
      width = 11, height = 7
    )
  }
  if ("te_regression" %in% selected && !is.null(fit) &&
    !is.null(fit$te_design)) {
    te_dir <- file.path(figures_dir, "te")
    dir.create(te_dir, recursive = TRUE, showWarnings = FALSE)
    te_figs <- te_report_figures(fit)
    for (nm in names(te_figs)) {
      save_plot(te_figs[[nm]], file.path(te_dir, paste0(nm, ".png")),
        width = 10, height = 6
      )
    }
    write_table_csv(
      countsynth_te_coef_table(fit), file.path(te_dir, "te_coefficients.csv")
    )
  }
  if ("ppc" %in% selected) {
    if (progress) cli::cli_progress_step("Running posterior predictive checks")
    ppc_dir <- diagnostics_dir
    dir.create(ppc_dir, recursive = TRUE, showWarnings = FALSE)
    ppc_source <- ppc_draws %||% reporting
    if (!is.null(ppc_draws)) {
      if (!is.null(aggregate_units)) {
        ppc_source <- add_aggregate_units(ppc_source, aggregate_units)
      }
      if (!is.null(aggregate_groups)) {
        ppc_source <- add_aggregate_groups(ppc_source, aggregate_groups)
      }
    }
    ppc <- countsynth_ppc_plots(
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
    # The figure's caption names only the first few; the full screen-out list
    # is written beside the p-values, and written even when it is empty so a
    # reader can tell "nothing was dropped" from "this run predates the file".
    if (!is.null(ppc$unit_corr_dropped)) {
      write_table_csv(
        ppc$unit_corr_dropped,
        file.path(ppc_dir, "ppc_unit_corr_dropped.csv")
      )
    }
  }
  # Rank-shrinkage diagnostics go beside the PPC figures deliberately: they
  # answer the question ("is this rank adequate?") that a PPC-driven rank
  # sweep was answering, and are read together with it. Skipped, not an
  # error, when the fit did not use shrinkage -- the same way te_regression
  # is skipped without a design.
  if (!is.null(fit) && !is.null(rank_shrinkage_source(fit))) {
    dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
    # The CSVs follow the always-write rule the other tables do: eff_rank is
    # how you learn whether the truncation was generous enough, which is a
    # property of the fit rather than a picture someone asked for. Only the
    # rendering is gated, since that is what the `figures` list is for.
    write_table_csv(
      countsynth_component_weight_summary(fit),
      file.path(diagnostics_dir, "rank_component_weight.csv")
    )
    write_table_csv(
      countsynth_eff_rank_summary(fit),
      file.path(diagnostics_dir, "rank_eff_rank.csv")
    )
    if ("rank_shrinkage" %in% selected) {
      rank_figs <- rank_report_figures(fit)
      for (nm in names(rank_figs)) {
        save_plot(rank_figs[[nm]], file.path(diagnostics_dir, paste0(nm, ".png")))
      }
    }
  }

  # Tables: always written, never gated by `figures`.
  if (progress) cli::cli_progress_step("Writing tables")
  # Returned (and printed) but not written: its rows are the target unit's
  # slice of summary_table_by_unit.csv below.
  summary_tbl <- countsynth_summary_table(
    reporting, target_unit,
    rate_normalizer = rate_normalizer,
    denominator_label = denominator_label,
    denominator_time_unit = denominator_time_unit
  )
  by_unit_tbl <- countsynth_summary_table_by_unit(
    reporting, treated_units,
    rate_normalizer = rate_normalizer,
    denominator_label = denominator_label,
    denominator_time_unit = denominator_time_unit
  )
  write_table_csv(
    by_unit_tbl, file.path(tables_dir, "summary_table_by_unit.csv")
  )
  detail <- countsynth_expected_vs_observed(reporting, target_unit)
  write_table_csv(
    detail, file.path(tables_dir, "expected_vs_observed.csv")
  )
  per_unit <- countsynth_post_treatment_summary(reporting)
  write_table_csv(
    per_unit, file.path(tables_dir, "post_treatment_summary.csv")
  )
  if (html_tables) {
    write_gt_tables(
      reporting, target_unit, tables_dir,
      rate_normalizer = rate_normalizer,
      denominator_label = denominator_label,
      denominator_time_unit = denominator_time_unit,
      rate_label = rate_label,
      denominator_may_be_affected = denominator_may_be_affected
    )
  }

  if (progress) cli::cli_progress_done()

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
    # column-anchored footnote (see countsynth_gt_table()), so say in words what
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
          "assumption (the denominator cancels out of them). To estimate the",
          "causal effect on the count itself, fit the model with no",
          "denominator. Set {.field output.denominator_may_be_affected} to",
          "FALSE to silence this."
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
    target_unit = target_unit, figures_dir = figures_dir,
    tables_dir = tables_dir, diagnostics_dir = diagnostics_dir,
    treated_units = treated_units
  ))
}
