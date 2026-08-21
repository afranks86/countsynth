# Table math on the tidy draws frame. Port of tables.py; the summary-table
# columns are pre-formatted strings on purpose (parity with the Python CSVs
# and the JAMA supplement layout). Quantiles use stats::quantile type 7, the
# same linear interpolation as numpy's default.

#' Posterior-predictive quantiles per (unit, time, group)
#'
#' Port of `tables._compute_quantiles`: mean / median / equal-tailed 95%
#' interval of `ypred`. Rows whose grouping keys contain NA (e.g. a missing
#' observed outcome) are dropped, matching pandas' groupby behavior in the
#' Python implementation. `outcome_imputed`, when present on `draws` (set by
#' [add_aggregate_units()] for cells where a source unit's suppressed count
#' was filled in), is carried through so plots can mark those points.
#'
#' @param draws A `bpnmf_draws` frame.
#' @return A tibble keyed by unit, time, group, outcome, treatment, with an
#'   `outcome_imputed` column (`FALSE` when `draws` doesn't have one).
#' @export
compute_quantiles <- function(draws) {
  has_imputed <- "outcome_imputed" %in% names(draws)
  group_vars <- c("unit", "time", "group", "outcome", "treatment")
  if (has_imputed) group_vars <- c(group_vars, "outcome_imputed")
  out <- draws |>
    dplyr::filter(!is.na(.data$outcome), !is.na(.data$treatment)) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      ypred_mean = mean(.data$ypred),
      ypred_lower = stats::quantile(.data$ypred, 0.025, names = FALSE),
      ypred_upper = stats::quantile(.data$ypred, 0.975, names = FALSE),
      ypred_median = stats::median(.data$ypred),
      .groups = "drop"
    )
  if (!has_imputed) out$outcome_imputed <- FALSE
  out
}

#' Pick the unit with the most post-treatment observations
#' @param draws A `bpnmf_draws` frame.
#' @export
auto_detect_target <- function(draws) {
  treated <- draws[!is.na(draws$treatment) & draws$treatment == 1, ]
  if (nrow(treated) == 0) {
    return(NULL)
  }
  counts <- tapply(treated$time, treated$unit, function(t) length(unique(t)))
  names(which.max(counts))
}

years_per_row <- function(df) {
  if (all(c("start_date", "end_date") %in% names(df))) {
    as.numeric(as.Date(df$end_date) - as.Date(df$start_date)) / 365.25
  } else {
    rep(1.0, nrow(df))
  }
}

fmt_ci <- function(mean, lower, upper, digits = 2, suffix = "") {
  sprintf(
    paste0("%.", digits, "f%s (%.", digits, "f%s, %.", digits, "f%s)"),
    mean, suffix, lower, suffix, upper, suffix
  )
}

#' Headline observed-vs-expected summary table, one row per group
#'
#' Port of `tables.make_summary_table`. Person-year-weighted rates
#' (`denominator * years` summed over post-treatment periods), rate
#' difference and percent change with equal-tailed 95% intervals, and a
#' two-sided posterior p-value (`*` marks p < 0.05 on the group label). A
#' dagger marks a group whose post-treatment window includes any cell
#' imputed by [add_aggregate_units()] (see the `Imputed` column).
#'
#' @param draws A `bpnmf_draws` frame.
#' @param target_unit Unit to summarize.
#' @param rate_normalizer Rates are per this many person-years (default 1000).
#' @return A tibble with pre-formatted CI columns (parity with the Python
#'   CSV) plus a logical `Imputed` column, or an empty tibble when the unit
#'   has no post-treatment rows.
#' @export
bpnmf_summary_table <- function(draws, target_unit = NULL,
                                rate_normalizer = 1000) {
  target_unit <- target_unit %||% auto_detect_target(draws)
  if (is.null(target_unit)) {
    cli::cli_abort("No treated units in draws and no {.arg target_unit} given.")
  }
  df <- draws[draws$unit == target_unit &
    !is.na(draws$treatment) & draws$treatment == 1, ]
  if (nrow(df) == 0) {
    return(tibble::tibble())
  }
  df$years <- years_per_row(df)
  has_imputed <- "outcome_imputed" %in% names(df)

  draw_stats <- df |>
    dplyr::group_by(.data$group, .data$.draw) |>
    dplyr::summarise(
      ypred = sum(.data$ypred),
      outcome = sum(.data$outcome),
      treated = sum(exp(.data$mu_treated)),
      untreated = sum(exp(.data$mu)),
      denom_val = sum(.data$denominator * .data$years),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      treated_rate = .data$treated / .data$denom_val * rate_normalizer,
      untreated_rate = .data$untreated / .data$denom_val * rate_normalizer,
      outcome_rate = .data$outcome / .data$denom_val * rate_normalizer,
      outcome_diff = .data$treated - .data$untreated
    )

  imputed_by_group <- if (has_imputed) {
    df |>
      dplyr::distinct(.data$time, .data$group, .keep_all = TRUE) |>
      dplyr::group_by(.data$group) |>
      dplyr::summarise(imputed = any(.data$outcome_imputed), .groups = "drop")
  } else {
    NULL
  }

  q <- function(x, p) stats::quantile(x, p, names = FALSE)
  rows <- lapply(unique(draw_stats$group), function(grp) {
    gd <- draw_stats[draw_stats$group == grp, ]
    outcome_mean <- mean(gd$outcome)
    diff <- gd$outcome_diff
    rate_diff <- gd$treated_rate - gd$untreated_rate
    pct <- 100 * (gd$treated_rate / gd$untreated_rate - 1)
    pval <- 2 * min(mean(gd$untreated > gd$treated), mean(gd$untreated < gd$treated))
    sig <- if (pval < 0.05) "*" else ""
    imputed <- !is.null(imputed_by_group) &&
      isTRUE(imputed_by_group$imputed[imputed_by_group$group == grp])
    tibble::tibble(
      Group = paste0(grp, sig, if (imputed) " \u2020" else ""),
      Imputed = imputed,
      `Person-Years` = as.integer(mean(gd$denom_val)),
      Observed = as.integer(outcome_mean),
      Expected = as.integer(outcome_mean - mean(diff)),
      `Diff (95% CI)` = sprintf(
        "%d (%d, %d)",
        as.integer(mean(diff)), as.integer(q(diff, 0.025)), as.integer(q(diff, 0.975))
      ),
      `Obs Rate` = round(mean(gd$outcome_rate), 2),
      `Exp Rate` = round(mean(gd$outcome_rate) - mean(rate_diff), 2),
      `Rate Diff CI` = fmt_ci(mean(rate_diff), q(rate_diff, 0.025), q(rate_diff, 0.975)),
      `Pct Change CI` = fmt_ci(
        mean(pct), q(pct, 0.025), q(pct, 0.975),
        digits = 1, suffix = "%"
      )
    )
  })
  dplyr::bind_rows(rows)
}

#' Per-(unit, group) post-treatment totals, ranked by percent excess
#'
#' Port of `tables._compute_per_unit_post_treatment`; estimands match the
#' upstream R `make_fertility_table` / JAMA supplement:
#' `expected = sum(exp(mu))` (counterfactual), `treated = sum(exp(mu_treated))`,
#' `excess = treated - expected`, `excess_pct = 100 * (treated/expected - 1)`,
#' each with the draw-level equal-tailed 95% interval. `observed` is retained
#' for transparency only. `observed_imputed` flags a unit-group whose
#' post-treatment window includes a cell imputed by [add_aggregate_units()].
#'
#' @param draws A `bpnmf_draws` frame.
#' @export
bpnmf_post_treatment_summary <- function(draws) {
  post <- draws[!is.na(draws$treatment) & draws$treatment == 1, ]
  cols <- c(
    "unit", "group", "n_periods", "observed", "observed_imputed",
    "expected_mean", "expected_lower_95", "expected_upper_95",
    "excess_mean", "excess_lower_95", "excess_upper_95",
    "excess_pct_mean", "excess_pct_lower_95", "excess_pct_upper_95"
  )
  if (nrow(post) == 0) {
    empty <- tibble::as_tibble(
      stats::setNames(rep(list(numeric(0)), length(cols)), cols)
    )
    return(empty)
  }
  has_imputed <- "outcome_imputed" %in% names(post)

  draw_sums <- post |>
    dplyr::group_by(.data$unit, .data$group, .data$.draw) |>
    dplyr::summarise(
      expected = sum(exp(.data$mu)),
      treated = sum(exp(.data$mu_treated)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      excess = .data$treated - .data$expected,
      excess_pct = 100 * (.data$treated / .data$expected - 1)
    )

  observed_totals <- post |>
    dplyr::distinct(.data$unit, .data$group, .data$time, .keep_all = TRUE) |>
    dplyr::group_by(.data$unit, .data$group) |>
    dplyr::summarise(
      n_periods = length(unique(.data$time)),
      observed = sum(.data$outcome),
      observed_imputed = if (has_imputed) any(.data$outcome_imputed) else FALSE,
      .groups = "drop"
    )

  q <- function(x, p) stats::quantile(x, p, names = FALSE)
  stats_df <- draw_sums |>
    dplyr::group_by(.data$unit, .data$group) |>
    dplyr::summarise(
      expected_mean = mean(.data$expected),
      expected_lower_95 = q(.data$expected, 0.025),
      expected_upper_95 = q(.data$expected, 0.975),
      excess_mean = mean(.data$excess),
      excess_lower_95 = q(.data$excess, 0.025),
      excess_upper_95 = q(.data$excess, 0.975),
      excess_pct_mean = mean(.data$excess_pct),
      excess_pct_lower_95 = q(.data$excess_pct, 0.025),
      excess_pct_upper_95 = q(.data$excess_pct, 0.975),
      .groups = "drop"
    )

  observed_totals |>
    dplyr::left_join(stats_df, by = c("unit", "group")) |>
    dplyr::arrange(dplyr::desc(.data$excess_pct_mean)) |>
    dplyr::select(dplyr::all_of(cols))
}

#' Per-(unit, time, group) expected-vs-observed detail
#'
#' Port of the `expected_vs_observed.csv` construction in `reports.py`:
#' renamed posterior-predictive quantiles plus `gap = observed -
#' expected_mean` and `gap_pct`. `observed_imputed` flags a cell where a
#' source unit's suppressed count was filled in by [add_aggregate_units()].
#'
#' @param draws A `bpnmf_draws` frame.
#' @param target_unit Unit flagged in the `treated_unit` column
#'   (auto-detected when `NULL`).
#' @export
bpnmf_expected_vs_observed <- function(draws, target_unit = NULL) {
  target_unit <- target_unit %||% auto_detect_target(draws)
  detail <- compute_quantiles(draws) |>
    dplyr::rename(
      observed = "outcome",
      observed_imputed = "outcome_imputed",
      expected_mean = "ypred_mean",
      expected_median = "ypred_median",
      expected_lower_95 = "ypred_lower",
      expected_upper_95 = "ypred_upper"
    ) |>
    dplyr::mutate(
      treated_unit = .data$unit == (target_unit %||% ""),
      gap = .data$observed - .data$expected_mean,
      gap_pct = .data$gap / .data$expected_mean * 100
    ) |>
    dplyr::select(
      "unit", "time", "group", "treatment", "treated_unit", "observed",
      "observed_imputed",
      "expected_mean", "expected_median", "expected_lower_95",
      "expected_upper_95", "gap", "gap_pct"
    ) |>
    dplyr::arrange(.data$unit, .data$time)
  detail
}

#' Units with any post-treatment rows
#' @keywords internal
identify_treated_units <- function(draws) {
  unique(draws$unit[!is.na(draws$treatment) & draws$treatment == 1])
}
