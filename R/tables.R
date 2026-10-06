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
#' @param draws A `countsynth_draws` frame.
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

#' Pick the headline unit for reporting
#'
#' A synthetic aggregate unit wins when the frame has one (see
#' [add_aggregate_units()]): having gone to the trouble of defining a pooled
#' unit, that pooled effect is the headline, not whichever single unit happens
#' to have the longest exposure. Otherwise the unit with the most
#' post-treatment periods wins, ties broken by unit order.
#'
#' @param draws A `countsynth_draws` frame.
#' @export
auto_detect_target <- function(draws) {
  treated <- draws[!is.na(draws$treatment) & draws$treatment == 1, ]
  if (nrow(treated) == 0) {
    return(NULL)
  }
  # A spec can resolve to nothing and be skipped, so go by the rows that are
  # actually there rather than trusting the attribute alone.
  aggregates <- intersect(aggregate_unit_names(draws), unique(treated$unit))
  if (length(aggregates) > 0) {
    return(aggregates[1])
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
#' Port of `tables.make_summary_table`. Exposure-weighted rates
#' (`denominator * years` summed over post-treatment periods), rate
#' difference and percent change with equal-tailed 95% intervals, and a
#' two-sided posterior p-value (`*` marks p < 0.05 on the group label). A
#' dagger marks a group whose post-treatment window includes any cell
#' imputed by [add_aggregate_units()] (see the `Imputed` column).
#'
#' `Pct Change CI` is `100 * (treated/untreated - 1)` on raw counts
#' (`treated = sum(exp(mu_treated))`, `untreated = sum(exp(mu))`); routing it
#' through `treated_rate`/`untreated_rate` first doesn't change that, since
#' both divide by the same `denom_val` and `rate_normalizer`, which cancel.
#' The *number* is therefore identical whether or not a denominator is
#' configured -- but what it *estimates* is not, and that's the point that
#' actually matters:
#'
#' - **With a denominator**, `mu_ctrl = log(rate) + log(denominator)`, where
#'   `log(denominator)` is the same fixed, observed offset on the treated
#'   and untreated side alike (see `stan_data_joint()`). It cancels not just
#'   algebraically but *causally*: the model holds the measured exposure
#'   fixed, so any difference between the two sides is attributed entirely
#'   to the rate term. `Pct Change` is therefore an estimate of the percent
#'   change in the underlying **rate** (events per unit of that measured,
#'   held-fixed exposure), and the count-level `Expected`/`Diff` columns are
#'   just that rate effect applied to the observed exposure.
#' - **Without a denominator** (equivalent to denominator `1` everywhere),
#'   there is no separately measured exposure for the model to hold fixed.
#'   `Pct Change` is then only an estimate of the percent change in the
#'   raw **count**, and that number is consistent with many different
#'   rate/exposure stories -- a higher per-unit rate, a change in some real
#'   but unmeasured denominator, or both. The model has no way to tell these
#'   apart, so without a denominator, read `Pct Change` as "the count
#'   changed by X%," not "the rate changed by X%."
#'
#' A related but distinct caveat, and one that lands on the **count**
#' columns only. The model parameterizes the rate directly, so the
#' counterfactual rate behind `Exp Rate` needs no assumption about what the
#' denominator would have been -- the denominator cancels out of it. Turning
#' that rate into a counterfactual *count* is what requires multiplying by
#' some denominator, and `Expected` uses the one actually observed. That's
#' the right estimand when the denominator is exogenous to treatment (e.g.
#' total population for a mortality rate). It understates the full causal
#' picture when treatment could plausibly change the denominator too -- e.g.
#' births as the denominator for an infant mortality rate, when the exposure
#' could also change the number of births: `Expected` then answers "how many
#' deaths, given the births that were actually observed under treatment"
#' rather than "how many deaths under a world with no treatment at all"
#' (which would need a counterfactual birth count too, and in general a
#' counterfactual denominator is not identified by this model at all).
#'
#' When the effect on the *count* is the target, fit with no denominator at
#' all: `mu_ctrl` is then the untreated log-count surface itself, so the
#' counterfactual count is extrapolated from the factor structure rather than
#' conditioned on an observed denominator treatment may have moved. The
#' trade is that nothing adjusts for exposure any more -- whatever the
#' denominator would have done has to be carried by the factors. See
#' `denominator_may_be_affected` in [countsynth_output_opts()] to surface this
#' caveat in the rendered report, where it is a footnote on `Expected`.
#'
#' @param draws A `countsynth_draws` frame.
#' @param target_unit Unit to summarize.
#' @param rate_normalizer Rates are per this many units of exposure
#'   (default 1000). A display scale only -- it cancels out of `Pct Change`.
#' @return A tibble with pre-formatted CI columns (parity with the Python
#'   CSV) plus a logical `Imputed` column, or an empty tibble when the unit
#'   has no post-treatment rows. The rate columns (the exposure column,
#'   `Obs Rate`, `Exp Rate`, `Rate Diff CI`) are present only when the run
#'   has a denominator -- without one there is no exposure to divide by, so
#'   they would be counts divided by the summed period lengths rather than
#'   rates. `Pct Change CI` is always present (the denominator cancels out of
#'   it), and the returned tibble carries a `has_denominator` attribute.
#' @export
countsynth_summary_table <- function(draws, target_unit = NULL,
                                rate_normalizer = 1000,
                                denominator_label = "denominator",
                                denominator_time_unit = "year") {
  target_unit <- target_unit %||% auto_detect_target(draws)
  if (is.null(target_unit)) {
    cli::cli_abort("No treated units in draws and no {.arg target_unit} given.")
  }
  df <- draws[draws$unit == target_unit &
    !is.na(draws$treatment) & draws$treatment == 1, ]
  if (nrow(df) == 0) {
    return(tibble::tibble())
  }
  df$years <- time_weight_per_row(df, denominator_time_unit)
  # countsynth_draws() always emits a denominator column (1 everywhere when none
  # was configured -- see build_cell_table()), but a hand-built frame may not.
  has_denom <- draws_has_denominator(draws)
  if (!"denominator" %in% names(df)) {
    df$denominator <- 1
  }
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

  exposure_col <- exposure_column_name(
    denominator_label, denominator_time_unit
  )
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
      !!exposure_col := as.integer(mean(gd$denom_val)),
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
  out <- dplyr::bind_rows(rows)
  # Without a configured denominator there is no exposure to divide by:
  # `denom_val` collapses to the summed period lengths (0.25 per quarterly
  # cell, say), so Exposure truncates to 0 and every rate is a count
  # divided by very nearly nothing. Those columns are not unstable so much as
  # meaningless, so drop them rather than print them -- Pct Change survives
  # because it is a ratio, and the shared denominator cancels out of it.
  if (!has_denom) {
    out[c(exposure_col, "Obs Rate", "Exp Rate", "Rate Diff CI")] <- NULL
  }
  attr(out, "has_denominator") <- has_denom
  attr(out, "exposure_col") <- exposure_col
  out
}

#' Headline summary table stacked over several units
#'
#' [countsynth_summary_table()] run per unit and row-bound, with a `Unit` column in
#' front. Units with no post-treatment rows drop out.
#'
#' @param draws A `countsynth_draws` frame.
#' @param units Units to include (`NULL` = every treated unit, aggregates
#'   included when the frame has them).
#' @param rate_normalizer Rates are per this many units of exposure.
#' @return A tibble, or an empty tibble when no unit has post-treatment rows.
#' @export
countsynth_summary_table_by_unit <- function(draws, units = NULL,
                                        rate_normalizer = 1000,
                                        denominator_label = "denominator",
                                        denominator_time_unit = "year") {
  units <- units %||% identify_treated_units(draws)
  dplyr::bind_rows(lapply(units, function(u) {
    tbl <- countsynth_summary_table(
      draws, u,
      rate_normalizer = rate_normalizer,
      denominator_label = denominator_label,
      denominator_time_unit = denominator_time_unit
    )
    if (nrow(tbl) == 0) {
      return(tbl)
    }
    tibble::add_column(tbl, Unit = u, .before = 1)
  }))
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
#' All raw counts here -- no denominator involved -- but `expected` still
#' assumes the observed denominator baked into `mu` is fixed/exogenous; see
#' [countsynth_summary_table()] for the caveat when treatment could also affect
#' the denominator.
#'
#' @param draws A `countsynth_draws` frame.
#' @export
countsynth_post_treatment_summary <- function(draws) {
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
#' @param draws A `countsynth_draws` frame.
#' @param target_unit Unit flagged in the `treated_unit` column
#'   (auto-detected when `NULL`).
#' @export
countsynth_expected_vs_observed <- function(draws, target_unit = NULL) {
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
