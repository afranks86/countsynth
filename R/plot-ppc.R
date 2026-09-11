# Posterior predictive checks (ports of the four PPC suites in plots.py).
# All residuals are computed against the counterfactual rate on control-period
# cells: obs_diff = outcome - exp(mu), pred_diff = ypred - exp(mu); every
# p-value is P(observed statistic < predicted statistic) over draws.

filter_ppc_units <- function(df, treated_units, ppc_units, ppc_exclude_units) {
  if (!is.null(ppc_units)) {
    treated_units <- ppc_units
  }
  if (length(treated_units) > 0) {
    df <- df[df$unit %in% treated_units, ]
  }
  if (!is.null(ppc_exclude_units)) {
    df <- df[!df$unit %in% ppc_exclude_units, ]
  }
  df
}

# A (unit, group) whose control-period outcome is missing everywhere -- a
# subgroup suppressed in every pre-treatment period, say -- has no residual to
# check. Left in, `max(abs(obs_diff), na.rm = TRUE)` over an all-NA vector
# warns once per draw and returns -Inf, which then beats every predicted
# statistic and reports a clean p = 1 for a cell holding no data at all. The
# RMSE check turns the same cell into NaN, and the histograms drop the rows
# with a "non-finite values" warning. Drop the cells up front instead, once,
# by name.
drop_empty_ppc_cells <- function(df) {
  if (nrow(df) == 0) {
    return(df)
  }
  cells <- df |>
    dplyr::group_by(.data$unit, .data$group) |>
    dplyr::summarise(has_obs = any(!is.na(.data$obs_diff)), .groups = "drop")
  empty <- cells[!cells$has_obs, c("unit", "group")]
  if (nrow(empty) == 0) {
    return(df)
  }
  labels <- sprintf("%s/%s", empty$unit, empty$group)
  cli::cli_warn(c(
    "Skipping {nrow(empty)} PPC cell{?s} with no observed control-period data:
     {.val {labels}}.",
    i = "Every residual there is missing, so the check has nothing to compare
         and would report a p-value built from no data."
  ))
  dplyr::anti_join(df, empty, by = c("unit", "group"))
}

#' Shared residual prep for the PPC checks
#' @keywords internal
prepare_ppc_residuals <- function(draws, categories = NULL, ppc_units = NULL,
                                  ppc_exclude_units = NULL,
                                  sort_by_time = FALSE) {
  df <- draws
  categories <- categories %||% unique(df$group)
  df <- df[df$group %in% categories, ]
  treated_units <- identify_treated_units(df)
  df <- df[!is.na(df$treatment) & df$treatment == 0, ]
  df <- filter_ppc_units(df, treated_units, ppc_units, ppc_exclude_units)
  df$pred_diff <- df$ypred - exp(df$mu)
  df$obs_diff <- df$outcome - exp(df$mu)
  df <- drop_empty_ppc_cells(df)
  if (sort_by_time) {
    df <- dplyr::arrange(df, .data$unit, .data$group, .data$.draw, .data$time)
  }
  list(df = df, categories = categories, units = unique(df$unit))
}

# A 0-row main data frame crashes ggplot2's facet_wrap with "Faceting
# variables must have at least one value" (combine_vars() requires the
# layout have something to lay out). This can happen for reasons besides "no
# data at all" -- e.g. bpnmf_ppc_acf's lag exceeding every cell's available
# control-period length turns every autocorrelation into NA, and those rows
# get filtered out entirely. Render an explanatory plot instead of crashing.
ppc_insufficient_data_plot <- function(label) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = label) +
    ggplot2::theme_void()
}

#' Autocorrelation at one lag (NaN if too short or zero-variance)
#' @keywords internal
autocorrelation_at_lag <- function(x, lag) {
  n <- length(x)
  if (n <= lag) {
    return(NA_real_)
  }
  xc <- x - mean(x, na.rm = TRUE)
  v <- sum(xc^2, na.rm = TRUE)
  if (is.na(v) || v == 0) {
    return(NA_real_)
  }
  sum(xc[seq_len(n - lag)] * xc[seq(lag + 1, n)], na.rm = TRUE) / v
}

ppc_histogram <- function(stats_df, pvals_df, x_col, title, xlabel,
                          facet_vars) {
  pvals_df$label <- sprintf("p = %.3f", pvals_df$pval)
  ggplot2::ggplot(stats_df, ggplot2::aes(x = .data[[x_col]])) +
    ggplot2::geom_histogram(
      bins = 30, fill = "steelblue", alpha = 0.5, color = "white"
    ) +
    ggplot2::geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    ggplot2::geom_text(
      data = pvals_df,
      ggplot2::aes(x = Inf, y = Inf, label = .data$label),
      hjust = 1.1, vjust = 1.5, color = "red", fontface = "bold", size = 3,
      inherit.aes = FALSE
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(!!!rlang::syms(facet_vars)),
      ncol = if (length(facet_vars) == 1) 2 else 3, scales = "free"
    ) +
    ggplot2::labs(title = title, x = xlabel, y = "Count") +
    theme_bpnmf()
}

#' PPC: max absolute residual per (unit, group, draw)
#' @param draws A `bpnmf_draws` frame.
#' @param categories,ppc_units,ppc_exclude_units Filters (see
#'   [bpnmf_output_opts()]).
#' @return List with `plot` (ggplot) and `pvals` (unit x group tibble).
#' @export
bpnmf_ppc_abs <- function(draws, categories = NULL, ppc_units = NULL,
                          ppc_exclude_units = NULL) {
  prep <- prepare_ppc_residuals(draws, categories, ppc_units, ppc_exclude_units)
  stats_df <- prep$df |>
    dplyr::group_by(.data$unit, .data$group, .data$.draw) |>
    dplyr::summarise(
      max_pred_diff = max(abs(.data$pred_diff), na.rm = TRUE),
      max_obs_diff = max(abs(.data$obs_diff), na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(diff_in_diff = .data$max_obs_diff - .data$max_pred_diff)
  pvals <- stats_df |>
    dplyr::group_by(.data$unit, .data$group) |>
    dplyr::summarise(pval = mean(.data$diff_in_diff < 0), .groups = "drop")
  if (nrow(stats_df) == 0) {
    return(list(
      plot = ppc_insufficient_data_plot("Insufficient data for abs-residual check"),
      pvals = pvals
    ))
  }
  list(
    plot = ppc_histogram(
      stats_df, pvals, "diff_in_diff",
      "Difference in Maximum Absolute Predicted Residual",
      "Observed - Predicted Max Residual", c("unit", "group")
    ),
    pvals = pvals
  )
}

#' PPC: residual autocorrelation at one lag
#' @inheritParams bpnmf_ppc_abs
#' @param lag Autocorrelation lag (default 1).
#' @export
bpnmf_ppc_acf <- function(draws, lag = 1, categories = NULL, ppc_units = NULL,
                          ppc_exclude_units = NULL) {
  prep <- prepare_ppc_residuals(
    draws, categories, ppc_units, ppc_exclude_units,
    sort_by_time = TRUE
  )
  stats_df <- prep$df |>
    dplyr::group_by(.data$unit, .data$group, .data$.draw) |>
    dplyr::summarise(
      obs_ac = autocorrelation_at_lag(.data$obs_diff, lag),
      pred_ac = autocorrelation_at_lag(.data$pred_diff, lag),
      .groups = "drop"
    ) |>
    dplyr::mutate(diff_in_ac = .data$obs_ac - .data$pred_ac) |>
    dplyr::filter(!is.na(.data$diff_in_ac))
  pvals <- stats_df |>
    dplyr::group_by(.data$unit, .data$group) |>
    dplyr::summarise(pval = mean(.data$diff_in_ac < 0), .groups = "drop")
  if (nrow(stats_df) == 0) {
    return(list(
      plot = ppc_insufficient_data_plot(sprintf(
        "Insufficient data for lag-%d autocorrelation check\n(no cell has more than %d control-period observations)",
        lag, lag
      )),
      pvals = pvals
    ))
  }
  list(
    plot = ppc_histogram(
      stats_df, pvals, "diff_in_ac",
      sprintf("Difference in Residual Autocorrelation (Lag %d)", lag),
      "Observed - Predicted Autocorrelation", c("unit", "group")
    ),
    pvals = pvals
  )
}

#' PPC: residual RMSE per (unit, group, draw)
#' @inheritParams bpnmf_ppc_abs
#' @export
bpnmf_ppc_rmse <- function(draws, categories = NULL, ppc_units = NULL,
                           ppc_exclude_units = NULL) {
  prep <- prepare_ppc_residuals(draws, categories, ppc_units, ppc_exclude_units)
  stats_df <- prep$df |>
    dplyr::group_by(.data$unit, .data$group, .data$.draw) |>
    dplyr::summarise(
      rmse_pred_diff = sqrt(mean(.data$pred_diff^2, na.rm = TRUE)),
      rmse_obs_diff = sqrt(mean(.data$obs_diff^2, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::mutate(diff_in_diff = .data$rmse_obs_diff - .data$rmse_pred_diff)
  pvals <- stats_df |>
    dplyr::group_by(.data$unit, .data$group) |>
    dplyr::summarise(pval = mean(.data$diff_in_diff < 0), .groups = "drop")
  if (nrow(stats_df) == 0) {
    return(list(
      plot = ppc_insufficient_data_plot("Insufficient data for RMSE check"),
      pvals = pvals
    ))
  }
  list(
    plot = ppc_histogram(
      stats_df, pvals, "diff_in_diff",
      "Difference in RMSE", "Observed - Predicted RMSE", c("unit", "group")
    ),
    pvals = pvals
  )
}

spectral_norm <- function(m) {
  keep <- colSums(!is.na(m)) > 0
  m <- m[, keep, drop = FALSE]
  if (ncol(m) < 2) {
    return(NA_real_)
  }
  cm <- suppressWarnings(stats::cor(m, use = "pairwise.complete.obs"))
  if (anyNA(cm)) {
    cm[is.na(cm)] <- 0
    diag(cm) <- 1
  }
  sqrt(max(max(eigen(cm, symmetric = TRUE, only.values = TRUE)$values), 0))
}

#' PPC: spectral norm of the cross-unit residual correlation matrix
#'
#' Tests whether the model captures cross-sectional dependence. Uses the
#' control period (or `time < max_treat_date` when given), every unit (this
#' check needs >= 2 units), capped at `ndraws` draws for speed; units with
#' more than 25% missing outcomes are dropped.
#'
#' @inheritParams bpnmf_ppc_abs
#' @param max_treat_date Optional `YYYY-MM-DD` cutoff replacing the
#'   control-period filter.
#' @param ndraws Cap on the number of draws used.
#' @export
bpnmf_ppc_unit_corr <- function(draws, max_treat_date = NULL,
                                categories = NULL, ndraws = 1000,
                                ppc_units = NULL, ppc_exclude_units = NULL) {
  df <- draws
  categories <- categories %||% unique(df$group)
  df <- df[df$group %in% categories, ]
  if (!is.null(ppc_units) || !is.null(ppc_exclude_units)) {
    df <- filter_ppc_units(df, unique(df$unit), ppc_units, ppc_exclude_units)
  }
  if (!is.null(max_treat_date)) {
    df <- df[df$time < as.Date(max_treat_date), ]
  } else {
    df <- df[!is.na(df$treatment) & df$treatment == 0, ]
  }
  draw_ids <- unique(df$.draw)
  if (length(draw_ids) > ndraws) {
    df <- df[df$.draw %in% draw_ids[seq_len(ndraws)], ]
  }
  df$obs_residual <- df$outcome - exp(df$mu)
  df$pred_residual <- df$ypred - exp(df$mu)
  # pred residuals are masked where the observation is missing
  df$pred_residual[is.na(df$obs_residual)] <- NA_real_

  na_frac <- df |>
    dplyr::group_by(.data$unit, .data$group) |>
    dplyr::summarise(na_frac = mean(is.na(.data$outcome)), .groups = "drop")
  valid <- na_frac[na_frac$na_frac < 0.25, c("unit", "group")]
  df <- dplyr::inner_join(df, valid, by = c("unit", "group"))

  results <- list()
  for (grp in categories) {
    gd <- df[df$group == grp, ]
    if (nrow(gd) == 0) {
      next
    }
    units <- unique(gd$unit)
    times <- sort(unique(gd$time))
    for (dr in unique(gd$.draw)) {
      dd <- gd[gd$.draw == dr, ]
      t_idx <- match(dd$time, times)
      u_idx <- match(dd$unit, units)
      obs_m <- matrix(NA_real_, length(times), length(units))
      pred_m <- matrix(NA_real_, length(times), length(units))
      obs_m[cbind(t_idx, u_idx)] <- dd$obs_residual
      pred_m[cbind(t_idx, u_idx)] <- dd$pred_residual
      results[[length(results) + 1]] <- tibble::tibble(
        group = grp, .draw = dr,
        obs_sval = spectral_norm(obs_m),
        pred_sval = spectral_norm(pred_m)
      )
    }
  }
  stats_df <- dplyr::bind_rows(results)
  if (nrow(stats_df) > 0) {
    stats_df$eval_diff <- stats_df$obs_sval - stats_df$pred_sval
    stats_df <- stats_df[!is.na(stats_df$eval_diff), ]
  }
  if (nrow(stats_df) == 0) {
    return(list(
      plot = ppc_insufficient_data_plot("Insufficient data for spectral norm computation"),
      pvals = tibble::tibble(group = character(), pval = numeric())
    ))
  }
  pvals <- stats_df |>
    dplyr::group_by(.data$group) |>
    dplyr::summarise(pval = mean(.data$eval_diff < 0), .groups = "drop")
  list(
    plot = ppc_histogram(
      stats_df, pvals, "eval_diff",
      "Difference in Unit Correlations",
      "Observed - Predicted Spectral Norm", "group"
    ),
    pvals = pvals
  )
}

#' Run the full PPC suite
#'
#' @inheritParams bpnmf_ppc_abs
#' @param checks Subset of `c("abs", "acf", "rmse", "unit_corr")`.
#' @param acf_lags Integer lags for the ACF check (one figure per lag).
#' @param max_treat_date Cutoff date for the unit-correlation check.
#' @return List with `plots` (named list of ggplots) and `pvals` (combined
#'   tibble with a `check_type` column).
#' @export
bpnmf_ppc_plots <- function(draws, checks = c("abs", "acf", "rmse", "unit_corr"),
                            acf_lags = 1, categories = NULL,
                            max_treat_date = NULL, ppc_units = NULL,
                            ppc_exclude_units = NULL) {
  checkmate::assert_subset(checks, c("abs", "acf", "rmse", "unit_corr"))
  plots <- list()
  pvals <- list()
  if ("abs" %in% checks) {
    r <- bpnmf_ppc_abs(draws, categories, ppc_units, ppc_exclude_units)
    plots$ppc_abs_residual <- r$plot
    pvals$abs <- dplyr::mutate(r$pvals, check_type = "abs")
  }
  if ("acf" %in% checks) {
    for (lag in acf_lags) {
      r <- bpnmf_ppc_acf(draws, lag, categories, ppc_units, ppc_exclude_units)
      plots[[sprintf("ppc_acf_lag%d", lag)]] <- r$plot
      pvals[[sprintf("acf_lag%d", lag)]] <-
        dplyr::mutate(r$pvals, check_type = sprintf("acf_lag%d", lag))
    }
  }
  if ("rmse" %in% checks) {
    r <- bpnmf_ppc_rmse(draws, categories, ppc_units, ppc_exclude_units)
    plots$ppc_rmse <- r$plot
    pvals$rmse <- dplyr::mutate(r$pvals, check_type = "rmse")
  }
  if ("unit_corr" %in% checks) {
    # Deliberately NOT filtered by ppc_units: the spectral norm needs >= 2
    # units (mirrors make_all_ppc_plots).
    r <- bpnmf_ppc_unit_corr(
      draws,
      max_treat_date = max_treat_date, categories = categories,
      ppc_exclude_units = ppc_exclude_units
    )
    plots$ppc_unit_corr <- r$plot
    if (nrow(r$pvals) > 0) {
      pvals$unit_corr <- dplyr::mutate(r$pvals, check_type = "unit_corr")
    }
  }
  list(plots = plots, pvals = dplyr::bind_rows(pvals))
}
