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
# data at all" -- e.g. countsynth_ppc_acf's lag exceeding every cell's available
# control-period length turns every autocorrelation into NA, and those rows
# get filtered out entirely. Render an explanatory plot instead of crashing.
ppc_insufficient_data_plot <- function(label) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = label) +
    ggplot2::theme_void()
}

# Left-aligned and small: the caption is a footnote about the data the figure
# could use, not a second title.
add_plot_caption <- function(p, caption) {
  if (is.null(caption)) {
    return(p)
  }
  p +
    ggplot2::labs(caption = caption) +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(
        hjust = 0, size = 8, colour = "grey30"
      ),
      plot.caption.position = "plot"
    )
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
    theme_countsynth()
}

#' PPC: max absolute residual per (unit, group, draw)
#' @param draws A `countsynth_draws` frame.
#' @param categories,ppc_units,ppc_exclude_units Filters (see
#'   [countsynth_output_opts()]).
#' @return List with `plot` (ggplot) and `pvals` (unit x group tibble).
#' @export
countsynth_ppc_abs <- function(draws, categories = NULL, ppc_units = NULL,
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
#' @inheritParams countsynth_ppc_abs
#' @param lag Autocorrelation lag (default 1).
#' @export
countsynth_ppc_acf <- function(draws, lag = 1, categories = NULL, ppc_units = NULL,
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
#' @inheritParams countsynth_ppc_abs
#' @export
countsynth_ppc_rmse <- function(draws, categories = NULL, ppc_units = NULL,
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

# A (unit, group) whose observed control-period outcome never moves -- a state
# with zero homicides in every period, say -- has residual `outcome - exp(mu)`
# equal to `-exp(mu)`: a smooth deterministic curve, not noise. Any two such
# units correlate near-perfectly, so a handful of them dominates the observed
# correlation matrix and inflates its spectral norm, while the predicted
# residuals keep their sampling variation. The check then reports a tiny
# p-value that is an artifact of the sparse units rather than evidence about
# cross-sectional dependence. Drop those cells, along with the mostly-missing
# and mostly-zero ones that are degenerate for the same reason.
#
# The drops are reported on the figure's caption rather than at the console:
# which cells a diagnostic plot could use belongs beside the plot, where it is
# still legible months later in `model_diagnostics/`.
drop_degenerate_corr_cells <- function(df, max_na_frac = 0.25,
                                       max_zero_frac = 0.5) {
  empty <- tibble::tibble(
    unit = character(), group = character(), reason = character()
  )
  if (nrow(df) == 0) {
    return(list(df = df, dropped = empty, n_cells = 0L))
  }
  # One row per (unit, group, time): the outcome is repeated across draws, and
  # the counts below are about the observed series, not the posterior.
  series <- unique(df[, c("unit", "group", "time", "outcome")])
  cells <- series |>
    dplyr::group_by(.data$unit, .data$group) |>
    dplyr::summarise(
      na_frac = mean(is.na(.data$outcome)),
      zero_frac = mean(.data$outcome == 0, na.rm = TRUE),
      flat = {
        obs <- .data$outcome[!is.na(.data$outcome)]
        length(obs) < 2 || stats::sd(obs) == 0
      },
      .groups = "drop"
    )
  cells$reason <- dplyr::case_when(
    cells$na_frac >= max_na_frac ~ "mostly missing",
    cells$flat ~ "no variation",
    cells$zero_frac > max_zero_frac ~ "mostly zeros",
    TRUE ~ NA_character_
  )
  bad <- cells[!is.na(cells$reason), c("unit", "group", "reason")]
  list(
    df = dplyr::anti_join(df, bad, by = c("unit", "group")),
    dropped = bad,
    n_cells = nrow(cells)
  )
}

# One or two lines for the figure's caption: the tally by reason, then as many
# cell names as fit before the list stops being readable.
dropped_cells_caption <- function(dropped, n_cells, max_named = 8) {
  if (nrow(dropped) == 0) {
    return(NULL)
  }
  by_reason <- table(dropped$reason)
  tally <- paste(sprintf("%d %s", as.integer(by_reason), names(by_reason)),
                 collapse = ", ")
  labels <- sprintf("%s/%s", dropped$unit, dropped$group)
  named <- paste(utils::head(labels, max_named), collapse = ", ")
  if (length(labels) > max_named) {
    named <- sprintf("%s, and %d more", named, length(labels) - max_named)
  }
  lines <- c(
    sprintf(
      paste("Dropped %d of %d unit-group cells whose observed series carries",
            "no usable variation (%s)."),
      nrow(dropped), n_cells, tally
    ),
    paste0(named, ".")
  )
  paste(strwrap(lines, width = 110), collapse = "\n")
}

#' PPC: spectral norm of the cross-unit residual correlation matrix
#'
#' Tests whether the model captures cross-sectional dependence. Uses the
#' control period (or `time < max_treat_date` when given), every unit (this
#' check needs >= 2 units), capped at `ndraws` draws for speed. Cells whose
#' observed series carries no usable variation are dropped first: too many
#' missing outcomes, too many zeros, or a series that never moves at all.
#'
#' @inheritParams countsynth_ppc_abs
#' @param max_treat_date Optional `YYYY-MM-DD` cutoff replacing the
#'   control-period filter.
#' @param ndraws Cap on the number of draws used.
#' @param max_na_frac Drop a (unit, group) whose control-period outcome is
#'   missing in at least this fraction of periods.
#' @param max_zero_frac Drop a (unit, group) whose observed control-period
#'   outcome is zero in more than this fraction of the periods it does have.
#' @export
countsynth_ppc_unit_corr <- function(draws, max_treat_date = NULL,
                                categories = NULL, ndraws = 1000,
                                ppc_units = NULL, ppc_exclude_units = NULL,
                                max_na_frac = 0.25, max_zero_frac = 0.5) {
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

  screened <- drop_degenerate_corr_cells(df, max_na_frac, max_zero_frac)
  df <- screened$df
  caption <- dropped_cells_caption(screened$dropped, screened$n_cells)

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
      plot = add_plot_caption(
        ppc_insufficient_data_plot("Insufficient data for spectral norm computation"),
        caption
      ),
      pvals = tibble::tibble(group = character(), pval = numeric()),
      dropped = screened$dropped
    ))
  }
  pvals <- stats_df |>
    dplyr::group_by(.data$group) |>
    dplyr::summarise(pval = mean(.data$eval_diff < 0), .groups = "drop")
  list(
    plot = add_plot_caption(
      ppc_histogram(
        stats_df, pvals, "eval_diff",
        "Difference in Unit Correlations",
        "Observed - Predicted Spectral Norm", "group"
      ),
      caption
    ),
    pvals = pvals,
    dropped = screened$dropped
  )
}

#' Run the full PPC suite
#'
#' @inheritParams countsynth_ppc_abs
#' @param checks Subset of `c("abs", "acf", "rmse", "unit_corr")`.
#' @param acf_lags Integer lags for the ACF check (one figure per lag).
#' @param max_treat_date Cutoff date for the unit-correlation check.
#' @return List with `plots` (named list of ggplots), `pvals` (combined
#'   tibble with a `check_type` column) and `unit_corr_dropped` (the cells the
#'   unit-correlation check screened out, `NULL` when that check did not run).
#' @export
countsynth_ppc_plots <- function(draws, checks = c("abs", "acf", "rmse", "unit_corr"),
                            acf_lags = 1, categories = NULL,
                            max_treat_date = NULL, ppc_units = NULL,
                            ppc_exclude_units = NULL) {
  checkmate::assert_subset(checks, c("abs", "acf", "rmse", "unit_corr"))
  plots <- list()
  pvals <- list()
  unit_corr_dropped <- NULL
  if ("abs" %in% checks) {
    r <- countsynth_ppc_abs(draws, categories, ppc_units, ppc_exclude_units)
    plots$ppc_abs_residual <- r$plot
    pvals$abs <- dplyr::mutate(r$pvals, check_type = "abs")
  }
  if ("acf" %in% checks) {
    for (lag in acf_lags) {
      r <- countsynth_ppc_acf(draws, lag, categories, ppc_units, ppc_exclude_units)
      plots[[sprintf("ppc_acf_lag%d", lag)]] <- r$plot
      pvals[[sprintf("acf_lag%d", lag)]] <-
        dplyr::mutate(r$pvals, check_type = sprintf("acf_lag%d", lag))
    }
  }
  if ("rmse" %in% checks) {
    r <- countsynth_ppc_rmse(draws, categories, ppc_units, ppc_exclude_units)
    plots$ppc_rmse <- r$plot
    pvals$rmse <- dplyr::mutate(r$pvals, check_type = "rmse")
  }
  if ("unit_corr" %in% checks) {
    # Deliberately NOT filtered by ppc_units: the spectral norm needs >= 2
    # units (mirrors make_all_ppc_plots).
    r <- countsynth_ppc_unit_corr(
      draws,
      max_treat_date = max_treat_date, categories = categories,
      ppc_exclude_units = ppc_exclude_units
    )
    plots$ppc_unit_corr <- r$plot
    unit_corr_dropped <- r$dropped
    if (nrow(r$pvals) > 0) {
      pvals$unit_corr <- dplyr::mutate(r$pvals, check_type = "unit_corr")
    }
  }
  list(
    plots = plots, pvals = dplyr::bind_rows(pvals),
    unit_corr_dropped = unit_corr_dropped
  )
}
