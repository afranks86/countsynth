# Descriptive and model-fit figures (ggplot2 ports of plots.py). All figure
# functions return ggplot objects; saving/filenames live in reports.R.

PALETTE_TREATED <- "#E41A1C"
PALETTE_CONTROL <- "#999999"
PALETTE_SEPARATE <- "#FF7F00"

#' countsynth ggplot theme (matches the Python matplotlib whitegrid look)
#' @param base_size Base font size.
#' @export
theme_countsynth <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      strip.background = ggplot2::element_rect(fill = "grey92")
    )
}

# Accept either a draws frame (has .draw/ypred) or a precomputed quantiles
# frame from compute_quantiles().
as_quantiles <- function(x) {
  if ("ypred_mean" %in% names(x)) x else compute_quantiles(x)
}

first_treated_time <- function(df, unit) {
  tt <- df$time[df$unit == unit & !is.na(df$treatment) & df$treatment == 1]
  if (length(tt) == 0) NULL else min(tt)
}

#' Observed vs predicted fit for one unit
#'
#' Observed counts as points, posterior-predictive mean as a line, with the
#' equal-tailed 95% interval ribbon and a dashed vertical line at the unit's
#' first treated period. Cells where a source unit's suppressed count was
#' filled in by [add_aggregate_units()] (e.g. a synthetic "all treated"
#' unit) are ringed in orange.
#'
#' @param x A `countsynth_draws` frame or [compute_quantiles()] output.
#' @param unit Unit to plot.
#' @param group Group to plot (default `"total"` if present, else the first).
#' @return A ggplot object.
#' @export
countsynth_unit_fit_plot <- function(x, unit, group = NULL) {
  q <- as_quantiles(x)
  group <- group %||% if ("total" %in% q$group) "total" else q$group[[1]]
  sub <- q[q$unit == unit & q$group == group, ]
  if (nrow(sub) == 0) {
    cli::cli_abort("No rows for unit {.val {unit}}, group {.val {group}}.")
  }
  if (!"outcome_imputed" %in% names(sub)) sub$outcome_imputed <- FALSE
  imputed <- sub[sub$outcome_imputed, ]
  vline <- first_treated_time(sub, unit)
  p <- ggplot2::ggplot(sub, ggplot2::aes(x = .data$time)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$ypred_lower, ymax = .data$ypred_upper),
      fill = "steelblue", alpha = 0.3
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$ypred_mean),
      color = "steelblue", linewidth = 0.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$outcome),
      size = 1.4, alpha = 0.8
    ) +
    ggplot2::geom_point(
      data = imputed, ggplot2::aes(y = .data$outcome),
      shape = 1, size = 3, color = "#E6550D", stroke = 1
    ) +
    ggplot2::labs(
      title = sprintf("Model fit: %s (%s)", unit, group),
      x = "Time", y = "Count",
      subtitle = if (nrow(imputed) > 0) {
        "Points: observed. Orange circles: a source unit's suppressed count was imputed from the model. Line/ribbon: posterior predictive mean and 95% CI."
      } else {
        "Points: observed. Line/ribbon: posterior predictive mean and 95% CI."
      }
    ) +
    theme_countsynth()
  if (!is.null(vline)) {
    p <- p + ggplot2::geom_vline(
      xintercept = vline, linetype = "dashed", color = "black"
    )
  }
  p
}

#' Relative gap (observed / predicted - 1) for one unit
#'
#' The ribbon is the inverted predictive interval (`obs/upper - 1` to
#' `obs/lower - 1`); pre- and post-treatment segments are drawn separately so
#' neither the line nor the ribbon bridges the treatment date. Rows with
#' non-positive predictions are dropped (the ratio is undefined there).
#'
#' @inheritParams countsynth_unit_fit_plot
#' @export
countsynth_unit_gap_plot <- function(x, unit, group = NULL) {
  q <- as_quantiles(x)
  group <- group %||% if ("total" %in% q$group) "total" else q$group[[1]]
  sub <- q[q$unit == unit & q$group == group, ]
  sub <- sub[!is.na(sub$ypred_mean) & sub$ypred_mean > 0 &
    sub$ypred_lower > 0 & sub$ypred_upper > 0, ]
  if (nrow(sub) == 0) {
    cli::cli_abort("No usable rows for unit {.val {unit}}, group {.val {group}}.")
  }
  if (!"outcome_imputed" %in% names(sub)) sub$outcome_imputed <- FALSE
  sub$gap <- sub$outcome / sub$ypred_mean - 1
  sub$gap_lower <- sub$outcome / sub$ypred_upper - 1
  sub$gap_upper <- sub$outcome / sub$ypred_lower - 1
  vline <- first_treated_time(sub, unit)
  sub$phase <- if (is.null(vline)) {
    "pre"
  } else {
    ifelse(sub$time < vline, "pre", "post")
  }
  imputed <- sub[sub$outcome_imputed, ]

  ggplot2::ggplot(sub, ggplot2::aes(x = .data$time, group = .data$phase)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$gap_lower, ymax = .data$gap_upper),
      fill = "#4C72B0", alpha = 0.25
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$gap),
      color = "#C44E52", linewidth = 0.8
    ) +
    ggplot2::geom_point(
      data = imputed, ggplot2::aes(y = .data$gap),
      shape = 1, size = 3, color = "#E6550D", stroke = 1
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    {
      if (!is.null(vline)) {
        ggplot2::geom_vline(
          xintercept = vline, linetype = "dashed", color = "black"
        )
      }
    } +
    ggplot2::scale_y_continuous(
      labels = function(v) sprintf("%+.0f%%", 100 * v)
    ) +
    ggplot2::labs(
      title = sprintf("Relative gap: %s (%s)", unit, group),
      x = "Time", y = "Observed / predicted - 1",
      subtitle = if (nrow(imputed) > 0) {
        "Orange circles: a source unit's suppressed count was imputed from the model."
      } else {
        NULL
      }
    ) +
    theme_countsynth()
}

# Deduplicate a draws frame to one row per (unit, time, group) of observed
# data (the observed columns are constant across draws).
observed_rows <- function(df) {
  dplyr::distinct(df, .data$unit, .data$time, .data$group, .keep_all = TRUE)
}

#' Raw outcome rate (or count) for treated vs control units
#'
#' @param df A `countsynth_draws` frame or standardized long data frame with
#'   `unit`, `time`, `group`, `outcome`, `denominator`, `treatment`.
#' @param group Restrict to one group (`NULL` pools all groups).
#' @param rate_multiplier Rates are per this many (default 1000).
#' @param treatment_dates Named character vector of dates for vertical marker
#'   lines (names are labels).
#' @param separate_unit Optional unit drawn as its own series.
#' @param smooth_window Centered rolling-mean window (`NULL` = raw values).
#' @param plot_type `"rate"` or `"count"`. `"rate"` falls back to `"count"`
#'   when the run has no denominator: there is no exposure to divide by, and
#'   plotting `count * rate_multiplier` under a "Rate per 1,000" axis would
#'   just be a mislabeled count.
#' @export
countsynth_raw_rate_plot <- function(df, group = NULL, rate_multiplier = 1000,
                                treatment_dates = NULL, separate_unit = NULL,
                                smooth_window = NULL, plot_type = "rate",
                                denominator_label = "denominator") {
  checkmate::assert_choice(plot_type, c("rate", "count"))
  if (!draws_has_denominator(df)) {
    plot_type <- "count"
  }
  df <- observed_rows(df)
  if (!"denominator" %in% names(df)) {
    df$denominator <- 1
  }
  if (!is.null(group)) {
    df <- df[df$group == group, ]
  }
  treated_units <- identify_treated_units(df)
  df$treatment_group <- ifelse(
    !is.null(separate_unit) & df$unit %in% (separate_unit %||% ""),
    separate_unit %||% "",
    ifelse(df$unit %in% treated_units, "Treated", "Control")
  )

  agg <- df |>
    dplyr::group_by(.data$treatment_group, .data$time) |>
    dplyr::summarise(
      outcome = sum(.data$outcome, na.rm = TRUE),
      denominator = sum(.data$denominator, na.rm = TRUE),
      .groups = "drop"
    )
  agg$y_value <- if (plot_type == "count") {
    agg$outcome
  } else {
    agg$outcome / agg$denominator * rate_multiplier
  }
  ylabel <- if (plot_type == "count") {
    "Count"
  } else {
    format_rate_label(rate_multiplier, denominator_label, time_unit = "none")
  }

  if (!is.null(smooth_window) && smooth_window > 1) {
    agg <- agg |>
      dplyr::arrange(.data$treatment_group, .data$time) |>
      dplyr::group_by(.data$treatment_group) |>
      dplyr::mutate(
        y_smooth = stats::filter(
          .data$y_value, rep(1 / smooth_window, smooth_window), sides = 2
        ) |> as.numeric()
      ) |>
      dplyr::ungroup()
    agg$y_smooth[is.na(agg$y_smooth)] <- agg$y_value[is.na(agg$y_smooth)]
  } else {
    agg$y_smooth <- agg$y_value
  }

  levels <- c("Control", "Treated", separate_unit)
  agg$treatment_group <- factor(agg$treatment_group, levels = levels)
  colors <- c(Treated = PALETTE_TREATED, Control = PALETTE_CONTROL)
  if (!is.null(separate_unit)) {
    colors[separate_unit] <- PALETTE_SEPARATE
  }

  p <- ggplot2::ggplot(
    agg,
    ggplot2::aes(
      x = .data$time, y = .data$y_smooth,
      color = .data$treatment_group, group = .data$treatment_group
    )
  ) +
    ggplot2::geom_line(
      ggplot2::aes(
        linewidth = .data$treatment_group != "Control",
        alpha = .data$treatment_group != "Control"
      )
    ) +
    ggplot2::scale_linewidth_manual(
      values = c(`TRUE` = 0.9, `FALSE` = 0.6), guide = "none"
    ) +
    ggplot2::scale_alpha_manual(
      values = c(`TRUE` = 1, `FALSE` = 0.7), guide = "none"
    ) +
    ggplot2::scale_color_manual(values = colors, name = NULL)
  if (!is.null(smooth_window) && smooth_window > 1) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(y = .data$y_value),
      alpha = 0.2, size = 0.8
    )
  }
  if (!is.null(treatment_dates)) {
    marker_colors <- c("#FF8C00", "#DC143C", "#9400D3", "#228B22")
    vlines <- data.frame(
      date = as.Date(unname(treatment_dates)),
      label = names(treatment_dates) %||% as.character(treatment_dates),
      color = marker_colors[(seq_along(treatment_dates) - 1) %% 4 + 1]
    )
    p <- p + ggplot2::geom_vline(
      data = vlines,
      ggplot2::aes(xintercept = .data$date),
      color = vlines$color, linetype = "dashed", alpha = 0.8
    )
  }
  group_label <- if (!is.null(group)) sprintf(" (%s)", group) else ""
  p +
    ggplot2::labs(
      title = sprintf(
        "%s by Treatment Group%s",
        if (plot_type == "rate") "Rate" else "Count", group_label
      ),
      x = "Time", y = paste0(ylabel, group_label)
    ) +
    theme_countsynth()
}

#' Treated vs control rate comparison, faceted by group
#'
#' @inheritParams countsynth_raw_rate_plot
#' @param groups Groups to facet over (`NULL` = all).
#' @export
countsynth_group_comparison_plot <- function(df, groups = NULL,
                                        rate_multiplier = 1000,
                                        treatment_dates = NULL,
                                        plot_type = "rate",
                                        denominator_label = "denominator") {
  checkmate::assert_choice(plot_type, c("rate", "count"))
  if (!draws_has_denominator(df)) {
    plot_type <- "count"
  }
  df <- observed_rows(df)
  if (!"denominator" %in% names(df)) {
    df$denominator <- 1
  }
  groups <- groups %||% unique(df$group)
  df <- df[df$group %in% groups, ]
  treated_units <- identify_treated_units(df)
  df$treatment_group <- ifelse(
    df$unit %in% treated_units, "Treated", "Control"
  )

  agg <- df |>
    dplyr::group_by(.data$group, .data$treatment_group, .data$time) |>
    dplyr::summarise(
      outcome = sum(.data$outcome, na.rm = TRUE),
      denominator = sum(.data$denominator, na.rm = TRUE),
      .groups = "drop"
    )
  agg$y_value <- if (plot_type == "count") {
    agg$outcome
  } else {
    agg$outcome / agg$denominator * rate_multiplier
  }
  ylabel <- if (plot_type == "count") {
    "Count"
  } else {
    format_rate_label(rate_multiplier, denominator_label, time_unit = "none")
  }

  p <- ggplot2::ggplot(
    agg,
    ggplot2::aes(
      x = .data$time, y = .data$y_value, color = .data$treatment_group
    )
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_color_manual(
      values = c(Treated = PALETTE_TREATED, Control = PALETTE_CONTROL),
      name = NULL
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$group),
      ncol = min(2, length(groups)), scales = "free_y"
    )
  if (!is.null(treatment_dates)) {
    p <- p + ggplot2::geom_vline(
      xintercept = as.Date(unname(treatment_dates)),
      linetype = "dashed", alpha = 0.8
    )
  }
  p +
    ggplot2::labs(
      title = sprintf(
        "%s Comparison by Group", if (plot_type == "rate") "Rate" else "Count"
      ),
      x = "Time", y = ylabel
    ) +
    theme_countsynth()
}
