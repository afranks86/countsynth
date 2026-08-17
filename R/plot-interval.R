# Causal-effect interval plot (port of plots.make_interval_plot). One row per
# unit with nested 67% + 95% credible intervals and a median point, sorted by
# median effect.

#' Per-draw causal effect for one (unit[, color group]) cell
#'
#' Port of `compute_draw_effect`: with `method = "mu"`,
#' `treated = sum(exp(mu_treated))`, `untreated = sum(exp(mu))`, and rates
#' divide by the summed person-time `sum(denominator * years)`; with
#' `method = "pred"`, observed vs posterior-predictive count totals.
#' @keywords internal
compute_draw_effects <- function(df, estimand, method, rate_normalizer,
                                 agg_cols) {
  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(agg_cols))) |>
    dplyr::summarise(
      causal_effect = if (method == "mu") {
        treated <- sum(exp(.data$mu_treated))
        untreated <- sum(exp(.data$mu))
        denom_val <- sum(.data$denominator * .data$years)
        treated_rate <- treated / denom_val * rate_normalizer
        untreated_rate <- untreated / denom_val * rate_normalizer
        if (estimand == "diff") {
          treated_rate - untreated_rate
        } else if (untreated_rate > 0) {
          100 * (treated_rate / untreated_rate - 1)
        } else {
          0
        }
      } else {
        outcome_rate <- (sum(.data$outcome) / mean(.data$years)) /
          (mean(.data$denominator) / rate_normalizer)
        ypred_rate <- (sum(.data$ypred) / mean(.data$years)) /
          (mean(.data$denominator) / rate_normalizer)
        if (estimand == "diff") {
          outcome_rate - ypred_rate
        } else if (ypred_rate > 0) {
          outcome_rate / ypred_rate
        } else {
          1.0
        }
      },
      .groups = "drop"
    )
}

#' Causal-effect interval plot
#'
#' Post-treatment effect per unit with nested credible intervals: a thin 95%
#' segment, a thick 67% segment, and a point at the posterior median; units
#' are ordered by median effect.
#'
#' @param draws A `bpnmf_draws` frame.
#' @param units Units to include (`NULL` = every treated unit).
#' @param categories Groups to include (`NULL` = all).
#' @param estimand `"diff"` (rate difference) or `"ratio"` (percent change
#'   under `method = "mu"`; rate ratio under `method = "pred"`).
#' @param method `"mu"` (model log-rate) or `"pred"` (posterior-predictive
#'   counts).
#' @param rate_normalizer Rates are per this many person-years.
#' @param color_group Optional column used to color/dodge points within a
#'   row (defaults to `"group"` when more than one group is present).
#' @return A ggplot object.
#' @export
bpnmf_interval_plot <- function(draws, units = NULL, categories = NULL,
                                estimand = c("ratio", "diff"),
                                method = c("mu", "pred"),
                                rate_normalizer = 1000,
                                color_group = NULL) {
  estimand <- match.arg(estimand)
  method <- match.arg(method)
  df <- draws
  if (!is.null(categories)) {
    df <- df[df$group %in% categories, ]
  }
  df <- df[!is.na(df$treatment) & df$treatment == 1, ]
  if (!is.null(units)) {
    df <- df[df$unit %in% units, ]
  }
  if (nrow(df) == 0) {
    cli::cli_abort("No post-treatment rows to plot.")
  }
  df$years <- years_per_row(df)

  multi_group <- length(unique(df$group)) > 1
  color_group <- color_group %||% if (multi_group) "group" else NULL
  agg_cols <- unique(c(".draw", "unit", color_group))

  effects <- compute_draw_effects(df, estimand, method, rate_normalizer, agg_cols)

  q <- function(x, p) stats::quantile(x, p, names = FALSE)
  by_cols <- unique(c("unit", color_group))
  plot_df <- effects |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by_cols))) |>
    dplyr::summarise(
      median = stats::median(.data$causal_effect),
      lower_95 = q(.data$causal_effect, 0.025),
      upper_95 = q(.data$causal_effect, 0.975),
      lower_67 = q(.data$causal_effect, 0.165),
      upper_67 = q(.data$causal_effect, 0.835),
      .groups = "drop"
    )

  unit_order <- plot_df |>
    dplyr::group_by(.data$unit) |>
    dplyr::summarise(m = stats::median(.data$median), .groups = "drop") |>
    dplyr::arrange(.data$m)
  plot_df$unit <- factor(plot_df$unit, levels = unit_order$unit)

  ref <- if (estimand == "ratio" && method == "pred") 1 else 0
  xlab <- if (estimand == "ratio") {
    if (method == "mu") "Percent Change (%)" else "Rate Ratio"
  } else {
    sprintf("Rate Difference (per %s person-years)", format(rate_normalizer, big.mark = ","))
  }

  aes_base <- if (!is.null(color_group)) {
    ggplot2::aes(y = .data$unit, color = .data[[color_group]])
  } else {
    ggplot2::aes(y = .data$unit)
  }
  dodge <- if (!is.null(color_group)) {
    ggplot2::position_dodge(width = 0.5)
  } else {
    ggplot2::position_identity()
  }

  p <- ggplot2::ggplot(plot_df, aes_base) +
    ggplot2::geom_vline(xintercept = ref, linetype = "dashed", color = "grey40") +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$lower_95, xend = .data$upper_95, yend = .data$unit),
      linewidth = 0.7, alpha = 0.4, position = dodge
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$lower_67, xend = .data$upper_67, yend = .data$unit),
      linewidth = 1.8, alpha = 0.9, position = dodge
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$median),
      shape = 21, fill = "white", size = 2.4, stroke = 0.9, position = dodge
    ) +
    ggplot2::labs(
      title = "Post-treatment effect by unit",
      subtitle = "Thick segment: 67% CI. Thin segment: 95% CI. Point: posterior median.",
      x = xlab, y = NULL
    ) +
    theme_bpnmf()
  if (is.null(color_group)) {
    p <- p + ggplot2::scale_color_manual(values = "#4C72B0", guide = "none")
  }
  p
}
