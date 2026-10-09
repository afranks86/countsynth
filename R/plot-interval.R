# Causal-effect interval plot (port of plots.make_interval_plot). One row per
# unit with nested 67% + 95% credible intervals and a median point, sorted by
# median effect. Aggregate/pooled units (see add_aggregate_units()) are not
# comparable to the individual units they pool over -- they share every draw
# with them -- so `separate_units` sets them into their own facet band,
# visually split off from the rest rather than sorted in among them.

#' Per-draw causal effect per unit (and optional color group)
#'
#' Port of `compute_draw_effect`: with `method = "mu"`,
#' `treated = sum(exp(mu_treated))`, `untreated = sum(exp(mu))`, and rates
#' divide by the summed person-time `sum(denominator * years)`; with
#' `method = "pred"`, observed vs posterior-predictive count totals.
#'
#' With `estimand = "ratio"` (the report default), `treated_rate` and
#' `untreated_rate` share the same `denom_val` and `rate_normalizer`, so both
#' cancel algebraically: the plotted percent change is numerically identical
#' to `100 * (treated/untreated - 1)` computed directly on raw counts, and is
#' the same *number* whether or not a denominator is configured at all. Only
#' `estimand = "diff"` actually depends on the denominator's value (a genuine
#' rate difference).
#'
#' That numeric equivalence is not an interpretive one, though -- see
#' [countsynth_summary_table()] for the full explanation, in short: with a
#' denominator, the offset it contributes to `mu`/`mu_treated` is identical
#' on both sides of the ratio, so the percent change is an estimate of the
#' percent change in the **rate** (holding the measured exposure fixed).
#' Without one, it is only an estimate of the percent change in the raw
#' **count**, which cannot be attributed to a rate change versus a change in
#' some unmeasured denominator. Either way, `mu`/`mu_treated` already have
#' `log(denominator)` baked in as a fixed, observed offset (see
#' `stan_data_joint()`), so `untreated`/`ypred_rate` are the counterfactual
#' conditional on that denominator holding -- not a counterfactual where the
#' denominator is also allowed to vary.
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
#' @param draws A `countsynth_draws` frame.
#' @param units Units to include (`NULL` = every treated unit).
#' @param categories Groups to include (`NULL` = all).
#' @param estimand `"diff"` (rate difference) or `"ratio"` (percent change
#'   under `method = "mu"`; rate ratio under `method = "pred"`).
#' @param method `"mu"` (model log-rate) or `"pred"` (posterior-predictive
#'   counts).
#' @param rate_normalizer Rates are per this many units of exposure.
#' @param color_group Optional column used to color/dodge points within a
#'   row (defaults to `"group"` when more than one group is present).
#' @param separate_units Units to split into their own band at the top of the
#'   plot, above a gap, instead of being sorted in with the rest. Defaults to
#'   the frame's own aggregate units (see [aggregate_unit_names()]); pass
#'   `character()` to rank everything together instead. Units not present in
#'   `draws` are ignored. Each band is still sorted by median effect.
#' @return A ggplot object.
#' @export
countsynth_interval_plot <- function(draws, units = NULL, categories = NULL,
                                estimand = c("ratio", "diff"),
                                method = c("mu", "pred"),
                                rate_normalizer = 1000,
                                denominator_label = "denominator",
                                denominator_time_unit = "year",
                                rate_label = NULL,
                                color_group = NULL,
                                separate_units = NULL) {
  estimand <- match.arg(estimand)
  method <- match.arg(method)
  checkmate::assert_character(
    separate_units,
    any.missing = FALSE, null.ok = TRUE
  )
  separate_units <- separate_units %||% aggregate_unit_names(draws)
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
  df$years <- time_weight_per_row(df, denominator_time_unit)
  # countsynth_draws() always emits a denominator column (1 everywhere when none
  # was configured -- see build_cell_table()), but a hand-built frame may not.
  if (!"denominator" %in% names(df)) {
    df$denominator <- 1
  }

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

  # Facet with free + proportional y so each band shows only its own units
  # and keeps one row's worth of height per unit.
  split_units <- intersect(separate_units, levels(plot_df$unit))
  faceted <- length(split_units) > 0
  if (faceted) {
    plot_df$.band <- factor(
      ifelse(plot_df$unit %in% split_units, "separate", "units"),
      levels = c("separate", "units")
    )
  }

  ref <- if (estimand == "ratio" && method == "pred") 1 else 0
  xlab <- if (estimand == "ratio") {
    if (method == "mu") "Percent Change (%)" else "Rate Ratio"
  } else {
    sprintf(
      "Rate Difference (%s)",
      format_rate_label(
        rate_normalizer, denominator_label, denominator_time_unit,
        prefix = "per", rate_label = rate_label
      )
    )
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
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = .data$lower_95, xmax = .data$upper_95),
      linewidth = 0.7, alpha = 0.4, position = dodge
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = .data$lower_67, xmax = .data$upper_67),
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
    theme_countsynth()
  if (faceted) {
    p <- p +
      ggplot2::facet_grid(
        rows = ggplot2::vars(.data$.band),
        scales = "free_y", space = "free_y"
      ) +
      ggplot2::theme(
        strip.text.y = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        panel.spacing.y = ggplot2::unit(0.5, "lines")
      )
  }
  if (is.null(color_group)) {
    p <- p + ggplot2::scale_color_manual(values = "#4C72B0", guide = "none")
  }
  p
}
