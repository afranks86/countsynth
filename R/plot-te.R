# Figures for the treatment-effect regression: the posterior of the fitted
# regression function over a covariate (curve + nested ribbons), and a forest
# plot of the coefficients. Both consume the tidy coefficient frame from
# countsynth_te_draws() plus the design's parsed formula, so joint and cut fits
# plot identically.

# Wide (n_draws x n_coefficient) matrix keyed by "term\rpredictor\rlevel".
te_coef_key <- function(term, predictor, level) {
  paste(term, predictor, ifelse(is.na(level), "", level), sep = "\r")
}

te_wide_matrix <- function(te_draws) {
  key <- te_coef_key(te_draws$term, te_draws$predictor, te_draws$level)
  ud <- sort(unique(te_draws$.draw))
  uk <- unique(key)
  m <- matrix(
    NA_real_, length(ud), length(uk),
    dimnames = list(NULL, uk)
  )
  m[cbind(match(te_draws$.draw, ud), match(key, uk))] <- te_draws$value
  m
}

# Design matrices for a raw (unstandardized) prediction frame, rebuilt with
# the same expressions and standardization as build_te_design().
te_matrices_for <- function(design, raw_newdata) {
  frame <- raw_newdata
  for (v in names(design$centers)) {
    frame[[v]] <- (frame[[v]] - design$centers[[v]]) / design$scales[[v]]
  }
  n <- nrow(frame)
  X <- if (is.null(design$fixed_expr)) {
    matrix(0, nrow = n, ncol = 0)
  } else {
    stats::model.matrix(
      stats::as.formula(call("~", design$fixed_expr), env = baseenv()), frame
    )
  }
  Z <- lapply(design$bars, function(b) {
    stats::model.matrix(
      stats::as.formula(call("~", b[[2]]), env = baseenv()), frame
    )
  })
  list(X = unname(X), Z = lapply(Z, unname))
}

# One-row frame with numerics at their median and factors at their modal
# level; `at` overrides individual variables.
te_reference_row <- function(design, at = NULL) {
  ref <- design$raw_frame[1, , drop = FALSE]
  for (nm in names(ref)) {
    col <- design$raw_frame[[nm]]
    if (is.factor(col)) {
      counts <- tabulate(as.integer(col), nlevels(col))
      ref[[nm]][1] <- levels(col)[which.max(counts)]
    } else if (is.numeric(col)) {
      ref[[nm]][1] <- stats::median(col)
    } else if (is.logical(col)) {
      ref[[nm]][1] <- stats::median(col) > 0.5
    }
  }
  for (nm in names(at)) {
    if (!nm %in% names(ref)) {
      cli::cli_abort("{.field at} names unknown variable {.val {nm}}.")
    }
    ref[[nm]][1] <- at[[nm]]
  }
  ref
}

te_predictor_grid <- function(col, n_grid) {
  u <- sort(unique(col))
  # Discrete-looking covariates (event_time, period counters) read better on
  # their own values than on a linspace through them.
  if (length(u) <= n_grid && all(u == round(u))) {
    return(u)
  }
  seq(min(col), max(col), length.out = n_grid)
}

te_default_predictor <- function(design) {
  vars <- unique(c(
    all.vars(design$fixed_expr %||% quote(0)),
    unlist(lapply(design$bars, function(b) all.vars(b[[2]])))
  ))
  is_num <- vapply(
    vars, function(v) is.numeric(design$raw_frame[[v]]), logical(1)
  )
  num <- vars[is_num]
  if (length(num) == 0) {
    cli::cli_abort(
      "No continuous predictor in the treatment-effect formula; pass
       {.field predictor} explicitly or use {.fn countsynth_te_coef_plot}."
    )
  }
  num[[1]]
}

te_design_of <- function(x) {
  d <- x$te_design
  if (is.null(d)) {
    stop_no_te_design()
  }
  d
}

#' Posterior of the fitted treatment-effect regression function
#'
#' Evaluates the regression surface over a grid of one covariate, holding the
#' other model variables at a reference value (median / modal level, override
#' with `at`), and summarizes the posterior as a median line with nested 67%
#' and 95% credible bands. With `by` naming a random-effect grouping factor,
#' one curve per level is drawn (the fixed surface plus that level's
#' deviation), letting subgroup regression surfaces be compared directly.
#'
#' The x axis is on the covariate's original scale even when the sampler ran
#' on standardized covariates.
#'
#' @param x A `countsynth_fit` or `countsynth_cut_fit` fit with a `treatment_effects`
#'   formula.
#' @param predictor Covariate on the x axis (defaults to the first continuous
#'   variable in the formula).
#' @param by Grouping label of a random-effect term (e.g. `"group"`,
#'   `"unit"`); `NULL` draws only the population (fixed-effect) surface.
#' @param levels_keep Optional subset of `by` levels to draw.
#' @param at Named list of values other model variables are held at.
#' @param n_grid Maximum number of grid points.
#' @param scale `"log"` plots the additive effect on the log-rate scale (the
#'   model's own scale); `"percent"` plots `100 * (exp(f) - 1)`, the percent
#'   change in the rate.
#' @param bands Draw the 67%/95% credible bands (`FALSE` gives median lines
#'   only, which reads better with many levels).
#' @return A ggplot object.
#' @export
countsynth_te_regression_plot <- function(x, predictor = NULL, by = NULL,
                                     levels_keep = NULL, at = NULL,
                                     n_grid = 60,
                                     scale = c("log", "percent"),
                                     bands = TRUE) {
  scale <- match.arg(scale)
  design <- te_design_of(x)
  te <- countsynth_te_draws(x)
  predictor <- predictor %||% te_default_predictor(design)
  if (!predictor %in% names(design$raw_frame)) {
    cli::cli_abort("{.val {predictor}} is not a variable in the design.")
  }
  if (!is.numeric(design$raw_frame[[predictor]])) {
    cli::cli_abort(
      "{.val {predictor}} is categorical; use {.fn countsynth_te_coef_plot} for
       categorical effects."
    )
  }
  curves <- te_curve_frame(
    design, te, predictor, by, levels_keep, at, n_grid, scale
  )
  facet_by <- if (is.null(by)) NULL else "level"

  ylab <- if (scale == "percent") {
    "Change in rate (%)"
  } else {
    "Treatment effect (log rate)"
  }
  n_lev <- length(unique(curves$level))
  aes_base <- if (is.null(by)) {
    ggplot2::aes(x = .data$x, y = .data$median)
  } else {
    ggplot2::aes(
      x = .data$x, y = .data$median, color = .data$level, fill = .data$level
    )
  }
  p <- ggplot2::ggplot(curves, aes_base) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40")
  if (bands) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lower_95, ymax = .data$upper_95),
        alpha = 0.15, color = NA
      ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lower_67, ymax = .data$upper_67),
        alpha = 0.3, color = NA
      )
  }
  p <- p +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::labs(
      title = "Posterior treatment-effect regression function",
      subtitle = if (is.null(by)) {
        "Population surface. Bands: 67% and 95% credible intervals."
      } else {
        sprintf(
          "By %s (fixed surface + level deviation). Bands: 67%% and 95%% CI.",
          by
        )
      },
      x = predictor, y = ylab, color = by, fill = by
    ) +
    theme_countsynth()
  # Many levels overplot badly; facet instead of relying on the legend.
  if (!is.null(facet_by) && bands && n_lev > 3) {
    p <- p +
      ggplot2::facet_wrap(facet_by) +
      ggplot2::guides(color = "none", fill = "none")
  }
  if (is.null(by)) {
    p <- p +
      ggplot2::scale_color_manual(values = "#4C72B0", guide = "none") +
      ggplot2::scale_fill_manual(values = "#4C72B0", guide = "none")
  }
  p
}

#' Posterior summary of the treatment-effect regression function
#'
#' The data behind [countsynth_te_regression_plot()]: one row per grid point (and
#' level), with the posterior median and the 67% / 95% credible bounds.
#'
#' @inheritParams countsynth_te_regression_plot
#' @return A tibble with `x`, `level`, `median`, `lower_67`, `upper_67`,
#'   `lower_95`, `upper_95`.
#' @export
countsynth_te_regression_summary <- function(x, predictor = NULL, by = NULL,
                                        levels_keep = NULL, at = NULL,
                                        n_grid = 60,
                                        scale = c("log", "percent")) {
  scale <- match.arg(scale)
  design <- te_design_of(x)
  te <- countsynth_te_draws(x)
  predictor <- predictor %||% te_default_predictor(design)
  te_curve_frame(design, te, predictor, by, levels_keep, at, n_grid, scale)
}

# Posterior curves over the predictor grid, for the population surface and
# (when `by` is given) each level of that grouping term.
te_curve_frame <- function(design, te, predictor, by, levels_keep, at,
                           n_grid, scale) {
  grid <- te_predictor_grid(design$raw_frame[[predictor]], n_grid)
  newdata <- te_reference_row(design, at)[rep(1, length(grid)), , drop = FALSE]
  newdata[[predictor]] <- grid

  term_labels <- vapply(design$terms, `[[`, character(1), "label")
  if (!is.null(by) && !by %in% term_labels) {
    cli::cli_abort(
      "{.val {by}} is not a random-effect grouping term; available:
       {.val {term_labels}}."
    )
  }
  W <- te_wide_matrix(te)
  n_draws <- nrow(W)
  mats <- te_matrices_for(design, newdata)

  # With no fixed part (~ 0 + (...|g)) the population surface is flat zero and
  # every level curve is its own deviation.
  fixed <- matrix(0, nrow = n_draws, ncol = length(grid))
  if (length(design$x_names) > 0) {
    beta <- W[, te_coef_key("(fixed)", design$x_names, NA), drop = FALSE]
    fixed <- beta %*% t(mats$X)
  }

  pieces <- list()
  if (is.null(by)) {
    pieces[["population"]] <- fixed
  } else {
    j <- match(by, term_labels)
    t_j <- design$terms[[j]]
    Zj <- mats$Z[[j]]
    lv <- levels_keep %||% t_j$levels
    unknown <- setdiff(lv, t_j$levels)
    if (length(unknown) > 0) {
      cli::cli_abort("Unknown {.val {by}} level{?s} {.val {unknown}}.")
    }
    for (l in lv) {
      u <- W[, te_coef_key(by, t_j$predictors, l), drop = FALSE]
      pieces[[l]] <- fixed + u %*% t(Zj)
    }
  }

  q <- function(m, p) apply(m, 2, stats::quantile, probs = p, names = FALSE)
  out <- dplyr::bind_rows(lapply(names(pieces), function(nm) {
    m <- pieces[[nm]]
    if (scale == "percent") {
      m <- 100 * (exp(m) - 1)
    }
    tibble::tibble(
      x = grid, level = nm,
      median = apply(m, 2, stats::median),
      lower_67 = q(m, 0.165), upper_67 = q(m, 0.835),
      lower_95 = q(m, 0.025), upper_95 = q(m, 0.975)
    )
  }))
  if (!is.null(by)) {
    out$level <- factor(out$level, levels = names(pieces))
  }
  out
}

# Standard figure set for the report: the population regression function for
# each continuous predictor, the same broken out by every grouping term, and
# one coefficient forest plot. Names become the PNG filenames.
te_report_figures <- function(x) {
  design <- te_design_of(x)
  vars <- unique(c(
    all.vars(design$fixed_expr %||% quote(0)),
    unlist(lapply(design$bars, function(b) all.vars(b[[2]])))
  ))
  continuous <- vars[vapply(
    vars, function(v) is.numeric(design$raw_frame[[v]]), logical(1)
  )]
  term_labels <- vapply(design$terms, `[[`, character(1), "label")

  out <- list()
  for (v in continuous) {
    out[[sprintf("te_regression_%s", v)]] <-
      countsynth_te_regression_plot(x, predictor = v)
    for (j in seq_along(design$terms)) {
      t_j <- design$terms[[j]]
      varies_with_v <- v %in% all.vars(design$bars[[j]][[2]])
      # An intercept-only term shifts the curve by level rather than bending
      # it, which is still worth seeing -- but only while the panel count
      # stays readable. Callers can always plot a wider term directly.
      if (!varies_with_v && length(t_j$levels) > 12) {
        next
      }
      slug <- gsub("[^a-z0-9]+", "_", tolower(t_j$label))
      out[[sprintf("te_regression_%s_by_%s", v, slug)]] <-
        countsynth_te_regression_plot(x, predictor = v, by = t_j$label)
    }
  }
  out[["te_coefficients"]] <- countsynth_te_coef_plot(
    x,
    terms = if (length(term_labels) > 0) "all" else NULL
  )
  out
}

#' Forest plot of treatment-effect regression coefficients
#'
#' Nested 67% / 95% credible intervals for the fixed-effect coefficients and,
#' when `terms` selects them, the per-level coefficients of the random-effect
#' terms. Level rows show the total coefficient (fixed surface + deviation)
#' by default, or the deviation alone with `deviation = TRUE`.
#'
#' @param x A `countsynth_fit` or `countsynth_cut_fit` fit with a `treatment_effects`
#'   formula.
#' @param terms Grouping labels to include beyond the fixed effects (`NULL`
#'   = fixed effects only, `"all"` = every term).
#' @param predictors Restrict to these design columns (`NULL` = all).
#' @param deviation Show level deviations from the fixed surface instead of
#'   the total per-level coefficient.
#' @return A ggplot object.
#' @export
countsynth_te_coef_plot <- function(x, terms = NULL, predictors = NULL,
                               deviation = FALSE) {
  design <- te_design_of(x)
  te <- countsynth_te_draws(x)
  term_labels <- vapply(design$terms, `[[`, character(1), "label")
  if (identical(terms, "all")) {
    terms <- term_labels
  }
  unknown <- setdiff(terms, term_labels)
  if (length(unknown) > 0) {
    cli::cli_abort(
      "Unknown random-effect term{?s} {.val {unknown}}; available:
       {.val {term_labels}}."
    )
  }
  W <- te_wide_matrix(te)

  rows <- list()
  add_row <- function(predictor, term, level, v) {
    rows[[length(rows) + 1]] <<- tibble::tibble(
      predictor = predictor, term = term, level = level,
      median = stats::median(v),
      lower_67 = stats::quantile(v, 0.165, names = FALSE),
      upper_67 = stats::quantile(v, 0.835, names = FALSE),
      lower_95 = stats::quantile(v, 0.025, names = FALSE),
      upper_95 = stats::quantile(v, 0.975, names = FALSE)
    )
  }
  keep_pred <- function(p) is.null(predictors) || p %in% predictors
  fixed_of <- function(p) {
    k <- te_coef_key("(fixed)", p, NA)
    if (k %in% colnames(W)) W[, k] else rep(0, nrow(W))
  }
  for (p in design$x_names) {
    if (keep_pred(p)) {
      add_row(p, "(fixed)", "(population)", fixed_of(p))
    }
  }
  for (tl in terms) {
    t_j <- design$terms[[match(tl, term_labels)]]
    for (p in t_j$predictors) {
      if (!keep_pred(p)) {
        next
      }
      base <- if (deviation) rep(0, nrow(W)) else fixed_of(p)
      for (l in t_j$levels) {
        add_row(p, tl, l, base + W[, te_coef_key(tl, p, l)])
      }
    }
  }
  if (length(rows) == 0) {
    cli::cli_abort("No coefficients selected.")
  }
  plot_df <- dplyr::bind_rows(rows)
  plot_df$label <- ifelse(
    plot_df$term == "(fixed)", "population",
    sprintf("%s: %s", plot_df$term, plot_df$level)
  )
  plot_df$label <- factor(plot_df$label, levels = rev(unique(plot_df$label)))

  ggplot2::ggplot(plot_df, ggplot2::aes(y = .data$label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$lower_95, xend = .data$upper_95, yend = .data$label
      ),
      linewidth = 0.7, alpha = 0.4
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$lower_67, xend = .data$upper_67, yend = .data$label
      ),
      linewidth = 1.8, alpha = 0.9
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$median),
      shape = 21, fill = "white", size = 2.4, stroke = 0.9
    ) +
    ggplot2::facet_wrap("predictor", scales = "free_x") +
    ggplot2::labs(
      title = if (deviation) {
        "Treatment-effect coefficient deviations"
      } else {
        "Treatment-effect regression coefficients"
      },
      subtitle = paste(
        "Thick: 67% CI. Thin: 95% CI.",
        if (design$standardize && length(design$centers) > 0) {
          "Continuous covariates are standardized (effect per SD)."
        } else {
          ""
        }
      ),
      x = "Coefficient (log rate)", y = NULL
    ) +
    theme_countsynth()
}
