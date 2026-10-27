# Mixed predictive check of the treatment-effect model (Marshall &
# Spiegelhalter 2007).
#
# The rest of the PPC suite checks the counterfactual on control-period
# cells. This one checks the treated side, which the reported effects are
# built from (`treated = sum(exp(mu_treated))`, not the observed counts), and
# in particular whether the treatment-effect hierarchy shrinks each unit
# toward the pooled surface by the right amount.
#
# A plain posterior predictive cannot answer that. A unit's own coefficients
# were fit to its exposed cells, so its replicates track the observation
# whatever the shrinkage, and the p-values bunch at 0.5. Here, per draw, every
# coefficient that belongs to one unit alone -- (1 | unit), (1 | group:unit),
# a random slope by unit, the iid exposed-cell term -- is redrawn from its
# population distribution at that draw's sampled scale. The fixed surface,
# terms whose levels span units, the scales themselves and mu_ctrl stay at
# their posterior draw. The question becomes: is this unit's observed
# post-period total plausible as a draw from the population the model pools
# it toward?
#
# The p-values are not exactly uniform under a correct model. The scales,
# the fixed surface and, in joint mode, mu_ctrl were still fit with the
# unit's data, so an outlying unit
# widens the spread it is judged against, and a high-leverage unit pulls the
# fixed slopes toward itself; both push its p-value toward 0.5. Excess mass
# near 0 and 1 is evidence of over-shrinkage; a clean calibration plot is
# weaker reassurance than it looks.
#
# Cost is one normal per redrawn coefficient and one count per observed
# exposed cell, per draw -- a fraction of the standard suite, which simulates
# every cell of the panel.

# The treatment-effect pieces the check redraws, one entry per sampled scale:
# `level` maps each exposed cell to its coefficient, `mult` is the design
# column the coefficient multiplies (NULL for the legacy intercepts), and
# `scale_var` / `z_cols` name the sampled scale and the non-centered
# coefficients in the Stan output. `scale_col` picks one element of a vector
# scale; NULL means the variable's only element, which is `x` in one model
# and `x[1]` in the other. A random-effect term is redrawn when each of its
# levels sits inside a single unit on the exposed cells the design was built
# on, so a categorical covariate that happens to be one-to-one with the
# treated units is redrawn too, and (1 | group) never is.
te_unit_components <- function(design, stan_data) {
  n_exp <- as.integer(stan_data$n_exposed)
  cell <- list(
    label = "exposed cell", level = seq_len(n_exp), n_levels = n_exp,
    mult = NULL, scale_var = "treatment_it_scale", scale_col = NULL,
    z_var = "treatment_kt_z",
    z_cols = sprintf("treatment_kt_z[%d]", seq_len(n_exp))
  )
  if (is.null(design)) {
    K <- as.integer(stan_data$K)
    D <- as.integer(stan_data$D)
    kd <- seq_len(K * D)
    return(list(
      list(
        label = "(1 | unit)", level = as.integer(stan_data$exp_d),
        n_levels = D, mult = NULL, scale_var = "treatment_unit_scale",
        scale_col = NULL, z_var = "unit_treatment_effect_z",
        z_cols = sprintf("unit_treatment_effect_z[%d]", seq_len(D))
      ),
      # exp_kd is (d - 1) * K + k, the column-major position of [k, d].
      list(
        label = "(1 | group:unit)", level = as.integer(stan_data$exp_kd),
        n_levels = K * D, mult = NULL, scale_var = "unit_category_scale",
        scale_col = NULL, z_var = "unit_category_te_z",
        z_cols = sprintf(
          "unit_category_te_z[%d,%d]", (kd - 1L) %% K + 1L, (kd - 1L) %/% K + 1L
        )
      ),
      cell
    ))
  }
  comps <- list()
  uo <- 0L
  qo <- 0L
  for (t in design$terms) {
    n_lev <- length(t$levels)
    units_per_level <- tapply(
      design$frame$unit, t$index, function(u) length(unique(u))
    )
    within_unit <- all(units_per_level == 1L)
    for (q in seq_along(t$predictors)) {
      if (within_unit) {
        pred <- t$predictors[[q]]
        comps[[length(comps) + 1L]] <- list(
          label = sprintf(
            "(%s | %s)", if (pred == "(Intercept)") "1" else pred, t$label
          ),
          level = t$index, n_levels = n_lev, mult = t$Z[, q],
          scale_var = "te_re_scale", scale_col = sprintf("te_re_scale[%d]", qo + q),
          z_var = "te_re_z", z_cols = sprintf("te_re_z[%d]", uo + seq_len(n_lev))
        )
      }
      uo <- uo + n_lev
    }
    qo <- qo + length(t$predictors)
  }
  c(comps, list(cell))
}

# Stan variables te_mixed_parts() reads, so a caller can fetch them in the
# same draws() call as anything else it needs: each new call on a fit whose
# CSV has not been read for those variables is another pass over the CSV.
te_mixed_vars <- function(components) {
  unique(c(
    "te",
    vapply(components, `[[`, character(1), "scale_var"),
    vapply(components, `[[`, character(1), "z_var")
  ))
}

# Per draw, the log mean each exposed cell keeps once the redrawn pieces are
# taken out (`base` = mu_ctrl + te minus their posterior contribution), and
# the sampled scale of each redrawn piece (`scale`, one column per
# component). `m` is a draws matrix holding te_mixed_vars(); `mu_exp` is
# mu_ctrl on the exposed cells, either a draws x cells matrix or one row
# recycled over draws (a cut component's frozen baseline).
te_mixed_parts <- function(m, components, mu_exp) {
  n_draws <- nrow(m)
  n_exp <- length(components[[length(components)]]$level)
  base <- m[, sprintf("te[%d]", seq_len(n_exp)), drop = FALSE]
  base <- if (is.matrix(mu_exp)) {
    base + mu_exp
  } else {
    base + rep(mu_exp, each = n_draws)
  }
  scale <- matrix(NA_real_, n_draws, length(components))
  for (i in seq_along(components)) {
    comp <- components[[i]]
    col <- comp$scale_col %||% intersect(
      c(comp$scale_var, paste0(comp$scale_var, "[1]")), colnames(m)
    )[[1]]
    s <- m[, col]
    # s recycles down the columns, so row i is scaled by draw i's scale.
    contrib <- s * m[, comp$z_cols, drop = FALSE][, comp$level, drop = FALSE]
    if (!is.null(comp$mult)) {
      contrib <- contrib * rep(comp$mult, each = n_draws)
    }
    base <- base - contrib
    scale[, i] <- s
  }
  dimnames(base) <- NULL
  list(base = base, scale = scale)
}

# The exposed cells' unit, group and observed outcome (NA where suppressed
# or missing), in the canonical flat order every exposed-cell vector uses.
te_mixed_cells <- function(data, stan_data) {
  cells <- build_cell_table(data)[as.integer(stan_data$exp_cell), ]
  list(
    unit = cells$unit, group = cells$group, outcome = cells$outcome,
    exp_d = as.integer(stan_data$exp_d)
  )
}

#' Inputs of the mixed predictive check
#'
#' Everything [countsynth_ppc_te_mixed()] needs, from either kind of fit:
#' `base` and `scale` (see `te_mixed_parts()`), `phi` (draws x units NB
#' concentration, `NULL` for Poisson), `components`, and the exposed cells'
#' `unit`, `group`, `outcome` and `exp_d`.
#' @param x A `countsynth_fit` or `countsynth_cut_fit`.
#' @keywords internal
te_mixed_source <- function(x) {
  UseMethod("te_mixed_source")
}

#' @export
te_mixed_source.countsynth_fit <- function(x) {
  sd <- x$stan_data
  if (!isTRUE(x$model_treated) || sd$n_exposed == 0) {
    cli::cli_abort(
      "The mixed predictive check needs a fit with {.code model_treated =
       TRUE} and exposed cells."
    )
  }
  components <- te_unit_components(x$te_design, sd)
  mu_names <- sprintf("mu_ctrl[%d]", as.integer(sd$exp_cell))
  is_nb <- isTRUE(sd$is_nb == 1)
  m <- extract_var_matrix(
    x$fit, c(mu_names, te_mixed_vars(components), if (is_nb) "phi_unit")
  )
  parts <- te_mixed_parts(m, components, m[, mu_names, drop = FALSE])
  parts$phi <- if (is_nb) {
    m[, sprintf("phi_unit[%d]", seq_len(sd$D)), drop = FALSE]
  } else {
    NULL
  }
  c(parts, list(components = components), te_mixed_cells(x$data, sd))
}

#' @export
te_mixed_source.countsynth_cut_fit <- function(x) {
  if (is.null(x$te_mixed)) {
    cli::cli_abort(c(
      "This cut fit carries no mixed predictive inputs.",
      i = "Cut fits save them at fit time; refit to run the check."
    ))
  }
  c(
    x$te_mixed[c("base", "scale", "phi", "components")],
    te_mixed_cells(x$data, x$stage1$stan_data)
  )
}

# Redraw, simulate and total over the exposed cells in `keep`. Returns the
# per-draw observed-minus-replicated totals by (unit, group) and their
# p-values, P(observed < replicated) with half weight on ties, the same
# direction as the rest of the suite.
te_mixed_check <- function(src, keep) {
  n_draws <- nrow(src$base)
  eta <- src$base[, keep, drop = FALSE]
  for (i in seq_along(src$components)) {
    comp <- src$components[[i]]
    u <- matrix(stats::rnorm(n_draws * comp$n_levels), n_draws, comp$n_levels)
    add <- src$scale[, i] * u[, comp$level[keep], drop = FALSE]
    if (!is.null(comp$mult)) {
      add <- add * rep(comp$mult[keep], each = n_draws)
    }
    eta <- eta + add
  }
  rate <- exp(pmin(eta, 20.79))
  yrep <- if (is.null(src$phi)) {
    stats::rpois(length(rate), rate)
  } else {
    stats::rnbinom(
      length(rate),
      size = src$phi[, src$exp_d[keep], drop = FALSE], mu = rate
    )
  }
  yrep <- matrix(yrep, n_draws, length(keep))

  cell_key <- paste(src$unit[keep], src$group[keep], sep = "\r")
  keys <- unique(cell_key)
  idx <- match(cell_key, keys)
  cells <- tibble::tibble(
    unit = src$unit[keep][!duplicated(idx)],
    group = src$group[keep][!duplicated(idx)]
  )
  # rowsum sorts by idx, which is 1..n_cells in `cells` order.
  t_rep <- rowsum(t(yrep), idx)
  t_obs <- rowsum(src$outcome[keep], idx)[, 1]
  diff <- t_obs - t_rep
  cells$pval <- rowMeans(diff < 0) + 0.5 * rowMeans(diff == 0)
  stats_df <- tibble::tibble(
    unit = rep(cells$unit, times = n_draws),
    group = rep(cells$group, times = n_draws),
    .draw = rep(seq_len(n_draws), each = nrow(cells)),
    diff_in_total = as.vector(diff)
  )
  list(stats = stats_df, pvals = cells)
}

te_mixed_caption <- function(redrawn, extra = NULL) {
  lines <- c(
    sprintf("Redrawn per unit from the population: %s.", toString(redrawn)),
    extra
  )
  paste(strwrap(lines, width = 110), collapse = "\n")
}

# Sorted p-values against uniform quantiles, one panel per group. With a
# dozen or so treated units a histogram of p-values is mostly empty bins;
# the quantile plot shows the same departure without binning. The units at
# either extreme are named so they can be looked up. The names run
# vertically, from the point into the panel: over-shrinkage stacks several
# units at p = 0 or 1, where horizontal names collide, but each sorted point
# has an x of its own, and nothing else sits directly above the bottom row
# or below the top one.
te_mixed_calibration_plot <- function(pvals, redrawn) {
  pv <- pvals |>
    dplyr::group_by(.data$group) |>
    dplyr::arrange(.data$pval, .by_group = TRUE) |>
    dplyr::mutate(uniform = (dplyr::row_number() - 0.5) / dplyr::n()) |>
    dplyr::ungroup()
  flagged <- pv[pv$pval < 0.05 | pv$pval > 0.95, ]
  flagged$hjust <- ifelse(flagged$pval < 0.5, -0.15, 1.15)
  p <- ggplot2::ggplot(pv, ggplot2::aes(x = .data$uniform, y = .data$pval)) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0, color = "red", linetype = "dashed"
    ) +
    ggplot2::geom_point(color = "steelblue", size = 2) +
    ggplot2::geom_text(
      data = flagged,
      ggplot2::aes(label = .data$unit, hjust = .data$hjust),
      angle = 90, size = 2.5, color = "grey20"
    ) +
    ggplot2::facet_wrap(ggplot2::vars(!!!rlang::syms("group")), ncol = 2) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Mixed Predictive P-values Against Uniform",
      x = "Uniform quantile", y = "Sorted p-value"
    ) +
    theme_countsynth()
  add_plot_caption(p, te_mixed_caption(redrawn, paste(
    "Points below the line on the left and above it on the right mean too",
    "many p-values near 0 and 1: unit effects shrunk harder than the data",
    "support. The check leans toward 0.5 even when the model is right, so a",
    "mild pull toward the middle is expected."
  )))
}

#' PPC: mixed predictive check of the treatment-effect shrinkage
#'
#' Checks the treated side of the model on the exposed cells: whether each
#' treated unit's observed post-period total is plausible as a draw from the
#' population of effects the model pools it toward. For every posterior
#' draw, the treatment-effect coefficients that belong to one unit alone --
#' random effects whose levels each sit inside a single unit, such as
#' `(1 | unit)` and `(1 | group:unit)`, and the iid exposed-cell effect --
#' are redrawn from their population distribution at that draw's sampled
#' scale. The fixed surface, random effects whose levels span units, and
#' `mu_ctrl` stay at their posterior draw. Counts are then simulated on the
#' observed exposed cells and totaled per (unit, group).
#'
#' Read the p-values together, on the calibration plot, rather than one at a
#' time. Over-shrinkage -- unit effects pulled toward the surface harder
#' than the data support -- puts too many of them near 0 and 1. The check
#' is conservative: the scales, the fixed surface and (in joint mode)
#' `mu_ctrl` still saw each unit's data, which pulls its p-value toward 0.5,
#' so a mild excess in the middle is expected even when the model is right.
#'
#' Unlike the rest of the suite this reads the fit rather than the draws
#' frame, since the frame does not carry the coefficients the check redraws.
#' A cut fit saves what the check needs at fit time.
#'
#' @param fit A `countsynth_fit` with `model_treated = TRUE`, or a
#'   `countsynth_cut_fit`.
#' @inheritParams countsynth_ppc_abs
#' @param seed Seed for the redrawn effects and simulated counts; defaults
#'   to the fit's `mcmc$seed`.
#' @return List with `plot` (per-cell histograms of observed minus
#'   replicated post-period totals), `calibration` (sorted p-values against
#'   uniform quantiles), `pvals` (unit x group tibble) and `redrawn` (labels
#'   of the redrawn components).
#' @references Marshall, E. C. and Spiegelhalter, D. J. (2007). Identifying
#'   outliers in Bayesian hierarchical models: a simulation-based approach.
#'   *Bayesian Analysis* 2(2), 409--444.
#' @export
countsynth_ppc_te_mixed <- function(fit, categories = NULL, ppc_units = NULL,
                                    ppc_exclude_units = NULL, seed = NULL) {
  src <- te_mixed_source(fit)
  seed <- seed %||% fit$config$mcmc$seed %||% 1L
  redrawn <- vapply(src$components, `[[`, character(1), "label")
  keep <- !is.na(src$outcome)
  if (!is.null(categories)) keep <- keep & src$group %in% categories
  if (!is.null(ppc_units)) keep <- keep & src$unit %in% ppc_units
  if (!is.null(ppc_exclude_units)) {
    keep <- keep & !src$unit %in% ppc_exclude_units
  }
  if (!any(keep)) {
    label <- "Insufficient data for mixed predictive check\n(no observed exposed cells)"
    return(list(
      plot = ppc_insufficient_data_plot(label),
      calibration = ppc_insufficient_data_plot(label),
      pvals = tibble::tibble(unit = character(), group = character(), pval = numeric()),
      redrawn = redrawn
    ))
  }
  res <- withr::with_seed(seed, te_mixed_check(src, which(keep)))
  list(
    plot = add_plot_caption(
      ppc_histogram(
        res$stats, res$pvals, "diff_in_total",
        "Mixed Predictive Check: Post-Period Total",
        "Observed - Replicated Total", c("unit", "group")
      ),
      te_mixed_caption(redrawn)
    ),
    calibration = te_mixed_calibration_plot(res$pvals, redrawn),
    pvals = res$pvals,
    redrawn = redrawn
  )
}
