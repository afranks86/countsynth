# Tidy posterior draws of the treatment-effect regression coefficients.
# Fixed-effect rows carry term = "(fixed)" and level = NA; random-effect rows
# carry the bar's grouping label (e.g. "group", "group:unit") and the level,
# with `value` holding the level's DEVIATION from the fixed surface (the
# non-centered z times its sampled scale). The full coefficient for a level
# is fixed + deviation; countsynth_te_regression_plot() composes that.

# Assemble the tidy coefficient frame from raw draw matrices. beta_mat is
# (n_draws x P) or NULL; scale_mat (n_draws x Qtot) and z_mat (n_draws x Utot)
# or NULL when the design has no random terms.
te_coef_tidy <- function(design, beta_mat, scale_mat, z_mat,
                         chain, iteration) {
  n_draws <- length(chain)
  base <- tibble::tibble(
    .draw = seq_len(n_draws),
    .chain = as.integer(chain),
    .iteration = as.integer(iteration)
  )
  blocks <- list()
  add_block <- function(term, predictor, level, values) {
    blocks[[length(blocks) + 1]] <<- tibble::tibble(
      base,
      term = term, predictor = predictor, level = level,
      value = as.numeric(values)
    )
  }
  for (p in seq_len(if (is.null(beta_mat)) 0 else ncol(beta_mat))) {
    add_block("(fixed)", design$x_names[[p]], NA_character_, beta_mat[, p])
  }
  uo <- 0L
  qo <- 0L
  for (t in design$terms) {
    n_lev <- length(t$levels)
    for (q in seq_along(t$predictors)) {
      sc <- scale_mat[, qo + q]
      for (l in seq_len(n_lev)) {
        add_block(
          t$label, t$predictors[[q]], t$levels[[l]], z_mat[, uo + l] * sc
        )
      }
      uo <- uo + n_lev
    }
    qo <- qo + length(t$predictors)
  }
  out <- dplyr::bind_rows(blocks)
  class(out) <- c("countsynth_te_draws", class(out))
  out
}

# Extract the regression draw matrices from a CmdStanMCMC fit; NULL entries
# where the design has no such block.
te_coef_matrices <- function(fit, design) {
  list(
    beta = if (ncol(design$X) > 0) extract_var_matrix(fit, "te_beta") else NULL,
    scale = if (length(design$terms) > 0) {
      extract_var_matrix(fit, "te_re_scale")
    } else {
      NULL
    },
    z = if (length(design$terms) > 0) {
      extract_var_matrix(fit, "te_re_z")
    } else {
      NULL
    }
  )
}

#' Posterior draws of the treatment-effect regression coefficients
#'
#' Available when the fit was configured with a
#' `treatment_effects` formula (see [countsynth_te_opts()]). Coefficients are on
#' the standardized covariate scale used in the sampler when
#' `standardize = TRUE`; [countsynth_te_regression_plot()] reports functions on
#' the original scale.
#'
#' @param x A `countsynth_fit` or `countsynth_cut_fit` object.
#' @param ... Unused.
#' @return A tibble with one row per draw and coefficient: `.draw`, `.chain`,
#'   `.iteration`, `term` (`"(fixed)"` or the random-effect grouping label),
#'   `predictor`, `level` (`NA` for fixed rows), and `value` (fixed rows:
#'   the coefficient; random rows: the level's deviation from the fixed
#'   surface). Cut fits add the provenance columns `cut_component`,
#'   `stage1_draw`, `stage1_chain`, `stage1_iteration`.
#' @export
countsynth_te_draws <- function(x, ...) {
  UseMethod("countsynth_te_draws")
}

stop_no_te_design <- function() {
  cli::cli_abort(
    "This fit has no treatment-effect regression; configure
     {.code countsynth_model_opts(treatment_effects = countsynth_te_opts(formula =
     ...))}."
  )
}

#' @export
countsynth_te_draws.countsynth_fit <- function(x, ...) {
  if (is.null(x$te_design)) {
    stop_no_te_design()
  }
  ci <- chain_iteration_vectors(x$fit)
  m <- te_coef_matrices(x$fit, x$te_design)
  te_coef_tidy(x$te_design, m$beta, m$scale, m$z, ci$chain, ci$iteration)
}

#' @export
countsynth_te_draws.countsynth_cut_fit <- function(x, ...) {
  if (is.null(x$te_draws)) {
    stop_no_te_design()
  }
  x$te_draws
}

#' Posterior summary table of the treatment-effect coefficients
#'
#' One row per coefficient with the posterior mean, median, equal-tailed 95%
#' interval, and the posterior probability that the coefficient is positive
#' (the usual two-sided read is `2 * min(p_positive, 1 - p_positive)`).
#' Random-effect rows report the level's deviation from the fixed surface.
#'
#' @param x A `countsynth_fit` or `countsynth_cut_fit` fit with a `treatment_effects`
#'   formula.
#' @return A tibble keyed by `term`, `predictor`, `level`.
#' @export
countsynth_te_coef_table <- function(x) {
  te <- countsynth_te_draws(x)
  te |>
    dplyr::group_by(.data$term, .data$predictor, .data$level) |>
    dplyr::summarise(
      mean = mean(.data$value),
      median = stats::median(.data$value),
      lower_95 = stats::quantile(.data$value, 0.025, names = FALSE),
      upper_95 = stats::quantile(.data$value, 0.975, names = FALSE),
      p_positive = mean(.data$value > 0),
      .groups = "drop"
    )
}

#' @export
print.countsynth_te_draws <- function(x, ...) {
  cli::cli_h1("countsynth treatment-effect regression draws")
  cli::cli_li("{length(unique(x$.draw))} draw{?s}")
  fixed <- unique(x$predictor[x$term == "(fixed)"])
  if (length(fixed) > 0) {
    cli::cli_li("fixed effects: {.val {fixed}}")
  }
  for (tm in setdiff(unique(x$term), "(fixed)")) {
    sub <- x[x$term == tm, ]
    cli::cli_li(
      "({toString(unique(sub$predictor))} | {tm}):
       {length(unique(sub$level))} level{?s}"
    )
  }
  NextMethod()
}
