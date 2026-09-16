# Tidy long draws frame -- the reporting contract every table and figure
# consumes. Parity with results.py: one row per draw x group x unit x time,
# 1-indexed .draw/.chain/.iteration (tidybayes-style), and columns
#   unit, time, group, outcome, denominator, treatment, ypred, mu, mu_treated
# where `mu` is the UNTREATED log-count surface (mu_ctrl) and
# `mu_treated = mu + te` (results.py:126-132). Cut mode adds provenance
# columns cut_component / stage1_draw / stage1_chain / stage1_iteration.

# Per-cell observed metadata in canonical flat-cell order (row-major k,d,n).
build_cell_table <- function(data) {
  K <- length(data$groups)
  D <- length(data$units)
  N <- length(data$times)
  sub <- kdn_from_flat(seq_len(K * D * N), D, N)
  cells <- tibble::tibble(
    group = data$groups[sub$k],
    unit = data$units[sub$d],
    time = data$times[sub$n]
  )
  obs_cols <- intersect(
    c("unit", "time", "group", "outcome", "denominator", "treatment",
      "start_date", "end_date"),
    names(data$df)
  )
  out <- dplyr::left_join(
    cells, data$df[obs_cols],
    by = c("unit", "time", "group")
  )
  # No denominator configured means denominator 1 everywhere (the same
  # convention build_model_arrays() uses for the Stan-facing array). Emit the
  # column either way so every downstream table and figure sees one schema --
  # dropping it instead made every rate-aware consumer fail with "Column
  # `denominator` not found" on an otherwise valid denominator-free run.
  # `has_denominator` (stamped in draws_frame_core()) is what distinguishes
  # the two cases for reporting; the column's presence no longer does.
  if (!"denominator" %in% names(out)) {
    out$denominator <- 1
  }
  out
}

# Core long-frame assembly. mu_mat / ypred_mat are (n_draws x KDN) matrices in
# flat-cell column order; te_mat is (n_draws x n_exposed) or NULL;
# chain/iteration are per-draw integer vectors.
draws_frame_core <- function(mu_mat, te_mat, ypred_mat, chain, iteration,
                             data, exposed_cell) {
  n_draws <- nrow(mu_mat)
  KDN <- ncol(mu_mat)
  cells <- build_cell_table(data)

  mu_treated_mat <- mu_mat
  if (!is.null(te_mat) && length(exposed_cell) > 0) {
    mu_treated_mat[, exposed_cell] <- mu_treated_mat[, exposed_cell] + te_mat
  }

  out <- tibble::tibble(
    .draw = rep(seq_len(n_draws), each = KDN),
    .chain = rep(as.integer(chain), each = KDN),
    .iteration = rep(as.integer(iteration), each = KDN),
    cells[rep(seq_len(KDN), times = n_draws), ],
    ypred = as.vector(t(ypred_mat)),
    mu = as.vector(t(mu_mat)),
    mu_treated = as.vector(t(mu_treated_mat))
  )
  attr(out, "groups") <- data$groups
  attr(out, "units") <- data$units
  attr(out, "times") <- data$times
  attr(out, "has_denominator") <- "denominator" %in% names(data$df)
  class(out) <- c("bpnmf_draws", class(out))
  out
}

# Was a real denominator configured, as opposed to the implicit 1 that
# build_cell_table() fills in? Reporting needs the distinction even though the
# column is always present: with a denominator, "Pct Change" estimates a
# change in the RATE (the measured exposure is held fixed on both sides of the
# ratio); without one it only estimates a change in the raw COUNT. Frames
# built by hand (tests, external callers) carry no flag, so fall back to
# column presence -- what the flag's absence used to mean.
draws_has_denominator <- function(draws) {
  flag <- attr(draws, "has_denominator")
  if (!is.null(flag)) {
    return(isTRUE(flag))
  }
  "denominator" %in% names(draws)
}

# Extract a variable from a CmdStanMCMC fit as a (n_draws x n_elements)
# matrix, rows ordered chain-major (chain 1 iterations, then chain 2, ...).
extract_var_matrix <- function(fit, var) {
  dm <- posterior::as_draws_matrix(fit$draws(variables = var))
  m <- unclass(dm)
  attr(m, "nchains") <- NULL
  m
}

chain_iteration_vectors <- function(fit) {
  # Read both counts from the draws themselves rather than from metadata:
  # Stan retains ceil(iter_sampling / thin), so dividing metadata values would
  # be wrong whenever thin does not divide iter_sampling evenly -- and a
  # variational fit has no `num_chains()` at all, just one stream of draws,
  # which `nchains()` correctly reports as a single chain.
  lp <- fit$draws(variables = "lp__")
  n_chains <- posterior::nchains(lp)
  per_chain <- posterior::niterations(lp)
  list(
    chain = rep(seq_len(n_chains), each = per_chain),
    iteration = rep(seq_len(per_chain), times = n_chains),
    per_chain = per_chain,
    n_chains = n_chains
  )
}

#' Build the tidy posterior draws frame
#'
#' @param x A `bpnmf_fit` (joint fit) or `bpnmf_cut_fit` object.
#' @param ... Unused.
#' @return A `bpnmf_draws` tibble; see the file header for the schema.
#' @export
bpnmf_draws <- function(x, ...) {
  UseMethod("bpnmf_draws")
}

#' @export
bpnmf_draws.bpnmf_fit <- function(x, ...) {
  ci <- chain_iteration_vectors(x$fit)
  mu_mat <- extract_var_matrix(x$fit, "mu_ctrl")
  ypred_mat <- extract_var_matrix(x$fit, "ypred")
  te_mat <- NULL
  exposed_cell <- integer()
  if (x$model_treated && x$stan_data$n_exposed > 0) {
    te_mat <- extract_var_matrix(x$fit, "te")
    exposed_cell <- as.integer(x$stan_data$exp_cell)
  }
  draws_frame_core(
    mu_mat, te_mat, ypred_mat, ci$chain, ci$iteration, x$data, exposed_cell
  )
}

#' @export
bpnmf_draws.bpnmf_cut_fit <- function(x, ...) {
  x$draws
}

#' Untreated posterior-predictive counts for given baseline surfaces
#'
#' Port of `cut.sample_untreated_predictions`: one count draw per (draw,
#' cell) at rate `exp(mu)`, NB2 when `phi_row` is supplied.
#'
#' @param mu_mat `(n_draws x KDN)` matrix of log-count surfaces.
#' @param phi_row Per-cell NB concentration vector of length KDN (recycled
#'   across draws), or `NULL` for Poisson.
#' @keywords internal
sample_untreated_predictions <- function(mu_mat, phi_row = NULL) {
  n <- length(mu_mat)
  rate <- exp(pmin(as.vector(t(mu_mat)), 20.79))
  if (is.null(phi_row)) {
    draws <- stats::rpois(n, rate)
  } else {
    size <- rep(phi_row, times = nrow(mu_mat))
    draws <- stats::rnbinom(n, size = size, mu = rate)
  }
  matrix(draws, nrow = nrow(mu_mat), byrow = TRUE)
}

#' @export
print.bpnmf_draws <- function(x, ...) {
  n_draws <- length(unique(x$.draw))
  cli::cli_h1("bpnmf draws")
  cli::cli_li("{n_draws} draw{?s} x {length(attr(x, 'groups'))} group{?s} x
               {length(attr(x, 'units'))} unit{?s} x {length(attr(x, 'times'))} time{?s}")
  if ("cut_component" %in% names(x)) {
    cli::cli_li("cut posterior: {length(unique(x$cut_component))} component{?s}")
  }
  NextMethod()
}
