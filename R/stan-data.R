# Construct the data lists passed to the Stan models. All cell masks are
# converted to 1-based integer index vectors in the canonical row-major
# (k, d, n) flat order (see flat_idx() in utils.R) so that Stan never touches
# booleans and the exposed-cell ordering matches numpy's reshape(-1) — that
# ordering defines the meaning of treatment_kt_z[e] and cut-mode provenance.

#' Build the Stan data list for the joint / baseline model
#'
#' @param data A `bpnmf_data` object.
#' @param rank Factorization rank.
#' @param model_treated Include the treatment block (joint model); `FALSE`
#'   gives the baseline (cut stage-1) model.
#' @param outcome_distribution `"NB"` or `"Poisson"`.
#' @param nb_disp Fixed NB dispersion.
#' @param sample_disp Sample per-unit dispersion.
#' @param adjust_for_missingness Integrate over censored small counts.
#' @param gen_ypred Emit the counterfactual posterior predictive in
#'   `generated quantities`.
#' @return A named list for `cmdstanr`'s `data` argument, plus attributes
#'   `dims` (K, D, N) and `exposed` (the exposed-cell subscripts).
#' @keywords internal
stan_data_joint <- function(data, rank, model_treated = TRUE,
                            outcome_distribution = "NB", nb_disp = 1e-4,
                            sample_disp = FALSE,
                            adjust_for_missingness = TRUE,
                            gen_ypred = TRUE) {
  K <- length(data$groups)
  D <- length(data$units)
  N <- length(data$times)
  KDN <- K * D * N

  y_flat <- flatten_kdn(data$Y)
  denom_flat <- flatten_kdn(data$denominators)
  control_flat <- flatten_kdn(data$control_idx_array)
  missing_flat <- flatten_kdn(data$missing_idx_array)

  # Direct-likelihood cells: not missing, and (all cells if model_treated,
  # else control cells only). Mirrors joint.py's mask logic.
  obs_mask <- !missing_flat & (model_treated | control_flat)
  obs_cell <- which(obs_mask)

  # Exposed cells in canonical flat order.
  exp_cell <- which(!control_flat)
  exp_sub <- kdn_from_flat(exp_cell, D, N)

  # Censoring adjustment set: all cells when model_treated (matching joint.py
  # passing an all-True mask), control-only otherwise.
  adj_flat <- if (model_treated) rep(TRUE, KDN) else control_flat
  cens_cell <- which(missing_flat & adj_flat)
  notcens_cell <- which(!missing_flat & adj_flat)

  unit_of <- function(cell) kdn_from_flat(cell, D, N)$d

  sd <- list(
    K = K, D = D, N = N, R = as.integer(rank), KDN = KDN,
    log_denom = log(denom_flat),
    cell_unit = as.array(unit_of(seq_len(KDN))),
    n_obs = length(obs_cell),
    obs_cell = as.array(obs_cell),
    y = as.array(as.integer(round(y_flat[obs_cell]))),
    obs_unit = as.array(unit_of(obs_cell)),
    n_exposed = length(exp_cell),
    exp_k = as.array(exp_sub$k),
    exp_d = as.array(exp_sub$d),
    exp_cell = as.array(exp_cell),
    n_cens = length(cens_cell),
    cens_cell = as.array(cens_cell),
    cens_unit = as.array(unit_of(cens_cell)),
    n_notcens = length(notcens_cell),
    notcens_cell = as.array(notcens_cell),
    notcens_unit = as.array(unit_of(notcens_cell)),
    model_treated = as.integer(model_treated),
    is_nb = as.integer(outcome_distribution == "NB"),
    sample_disp = as.integer(sample_disp),
    adjust_missing = as.integer(adjust_for_missingness),
    nb_disp = nb_disp,
    gen_ypred = as.integer(gen_ypred)
  )
  attr(sd, "dims") <- c(K = K, D = D, N = N)
  attr(sd, "exposed") <- exp_sub
  sd
}

#' Build the Stan data list for the cut stage-2 model
#'
#' @param data A `bpnmf_data` object.
#' @param mu_ctrl_flat One stage-1 draw of the baseline log-rate surface, as a
#'   flat vector in canonical row-major order (length `K*D*N`).
#' @param phi_unit Matched per-unit NB concentration for the same stage-1 draw
#'   (`NULL` for Poisson).
#' @param outcome_distribution,adjust_for_missingness See [stan_data_joint()].
#' @keywords internal
stan_data_stage2 <- function(data, mu_ctrl_flat, phi_unit = NULL,
                             outcome_distribution = "NB",
                             adjust_for_missingness = TRUE) {
  K <- length(data$groups)
  D <- length(data$units)
  N <- length(data$times)
  KDN <- K * D * N
  stopifnot(length(mu_ctrl_flat) == KDN)

  y_flat <- flatten_kdn(data$Y)
  control_flat <- flatten_kdn(data$control_idx_array)
  missing_flat <- flatten_kdn(data$missing_idx_array)

  exp_cell <- which(!control_flat)
  exp_sub <- kdn_from_flat(exp_cell, D, N)

  # Stage-2 likelihood/censoring subsets are positions WITHIN the exposed-cell
  # list (1..n_exposed), mirroring cut_treatment.py's exposed-only masks.
  exp_missing <- missing_flat[exp_cell]
  obs_e <- which(!exp_missing)
  cens_e <- which(exp_missing)

  is_nb <- outcome_distribution == "NB"
  if (is_nb && is.null(phi_unit)) {
    cli::cli_abort("NB stage-2 model requires a matched {.field phi_unit} vector.")
  }

  sd <- list(
    K = K, D = D, KDN = KDN,
    mu_ctrl = as.numeric(mu_ctrl_flat),
    n_exposed = length(exp_cell),
    exp_k = as.array(exp_sub$k),
    exp_d = as.array(exp_sub$d),
    exp_cell = as.array(exp_cell),
    n_obs = length(obs_e),
    obs_e = as.array(obs_e),
    y = as.array(as.integer(round(y_flat[exp_cell][obs_e]))),
    n_cens = length(cens_e),
    cens_e = as.array(cens_e),
    is_nb = as.integer(is_nb),
    adjust_missing = as.integer(adjust_for_missingness),
    phi_unit = if (is_nb) as.array(as.numeric(phi_unit)) else numeric(0)
  )
  attr(sd, "dims") <- c(K = K, D = D, N = N)
  attr(sd, "exposed") <- exp_sub
  sd
}

#' Error early when a cut run would have an empty likelihood in either stage
#' @keywords internal
validate_cut_data <- function(data) {
  control_flat <- flatten_kdn(data$control_idx_array)
  missing_flat <- flatten_kdn(data$missing_idx_array)
  if (!any(control_flat & !missing_flat)) {
    cli::cli_abort(
      "Cut stage 1 has an empty likelihood: no non-missing control cells."
    )
  }
  if (!any(!control_flat & !missing_flat)) {
    cli::cli_abort(
      "Cut stage 2 has an empty likelihood: no non-missing exposed cells."
    )
  }
  invisible(NULL)
}
