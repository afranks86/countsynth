# Joint / baseline model fitting.

#' Fit the joint (or baseline) bpnmf model
#'
#' Runs NUTS on the joint model via cmdstanr. With
#' `config$model$model_treated = FALSE` this fits the untreated baseline
#' model (the same model cut mode uses for stage 1).
#'
#' @param data A `bpnmf_data` object from [bpnmf_data()].
#' @param rank Factorization rank. Defaults to the first entry of the type's
#'   `ranks_to_test`.
#' @param config The [bpnmf_config()] object (MCMC and model options are read
#'   from it).
#' @param model_treated Override `config$model$model_treated`.
#' @param gen_ypred Emit the counterfactual posterior predictive (default
#'   `TRUE`).
#' @param init Optional cmdstanr `init` argument passthrough (escape hatch for
#'   difficult initializations).
#' @param ... Additional arguments passed to `CmdStanModel$sample()`.
#' @return A `bpnmf_fit` object.
#' @export
bpnmf_fit <- function(data, rank = NULL, config, model_treated = NULL,
                      gen_ypred = TRUE, init = NULL, ...) {
  checkmate::assert_class(data, "bpnmf_data")
  checkmate::assert_class(config, "bpnmf_config")
  type_spec <- config$model$types[[data$type]]
  rank <- rank %||% type_spec$ranks_to_test[[1]]
  checkmate::assert_int(rank, lower = 1)
  model_treated <- model_treated %||% config$model$model_treated

  sd <- stan_data_joint(
    data,
    rank = rank,
    model_treated = model_treated,
    outcome_distribution = config$model$outcome_distribution,
    nb_disp = config$model$nb_disp,
    sample_disp = config$model$sample_disp,
    adjust_for_missingness = config$model$adjust_for_missingness,
    gen_ypred = gen_ypred
  )
  if (model_treated && sd$n_exposed == 0) {
    cli::cli_abort(
      "model_treated = TRUE but the panel has no exposed (treated) cells."
    )
  }

  mcmc <- config$mcmc
  ch <- resolve_chains(mcmc)
  model <- bpnmf_stan_model("joint")

  args <- list(
    data = sd,
    chains = ch$chains,
    parallel_chains = ch$parallel_chains,
    iter_warmup = mcmc$iter_warmup,
    iter_sampling = mcmc$iter_sampling,
    thin = mcmc$thin,
    adapt_delta = mcmc$adapt_delta,
    seed = mcmc$seed,
    refresh = if (mcmc$progress) NULL else 0,
    show_messages = mcmc$progress
  )
  if (!is.null(init)) args$init <- init
  args <- c(args[!vapply(args, is.null, logical(1))], list(...))
  fit <- do.call(model$sample, args)

  new_bpnmf_class(
    list(
      fit = fit,
      data = data,
      config = config,
      stan_data = sd,
      rank = as.integer(rank),
      type = data$type,
      model_treated = model_treated,
      inference_mode = "joint"
    ),
    "bpnmf_fit"
  )
}
