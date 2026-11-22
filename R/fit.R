# Joint / baseline model fitting. `method = "sample"` runs NUTS; `method =
# "variational"` runs Stan's ADVI, which cut mode can use for stage 1 to trade
# an exact posterior for a fast approximate one.

#' Fit the joint (or baseline) countsynth model
#'
#' Runs NUTS (or, with `method = "variational"`, Stan's ADVI) on the joint
#' model via cmdstanr. With `config$model$model_treated = FALSE` this fits the
#' untreated baseline model (the same model cut mode uses for stage 1).
#'
#' @param data A `countsynth_data` object from [countsynth_data()].
#' @param rank Factorization rank. Defaults to the first entry of the type's
#'   `ranks_to_test`.
#' @param config The [countsynth_config()] object (MCMC and model options are read
#'   from it).
#' @param model_treated Override `config$model$model_treated`.
#' @param gen_ypred Emit the counterfactual posterior predictive (default
#'   `TRUE`).
#' @param init Optional cmdstanr `init` argument passthrough (escape hatch for
#'   difficult initializations).
#' @param show_exceptions Print the sampler's informational exception
#'   messages (default `FALSE`). During early warmup the positivity
#'   constraints (e.g. on `time_fac`) can overflow on a proposed leapfrog
#'   step; Stan rejects the proposal, prints an informational note, and
#'   adapts -- the same rejection NumPyro performs silently. Sporadic
#'   occurrences are harmless, so the notes are suppressed by default; the
#'   convergence gate (R-hat / ESS / divergences) is the real health signal.
#'   Set to `TRUE` to see them, or read them later via `fit$fit$output()`.
#' @param method `"sample"` (NUTS, the default) or `"variational"` (Stan's
#'   ADVI). ADVI returns an approximation, not a posterior sample: it
#'   understates variance, ignores posterior correlation, and has no R-hat /
#'   ESS / divergence diagnostics, so [convergence_gate()] reports
#'   `converged = NA` and [parameter_diagnostics()] errors. It is meant for
#'   fast iteration -- most usefully as cut stage 1, via
#'   `cut$stage1_method = "variational"`.
#' @param variational Named list of ADVI tuning arguments forwarded to
#'   `CmdStanModel$variational()` (`r toString(VARIATIONAL_KEYS)`).
#'   Ignored when `method = "sample"`.
#' @param ... Additional arguments passed to `CmdStanModel$sample()` (or
#'   `$variational()`).
#' @return A `countsynth_fit` object; `fit_method` records which algorithm ran.
#' @export
countsynth_fit <- function(data, rank = NULL, config, model_treated = NULL,
                      gen_ypred = TRUE, init = NULL,
                      show_exceptions = FALSE,
                      method = c("sample", "variational"),
                      variational = NULL, ...) {
  checkmate::assert_class(data, "countsynth_data")
  checkmate::assert_class(config, "countsynth_config")
  method <- match.arg(method)
  type_spec <- config$model$types[[data$type]]
  rank <- rank %||% type_spec$ranks_to_test[[1]]
  checkmate::assert_int(rank, lower = 1)
  model_treated <- model_treated %||% config$model$model_treated

  # The covariate regression only exists in the treated model; the baseline
  # (cut stage-1) fit ignores it.
  te_design <- if (model_treated && !is.null(config$model$treatment_effects)) {
    build_te_design(config$model$treatment_effects, data)
  } else {
    NULL
  }

  sd <- stan_data_joint(
    data,
    rank = rank,
    model_treated = model_treated,
    outcome_distribution = config$model$outcome_distribution,
    nb_disp = config$model$nb_disp,
    sample_disp = config$model$sample_disp,
    # A config built before this option existed gets the new default.
    disp_by = config$model[["disp_by"]] %||% "group",
    disp_prior_rate = config$model[["disp_prior_rate"]] %||% DEFAULT_DISP_PRIOR_RATE,
    adjust_for_missingness = config$model$adjust_for_missingness,
    gen_ypred = gen_ypred,
    te_design = te_design,
    time_fac_shape = gamma_shape_from_pct(
      config$model$factor_variation_pct, DEFAULT_TIME_FAC_SHAPE
    ),
    time_fe_shape = gamma_shape_from_pct(
      config$model$time_level_variation_pct, DEFAULT_TIME_FE_SHAPE
    ),
    rank_shrinkage = config$model$rank_shrinkage,
    # A config built before this option existed has no field: it gets the
    # new default, like any config that does not mention it.
    # [[ ]], not $: on a config built before this field existed, `$` would
    # partial-match time_level_variation_pct.
    time_level = config$model[["time_level"]] %||% "centered"
  )
  if (model_treated && sd$n_exposed == 0) {
    cli::cli_abort(
      "model_treated = TRUE but the panel has no exposed (treated) cells."
    )
  }

  mcmc <- config$mcmc
  # [[ ]]: a config built before `backend` existed has no field.
  model <- countsynth_stan_model(
    "joint", backend = mcmc[["backend"]] %||% "cmdstan"
  )

  if (method == "variational") {
    # ADVI has no chains, warmup, thinning or adapt_delta; the only MCMC
    # option it shares is the seed, so stage 1 runs at the same seed either
    # way and `stage1_variational` supplies the rest.
    args <- list(
      data = sd,
      seed = mcmc$seed,
      refresh = if (mcmc$progress) NULL else 0,
      show_messages = mcmc$progress,
      show_exceptions = show_exceptions
    )
    args <- c(args, check_variational_args(variational))
    run <- model$variational
  } else {
    ch <- resolve_chains(mcmc)
    args <- list(
      data = sd,
      chains = ch$chains,
      parallel_chains = ch$parallel_chains,
      iter_warmup = mcmc$iter_warmup,
      iter_sampling = mcmc$iter_sampling,
      thin = mcmc$thin,
      adapt_delta = mcmc$adapt_delta,
      max_treedepth = mcmc$max_treedepth %||% 10L,
      seed = mcmc$seed,
      refresh = if (mcmc$progress) NULL else 0,
      show_messages = mcmc$progress,
      show_exceptions = show_exceptions,
      # Divergences and treedepth are narrated against the gate's thresholds
      # by diagnostic_context_notes() / gate_failure_bullets(); leaving
      # cmdstanr's own unglossed warnings on reports each one twice. E-BFMI
      # has no countsynth equivalent, so it stays.
      diagnostics = "ebfmi"
    )
    run <- model$sample
  }
  if (is.null(init)) {
    init <- disp_init(sd, mcmc$seed)
    # The start sets disp alone, by design; cmdstanr would say so per chain.
    if (!is.null(init)) withr::local_options(cmdstanr_warn_inits = FALSE)
  }
  if (!is.null(init)) args$init <- init
  args <- args[!vapply(args, is.null, logical(1))]
  # Explicit `...` arguments win over the config-derived defaults.
  dots <- list(...)
  named <- names(dots)
  if (!is.null(named)) args <- args[!(names(args) %in% named[nzchar(named)])]
  args <- c(args, dots)
  fit <- do.call(run, args)
  if (method == "variational") warn_advi_quality(fit)

  new_countsynth_class(
    list(
      fit = fit,
      data = data,
      config = config,
      stan_data = sd,
      rank = as.integer(rank),
      type = data$type,
      model_treated = model_treated,
      te_design = te_design,
      inference_mode = "joint",
      fit_method = method
    ),
    "countsynth_fit"
  )
}

#' Starting values for a sampled NB dispersion
#'
#' Each chain's `disp` starts near the Poisson, at the prior's own scale
#' `(1 / rate)^2` (the prior mean of `sqrt(disp)` is `1 / rate`), jittered by
#' up to a factor of e so the chains still start apart, and moves away only
#' as far as the data push it. Every other parameter keeps Stan's random
#' init. Left to Stan, `disp` would start uniform on (-2, 2) in logit space,
#' i.e. at 0.12-0.88, where a large dispersion can absorb what the factors
#' should explain: on the infant-mortality panel, two of four chains from such
#' starts settled with every state's loadings collapsed onto the shared
#' profile (`unit_weight_sd` ~1e-4) and a dispersion of ~0.04, an lp
#' 1600-2500 below the other two, and never left.
#'
#' @param stan_data The list from [stan_data_joint()].
#' @param seed The MCMC seed; the jitter for chain `c` is drawn at `seed + c`.
#' @return `NULL` unless the dispersion is sampled, else a cmdstanr `init`
#'   function of `chain_id`.
#' @keywords internal
disp_init <- function(stan_data, seed) {
  if (!(stan_data$is_nb == 1 && stan_data$sample_disp == 1)) {
    return(NULL)
  }
  n <- stan_data$n_disp
  scale <- 1 / stan_data$disp_prior_rate^2
  function(chain_id) {
    jitter <- withr::with_seed((seed %||% 0L) + chain_id, stats::runif(n, -1, 1))
    # disp is bounded by 1, which only a rate below ~1.7 could reach, and by
    # disp_floor; starts stay 10x above the floor, which binds only for a
    # rate above ~2e4.
    start <- pmin(pmax(scale * exp(jitter), 10 * stan_data$disp_floor), 0.5)
    list(disp = as.array(start))
  }
}

#' Validate the ADVI tuning list against the settable knobs
#'
#' The seed is deliberately not settable: stage 1 runs at `mcmc$seed`
#' whichever algorithm it uses, so a run stays reproducible from one field.
#' @keywords internal
check_variational_args <- function(variational) {
  if (is.null(variational)) {
    return(list())
  }
  checkmate::assert_list(variational, names = "unique")
  unknown <- setdiff(names(variational), VARIATIONAL_KEYS)
  if (length(unknown) > 0) {
    cli::cli_abort(c(
      "Unknown ADVI {cli::qty(length(unknown))}argument{?s} {.val {unknown}}.",
      i = "Valid arguments: {.val {VARIATIONAL_KEYS}}.",
      i = if ("seed" %in% unknown) {
        "The seed comes from {.field mcmc.seed}."
      }
    ))
  }
  variational
}

# CmdStan reports a poor ADVI approximation on stdout rather than through the
# exit status, so the console output is the only place the warning lives.
# These are its run-specific complaints only -- deliberately NOT the
# "EXPERIMENTAL ALGORITHM ... may be unstable" banner, which CmdStan prints on
# every ADVI run and which would make this warning cry wolf.
ADVI_TROUBLE <- paste(
  "may be poor", "may not have converged", "MAY BE DIVERGING",
  "Maximum number of iterations",
  sep = "|"
)

#' Surface CmdStan's own complaints about an ADVI fit
#' @keywords internal
warn_advi_quality <- function(fit) {
  out <- tryCatch(
    utils::capture.output(fit$output()),
    error = function(e) character()
  )
  hits <- unique(trimws(grep(ADVI_TROUBLE, out, value = TRUE)))
  # CmdStan's text is data, not a glue template.
  hits <- gsub("}", "}}", gsub("{", "{{", hits, fixed = TRUE), fixed = TRUE)
  if (length(hits) > 0) {
    cli::cli_warn(c(
      "CmdStan flagged the ADVI approximation:",
      stats::setNames(hits, rep("*", length(hits))),
      i = "Re-fit with {.code method = \"sample\"} before trusting these
           estimates."
    ))
  }
  invisible(fit)
}
