# Config constructors + validation. Mirrors the pydantic v2 schema in the
# Python package (config.py): every constructor validates immediately, unknown
# fields are impossible by construction here and rejected by the YAML loader
# (config-yaml.R), and cross-field rules (outcomes XOR prefixes, cut requires
# model_treated, stage2 seed authority, aggregate-unit single selector) are
# enforced at build time so a bad config never reaches the sampler.

# Canonical figure names output$figures may select. Must stay in sync with the
# figure functions wired up in reports.R (the Python side derives this from
# PLOT_REGISTRY for the same no-drift reason). "te_regression" is the one
# exception: it needs the fit object (not just draws), so bpnmf_run() renders
# it and bpnmf_report() ignores it.
FIGURE_NAMES <- c(
  "unit_fit", "unit_gap", "raw_rate", "interval", "group_comparison", "ppc",
  "te_regression", "rank_shrinkage"
)

# Historical shapes of the Gamma(shape, shape) priors, carried over from the
# fertility application. Kept exact as defaults so a config that does not
# mention either prior reproduces earlier runs.
DEFAULT_TIME_FAC_SHAPE <- 20   # unit-specific temporal basis
DEFAULT_TIME_FE_SHAPE <- 1     # common time level (Gamma(1, 1))

#' Gamma shape implied by an expected multiplicative swing
#'
#' Both temporal priors are `Gamma(shape, shape)`: mean 1, entering the
#' log-rate as `log(x)`, with `Var[log x] = trigamma(shape)` exactly. Asking
#' for a swing of `pct` percent means `sd(log x) = log1p(pct / 100)`, so the
#' shape solves `trigamma(shape) = log1p(pct / 100)^2`.
#'
#' Inverted numerically rather than via the usual `1 / sd^2` approximation:
#' that approximation is good above shape 10 or so (20.1 vs 20.6 at 25
#' percent) but is 28 percent off at shape 1, which is exactly where the
#' `time_fe` prior sits.
#'
#' @param pct Expected multiplicative variation, in percent. `NULL` returns
#'   `default_shape` unchanged.
#' @param default_shape Shape to return when `pct` is `NULL`.
#' @return The Gamma shape (and rate; they are equal).
#' @export
gamma_shape_from_pct <- function(pct, default_shape = DEFAULT_TIME_FAC_SHAPE) {
  if (is.null(pct)) {
    return(default_shape)
  }
  checkmate::assert_number(pct, lower = 1e-6, finite = TRUE)
  target <- log1p(pct / 100)^2
  # trigamma decreases from Inf at 0 to 0 at Inf, so the root is unique.
  stats::uniroot(
    function(s) trigamma(s) - target,
    interval = c(1e-8, 1e9), tol = .Machine$double.eps^0.5
  )$root
}

AGGREGATION_PERIODS <- c("monthly", "bimonthly", "quarterly", "yearly")

# Tuning knobs a variational (ADVI) stage-1 fit may set. cmdstanr's
# `variational()` also takes data / seed / refresh / message flags, but those
# are owned by the fitter (the seed comes from mcmc.seed) so they are not
# settable here.
VARIATIONAL_KEYS <- c(
  "algorithm", "iter", "grad_samples", "elbo_samples", "eta",
  "adapt_engaged", "adapt_iter", "tol_rel_obj", "eval_elbo",
  "output_samples", "draws"
)

# Every option group is its own constructor, which keeps 80-odd arguments from
# collapsing into one unreadable signature -- but it also meant an R script had
# to build objects the YAML loader lets you write as plain nested maps. These
# coercions close that gap: a named list is passed through the very same
# constructor, so nothing skips validation, and `bpnmf_type(...)` and
# `list(...)` are interchangeable wherever one is expected.
coerce_bpnmf <- function(x, ctor_name, class, field) {
  if (is.null(x) || inherits(x, class)) {
    return(x)
  }
  ctor <- get(ctor_name, envir = asNamespace("bpnmf"))
  if (!is.list(x) || (length(x) > 0 && is.null(names(x)))) {
    cli::cli_abort(
      "{.field {field}} must be a {.fn {ctor_name}} object or a named list of
       its arguments."
    )
  }
  unknown <- setdiff(names(x), names(formals(ctor)))
  if (length(unknown) > 0) {
    cli::cli_abort(c(
      "Unknown {cli::qty(length(unknown))}name{?s} {.val {unknown}} in
       {.field {field}}.",
      i = "Valid arguments: {.val {setdiff(names(formals(ctor)), '...')}}."
    ))
  }
  do.call(ctor, x)
}

# The same, for a list of them (model types, outcomes, aggregate units).
coerce_bpnmf_list <- function(x, ctor_name, class, field) {
  if (is.null(x)) {
    return(x)
  }
  if (!is.list(x)) {
    cli::cli_abort("{.field {field}} must be a list.")
  }
  nms <- names(x) %||% rep("", length(x))
  out <- lapply(seq_along(x), function(i) {
    label <- if (nzchar(nms[i])) {
      sprintf("%s$%s", field, nms[i])
    } else {
      sprintf("%s[[%d]]", field, i)
    }
    coerce_bpnmf(x[[i]], ctor_name, class, label)
  })
  names(out) <- names(x)
  out
}

new_bpnmf_class <- function(x, class) {
  structure(x, class = c(class, "list"))
}

#' Define one explicit outcome column
#'
#' @param outcome_col Name of the count column in the input CSV.
#' @param label Group label the outcome is reported under.
#' @param denominator_col Optional exposure/denominator column.
#' @return A `bpnmf_outcome` spec.
#' @export
bpnmf_outcome <- function(outcome_col, label, denominator_col = NULL) {
  checkmate::assert_string(outcome_col, min.chars = 1)
  checkmate::assert_string(label, min.chars = 1)
  checkmate::assert_string(denominator_col, min.chars = 1, null.ok = TRUE)
  new_bpnmf_class(
    list(
      outcome_col = outcome_col, label = label,
      denominator_col = denominator_col
    ),
    "bpnmf_outcome"
  )
}

#' Derive outcomes from shared column prefixes
#'
#' @param outcome_prefix Prefix of the count columns (e.g. `"births_"`);
#'   labels are the column names with the prefix stripped.
#' @param denominator_prefix Optional prefix of matching denominator columns.
#' @param include Optional character vector restricting which labels are kept.
#' @return A `bpnmf_prefixes` spec.
#' @export
bpnmf_prefixes <- function(outcome_prefix, denominator_prefix = NULL,
                           include = NULL) {
  checkmate::assert_string(outcome_prefix, min.chars = 1)
  checkmate::assert_string(denominator_prefix, min.chars = 1, null.ok = TRUE)
  checkmate::assert_character(
    include,
    min.len = 1, any.missing = FALSE, null.ok = TRUE
  )
  new_bpnmf_class(
    list(
      outcome_prefix = outcome_prefix,
      denominator_prefix = denominator_prefix,
      include = include
    ),
    "bpnmf_prefixes"
  )
}

#' Column roles and outcome definitions for the input data
#'
#' Exactly one of `outcomes` or `outcomes_from_prefixes` must be given.
#'
#' @param unit_col,time_col,treatment_col Column names for the panel unit,
#'   time period, and 0/1 treatment indicator.
#' @param outcomes A list of [bpnmf_outcome()] specs.
#' @param outcomes_from_prefixes A [bpnmf_prefixes()] spec.
#' @return A `bpnmf_schema` object.
#' @export
bpnmf_schema <- function(unit_col, time_col, treatment_col,
                         outcomes = NULL, outcomes_from_prefixes = NULL) {
  checkmate::assert_string(unit_col, min.chars = 1)
  checkmate::assert_string(time_col, min.chars = 1)
  checkmate::assert_string(treatment_col, min.chars = 1)
  if (is.null(outcomes) == is.null(outcomes_from_prefixes)) {
    cli::cli_abort(
      "schema must define exactly one of {.field outcomes} or {.field outcomes_from_prefixes}."
    )
  }
  if (!is.null(outcomes)) {
    outcomes <- coerce_bpnmf_list(outcomes, "bpnmf_outcome", "bpnmf_outcome", "outcomes")
    checkmate::assert_list(outcomes, min.len = 1, types = "bpnmf_outcome")
  }
  if (!is.null(outcomes_from_prefixes)) {
    outcomes_from_prefixes <- coerce_bpnmf(
      outcomes_from_prefixes, "bpnmf_prefixes", "bpnmf_prefixes",
      "outcomes_from_prefixes"
    )
    checkmate::assert_class(outcomes_from_prefixes, "bpnmf_prefixes")
  }
  new_bpnmf_class(
    list(
      unit_col = unit_col, time_col = time_col, treatment_col = treatment_col,
      outcomes = outcomes, outcomes_from_prefixes = outcomes_from_prefixes
    ),
    "bpnmf_schema"
  )
}

#' Temporal aggregation settings
#'
#' Two ways to coarsen the time axis, one of which must be chosen when
#' `enabled` is `TRUE`:
#'
#' * `period` bins by the **calendar**: rows are grouped by calendar month,
#'   two-month block, quarter, or year, so bins land on calendar boundaries
#'   (Q1 is always Jan-Mar) no matter when the panel starts. This assumes the
#'   input is monthly or finer.
#' * `n_periods` bins by **position**: every `n_periods` consecutive time
#'   points in the panel are combined, whatever the underlying resolution --
#'   daily, weekly, monthly, or an arbitrary integer index. Blocks are cut
#'   from the panel's sorted distinct times, so all units stay aligned. Use
#'   this when the data is not monthly, or when you want blocks anchored to
#'   the panel's own start rather than to the calendar.
#'
#' @param enabled Aggregate input rows into coarser periods?
#' @param period Calendar bin: one of `"monthly"`, `"bimonthly"`,
#'   `"quarterly"`, `"yearly"`. Mutually exclusive with `n_periods`.
#' @param n_periods Number of consecutive input periods to combine into one.
#'   Mutually exclusive with `period`. A trailing block with fewer than
#'   `n_periods` points is kept (and warned about); its shorter exposure is
#'   carried in `start_date`/`end_date`, so exposure-weighted rates stay
#'   correct.
#' @export
bpnmf_time_aggregation <- function(enabled = FALSE, period = NULL,
                                   n_periods = NULL) {
  checkmate::assert_flag(enabled)
  checkmate::assert_choice(period, AGGREGATION_PERIODS, null.ok = TRUE)
  checkmate::assert_int(n_periods, lower = 1, null.ok = TRUE)
  if (!is.null(period) && !is.null(n_periods)) {
    cli::cli_abort(
      "Set only one of {.field period} (calendar bins) and
       {.field n_periods} (combine N consecutive periods)."
    )
  }
  # Keeping the historical default here rather than in the signature means
  # `period`/`n_periods` can stay NULL-by-default and still be exclusive.
  if (enabled && is.null(period) && is.null(n_periods)) {
    period <- "bimonthly"
  }
  new_bpnmf_class(
    list(enabled = enabled, period = period, n_periods = n_periods),
    "bpnmf_time_aggregation"
  )
}

#' One model type: a set of groups fit independently at one or more ranks
#'
#' @param groups Character vector of group labels this type models. May
#'   include the synthetic `"total"` group (see `total_from` / `total_all`).
#' @param ranks_to_test Integer vector of factorization ranks to fit.
#' @param total_from Labels summed to build the synthetic `"total"` group.
#' @param total_all If `TRUE`, `"total"` sums every resolved outcome label.
#' @param exclude_units Units dropped before fitting this type.
#' @export
bpnmf_type <- function(groups, ranks_to_test, total_from = NULL,
                       total_all = FALSE, exclude_units = NULL) {
  checkmate::assert_character(groups, min.len = 1, any.missing = FALSE)
  checkmate::assert_integerish(
    ranks_to_test,
    lower = 1, min.len = 1, any.missing = FALSE
  )
  checkmate::assert_character(
    total_from,
    min.len = 1, any.missing = FALSE, null.ok = TRUE
  )
  checkmate::assert_flag(total_all)
  checkmate::assert_character(exclude_units, any.missing = FALSE, null.ok = TRUE)
  new_bpnmf_class(
    list(
      groups = groups, ranks_to_test = as.integer(ranks_to_test),
      total_from = total_from, total_all = total_all,
      exclude_units = exclude_units
    ),
    "bpnmf_type"
  )
}

#' Rank-shrinkage options: shrink unused factor components
#'
#' Turns the flat per-unit component weights into a two-level hierarchy over
#' the `R` candidate temporal curves -- the finite-dimensional form of a
#' hierarchical Dirichlet process. Each group gets a shared
#' component-popularity profile `group_weight[k]`, built by truncated
#' stick-breaking (`stick ~ Beta(1, group_weight_mass)`), and each unit's
#' weights are that profile perturbed multiplicatively and renormalized:
#' `softmax(log group_weight[k] + unit_weight_sd * z)` with `z` standard
#' normal.
#'
#' Two things follow. Components past the rank a group actually uses get
#' near-zero profile weight for *every* unit in it, so an over-specified rank
#' stops changing the fit -- pick a generous `R` once instead of sweeping.
#' And the prior amplitude of the low-rank term stops depending on `R`: under
#' the default `Dirichlet(1, ..., 1)`,
#' `E[sum_r unit_weight^2] = 2 / (R + 1)`, so the variance of the low-rank
#' term falls off like `1 / R` and [bpnmf_model_opts()]'s
#' `factor_variation_pct` only means what it says at a fixed rank.
#'
#' The fit reports `eff_rank[k] = 1 / sum_r group_weight[k, r]^2` -- the
#' effective number of components group `k` occupies. If its posterior sits
#' well below `R` the truncation was generous enough; if it presses against
#' `R`, raise the rank and refit. See [bpnmf_component_weight_plot()] and
#' [bpnmf_rank_hyper_plot()].
#'
#' `NULL` in [bpnmf_model_opts()] leaves the weights iid uniform on the
#' simplex; the Stan parameters this adds are then zero-size, so a fit is
#' unchanged draw for draw.
#'
#' @param group_mass_prior `Gamma(shape, rate)` prior on
#'   `group_weight_mass`, the DP mass behind the stick-breaking profile
#'   (`stick ~ Beta(1, mass)`). Smaller mass means a sparser profile: the
#'   implied effective rank is about `1 + mass`, so the default `c(2, 1)`
#'   (mean 2) centres it near three components and reaches about seven.
#' @param unit_sd_prior Half-normal prior scale for `unit_weight_sd`, the
#'   multiplicative spread of a unit's component loadings about its group's
#'   profile (a standard deviation on the log scale, so `0.5` is roughly a
#'   1.6-fold departure and `1` a 2.7-fold one). Near zero pins every unit to
#'   the shared profile, large lets each unit load on its own few components;
#'   it is sampled rather than fixed because how much cross-unit sharing the
#'   panel supports is exactly what is not known in advance. Default `1`.
#'
#'   Logistic-normal rather than the textbook HDP's
#'   `Dirichlet(conc * profile)` for a computational reason, not a modelling
#'   one: the Dirichlet form measured 5.4x slower on the bundled 51-unit
#'   fertility panel (1193s against 223s for 300 iterations at rank 8), at
#'   the same leapfrog count and with no divergences either way, because a
#'   Dirichlet whose concentration is itself sampled costs an `lgamma` and a
#'   `digamma` per component per unit on every gradient. It also needs a
#'   concentration floor to keep unused components off the simplex boundary,
#'   and no floor is both large enough to do that and small enough to keep
#'   the rank-invariance. A concentration `c` corresponds to `unit_sd_prior`
#'   near `sqrt(2 / c)`.
#' @return A `bpnmf_rank_shrinkage_opts` object.
#' @export
bpnmf_rank_shrinkage_opts <- function(group_mass_prior = c(2, 1),
                                      unit_sd_prior = 1) {
  checkmate::assert_numeric(
    group_mass_prior,
    lower = .Machine$double.xmin, finite = TRUE, any.missing = FALSE, len = 2,
    .var.name = "group_mass_prior"
  )
  checkmate::assert_number(
    unit_sd_prior, lower = .Machine$double.xmin, finite = TRUE
  )
  new_bpnmf_class(
    list(
      group_mass_prior = as.numeric(group_mass_prior),
      unit_sd_prior = as.numeric(unit_sd_prior)
    ),
    "bpnmf_rank_shrinkage_opts"
  )
}

#' Model options: likelihood family, model types, dispersion, treatment
#'
#' @param outcome_distribution `"NB"` or `"Poisson"`.
#' @param types Named list of [bpnmf_type()] specs; names are the model-type
#'   labels (e.g. `"total"`, `"age"`).
#' @param nb_disp Fixed NB dispersion (concentration is `1/nb_disp`).
#' @param sample_disp Sample per-unit dispersion instead of fixing it.
#' @param adjust_for_missingness Integrate the likelihood over censored small
#'   counts (values 1-9 suppressed in the source data).
#' @param model_treated Include the treatment-effect block.
#' @param inference_mode `NULL` (default joint), `"joint"`, or `"cut"`.
#' @param treatment_effects Optional [bpnmf_te_opts()] object replacing the
#'   default treatment-effect hierarchy with a covariate regression surface
#'   (applies to both joint inference and cut stage 2). `NULL` keeps the
#'   legacy `(1 | group) + (1 | unit) + (1 | group:unit)` model.
#' @param factor_variation_pct Expected multiplicative variation of the
#'   low-rank temporal basis, in percent, setting the `Gamma(shape, shape)`
#'   prior on `time_fac` via [gamma_shape_from_pct()]. `NULL` (default)
#'   keeps the historical shape of 20, which corresponds to about 25 percent.
#'   Raise it for panels whose within-unit seasonality or volatility is larger
#'   than that; lower it for smooth series.
#' @param rank_shrinkage Optional [bpnmf_rank_shrinkage_opts()] object
#'   shrinking unused factor components, so that an over-specified rank stops
#'   changing the fit. `NULL` (default) keeps the flat
#'   `Dirichlet(1, ..., 1)` weights, in which case results do depend on the
#'   rank chosen and should be checked by sweeping `ranks_to_test`.
#' @param time_level_variation_pct Expected multiplicative variation of the
#'   *common* time level (`time_fe`, shared by all units within a group), in
#'   percent. `NULL` (default) keeps the historical `Gamma(1, 1)`, which is
#'   very diffuse: `sd(log time_fe)` is 1.28, a central 95% range of x0.03 to
#'   x3.7. Because `unit_fe_mu` has a flat prior and the likelihood sees only
#'   `log(time_fe) + unit_fe_mu`, this prior is the only thing separating the
#'   two, so tightening it improves identification as well as encoding a
#'   belief.
#' @export
bpnmf_model_opts <- function(outcome_distribution = "NB", types = list(),
                             nb_disp = 1e-4, sample_disp = FALSE,
                             adjust_for_missingness = TRUE,
                             model_treated = TRUE, inference_mode = NULL,
                             treatment_effects = NULL,
                             factor_variation_pct = NULL,
                             time_level_variation_pct = NULL,
                             rank_shrinkage = NULL) {
  checkmate::assert_choice(outcome_distribution, c("NB", "Poisson"))
  types <- coerce_bpnmf_list(types, "bpnmf_type", "bpnmf_type", "types")
  checkmate::assert_list(types, types = "bpnmf_type")
  if (length(types) > 0) {
    checkmate::assert_names(names(types), type = "unique")
  }
  checkmate::assert_number(nb_disp, lower = .Machine$double.xmin)
  checkmate::assert_flag(sample_disp)
  checkmate::assert_flag(adjust_for_missingness)
  checkmate::assert_flag(model_treated)
  checkmate::assert_choice(inference_mode, c("joint", "cut"), null.ok = TRUE)
  if (identical(inference_mode, "cut") && !model_treated) {
    cli::cli_abort(
      "{.field model.inference_mode}='cut' requires {.field model_treated}=TRUE
       (the cut model estimates treatment effects)."
    )
  }
  if (sample_disp && outcome_distribution == "Poisson") {
    cli::cli_abort(
      "{.field model.sample_disp}=TRUE requires {.field outcome_distribution}='NB'."
    )
  }
  checkmate::assert_number(
    factor_variation_pct,
    lower = 1e-6, finite = TRUE, null.ok = TRUE
  )
  checkmate::assert_number(
    time_level_variation_pct,
    lower = 1e-6, finite = TRUE, null.ok = TRUE
  )
  # `TRUE` / `FALSE` are the natural shorthand for "on at the defaults" /
  # "off", and the YAML loader accepts the same, so honour them here too.
  if (is.logical(rank_shrinkage) && length(rank_shrinkage) == 1 &&
    !is.na(rank_shrinkage)) {
    rank_shrinkage <- if (rank_shrinkage) bpnmf_rank_shrinkage_opts() else NULL
  }
  rank_shrinkage <- coerce_bpnmf(
    rank_shrinkage, "bpnmf_rank_shrinkage_opts", "bpnmf_rank_shrinkage_opts",
    "rank_shrinkage"
  )
  checkmate::assert_class(treatment_effects, "bpnmf_te_opts", null.ok = TRUE)
  if (!is.null(treatment_effects) && !is.null(treatment_effects$formula) &&
    !model_treated) {
    cli::cli_abort(
      "{.field model.treatment_effects} requires {.field model_treated}=TRUE."
    )
  }
  new_bpnmf_class(
    list(
      outcome_distribution = outcome_distribution, types = types,
      nb_disp = nb_disp, sample_disp = sample_disp,
      adjust_for_missingness = adjust_for_missingness,
      model_treated = model_treated, inference_mode = inference_mode,
      treatment_effects = treatment_effects,
      factor_variation_pct = factor_variation_pct,
      time_level_variation_pct = time_level_variation_pct,
      rank_shrinkage = rank_shrinkage
    ),
    "bpnmf_model_opts"
  )
}

#' Convergence-gate thresholds
#'
#' A parameter PASSes below `rhat_warn` with ESS at/above `ess_min`; FAILs at/
#' above `rhat_fail` or below `ess_min * ess_fail_fraction`; WARNs in between.
#' ESS is `min(bulk, tail)` per parameter.
#'
#' Divergent transitions are gated separately, as a rate over the retained
#' draws rather than a raw count, so the threshold does not drift as the
#' sampling length changes. `divergence_fail_fraction = 0` restores the
#' zero-divergence rule.
#' @param rhat_warn,rhat_fail,ess_min,ess_fail_fraction Gate thresholds.
#' @param divergence_fail_fraction Largest share of retained transitions that
#'   may be divergent while still counting as converged.
#' @export
bpnmf_convergence <- function(rhat_warn = 1.01, rhat_fail = 1.05,
                              ess_min = 400, ess_fail_fraction = 0.25,
                              divergence_fail_fraction = 0.01) {
  checkmate::assert_number(rhat_warn, lower = 1)
  checkmate::assert_number(rhat_fail, lower = 1)
  checkmate::assert_number(ess_min, lower = 0)
  checkmate::assert_number(ess_fail_fraction, lower = 0, upper = 1)
  checkmate::assert_number(divergence_fail_fraction, lower = 0, upper = 1)
  new_bpnmf_class(
    list(
      rhat_warn = rhat_warn, rhat_fail = rhat_fail,
      ess_min = ess_min, ess_fail_fraction = ess_fail_fraction,
      divergence_fail_fraction = divergence_fail_fraction
    ),
    "bpnmf_convergence"
  )
}

#' MCMC options
#'
#' Names follow cmdstanr conventions; the YAML loader translates the Python
#' names (`num_warmup` -> `iter_warmup`, `num_samples` -> `iter_sampling`,
#' `thinning` -> `thin`, `target_accept` -> `adapt_delta`,
#' `random_seed` -> `seed`).
#'
#' @param auto_parallelism Choose chains/parallelism from available cores.
#' @param max_chains Chain cap under `auto_parallelism`.
#' @param chains,parallel_chains Manual chain settings (used when
#'   `auto_parallelism = FALSE`; both default to 4 / sequential).
#' @param iter_warmup,iter_sampling,thin Warmup / sampling iterations and
#'   retention stride (both implementations run `iter_sampling` post-warmup
#'   iterations and retain `iter_sampling / thin`).
#' @param adapt_delta NUTS target acceptance probability in (0, 1).
#' @param seed Base RNG seed (cut-mode stage seeds derive from it).
#' @param progress Show sampler progress.
#' @param gate_params Parameter-name prefixes the convergence gate is
#'   restricted to (divergences always count run-wide). Defaults to
#'   `c("mu_ctrl", "te")` -- the counterfactual surface and the treatment
#'   effect, which is what the reported estimands are built from. Pass a wider
#'   vector to gate on more, or `"all"` to gate on every diagnosable variable
#'   (the behaviour before this became a default). [parameter_diagnostics()]
#'   always reports every variable regardless.
#' @param convergence A [bpnmf_convergence()] object.
#' @export
bpnmf_mcmc_opts <- function(auto_parallelism = TRUE, max_chains = 4,
                            chains = NULL, parallel_chains = NULL,
                            iter_warmup = 1000, iter_sampling = 2500,
                            thin = 10, adapt_delta = 0.8, seed = 8675309,
                            progress = TRUE,
                            gate_params = DEFAULT_GATE_PARAMS,
                            convergence = bpnmf_convergence()) {
  checkmate::assert_flag(auto_parallelism)
  checkmate::assert_int(max_chains, lower = 1)
  checkmate::assert_int(chains, lower = 1, null.ok = TRUE)
  checkmate::assert_int(parallel_chains, lower = 1, null.ok = TRUE)
  checkmate::assert_int(iter_warmup, lower = 1)
  checkmate::assert_int(iter_sampling, lower = 1)
  checkmate::assert_int(thin, lower = 1)
  checkmate::assert_number(adapt_delta)
  if (adapt_delta <= 0 || adapt_delta >= 1) {
    cli::cli_abort("{.field adapt_delta} must be strictly between 0 and 1.")
  }
  checkmate::assert_int(seed)
  checkmate::assert_flag(progress)
  # Resolve here rather than relying on the signature default: the YAML
  # loader passes gate_params = NULL explicitly when the key is absent, and an
  # explicit NULL skips an R default.
  gate_params <- gate_params %||% DEFAULT_GATE_PARAMS
  checkmate::assert_character(gate_params, min.len = 1, any.missing = FALSE)
  convergence <- coerce_bpnmf(
    convergence, "bpnmf_convergence", "bpnmf_convergence", "convergence"
  )
  checkmate::assert_class(convergence, "bpnmf_convergence")
  new_bpnmf_class(
    list(
      auto_parallelism = auto_parallelism, max_chains = as.integer(max_chains),
      chains = if (is.null(chains)) NULL else as.integer(chains),
      parallel_chains = if (is.null(parallel_chains)) NULL else as.integer(parallel_chains),
      iter_warmup = as.integer(iter_warmup),
      iter_sampling = as.integer(iter_sampling), thin = as.integer(thin),
      adapt_delta = adapt_delta, seed = as.integer(seed), progress = progress,
      gate_params = gate_params, convergence = convergence
    ),
    "bpnmf_mcmc_opts"
  )
}

#' One synthetic reporting-only aggregate unit
#'
#' Exactly one of `include_treated_units`, `include_all_units`, or
#' `include_units` must select the member units.
#'
#' @param unit Name of the synthetic unit.
#' @param include_treated_units,include_all_units Selector flags.
#' @param include_units Explicit member unit names.
#' @param exclude_units Units removed after selection.
#' @param strict Error (instead of warn) when a listed unit is absent.
#' @param overwrite Replace an existing unit of the same name.
#' @export
bpnmf_aggregate_unit <- function(unit, include_treated_units = FALSE,
                                 include_all_units = FALSE,
                                 include_units = NULL, exclude_units = NULL,
                                 strict = FALSE, overwrite = FALSE) {
  checkmate::assert_string(unit, min.chars = 1)
  checkmate::assert_flag(include_treated_units)
  checkmate::assert_flag(include_all_units)
  checkmate::assert_character(
    include_units,
    min.len = 1, any.missing = FALSE, null.ok = TRUE
  )
  checkmate::assert_character(exclude_units, any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_flag(strict)
  checkmate::assert_flag(overwrite)
  active <- c(
    include_treated_units = include_treated_units,
    include_all_units = include_all_units,
    include_units = !is.null(include_units)
  )
  if (sum(active) != 1) {
    cli::cli_abort(
      "aggregate unit {.val {unit}} must have exactly one include selector;
       got {.val {names(active)[active]}}."
    )
  }
  new_bpnmf_class(
    list(
      unit = unit, include_treated_units = include_treated_units,
      include_all_units = include_all_units, include_units = include_units,
      exclude_units = exclude_units, strict = strict, overwrite = overwrite
    ),
    "bpnmf_aggregate_unit"
  )
}

normalize_figures <- function(v) {
  if (is.logical(v) && length(v) == 1 && !is.na(v)) {
    return(if (v) sort(FIGURE_NAMES) else character())
  }
  if (is.character(v) && length(v) == 1 && v %in% c("all", "none")) {
    return(if (v == "all") sort(FIGURE_NAMES) else character())
  }
  if (is.character(v) || (is.list(v) && all(vapply(v, is.character, logical(1))))) {
    v <- as.character(v)
    unknown <- setdiff(v, FIGURE_NAMES)
    if (length(unknown) > 0) {
      cli::cli_abort(c(
        "output.figures contains unknown figure name{?s} {.val {unknown}}.",
        i = "Valid names: {.val {sort(FIGURE_NAMES)}}."
      ))
    }
    return(v)
  }
  cli::cli_abort(
    "output.figures must be TRUE/FALSE, 'all', 'none', or a character vector
     of figure names."
  )
}

#' Output options: figures, tables, reporting/PPC filters
#'
#' @param figures `TRUE`/`FALSE`, `"all"`/`"none"`, or a character vector of
#'   figure names from `r toString(FIGURE_NAMES)`. Tables always render.
#' @param clean Remove the per-type output dir before writing.
#' @param save_traces Save full posterior draws (RDS via \pkg{posterior}).
#' @param target_unit Unit highlighted in tables and per-unit figures.
#' @param report_groups Restrict per-unit figures to these groups.
#' @param fit_gap_per_unit Also render fit/gap for every treated unit
#'   (restricted to the `"total"` group).
#' @param print_tables,print_target_table Terminal-table switches.
#'   `print_target_table` adds the target unit's own table above the by-unit
#'   one; its rows are already in the by-unit table, so it defaults to off.
#' @param html_tables Also write the `gt` HTML summary tables.
#' @param aggregate_units List of [bpnmf_aggregate_unit()] specs.
#' @param ppc_units,ppc_exclude_units Unit filters for the PPC suite.
#' @param ppc_acf_lags Integer lags for the ACF check (default 1).
#' @param ppc_unit_corr_max_time Cutoff date for the unit-correlation check.
#' @param draws_format `"csv"` or `"parquet"` for the draws artifact.
#' @param rate_normalizer Reported rates are per this many units of exposure
#'   (default 1000). Purely a display scale: it cancels out of every
#'   percent-change figure and only sets the scale of the rate columns and
#'   the rate-difference axis.
#' @param denominator_label Singular noun for one unit of the denominator --
#'   `"person"`, `"birth"`, `"vehicle-mile"`, whatever your denominator
#'   column counts. Combined with `denominator_time_unit` to name the
#'   exposure ("person-years"), which labels the rate columns and axes.
#'   Defaults to the neutral `"denominator"`, since the package has no way to
#'   know what yours measures.
#' @param denominator_time_unit Time unit the denominator is weighted by when
#'   forming exposure: one of `"year"` (default), `"month"`, `"week"`,
#'   `"day"`, or `"none"`. This is part of the quantity, not decoration --
#'   rates are counts per `denominator * time`, so changing it changes the
#'   numbers. Use `"none"` when the denominator is itself a per-period flow
#'   (births *during* each quarter) rather than a level standing through the
#'   period (population), since time-weighting a flow double-counts the
#'   period length.
#' @param denominator_may_be_affected The model treats the denominator as a
#'   fixed, treatment-unaffected exposure: `mu_ctrl` bakes in
#'   `log(denominator)` as a constant offset, so "Expected"/"Pct Change" are
#'   the counterfactual count conditional on the *observed* denominator, not
#'   a counterfactual where the denominator itself is also allowed to vary.
#'   That's usually fine (e.g. population as the denominator for a mortality
#'   rate), but not when treatment could plausibly change the denominator
#'   too (e.g. births as the denominator for an infant mortality rate, when
#'   the exposure could also change the number of births). `TRUE` (the
#'   default) surfaces that caveat in the printed report and the `gt`
#'   tables; it changes no computation. Set `FALSE` only once you've
#'   confirmed the denominator is exogenous to treatment, to quiet the
#'   caveat -- the default errs toward showing it, since a user who doesn't
#'   already know about the issue has no reason to go looking for a flag to
#'   turn it on.
#' @export
bpnmf_output_opts <- function(figures = FALSE, clean = FALSE,
                              save_traces = FALSE, target_unit = NULL,
                              report_groups = NULL, fit_gap_per_unit = FALSE,
                              interval_aggregates = TRUE,
                              print_tables = TRUE, print_target_table = FALSE,
                              html_tables = TRUE,
                              aggregate_units = NULL, ppc_units = NULL,
                              ppc_exclude_units = NULL, ppc_acf_lags = NULL,
                              ppc_unit_corr_max_time = NULL,
                              draws_format = "csv",
                              rate_normalizer = 1000,
                              denominator_label = "denominator",
                              denominator_time_unit = "year",
                              denominator_may_be_affected = TRUE) {
  figures <- normalize_figures(figures)
  checkmate::assert_flag(clean)
  checkmate::assert_flag(save_traces)
  checkmate::assert_string(target_unit, min.chars = 1, null.ok = TRUE)
  checkmate::assert_character(
    report_groups,
    min.len = 1, any.missing = FALSE, null.ok = TRUE
  )
  checkmate::assert_flag(fit_gap_per_unit)
  checkmate::assert_flag(interval_aggregates)
  checkmate::assert_flag(print_tables)
  checkmate::assert_flag(print_target_table)
  checkmate::assert_flag(html_tables)
  checkmate::assert_number(rate_normalizer, lower = .Machine$double.eps)
  checkmate::assert_string(denominator_label, min.chars = 1)
  checkmate::assert_choice(denominator_time_unit, DENOMINATOR_TIME_UNITS)
  checkmate::assert_flag(denominator_may_be_affected)
  aggregate_units <- coerce_bpnmf_list(
    aggregate_units, "bpnmf_aggregate_unit", "bpnmf_aggregate_unit",
    "aggregate_units"
  )
  checkmate::assert_list(
    aggregate_units,
    types = "bpnmf_aggregate_unit", null.ok = TRUE
  )
  checkmate::assert_character(ppc_units, any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_character(
    ppc_exclude_units,
    any.missing = FALSE, null.ok = TRUE
  )
  checkmate::assert_integerish(
    ppc_acf_lags,
    lower = 1, any.missing = FALSE, null.ok = TRUE
  )
  checkmate::assert_string(ppc_unit_corr_max_time, null.ok = TRUE)
  checkmate::assert_choice(draws_format, c("csv", "parquet"))
  new_bpnmf_class(
    list(
      figures = figures, clean = clean, save_traces = save_traces,
      target_unit = target_unit, report_groups = report_groups,
      fit_gap_per_unit = fit_gap_per_unit,
      interval_aggregates = interval_aggregates, print_tables = print_tables,
      print_target_table = print_target_table, html_tables = html_tables,
      aggregate_units = aggregate_units, ppc_units = ppc_units,
      ppc_exclude_units = ppc_exclude_units,
      ppc_acf_lags = if (is.null(ppc_acf_lags)) NULL else as.integer(ppc_acf_lags),
      ppc_unit_corr_max_time = ppc_unit_corr_max_time,
      draws_format = draws_format,
      rate_normalizer = rate_normalizer,
      denominator_label = denominator_label,
      denominator_time_unit = denominator_time_unit,
      denominator_may_be_affected = denominator_may_be_affected
    ),
    "bpnmf_output_opts"
  )
}

#' Two-stage cut-posterior settings
#'
#' @param num_stage1_draws Stage-1 posterior draws promoted to cut components.
#' @param stage2_draws_per_component Retained draws per component after evenly
#'   strided output thinning (`NULL` keeps all; counts must be equal across
#'   components so pooling weights them equally).
#' @param selection_seed Seed for the chain-stratified stage-1 draw selection
#'   (default `seed + 2`).
#' @param stage2_seed Base seed for stage-2 sampling (default `seed + 3`;
#'   component `i` runs at `stage2_seed + i`).
#' @param stage2_mcmc Named list shallow-merged over the top-level MCMC options
#'   for stage 2 (accepts both R and Python key names; must not set the seed --
#'   `stage2_seed` is the authority).
#' @param stage1_method How stage 1 is fit: `"sample"` (NUTS, the default) or
#'   `"variational"` (Stan's ADVI). ADVI turns a stage-1 fit that takes hours
#'   into one that takes minutes, at the cost of an approximation: mean-field
#'   ADVI understates posterior variance and ignores posterior correlations,
#'   so the stage-1 components it hands to stage 2 are drawn from too narrow a
#'   spread and the pooled cut posterior for `te` comes out over-confident.
#'   There is also no R-hat / ESS / divergence gate for a variational fit, so
#'   `manifest$stage1$converged` is `NA` rather than `TRUE`/`FALSE`. Use it to
#'   iterate; re-run with `"sample"` for results you intend to report.
#' @param stage1_variational Named list of ADVI tuning arguments passed to
#'   cmdstanr's `variational()` -- any of
#'   `r toString(VARIATIONAL_KEYS)`. The seed is not settable here (stage 1
#'   runs at `mcmc$seed`, as it does under `"sample"`). Ignored when
#'   `stage1_method = "sample"`.
#' @export
bpnmf_cut_opts <- function(num_stage1_draws = 25,
                           stage2_draws_per_component = 100,
                           selection_seed = NULL, stage2_seed = NULL,
                           stage2_mcmc = NULL, stage1_method = "sample",
                           stage1_variational = NULL) {
  checkmate::assert_int(num_stage1_draws, lower = 1)
  checkmate::assert_int(stage2_draws_per_component, lower = 1, null.ok = TRUE)
  checkmate::assert_int(selection_seed, null.ok = TRUE)
  checkmate::assert_int(stage2_seed, null.ok = TRUE)
  checkmate::assert_list(stage2_mcmc, names = "unique", null.ok = TRUE)
  checkmate::assert_choice(stage1_method, c("sample", "variational"))
  checkmate::assert_list(stage1_variational, names = "unique", null.ok = TRUE)
  if (!is.null(stage2_mcmc) &&
    any(c("random_seed", "seed") %in% names(stage2_mcmc))) {
    cli::cli_abort(
      "cut.stage2_mcmc must not set the seed; {.field cut.stage2_seed} is the
       Stage-2 seed authority."
    )
  }
  if (!is.null(stage1_variational)) {
    unknown <- setdiff(names(stage1_variational), VARIATIONAL_KEYS)
    if (length(unknown) > 0) {
      cli::cli_abort(c(
        "Unknown {cli::qty(length(unknown))}key{?s} {.val {unknown}} in
         {.field cut.stage1_variational}.",
        i = "Valid keys: {.val {VARIATIONAL_KEYS}}."
      ))
    }
    if (identical(stage1_method, "sample")) {
      cli::cli_warn(
        "{.field cut.stage1_variational} is ignored because
         {.field cut.stage1_method} is {.val sample}."
      )
    }
  }
  new_bpnmf_class(
    list(
      num_stage1_draws = as.integer(num_stage1_draws),
      stage2_draws_per_component =
        if (is.null(stage2_draws_per_component)) {
          NULL
        } else {
          as.integer(stage2_draws_per_component)
        },
      selection_seed = if (is.null(selection_seed)) NULL else as.integer(selection_seed),
      stage2_seed = if (is.null(stage2_seed)) NULL else as.integer(stage2_seed),
      stage2_mcmc = stage2_mcmc,
      stage1_method = stage1_method,
      stage1_variational = stage1_variational
    ),
    "bpnmf_cut_opts"
  )
}

#' Top-level bpnmf configuration
#'
#' Anywhere a `bpnmf_*` options object is expected -- here and in
#' [bpnmf_schema()], [bpnmf_model_opts()], [bpnmf_mcmc_opts()] and
#' [bpnmf_output_opts()] -- a plain **named list** of that constructor's
#' arguments is accepted and passed through the constructor itself. Validation
#' and defaults are identical either way, and an unrecognized name is an error
#' rather than a silently ignored option, so the list form is a shorthand
#' rather than a way around the checks. Use the constructors for argument
#' completion, or lists for a single call shaped like the YAML config.
#'
#' @param input_file Path to the input CSV.
#' @param output_dir Directory run artifacts are written to.
#' @param schema A [bpnmf_schema()] object, or a named list of its arguments.
#' @param model A [bpnmf_model_opts()] object, or a named list of its arguments.
#' @param mcmc A [bpnmf_mcmc_opts()] object, or a named list of its arguments.
#' @param output A [bpnmf_output_opts()] object, or a named list of its arguments.
#' @param cut A [bpnmf_cut_opts()] object, a named list, or `NULL`.
#' @param date_format `"auto"` or a `strptime` format for the time column.
#' @param start_date,end_date Optional date filter; `start_date` inclusive,
#'   `end_date` **exclusive**.
#' @param time_aggregation A [bpnmf_time_aggregation()] object, or a named
#'   list of its arguments.
#' @param allow_unbalanced_panel Treat structurally absent (unit, time) cells
#'   as missing instead of erroring.
#' @param outcome Optional label used in draws filenames (falls back to the
#'   outcome prefix, then `"outcome"`).
#' @return A validated `bpnmf_config` object.
#' @export
bpnmf_config <- function(input_file, output_dir, schema,
                         model = bpnmf_model_opts(),
                         mcmc = bpnmf_mcmc_opts(),
                         output = bpnmf_output_opts(),
                         cut = NULL,
                         date_format = "auto", start_date = NULL,
                         end_date = NULL,
                         time_aggregation = bpnmf_time_aggregation(),
                         allow_unbalanced_panel = FALSE, outcome = NULL) {
  checkmate::assert_string(input_file, min.chars = 1)
  checkmate::assert_string(output_dir, min.chars = 1)
  schema <- coerce_bpnmf(schema, "bpnmf_schema", "bpnmf_schema", "schema")
  model <- coerce_bpnmf(model, "bpnmf_model_opts", "bpnmf_model_opts", "model")
  mcmc <- coerce_bpnmf(mcmc, "bpnmf_mcmc_opts", "bpnmf_mcmc_opts", "mcmc")
  output <- coerce_bpnmf(output, "bpnmf_output_opts", "bpnmf_output_opts", "output")
  cut <- coerce_bpnmf(cut, "bpnmf_cut_opts", "bpnmf_cut_opts", "cut")
  checkmate::assert_class(schema, "bpnmf_schema")
  checkmate::assert_class(model, "bpnmf_model_opts")
  checkmate::assert_class(mcmc, "bpnmf_mcmc_opts")
  checkmate::assert_class(output, "bpnmf_output_opts")
  checkmate::assert_class(cut, "bpnmf_cut_opts", null.ok = TRUE)
  checkmate::assert_string(date_format, min.chars = 1)
  checkmate::assert_string(start_date, null.ok = TRUE)
  checkmate::assert_string(end_date, null.ok = TRUE)
  time_aggregation <- coerce_bpnmf(
    time_aggregation, "bpnmf_time_aggregation", "bpnmf_time_aggregation",
    "time_aggregation"
  )
  checkmate::assert_class(time_aggregation, "bpnmf_time_aggregation")
  checkmate::assert_flag(allow_unbalanced_panel)
  checkmate::assert_string(outcome, min.chars = 1, null.ok = TRUE)
  new_bpnmf_class(
    list(
      input_file = input_file, output_dir = output_dir, schema = schema,
      model = model, mcmc = mcmc, output = output, cut = cut,
      date_format = date_format, start_date = start_date, end_date = end_date,
      time_aggregation = time_aggregation,
      allow_unbalanced_panel = allow_unbalanced_panel, outcome = outcome
    ),
    "bpnmf_config"
  )
}

#' @export
print.bpnmf_config <- function(x, ...) {
  cli::cli_h1("bpnmf config")
  cli::cli_li("input: {.file {x$input_file}}")
  cli::cli_li("output dir: {.file {x$output_dir}}")
  cli::cli_li("distribution: {x$model$outcome_distribution}")
  mode <- x$model$inference_mode %||% "joint"
  cli::cli_li("inference mode: {mode}")
  if (identical(mode, "cut")) {
    cli::cli_li("cut stage 1: {fit_method_label(x$cut$stage1_method %||% 'sample')}")
  }
  if (length(x$model$types) > 0) {
    for (nm in names(x$model$types)) {
      tp <- x$model$types[[nm]]
      cli::cli_li(
        "type {.strong {nm}}: {length(tp$groups)} group{?s}, rank{?s} {.val {tp$ranks_to_test}}"
      )
    }
  }
  cli::cli_li(
    "mcmc: warmup {x$mcmc$iter_warmup}, sampling {x$mcmc$iter_sampling}, thin {x$mcmc$thin}, adapt_delta {x$mcmc$adapt_delta}, seed {x$mcmc$seed}"
  )
  invisible(x)
}
