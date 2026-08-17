# Config constructors + validation. Mirrors the pydantic v2 schema in the
# Python package (config.py): every constructor validates immediately, unknown
# fields are impossible by construction here and rejected by the YAML loader
# (config-yaml.R), and cross-field rules (outcomes XOR prefixes, cut requires
# model_treated, stage2 seed authority, aggregate-unit single selector) are
# enforced at build time so a bad config never reaches the sampler.

# Canonical figure names output$figures may select. Must stay in sync with the
# figure functions wired up in reports.R (the Python side derives this from
# PLOT_REGISTRY for the same no-drift reason).
FIGURE_NAMES <- c(
  "unit_fit", "unit_gap", "raw_rate", "interval", "group_comparison", "ppc"
)

AGGREGATION_PERIODS <- c("monthly", "bimonthly", "quarterly", "yearly")

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
    checkmate::assert_list(outcomes, min.len = 1, types = "bpnmf_outcome")
  }
  if (!is.null(outcomes_from_prefixes)) {
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
#' @param enabled Aggregate input rows into coarser periods?
#' @param period One of `"monthly"`, `"bimonthly"`, `"quarterly"`, `"yearly"`.
#' @export
bpnmf_aggregation <- function(enabled = FALSE, period = "bimonthly") {
  checkmate::assert_flag(enabled)
  checkmate::assert_choice(period, AGGREGATION_PERIODS)
  new_bpnmf_class(
    list(enabled = enabled, period = period),
    "bpnmf_aggregation"
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
#' @export
bpnmf_model_opts <- function(outcome_distribution = "NB", types = list(),
                             nb_disp = 1e-4, sample_disp = FALSE,
                             adjust_for_missingness = TRUE,
                             model_treated = TRUE, inference_mode = NULL) {
  checkmate::assert_choice(outcome_distribution, c("NB", "Poisson"))
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
  new_bpnmf_class(
    list(
      outcome_distribution = outcome_distribution, types = types,
      nb_disp = nb_disp, sample_disp = sample_disp,
      adjust_for_missingness = adjust_for_missingness,
      model_treated = model_treated, inference_mode = inference_mode
    ),
    "bpnmf_model_opts"
  )
}

#' Convergence-gate thresholds
#'
#' A parameter PASSes below `rhat_warn` with ESS at/above `ess_min`; FAILs at/
#' above `rhat_fail` or below `ess_min * ess_fail_fraction`; WARNs in between.
#' ESS is `min(bulk, tail)` per parameter.
#' @param rhat_warn,rhat_fail,ess_min,ess_fail_fraction Gate thresholds.
#' @export
bpnmf_convergence <- function(rhat_warn = 1.01, rhat_fail = 1.05,
                              ess_min = 400, ess_fail_fraction = 0.25) {
  checkmate::assert_number(rhat_warn, lower = 1)
  checkmate::assert_number(rhat_fail, lower = 1)
  checkmate::assert_number(ess_min, lower = 0)
  checkmate::assert_number(ess_fail_fraction, lower = 0, upper = 1)
  new_bpnmf_class(
    list(
      rhat_warn = rhat_warn, rhat_fail = rhat_fail,
      ess_min = ess_min, ess_fail_fraction = ess_fail_fraction
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
#' @param gate_params Optional character vector of parameter-name prefixes the
#'   convergence gate is restricted to (divergences always count run-wide).
#' @param convergence A [bpnmf_convergence()] object.
#' @export
bpnmf_mcmc_opts <- function(auto_parallelism = TRUE, max_chains = 4,
                            chains = NULL, parallel_chains = NULL,
                            iter_warmup = 1000, iter_sampling = 2500,
                            thin = 10, adapt_delta = 0.8, seed = 8675309,
                            progress = TRUE, gate_params = NULL,
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
  checkmate::assert_character(
    gate_params,
    min.len = 1, any.missing = FALSE, null.ok = TRUE
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
#' @param aggregate_units List of [bpnmf_aggregate_unit()] specs.
#' @param ppc_units,ppc_exclude_units Unit filters for the PPC suite.
#' @param ppc_acf_lags Integer lags for the ACF check (default 6).
#' @param ppc_unit_corr_max_time Cutoff date for the unit-correlation check.
#' @param draws_format `"csv"` or `"parquet"` for the draws artifact.
#' @export
bpnmf_output_opts <- function(figures = FALSE, clean = FALSE,
                              save_traces = FALSE, target_unit = NULL,
                              report_groups = NULL, fit_gap_per_unit = FALSE,
                              print_tables = TRUE, print_target_table = TRUE,
                              aggregate_units = NULL, ppc_units = NULL,
                              ppc_exclude_units = NULL, ppc_acf_lags = NULL,
                              ppc_unit_corr_max_time = NULL,
                              draws_format = "csv") {
  figures <- normalize_figures(figures)
  checkmate::assert_flag(clean)
  checkmate::assert_flag(save_traces)
  checkmate::assert_string(target_unit, min.chars = 1, null.ok = TRUE)
  checkmate::assert_character(
    report_groups,
    min.len = 1, any.missing = FALSE, null.ok = TRUE
  )
  checkmate::assert_flag(fit_gap_per_unit)
  checkmate::assert_flag(print_tables)
  checkmate::assert_flag(print_target_table)
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
      fit_gap_per_unit = fit_gap_per_unit, print_tables = print_tables,
      print_target_table = print_target_table,
      aggregate_units = aggregate_units, ppc_units = ppc_units,
      ppc_exclude_units = ppc_exclude_units,
      ppc_acf_lags = if (is.null(ppc_acf_lags)) NULL else as.integer(ppc_acf_lags),
      ppc_unit_corr_max_time = ppc_unit_corr_max_time,
      draws_format = draws_format
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
#' @export
bpnmf_cut_opts <- function(num_stage1_draws = 25,
                           stage2_draws_per_component = 100,
                           selection_seed = NULL, stage2_seed = NULL,
                           stage2_mcmc = NULL) {
  checkmate::assert_int(num_stage1_draws, lower = 1)
  checkmate::assert_int(stage2_draws_per_component, lower = 1, null.ok = TRUE)
  checkmate::assert_int(selection_seed, null.ok = TRUE)
  checkmate::assert_int(stage2_seed, null.ok = TRUE)
  checkmate::assert_list(stage2_mcmc, names = "unique", null.ok = TRUE)
  if (!is.null(stage2_mcmc) &&
    any(c("random_seed", "seed") %in% names(stage2_mcmc))) {
    cli::cli_abort(
      "cut.stage2_mcmc must not set the seed; {.field cut.stage2_seed} is the
       Stage-2 seed authority."
    )
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
      stage2_mcmc = stage2_mcmc
    ),
    "bpnmf_cut_opts"
  )
}

#' Top-level bpnmf configuration
#'
#' @param input_file Path to the input CSV.
#' @param output_dir Directory run artifacts are written to.
#' @param schema A [bpnmf_schema()] object.
#' @param model A [bpnmf_model_opts()] object.
#' @param mcmc A [bpnmf_mcmc_opts()] object.
#' @param output A [bpnmf_output_opts()] object.
#' @param cut A [bpnmf_cut_opts()] object, or `NULL`.
#' @param date_format `"auto"` or a `strptime` format for the time column.
#' @param start_date,end_date Optional date filter; `start_date` inclusive,
#'   `end_date` **exclusive**.
#' @param aggregation A [bpnmf_aggregation()] object.
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
                         aggregation = bpnmf_aggregation(),
                         allow_unbalanced_panel = FALSE, outcome = NULL) {
  checkmate::assert_string(input_file, min.chars = 1)
  checkmate::assert_string(output_dir, min.chars = 1)
  checkmate::assert_class(schema, "bpnmf_schema")
  checkmate::assert_class(model, "bpnmf_model_opts")
  checkmate::assert_class(mcmc, "bpnmf_mcmc_opts")
  checkmate::assert_class(output, "bpnmf_output_opts")
  checkmate::assert_class(cut, "bpnmf_cut_opts", null.ok = TRUE)
  checkmate::assert_string(date_format, min.chars = 1)
  checkmate::assert_string(start_date, null.ok = TRUE)
  checkmate::assert_string(end_date, null.ok = TRUE)
  checkmate::assert_class(aggregation, "bpnmf_aggregation")
  checkmate::assert_flag(allow_unbalanced_panel)
  checkmate::assert_string(outcome, min.chars = 1, null.ok = TRUE)
  new_bpnmf_class(
    list(
      input_file = input_file, output_dir = output_dir, schema = schema,
      model = model, mcmc = mcmc, output = output, cut = cut,
      date_format = date_format, start_date = start_date, end_date = end_date,
      aggregation = aggregation,
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
