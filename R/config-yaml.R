# YAML config loader. Accepts the Python package's config schema verbatim
# (same section names, same key names) and translates to the R-native
# countsynth_config object. Replicates pydantic's extra="forbid" (unknown keys at
# any nesting level are an error, with the YAML path in the message) and
# StrictBool (a quoted "false" is rejected rather than treated as truthy).

check_known_keys <- function(x, known, path) {
  unknown <- setdiff(names(x), known)
  if (length(unknown) > 0) {
    cli::cli_abort(c(
      "Unknown config {cli::qty(length(unknown))}key{?s} {.val {unknown}}
       under {.field {path}}.",
      i = "Valid keys: {.val {sort(known)}}."
    ))
  }
}

yaml_flag <- function(x, path, default = NULL) {
  if (is.null(x)) {
    return(default)
  }
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    cli::cli_abort(
      "{.field {path}} must be a YAML boolean (true/false), got {.val {x}}."
    )
  }
  x
}

yaml_chr <- function(x) {
  if (is.null(x)) NULL else as.character(unlist(x, use.names = FALSE))
}

parse_yaml_schema <- function(x, path = "data.schema") {
  check_known_keys(
    x,
    c("unit_col", "time_col", "treatment_col", "outcomes",
      "outcomes_from_prefixes"),
    path
  )
  # NB: [[ ]] everywhere below -- `$` partial matching would let `x$outcomes`
  # silently resolve to `outcomes_from_prefixes`.
  outcomes <- NULL
  if (!is.null(x[["outcomes"]])) {
    outcomes <- lapply(seq_along(x[["outcomes"]]), function(i) {
      o <- x[["outcomes"]][[i]]
      check_known_keys(
        o, c("outcome_col", "label", "denominator_col"),
        glue::glue("{path}.outcomes[{i}]")
      )
      countsynth_outcome(o$outcome_col, o$label, o$denominator_col)
    })
  }
  prefixes <- NULL
  if (!is.null(x[["outcomes_from_prefixes"]])) {
    p <- x[["outcomes_from_prefixes"]]
    check_known_keys(
      p, c("outcome_prefix", "denominator_prefix", "include"),
      glue::glue("{path}.outcomes_from_prefixes")
    )
    prefixes <- countsynth_prefixes(
      p$outcome_prefix, p$denominator_prefix, yaml_chr(p$include)
    )
  }
  countsynth_schema(
    unit_col = x$unit_col, time_col = x$time_col,
    treatment_col = x$treatment_col,
    outcomes = outcomes, outcomes_from_prefixes = prefixes
  )
}

parse_yaml_model <- function(x, path = "model") {
  if (is.null(x)) {
    return(countsynth_model_opts())
  }
  check_known_keys(
    x,
    c("outcome_distribution", "types", "nb_disp", "sample_disp",
      "adjust_for_missingness", "model_treated", "inference_mode",
      "treatment_effects", "factor_variation_pct",
      "time_level_variation_pct", "rank_shrinkage", "shared_curves",
      "time_level"),
    path
  )
  rank_shrinkage <- parse_yaml_rank_shrinkage(
    x$rank_shrinkage, glue::glue("{path}.rank_shrinkage")
  )
  treatment_effects <- NULL
  if (!is.null(x$treatment_effects)) {
    te <- x$treatment_effects
    tpath <- glue::glue("{path}.treatment_effects")
    check_known_keys(
      te,
      c("formula", "covariates_file", "standardize", "coef_prior_scale",
        "re_prior_scale"),
      tpath
    )
    if (is.null(te$formula)) {
      cli::cli_abort("{.field {tpath}} must define a {.field formula}.")
    }
    treatment_effects <- countsynth_te_opts(
      formula = te$formula,
      covariates = te$covariates_file,
      standardize = yaml_flag(te$standardize, glue::glue("{tpath}.standardize"), TRUE),
      coef_prior_scale = unlist(te$coef_prior_scale %||% 1, use.names = FALSE),
      re_prior_scale = unlist(te$re_prior_scale %||% 1, use.names = FALSE)
    )
  }
  types <- list()
  if (!is.null(x$types)) {
    types <- lapply(names(x$types), function(nm) {
      tp <- x$types[[nm]]
      check_known_keys(
        tp,
        c("groups", "ranks_to_test", "total_from", "total_all",
          "exclude_units"),
        glue::glue("{path}.types.{nm}")
      )
      countsynth_type(
        groups = yaml_chr(tp$groups),
        ranks_to_test = unlist(tp$ranks_to_test, use.names = FALSE),
        total_from = yaml_chr(tp$total_from),
        total_all = yaml_flag(tp$total_all, glue::glue("{path}.types.{nm}.total_all"), FALSE),
        exclude_units = yaml_chr(tp$exclude_units)
      )
    })
    names(types) <- names(x$types)
  }
  countsynth_model_opts(
    outcome_distribution = x$outcome_distribution %||% "NB",
    types = types,
    nb_disp = x$nb_disp %||% 1e-4,
    sample_disp = yaml_flag(x$sample_disp, glue::glue("{path}.sample_disp"), FALSE),
    adjust_for_missingness =
      yaml_flag(x$adjust_for_missingness, glue::glue("{path}.adjust_for_missingness"), TRUE),
    model_treated = yaml_flag(x$model_treated, glue::glue("{path}.model_treated"), TRUE),
    inference_mode = x$inference_mode,
    treatment_effects = treatment_effects,
    factor_variation_pct = x$factor_variation_pct,
    time_level_variation_pct = x$time_level_variation_pct,
    rank_shrinkage = rank_shrinkage,
    shared_curves = yaml_flag(
      x$shared_curves, glue::glue("{path}.shared_curves"), FALSE
    ),
    # [[ ]], not $: `$` partial-matches, and with no `time_level` key it
    # would return time_level_variation_pct's value instead.
    time_level = x[["time_level"]] %||% "centered"
  )
}

# `rank_shrinkage: true` is the shorthand for the defaults; a mapping
# overrides individual hyperparameters. The two Gamma priors are written as
# two-element `[shape, rate]` sequences, which YAML gives us as a list.
parse_yaml_rank_shrinkage <- function(x, path) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.logical(x) && length(x) == 1 && !is.na(x)) {
    return(if (x) countsynth_rank_shrinkage_opts() else NULL)
  }
  check_known_keys(
    x, c("group_mass_prior", "unit_sd_prior", "group_sd_prior",
         "shared_unit_sd_prior"), path
  )
  num <- function(key, default) {
    if (is.null(x[[key]])) {
      return(default)
    }
    v <- unlist(x[[key]], use.names = FALSE)
    if (!is.numeric(v)) {
      cli::cli_abort("{.field {path}.{key}} must be numeric.")
    }
    v
  }
  countsynth_rank_shrinkage_opts(
    group_mass_prior = num("group_mass_prior", c(2, 1)),
    unit_sd_prior = num("unit_sd_prior", 1),
    group_sd_prior = num("group_sd_prior", 0.5),
    shared_unit_sd_prior = num("shared_unit_sd_prior", 1)
  )
}

parse_yaml_mcmc <- function(x, path = "mcmc") {
  if (is.null(x)) {
    return(countsynth_mcmc_opts())
  }
  check_known_keys(
    x,
    c("auto_parallelism", "max_chains", "num_chains", "chain_method",
      "num_warmup", "num_samples", "thinning", "target_accept", "random_seed",
      "progress_bar", "gate_params", "convergence", "max_treedepth"),
    path
  )
  convergence <- countsynth_convergence()
  if (!is.null(x$convergence)) {
    cv <- x$convergence
    check_known_keys(
      cv, c("rhat_warn", "rhat_fail", "ess_min", "ess_fail_fraction",
            "divergence_fail_fraction"),
      glue::glue("{path}.convergence")
    )
    convergence <- countsynth_convergence(
      rhat_warn = cv$rhat_warn %||% 1.01,
      rhat_fail = cv$rhat_fail %||% 1.05,
      ess_min = cv$ess_min %||% 400,
      ess_fail_fraction = cv$ess_fail_fraction %||% 0.25,
      divergence_fail_fraction = cv$divergence_fail_fraction %||% 0.01
    )
  }
  chains <- x$num_chains
  parallel_chains <- NULL
  if (!is.null(x$chain_method)) {
    checkmate::assert_choice(
      x$chain_method, c("sequential", "parallel", "vectorized")
    )
    if (x$chain_method == "vectorized") {
      cli::cli_warn(
        "{.field {path}.chain_method}='vectorized' has no Stan analogue;
         running chains in parallel instead."
      )
    }
    parallel_chains <- switch(x$chain_method,
      sequential = 1L,
      parallel = chains,
      vectorized = chains
    )
  }
  countsynth_mcmc_opts(
    auto_parallelism = yaml_flag(x$auto_parallelism, glue::glue("{path}.auto_parallelism"), TRUE),
    max_chains = x$max_chains %||% 4L,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = x$num_warmup %||% 1000L,
    iter_sampling = x$num_samples %||% 2500L,
    thin = x$thinning %||% 10L,
    adapt_delta = x$target_accept %||% 0.8,
    seed = x$random_seed %||% 8675309L,
    progress = yaml_flag(x$progress_bar, glue::glue("{path}.progress_bar"), TRUE),
    # cmdstanr's name: the Python package has no tree-depth setting to mirror.
    max_treedepth = x$max_treedepth %||% 10L,
    gate_params = yaml_chr(x$gate_params),
    convergence = convergence
  )
}

parse_yaml_output <- function(x, path = "output") {
  if (is.null(x)) {
    return(countsynth_output_opts())
  }
  check_known_keys(
    x,
    c("figures", "clean", "save_traces", "target_unit", "report_groups",
      "fit_gap_per_unit", "interval_aggregates",
      "print_tables", "print_target_table", "html_tables",
      "aggregate_units", "ppc_units", "ppc_exclude_units", "ppc_acf_lags",
      "ppc_unit_corr_max_time", "draws_format", "rate_normalizer",
      "denominator_label", "denominator_time_unit",
      "denominator_may_be_affected"),
    path
  )
  figures <- x$figures %||% FALSE
  if (is.character(figures) && length(figures) == 1 &&
    !figures %in% c("all", "none", FIGURE_NAMES)) {
    cli::cli_abort(
      "{.field {path}.figures} string value must be 'all' or 'none', got {.val {figures}}."
    )
  }
  if (is.list(figures)) figures <- yaml_chr(figures)
  aggregate_units <- NULL
  if (!is.null(x$aggregate_units)) {
    aggregate_units <- lapply(seq_along(x$aggregate_units), function(i) {
      a <- x$aggregate_units[[i]]
      apath <- glue::glue("{path}.aggregate_units[{i}]")
      check_known_keys(
        a,
        c("unit", "include_treated_units", "include_all_units",
          "include_units", "exclude_units", "strict", "overwrite"),
        apath
      )
      countsynth_aggregate_unit(
        unit = a$unit,
        include_treated_units =
          yaml_flag(a$include_treated_units, glue::glue("{apath}.include_treated_units"), FALSE),
        include_all_units =
          yaml_flag(a$include_all_units, glue::glue("{apath}.include_all_units"), FALSE),
        include_units = yaml_chr(a$include_units),
        exclude_units = yaml_chr(a$exclude_units),
        strict = yaml_flag(a$strict, glue::glue("{apath}.strict"), FALSE),
        overwrite = yaml_flag(a$overwrite, glue::glue("{apath}.overwrite"), FALSE)
      )
    })
  }
  countsynth_output_opts(
    figures = figures,
    clean = yaml_flag(x$clean, glue::glue("{path}.clean"), FALSE),
    save_traces = yaml_flag(x$save_traces, glue::glue("{path}.save_traces"), FALSE),
    target_unit = x$target_unit,
    report_groups = yaml_chr(x$report_groups),
    fit_gap_per_unit =
      yaml_flag(x$fit_gap_per_unit, glue::glue("{path}.fit_gap_per_unit"), FALSE),
    interval_aggregates =
      yaml_flag(x$interval_aggregates, glue::glue("{path}.interval_aggregates"), TRUE),
    print_tables = yaml_flag(x$print_tables, glue::glue("{path}.print_tables"), TRUE),
    print_target_table =
      yaml_flag(x$print_target_table, glue::glue("{path}.print_target_table"), FALSE),
    html_tables = yaml_flag(x$html_tables, glue::glue("{path}.html_tables"), TRUE),
    aggregate_units = aggregate_units,
    ppc_units = yaml_chr(x$ppc_units),
    ppc_exclude_units = yaml_chr(x$ppc_exclude_units),
    ppc_acf_lags =
      if (is.null(x$ppc_acf_lags)) NULL else unlist(x$ppc_acf_lags, use.names = FALSE),
    ppc_unit_corr_max_time = x$ppc_unit_corr_max_time,
    draws_format = x$draws_format %||% "csv",
    rate_normalizer = x$rate_normalizer %||% 1000,
    denominator_label = x$denominator_label %||% "denominator",
    denominator_time_unit = x$denominator_time_unit %||% "year",
    denominator_may_be_affected = yaml_flag(
      x$denominator_may_be_affected,
      glue::glue("{path}.denominator_may_be_affected"), TRUE
    )
  )
}

parse_yaml_cut <- function(x, path = "cut") {
  if (is.null(x)) {
    return(NULL)
  }
  check_known_keys(
    x,
    c("num_stage1_draws", "stage2_draws_per_component", "selection_seed",
      "stage2_seed", "stage2_mcmc", "stage1_method", "stage1_variational"),
    path
  )
  if (!is.null(x$stage2_mcmc)) {
    check_known_keys(
      x$stage2_mcmc,
      c("num_warmup", "num_samples", "thinning", "target_accept",
        "progress_bar", "random_seed"),
      glue::glue("{path}.stage2_mcmc")
    )
  }
  # [[ ]], per the note in parse_yaml_schema: `$` partial matching would let
  # a key resolve to a longer neighbour.
  advi <- x[["stage1_variational"]]
  if (!is.null(advi)) {
    check_known_keys(advi, VARIATIONAL_KEYS, glue::glue("{path}.stage1_variational"))
    if (!is.null(advi[["adapt_engaged"]])) {
      advi[["adapt_engaged"]] <- yaml_flag(
        advi[["adapt_engaged"]],
        glue::glue("{path}.stage1_variational.adapt_engaged")
      )
    }
  }
  countsynth_cut_opts(
    num_stage1_draws = x$num_stage1_draws %||% 25L,
    stage2_draws_per_component =
      if ("stage2_draws_per_component" %in% names(x)) {
        x$stage2_draws_per_component
      } else {
        100L
      },
    selection_seed = x$selection_seed,
    stage2_seed = x$stage2_seed,
    stage2_mcmc = x$stage2_mcmc,
    stage1_method = x[["stage1_method"]] %||% "sample",
    stage1_variational = advi
  )
}

#' Read a countsynth YAML config
#'
#' Accepts the Python `bayesian_panel_nmf` config schema verbatim (Python key
#' names like `num_warmup` / `target_accept` / `random_seed` are translated to
#' their cmdstanr equivalents). Unknown keys at any level are an error, and
#' quoted string booleans (`"false"`) are rejected where a boolean is
#' required, matching the Python loader's strictness.
#'
#' Relative `data.input_file` and `data.output_dir` paths are resolved
#' relative to the current working directory (same as the Python CLI).
#'
#' @param path Path to a YAML config file.
#' @return A validated [countsynth_config()] object.
#' @export
read_countsynth_config <- function(path) {
  checkmate::assert_file_exists(path, extension = c("yaml", "yml"))
  raw <- yaml::read_yaml(path)
  check_known_keys(raw, c("data", "model", "mcmc", "output", "cut"), "<top level>")
  if (is.null(raw$data)) {
    cli::cli_abort("Config must have a {.field data} section.")
  }
  d <- raw$data
  check_known_keys(
    d,
    c("input_file", "output_dir", "schema", "date_format", "start_date",
      "end_date", "time_aggregation", "aggregation",
      "allow_unbalanced_panel", "outcome"),
    "data"
  )
  if (is.null(d$schema)) {
    cli::cli_abort("Config must define {.field data.schema}.")
  }
  # `aggregation` was the key's name before it was made specific; still read
  # so existing configs keep loading, but only one of the two may be set.
  if (!is.null(d$aggregation)) {
    if (!is.null(d$time_aggregation)) {
      cli::cli_abort(
        "Config sets both {.field data.time_aggregation} and the deprecated
         {.field data.aggregation}; keep only {.field data.time_aggregation}."
      )
    }
    cli::cli_warn(
      "{.field data.aggregation} is deprecated; rename it to
       {.field data.time_aggregation}."
    )
    d$time_aggregation <- d$aggregation
  }
  time_aggregation <- countsynth_time_aggregation()
  if (!is.null(d$time_aggregation)) {
    ta <- d$time_aggregation
    check_known_keys(ta, c("enabled", "period", "n_periods"), "data.time_aggregation")
    time_aggregation <- countsynth_time_aggregation(
      enabled = yaml_flag(ta$enabled, "data.time_aggregation.enabled", FALSE),
      period = ta$period,
      n_periods = ta$n_periods
    )
  }
  countsynth_config(
    input_file = d$input_file,
    output_dir = d$output_dir,
    schema = parse_yaml_schema(d$schema),
    model = parse_yaml_model(raw$model),
    mcmc = parse_yaml_mcmc(raw$mcmc),
    output = parse_yaml_output(raw$output),
    cut = parse_yaml_cut(raw$cut),
    date_format = d$date_format %||% "auto",
    start_date = d$start_date,
    end_date = d$end_date,
    time_aggregation = time_aggregation,
    allow_unbalanced_panel =
      yaml_flag(d$allow_unbalanced_panel, "data.allow_unbalanced_panel", FALSE),
    outcome = d$outcome
  )
}
