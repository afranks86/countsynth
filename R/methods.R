# S3 methods for fit objects.

#' @export
print.bpnmf_fit <- function(x, ...) {
  cli::cli_h1("bpnmf fit ({x$type}, rank {x$rank})")
  cli::cli_li("model: {if (x$model_treated) 'joint (treated)' else 'baseline (untreated)'}")
  cli::cli_li(
    "distribution: {x$config$model$outcome_distribution}, censoring adjustment:
     {x$config$model$adjust_for_missingness}"
  )
  cli::cli_li(
    "{length(x$data$groups)} group{?s} x {length(x$data$units)} unit{?s} x
     {length(x$data$times)} time{?s}; {x$stan_data$n_exposed} exposed cell{?s}"
  )
  invisible(x)
}

#' @export
summary.bpnmf_fit <- function(object, ...) {
  print(object)
  gate <- convergence_gate(object)
  status <- if (isTRUE(gate$converged)) "PASS" else "FAIL"
  cli::cli_li(
    "gate: {status} (max R-hat {round(gate$rhat_max, 4)}, min bulk ESS
     {round(gate$ess_bulk_min)}, {gate$divergences} divergence{?s})"
  )
  invisible(gate)
}

#' @export
print.bpnmf_cut_fit <- function(x, ...) {
  cli::cli_h1("bpnmf cut fit ({x$type}, rank {x$rank})")
  n_comp <- length(x$component_records)
  n_pass <- sum(vapply(
    x$component_records, function(r) isTRUE(r$converged), logical(1)
  ))
  cli::cli_li("{n_comp} stage-2 component{?s}, {n_pass} converged")
  cli::cli_li(
    "stage 1: {if (isTRUE(x$manifest$stage1$converged)) 'PASS' else 'FAIL'}
     (max R-hat {round(x$manifest$stage1$rhat_max, 4)})"
  )
  cli::cli_li("{nrow(x$draws)} pooled draw row{?s}")
  invisible(x)
}

#' @export
summary.bpnmf_cut_fit <- function(object, ...) {
  print(object)
  cut_component_table(object)
  invisible(object$manifest)
}

#' Plot method for bpnmf fits
#'
#' @param x A `bpnmf_fit` or `bpnmf_cut_fit`.
#' @param which One of `"unit_fit"`, `"unit_gap"`, `"raw_rate"`,
#'   `"group_comparison"`, `"interval"`, `"trace"`, `"te_regression"`, or
#'   `"te_coef"` (the last two require a `treatment_effects` formula).
#' @param unit Unit for the per-unit figures (auto-detected target when
#'   `NULL`).
#' @param group Group for the per-unit figures.
#' @param ... Passed to the underlying figure function.
#' @export
plot.bpnmf_fit <- function(x, which = "unit_fit", unit = NULL, group = NULL,
                           ...) {
  checkmate::assert_choice(
    which,
    c("unit_fit", "unit_gap", "raw_rate", "group_comparison", "interval",
      "trace", "te_regression", "te_coef")
  )
  if (which == "trace") {
    return(bpnmf_trace_plot(x, ...))
  }
  # The treatment-effect figures read the design and coefficient draws off
  # the fit, not the tidy draws frame.
  if (which == "te_regression") {
    return(bpnmf_te_regression_plot(x, ...))
  }
  if (which == "te_coef") {
    return(bpnmf_te_coef_plot(x, ...))
  }
  draws <- bpnmf_draws(x)
  unit <- unit %||% auto_detect_target(draws)
  switch(which,
    unit_fit = bpnmf_unit_fit_plot(draws, unit, group),
    unit_gap = bpnmf_unit_gap_plot(draws, unit, group),
    raw_rate = bpnmf_raw_rate_plot(draws, group = group, separate_unit = unit, ...),
    group_comparison = bpnmf_group_comparison_plot(draws, ...),
    interval = bpnmf_interval_plot(draws, ...)
  )
}

#' @export
plot.bpnmf_cut_fit <- plot.bpnmf_fit
