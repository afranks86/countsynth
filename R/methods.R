# S3 methods for fit objects.

#' @export
print.countsynth_fit <- function(x, ...) {
  cli::cli_h1("countsynth fit ({x$type}, rank {x$rank})")
  cli::cli_li("model: {if (x$model_treated) 'joint (treated)' else 'baseline (untreated)'}")
  cli::cli_li("method: {fit_method_label(x$fit_method %||% 'sample')}")
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
summary.countsynth_fit <- function(object, ...) {
  print(object)
  gate <- convergence_gate(object)
  if (is.na(gate$converged)) {
    cli::cli_li("gate: not applicable -- ADVI has no R-hat, ESS or divergences")
    return(invisible(gate))
  }
  status <- if (isTRUE(gate$converged)) "PASS" else "FAIL"
  cli::cli_li(
    "gate: {status} (max R-hat {round(gate$rhat_max, 4)}, min bulk ESS
     {round(gate$ess_bulk_min)}, {gate$divergences} divergence{?s} =
     {sprintf('%.2f%%', 100 * gate$divergence_fraction)})"
  )
  invisible(gate)
}

#' @export
print.countsynth_cut_fit <- function(x, ...) {
  cli::cli_h1("countsynth cut fit ({x$type}, rank {x$rank})")
  n_comp <- length(x$component_records)
  n_pass <- sum(vapply(
    x$component_records, function(r) isTRUE(r$converged), logical(1)
  ))
  cli::cli_li("{n_comp} stage-2 component{?s}, {n_pass} converged")
  s1 <- x$manifest$stage1
  if (is.na(s1$converged)) {
    cli::cli_li("stage 1: {fit_method_label('variational')}, not gated")
  } else {
    cli::cli_li(
      "stage 1: {if (isTRUE(s1$converged)) 'PASS' else 'FAIL'}
       (max R-hat {round(s1$rhat_max, 4)})"
    )
  }
  cli::cli_li("{nrow(x$draws)} pooled draw row{?s}")
  invisible(x)
}

#' @export
summary.countsynth_cut_fit <- function(object, ...) {
  print(object)
  cut_component_table(object)
  invisible(object$manifest)
}

#' Plot method for countsynth fits
#'
#' @param x A `countsynth_fit` or `countsynth_cut_fit`.
#' @param which One of `"unit_fit"`, `"unit_gap"`, `"raw_rate"`,
#'   `"group_comparison"`, `"interval"`, `"trace"`, `"te_regression"`, or
#'   `"te_coef"` (the last two require a `treatment_effects` formula).
#' @param unit Unit for the per-unit figures (auto-detected target when
#'   `NULL`).
#' @param group Group for the per-unit figures.
#' @param ... Passed to the underlying figure function.
#' @export
plot.countsynth_fit <- function(x, which = "unit_fit", unit = NULL, group = NULL,
                           ...) {
  checkmate::assert_choice(
    which,
    c("unit_fit", "unit_gap", "raw_rate", "group_comparison", "interval",
      "trace", "te_regression", "te_coef")
  )
  if (which == "trace") {
    return(countsynth_trace_plot(x, ...))
  }
  # The treatment-effect figures read the design and coefficient draws off
  # the fit, not the tidy draws frame.
  if (which == "te_regression") {
    return(countsynth_te_regression_plot(x, ...))
  }
  if (which == "te_coef") {
    return(countsynth_te_coef_plot(x, ...))
  }
  draws <- countsynth_draws(x)
  unit <- unit %||% auto_detect_target(draws)
  switch(which,
    unit_fit = countsynth_unit_fit_plot(draws, unit, group),
    unit_gap = countsynth_unit_gap_plot(draws, unit, group),
    raw_rate = countsynth_raw_rate_plot(draws, group = group, separate_unit = unit, ...),
    group_comparison = countsynth_group_comparison_plot(draws, ...),
    interval = countsynth_interval_plot(draws, ...)
  )
}

#' @export
plot.countsynth_cut_fit <- plot.countsynth_fit
