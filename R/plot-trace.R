# Trace plots for MCMC diagnostics.

#' Trace plot for selected parameters
#'
#' Uses \pkg{bayesplot} when installed; otherwise a plain ggplot of chains
#' over iterations, faceted by parameter. Elements are capped at
#' `max_elements` (highest-variance elements kept, mirroring the Python
#' trace CLI's informative-coordinate subsampling).
#'
#' @param fit A `bpnmf_fit` / `bpnmf_cut_fit`, or a raw `CmdStanMCMC`.
#' @param variables Variable names (prefixes match whole variables); defaults
#'   to the scalar treatment/scale parameters.
#' @param max_elements Cap on the number of parameter elements plotted.
#' @return A ggplot object.
#' @export
bpnmf_trace_plot <- function(fit, variables = NULL, max_elements = 20) {
  if (inherits(fit, "bpnmf_fit") || inherits(fit, "bpnmf_cut_fit")) {
    fit <- fit$fit
  }
  if (is_variational_fit(fit)) {
    cli::cli_abort(c(
      "Trace plots need MCMC chains.",
      i = "This fit came from ADVI ({.code method = \"variational\"}); its
           draws are independent samples from the fitted approximation, so a
           trace over them shows nothing about convergence."
    ))
  }
  all_vars <- setdiff(fit$metadata()$stan_variables, DIAG_EXCLUDE)
  defaults <- intersect(
    c("disp", "treatment_it_scale", "treatment_state_scale",
      "treatment_category_scale", "state_category_scale", "state_fe_mu",
      "state_fe_sigma", "category_treatment_effect"),
    all_vars
  )
  variables <- variables %||% defaults
  variables <- intersect(variables, all_vars)
  if (length(variables) == 0) {
    cli::cli_abort("No matching posterior variables to plot.")
  }

  draws <- fit$draws(variables = variables)
  # Keep the max_elements highest-variance elements: structurally-zero or
  # fixed cells produce empty traces.
  variances <- apply(posterior::as_draws_matrix(draws), 2, stats::var)
  keep <- names(sort(variances, decreasing = TRUE))
  keep <- utils::head(keep[variances[keep] > 1e-12], max_elements)
  if (length(keep) == 0) {
    keep <- utils::head(names(variances), max_elements)
  }
  draws <- posterior::subset_draws(draws, variable = keep)

  if (requireNamespace("bayesplot", quietly = TRUE)) {
    return(bayesplot::mcmc_trace(draws) + theme_bpnmf())
  }
  long <- posterior::as_draws_df(draws) |>
    tidyr::pivot_longer(
      cols = -dplyr::all_of(c(".chain", ".iteration", ".draw")),
      names_to = "parameter", values_to = "value"
    )
  ggplot2::ggplot(
    long,
    ggplot2::aes(
      x = .data$.iteration, y = .data$value,
      color = factor(.data$.chain)
    )
  ) +
    ggplot2::geom_line(linewidth = 0.3, alpha = 0.8) +
    ggplot2::facet_wrap(ggplot2::vars(.data$parameter), scales = "free_y") +
    ggplot2::labs(x = "Iteration", y = NULL, color = "Chain") +
    theme_bpnmf()
}
