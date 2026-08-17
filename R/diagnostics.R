# Convergence diagnostics: R-hat / ESS gate + per-parameter table. Port of
# diagnostics.py. R-hat and ESS come from the posterior package (the same
# rank-normalized Vehtari et al. 2021 definitions ArviZ uses); divergences
# are always counted over the full run regardless of gate_params.

#' PASS / WARN / FAIL for one (rhat, ess) pair
#'
#' FAIL at/above `rhat_fail` or below `ess_min * ess_fail_fraction`; WARN
#' at/above `rhat_warn` or below `ess_min`; PASS otherwise.
#' @keywords internal
convergence_status <- function(rhat, ess, thresholds) {
  if (rhat >= thresholds$rhat_fail ||
    ess < thresholds$ess_min * thresholds$ess_fail_fraction) {
    return("FAIL")
  }
  if (rhat >= thresholds$rhat_warn || ess < thresholds$ess_min) {
    return("WARN")
  }
  "PASS"
}

# Variables never gated: lp__ is the sampler's own quantity, ypred is a
# generated-quantities RNG draw, not a posterior parameter.
DIAG_EXCLUDE <- c("lp__", "ypred")

diag_variables <- function(fit) {
  vars <- fit$metadata()$stan_variables
  setdiff(vars, DIAG_EXCLUDE)
}

# Per-element diagnostics for the selected variables, with base-name grouping
# and empirical fixed-site detection (zero span across all chains/draws).
element_diagnostics <- function(fit, variables) {
  draws <- fit$draws(variables = variables)
  summ <- posterior::summarise_draws(
    draws,
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail,
    span = \(x) max(x) - min(x)
  )
  summ$base <- sub("\\[.*$", "", summ$variable)
  summ
}

#' Match gate_params prefixes against base variable names
#' @keywords internal
match_gate_params <- function(base_names, gate_params) {
  if (is.null(gate_params)) {
    return(rep(TRUE, length(base_names)))
  }
  gated <- Reduce(
    `|`,
    lapply(gate_params, function(p) startsWith(base_names, p))
  )
  if (!any(gated)) {
    cli::cli_abort(c(
      "mcmc.gate_params {.val {gate_params}} matched no posterior variables.",
      i = "Available variables: {.val {sort(unique(base_names))}}."
    ))
  }
  gated
}

count_divergences <- function(fit) {
  sd <- fit$sampler_diagnostics()
  sum(posterior::as_draws_matrix(
    posterior::subset_draws(sd, variable = "divergent__")
  ))
}

#' Run-level convergence gate
#'
#' Port of `diagnostics.convergence_summary`: worst R-hat and smallest
#' bulk/tail ESS over the gated parameters, plus the run-wide divergence
#' count. `converged` requires PASS-band status and zero divergences.
#'
#' @param fit A `bpnmf_fit` / `bpnmf_cut_fit`, or a raw `CmdStanMCMC`.
#' @param gate_params Optional character vector of variable-name prefixes the
#'   R-hat/ESS gate is restricted to (defaults to the config's
#'   `mcmc$gate_params` when `fit` is a `bpnmf_fit`).
#' @param thresholds A [bpnmf_convergence()] object.
#' @return A list: `rhat_max`, `ess_bulk_min`, `ess_tail_min`, `divergences`,
#'   `converged` (+ `gate_params` when set).
#' @export
convergence_gate <- function(fit, gate_params = NULL, thresholds = NULL) {
  if (inherits(fit, "bpnmf_fit") || inherits(fit, "bpnmf_cut_fit")) {
    gate_params <- gate_params %||% fit$config$mcmc$gate_params
    thresholds <- thresholds %||% fit$config$mcmc$convergence
    fit <- fit$fit
  }
  thresholds <- thresholds %||% bpnmf_convergence()

  summ <- element_diagnostics(fit, diag_variables(fit))
  # Fixed elements (zero posterior span) carry no convergence information.
  summ <- summ[summ$span > 0 | is.na(summ$span), ]
  gated <- match_gate_params(summ$base, gate_params)
  gated_summ <- summ[gated, ]

  rhat_max <- max(gated_summ$rhat, na.rm = TRUE)
  ess_bulk_min <- min(gated_summ$ess_bulk, na.rm = TRUE)
  ess_tail_min <- min(gated_summ$ess_tail, na.rm = TRUE)
  divergences <- count_divergences(fit)
  status <- convergence_status(
    rhat_max, min(ess_bulk_min, ess_tail_min), thresholds
  )

  out <- list(
    rhat_max = rhat_max,
    ess_bulk_min = ess_bulk_min,
    ess_tail_min = ess_tail_min,
    divergences = divergences,
    converged = status == "PASS" && divergences == 0
  )
  if (!is.null(gate_params)) {
    out$gate_params <- gate_params
  }
  out
}

#' Per-parameter diagnostics table
#'
#' One row per base variable: worst R-hat, smallest `min(ess_bulk, ess_tail)`
#' across its elements, PASS/WARN/FAIL status, and whether the variable is
#' gated. Variables whose every element has zero posterior span are tagged
#' `"fixed"` (no R-hat/ESS, never FAIL) — this catches a fixed dispersion
#' empirically, with no model-specific knowledge. Sorted worst R-hat first.
#'
#' @inheritParams convergence_gate
#' @return A tibble with columns `parameter`, `rhat`, `ess`, `status`,
#'   `gated`.
#' @export
parameter_diagnostics <- function(fit, gate_params = NULL, thresholds = NULL) {
  if (inherits(fit, "bpnmf_fit") || inherits(fit, "bpnmf_cut_fit")) {
    gate_params <- gate_params %||% fit$config$mcmc$gate_params
    thresholds <- thresholds %||% fit$config$mcmc$convergence
    fit <- fit$fit
  }
  thresholds <- thresholds %||% bpnmf_convergence()

  summ <- element_diagnostics(fit, diag_variables(fit))
  gated_by_base <- match_gate_params(summ$base, gate_params)
  summ$gated <- gated_by_base

  rows <- summ |>
    dplyr::group_by(.data$base) |>
    dplyr::summarise(
      fixed = all(.data$span == 0, na.rm = TRUE),
      rhat = if (all(.data$span == 0, na.rm = TRUE)) {
        NA_real_
      } else {
        max(.data$rhat[.data$span > 0], na.rm = TRUE)
      },
      ess = if (all(.data$span == 0, na.rm = TRUE)) {
        NA_real_
      } else {
        min(
          pmin(
            .data$ess_bulk[.data$span > 0],
            .data$ess_tail[.data$span > 0]
          ),
          na.rm = TRUE
        )
      },
      gated = any(.data$gated),
      .groups = "drop"
    ) |>
    dplyr::rename(parameter = "base")

  rows$status <- vapply(
    seq_len(nrow(rows)),
    function(i) {
      if (rows$fixed[i]) {
        "fixed"
      } else {
        convergence_status(rows$rhat[i], rows$ess[i], thresholds)
      }
    },
    character(1)
  )
  rows$fixed <- NULL
  rows <- rows[order(-ifelse(is.na(rows$rhat), -Inf, rows$rhat)), ]
  class(rows) <- c("bpnmf_diagnostics", class(rows))
  attr(rows, "thresholds") <- thresholds
  rows
}

#' @export
print.bpnmf_diagnostics <- function(x, ...) {
  style <- c(
    PASS = "green", WARN = "yellow", FAIL = "red", fixed = "silver"
  )
  cli::cli_h1("Parameter diagnostics")
  fmt <- "%-28s %9s %9s  %-5s %s"
  cli::cli_verbatim(sprintf(fmt, "parameter", "max R-hat", "min ESS", "status", "gate"))
  for (i in seq_len(nrow(x))) {
    line <- sprintf(
      fmt,
      x$parameter[i],
      ifelse(is.na(x$rhat[i]), "-", sprintf("%.4f", x$rhat[i])),
      ifelse(is.na(x$ess[i]), "-", sprintf("%.0f", x$ess[i])),
      x$status[i],
      if (x$gated[i]) "✓" else ""
    )
    switch(x$status[i],
      FAIL = cli::cli_verbatim(cli::col_red(line)),
      WARN = cli::cli_verbatim(cli::col_yellow(line)),
      fixed = cli::cli_verbatim(cli::col_silver(line)),
      cli::cli_verbatim(line)
    )
  }
  invisible(x)
}

#' Write a convergence gate as JSON (artifact parity with Python)
#' @param gate A list from [convergence_gate()] or a cut manifest.
#' @param path Output file path.
#' @export
write_convergence_json <- function(gate, path) {
  jsonlite::write_json(gate, path, auto_unbox = TRUE, pretty = TRUE, digits = 8)
  invisible(path)
}
