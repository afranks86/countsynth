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

# Gated by default: the counterfactual log-rate surface and the treatment
# effect -- the two quantities every reported estimand is built from. The
# factor parameters (time_fac, unit_weight, unit_fe_*, and under rank
# shrinkage stick / group_weight / unit_weight_z, plus group_profile_z /
# unit_shared_z with shared curves) are deliberately left out:
# an NMF is invariant to permuting and rescaling its factors, so chains that
# settle on different labelings give those parameters an enormous R-hat while
# mu_ctrl and te are converged. Gating on them measures label disagreement,
# not convergence of anything reported. (unit_weight_z carries a genuinely
# flat direction on top of that -- see the softmax note in joint.stan -- and
# the rank diagnostics are built only from permutation-invariant functions:
# eff_rank, and the profile's order statistics.)
#
# Matching is by prefix over base names, and neither prefix over-matches:
# "te" does not reach `treatment_*` (which begins "tr"), and "mu_ctrl" is
# exact. parameter_diagnostics() still reports every variable, so the
# ungated ones stay visible.
DEFAULT_GATE_PARAMS <- c("mu_ctrl", "te")

#' Was this fit produced by ADVI rather than NUTS?
#'
#' R-hat, ESS and divergences are all chain-based, and ADVI produces one
#' stream of independent draws from an approximating family -- so none of them
#' exist for a variational fit, and computing them anyway (split-R-hat over a
#' single chain) would report a number that means nothing.
#' @param fit A `countsynth_fit` / `countsynth_cut_fit`, or a raw CmdStan fit object.
#' @keywords internal
is_variational_fit <- function(fit) {
  if (inherits(fit, "countsynth_fit") || inherits(fit, "countsynth_cut_fit")) {
    fit <- fit$fit
  }
  inherits(fit, "CmdStanVB")
}

#' Human-readable name for a fit method
#' @keywords internal
fit_method_label <- function(method) {
  if (identical(method, "variational")) "ADVI (variational)" else "NUTS (MCMC)"
}

# Gate-shaped result for a variational fit. The fields are kept so every
# consumer -- manifests, the convergence JSON, printing -- sees one shape;
# `converged = NA` means "not gated", which is distinct from FALSE ("gated
# and failed") and must not be collapsed into it.
variational_gate <- function() {
  list(
    method = "variational",
    rhat_max = NA_real_,
    ess_bulk_min = NA_real_,
    ess_tail_min = NA_real_,
    divergences = NA_integer_,
    divergence_fraction = NA_real_,
    treedepth_hits = NA_integer_,
    treedepth_fraction = NA_real_,
    max_treedepth = NA_integer_,
    converged = NA
  )
}

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
  # NULL for a direct call with no config; "all" is the opt-out that restores
  # gating on every diagnosable variable.
  if (is.null(gate_params) || identical(gate_params, "all")) {
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

# Divergent transitions over the retained draws, with the denominator needed
# to turn the count into a rate. Counting retained rather than all post-warmup
# transitions keeps the rate comparable to the ESS figures, which are also
# computed on the retained draws.
divergence_summary <- function(fit) {
  sd <- posterior::as_draws_matrix(
    posterior::subset_draws(fit$sampler_diagnostics(), variable = "divergent__")
  )
  transitions <- nrow(sd)
  list(
    count = sum(sd),
    transitions = transitions,
    fraction = if (transitions > 0) sum(sd) / transitions else 0
  )
}

# Treedepth saturation over the retained draws. A trajectory cut off at the
# depth ceiling is still a valid Hamiltonian proposal -- it costs mixing
# efficiency, not correctness -- so unlike a divergence this is never gated
# (see convergence_gate()) and never fails a run. It is tracked here purely
# so countsynth can report it with that context instead of leaving cmdstanr's own
# unglossed "N transitions hit the maximum treedepth" warning as the only
# signal a user sees.
treedepth_summary <- function(fit) {
  sd <- posterior::as_draws_matrix(
    posterior::subset_draws(fit$sampler_diagnostics(), variable = "treedepth__")
  )
  transitions <- nrow(sd)
  # cmdstan's own default when a fit doesn't expose it back (older cmdstanr).
  max_treedepth <- fit$metadata()$max_treedepth %||% 10L
  hits <- sum(sd >= max_treedepth)
  list(
    hits = hits,
    transitions = transitions,
    fraction = if (transitions > 0) hits / transitions else 0,
    max_treedepth = max_treedepth
  )
}

#' Run-level convergence gate
#'
#' Port of `diagnostics.convergence_summary`: worst R-hat and smallest
#' bulk/tail ESS over the gated parameters, plus the run-wide divergence
#' count. `converged` requires PASS-band status and a divergence rate at or
#' below `divergence_fail_fraction` (0 for that threshold restores the older
#' zero-divergence rule).
#'
#' @param fit A `countsynth_fit` / `countsynth_cut_fit`, or a raw `CmdStanMCMC`.
#' @param gate_params Optional character vector of variable-name prefixes the
#'   R-hat/ESS gate is restricted to (defaults to the config's
#'   `mcmc$gate_params` when `fit` is a `countsynth_fit`).
#' @param thresholds A [countsynth_convergence()] object.
#' @return A list: `rhat_max`, `ess_bulk_min`, `ess_tail_min`, `divergences`,
#'   `divergence_fraction`, `treedepth_hits`, `treedepth_fraction`,
#'   `max_treedepth`, `converged` (+ `gate_params` when set). Treedepth is
#'   informational only -- it is never part of `converged` (see
#'   [treedepth_summary()]). For a variational (ADVI) fit none of those
#'   quantities exist, so they are `NA` and `converged` is `NA` -- "not
#'   gated", not "failed".
#' @export
convergence_gate <- function(fit, gate_params = NULL, thresholds = NULL) {
  if (inherits(fit, "countsynth_fit") || inherits(fit, "countsynth_cut_fit")) {
    gate_params <- gate_params %||% fit$config$mcmc$gate_params
    thresholds <- thresholds %||% fit$config$mcmc$convergence
    fit <- fit$fit
  }
  if (is_variational_fit(fit)) {
    return(variational_gate())
  }
  thresholds <- thresholds %||% countsynth_convergence()
  gate_params <- gate_params %||% DEFAULT_GATE_PARAMS

  summ <- element_diagnostics(fit, diag_variables(fit))
  # Fixed elements (zero posterior span) carry no convergence information.
  summ <- summ[summ$span > 0 | is.na(summ$span), ]
  gated <- match_gate_params(summ$base, gate_params)
  gated_summ <- summ[gated, ]

  rhat_max <- max(gated_summ$rhat, na.rm = TRUE)
  ess_bulk_min <- min(gated_summ$ess_bulk, na.rm = TRUE)
  ess_tail_min <- min(gated_summ$ess_tail, na.rm = TRUE)
  div <- divergence_summary(fit)
  td <- treedepth_summary(fit)
  status <- convergence_status(
    rhat_max, min(ess_bulk_min, ess_tail_min), thresholds
  )

  out <- list(
    rhat_max = rhat_max,
    ess_bulk_min = ess_bulk_min,
    ess_tail_min = ess_tail_min,
    divergences = div$count,
    divergence_fraction = div$fraction,
    treedepth_hits = td$hits,
    treedepth_fraction = td$fraction,
    max_treedepth = td$max_treedepth,
    converged = status == "PASS" &&
      div$fraction <= thresholds$divergence_fail_fraction
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
#' `"fixed"` (no R-hat/ESS, never FAIL); this catches a fixed dispersion
#' empirically, with no model-specific knowledge. Sorted worst R-hat first.
#'
#' @inheritParams convergence_gate
#' @return A tibble with columns `parameter`, `rhat`, `ess`, `status`,
#'   `gated`.
#' @export
parameter_diagnostics <- function(fit, gate_params = NULL, thresholds = NULL) {
  if (inherits(fit, "countsynth_fit") || inherits(fit, "countsynth_cut_fit")) {
    gate_params <- gate_params %||% fit$config$mcmc$gate_params
    thresholds <- thresholds %||% fit$config$mcmc$convergence
    fit <- fit$fit
  }
  if (is_variational_fit(fit)) {
    cli::cli_abort(c(
      "Per-parameter R-hat / ESS diagnostics need multiple MCMC chains.",
      i = "This fit came from ADVI ({.code method = \"variational\"}), which
           has neither.",
      i = "Re-fit with {.code method = \"sample\"} to diagnose convergence."
    ))
  }
  thresholds <- thresholds %||% countsynth_convergence()
  gate_params <- gate_params %||% DEFAULT_GATE_PARAMS

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
  class(rows) <- c("countsynth_diagnostics", class(rows))
  attr(rows, "thresholds") <- thresholds
  rows
}

#' @export
print.countsynth_diagnostics <- function(x, ...) {
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
      if (x$gated[i]) "\u2713" else ""
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

# Why a gate came back FALSE, as cli bullets. A bare "FAILED" leaves the user
# to go dig the numbers out of the JSON; naming the tripped criterion, the
# threshold it tripped, and the worst offending parameters answers the
# question on the spot.
#
# The gate demands a clean PASS (see convergence_gate()), so a WARN-level
# R-hat or ESS fails it just as a FAIL-level one does. Reporting only the
# FAIL thresholds would then print "FAILED" next to an R-hat of 1.05 and
# explain nothing, so each criterion reports the level it actually reached.
gate_failure_bullets <- function(gate, thresholds = NULL, fit = NULL,
                                 max_params = 3L) {
  thresholds <- thresholds %||% countsynth_convergence()
  bullets <- character()
  warn_only <- TRUE

  rhat <- gate$rhat_max
  if (isTRUE(rhat >= thresholds$rhat_fail)) {
    warn_only <- FALSE
    bullets <- c(bullets, sprintf(
      "max R-hat %.3g, at or above the fail threshold %.3g",
      rhat, thresholds$rhat_fail
    ))
  } else if (isTRUE(rhat >= thresholds$rhat_warn)) {
    bullets <- c(bullets, sprintf(
      "max R-hat %.3g, at or above the warn threshold %.3g",
      rhat, thresholds$rhat_warn
    ))
  }

  ess <- suppressWarnings(min(gate$ess_bulk_min, gate$ess_tail_min, na.rm = TRUE))
  ess_floor <- thresholds$ess_min * thresholds$ess_fail_fraction
  if (isTRUE(ess < ess_floor)) {
    warn_only <- FALSE
    bullets <- c(bullets, sprintf(
      "min ESS %.3g, below the fail floor %.3g (ess_min %.3g x ess_fail_fraction %.3g)",
      ess, ess_floor, thresholds$ess_min, thresholds$ess_fail_fraction
    ))
  } else if (isTRUE(ess < thresholds$ess_min)) {
    bullets <- c(bullets, sprintf(
      "min ESS %.3g, below ess_min %.3g", ess, thresholds$ess_min
    ))
  }

  if (isTRUE(gate$divergence_fraction > thresholds$divergence_fail_fraction)) {
    warn_only <- FALSE
    bullets <- c(bullets, sprintf(
      "%d divergence%s = %.2f%%, above divergence_fail_fraction %.2f%% --
       a fraction this high usually means some region of the posterior is
       poorly explored, so treat estimates cautiously; raising
       {.field mcmc.target_accept} (adapt_delta) toward 0.95-0.99 or
       simplifying the model often resolves it",
      gate$divergences, if (isTRUE(gate$divergences == 1)) "" else "s",
      100 * gate$divergence_fraction,
      100 * thresholds$divergence_fail_fraction
    ))
  } else if (isTRUE(gate$divergences > 0)) {
    bullets <- c(bullets, sprintf(
      "%d divergence%s = %.2f%%, within the %.2f%% allowance -- a small,
       isolated count like this is common and usually not a validity
       concern on its own",
      gate$divergences, if (isTRUE(gate$divergences == 1)) "" else "s",
      100 * gate$divergence_fraction,
      100 * thresholds$divergence_fail_fraction
    ))
  }

  # A cut manifest has no top-level R-hat/ESS, so nothing above fires.
  if (length(bullets) == 0) {
    return("the gate reported no passing status; see the convergence JSON")
  }
  if (warn_only) {
    bullets <- c(bullets, paste(
      "nothing reached a fail threshold -- the gate requires a clean PASS,",
      "so a WARN-level criterion fails it"
    ))
  }

  worst <- gate_worst_parameters(fit, max_params)
  if (length(worst) > 0) {
    bullets <- c(bullets, sprintf(
      "worst gated parameter%s: %s",
      if (length(worst) == 1) "" else "s", paste(worst, collapse = ", ")
    ))
  }
  bullets
}

# Top offenders among the gated parameters, formatted "name (R-hat x, ESS y)".
gate_worst_parameters <- function(fit, max_params = 3L) {
  if (is.null(fit)) {
    return(character())
  }
  diag <- tryCatch(parameter_diagnostics(fit), error = function(e) NULL)
  if (is.null(diag) || nrow(diag) == 0) {
    return(character())
  }
  bad <- diag[diag$gated & diag$status == "FAIL" & !is.na(diag$rhat), ]
  if (nrow(bad) == 0) {
    return(character())
  }
  bad <- utils::head(bad, max_params)
  sprintf("%s (R-hat %.3g, ESS %.3g)", bad$parameter, bad$rhat, bad$ess)
}

# Concrete next steps for a failed gate, keyed off which criterion actually
# tripped -- gate_failure_bullets() says what went wrong; this says what to
# try. Silently returns nothing for a gate missing the relevant fields (e.g.
# a cut manifest, which carries no top-level R-hat/ESS/divergence rate).
gate_failure_advice <- function(gate, thresholds = NULL) {
  thresholds <- thresholds %||% countsynth_convergence()
  advice <- character()

  ess <- suppressWarnings(min(gate$ess_bulk_min, gate$ess_tail_min, na.rm = TRUE))
  if (isTRUE(gate$rhat_max >= thresholds$rhat_warn) ||
    isTRUE(ess < thresholds$ess_min)) {
    advice <- c(advice, paste(
      "R-hat/ESS: try more warmup/sampling iterations first",
      "({.field mcmc.num_warmup} / {.field mcmc.num_samples}); if it persists",
      "at a given rank, a lower rank often mixes faster (an overly high rank",
      "can leave components weakly identified), and",
      "{.code parameter_diagnostics(fit)} shows which parameter is the",
      "bottleneck."
    ))
  }
  if (isTRUE(gate$divergence_fraction > thresholds$divergence_fail_fraction)) {
    advice <- c(advice, paste(
      "Divergences: try raising {.field mcmc.target_accept} (adapt_delta)",
      "toward 0.95-0.99; this trades sampling speed for a smaller step size",
      "and usually clears a high divergence rate."
    ))
  }
  advice
}

# Non-gating diagnostic context, meant to be shown alongside a PASSING gate.
# gate_failure_bullets() already explains any divergences on a FAILING one,
# so the divergence note here only fires when the gate passed -- otherwise a
# harmless divergence count would be narrated twice. Treedepth is never part
# of gate_failure_bullets() at all (it is never gated -- see
# convergence_gate()), so it is reported here unconditionally: this is meant
# to replace cmdstanr's own unglossed "N transitions hit the maximum
# treedepth" warning as the thing a user actually reads, on a pass or a fail.
diagnostic_context_notes <- function(gate, thresholds = NULL) {
  thresholds <- thresholds %||% countsynth_convergence()
  notes <- character()

  if (isTRUE(gate$converged) && isTRUE(gate$divergences > 0)) {
    notes <- c(notes, sprintf(
      "%d divergent transition%s (%.2f%% of retained draws): a small,
       isolated count like this is common and does not usually affect
       validity -- it's only worth a closer look if the count keeps growing
       with more draws, or the fraction reaches several percent or more.",
      gate$divergences, if (isTRUE(gate$divergences == 1)) "" else "s",
      100 * gate$divergence_fraction
    ))
  }

  if (isTRUE(gate$treedepth_hits > 0)) {
    # Not target_accept: raising it shrinks the step size, so trajectories
    # need more steps and hit the cap more often. That is the remedy for
    # divergences, and it makes this worse.
    action <- if (isTRUE(gate$treedepth_fraction >= 0.1)) {
      "This share is large enough to be slowing the run down substantially.
       Raising {.field mcmc.max_treedepth} (e.g. to 12) lets those
       trajectories finish, at up to twice the cost per iteration for each
       level added; reparameterizing the model is the durable fix. Raising
       {.field mcmc.target_accept} does not help here -- it shrinks the step
       size, so trajectories need more steps and hit the cap more often."
    } else {
      "A share this small is common and rarely worth acting on."
    }
    notes <- c(notes, sprintf(
      "%d transition%s (%.2f%%) hit the maximum treedepth (%d): this only
       means those trajectories were cut short for speed and slows effective
       mixing -- unlike a divergence it does not bias the posterior, so it
       is not part of the convergence gate above. %s",
      gate$treedepth_hits, if (isTRUE(gate$treedepth_hits == 1)) "" else "s",
      100 * gate$treedepth_fraction, gate$max_treedepth, action
    ))
  }
  notes
}

#' Write a convergence gate as JSON (artifact parity with Python)
#' @param gate A list from [convergence_gate()] or a cut manifest.
#' @param path Output file path.
#' @export
write_convergence_json <- function(gate, path) {
  # `na = "null"`: an ungated variational stage 1 carries NA fields, and JSON
  # null reads as absent rather than as the string "NA".
  jsonlite::write_json(
    gate, path,
    auto_unbox = TRUE, pretty = TRUE, digits = 8, na = "null"
  )
  invisible(path)
}
