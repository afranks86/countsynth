# Two-stage cut (modular) inference. Port of cut.py:
#   p_cut(phi, theta | Z, Y) = p(phi | Z) * p(theta | Y, phi)
# Stage 1 fits the baseline on control cells; a chain-stratified, seeded
# subset of its draws becomes the cut components; stage 2 runs one small MCMC
# per component with that draw's mu_ctrl frozen as data, so exposed outcomes
# can never feed back into the baseline.
#
# Stage 1 is the expensive half (the full factorization over every control
# cell), so `cut$stage1_method = "variational"` swaps its NUTS run for ADVI.
# Everything downstream is unchanged: ADVI still yields a set of draws of
# mu_ctrl, and the same seeded selection promotes some of them to components.
# What changes is what those components mean -- see the warning in
# `warn_variational_stage1()`.

resolve_cut_settings <- function(config) {
  cut <- config$cut %||% countsynth_cut_opts()
  seed <- config$mcmc$seed
  selection_seed <- cut$selection_seed %||% (seed + 2L)
  stage2_seed <- cut$stage2_seed %||% (seed + 3L)
  reserved <- c(inference = seed, ppc = seed + 1L)
  for (nm in names(reserved)) {
    if (selection_seed == reserved[[nm]]) {
      cli::cli_warn(
        "cut selection_seed ({selection_seed}) collides with the {nm} seed;
         selection should be independent of the draws being selected."
      )
    }
    if (stage2_seed == reserved[[nm]]) {
      cli::cli_warn(
        "cut stage2_seed ({stage2_seed}) collides with the {nm} seed."
      )
    }
  }

  # Shallow overlay of stage2_mcmc over the top-level MCMC options; accepts
  # both Python and cmdstanr key names.
  overlay <- cut$stage2_mcmc %||% list()
  key_map <- c(
    num_warmup = "iter_warmup", num_samples = "iter_sampling",
    thinning = "thin", target_accept = "adapt_delta",
    progress_bar = "progress"
  )
  for (py in names(key_map)) {
    if (py %in% names(overlay)) {
      overlay[[key_map[[py]]]] <- overlay[[py]]
      overlay[[py]] <- NULL
    }
  }
  stage2_mcmc <- config$mcmc
  for (nm in names(overlay)) {
    stage2_mcmc[[nm]] <- overlay[[nm]]
  }

  list(
    num_stage1_draws = cut$num_stage1_draws,
    stage2_draws_per_component = cut$stage2_draws_per_component,
    selection_seed = as.integer(selection_seed),
    stage2_seed = as.integer(stage2_seed),
    stage2_mcmc = stage2_mcmc,
    stage1_method = cut$stage1_method %||% "sample",
    stage1_variational = cut$stage1_variational
  )
}

#' Warn that an ADVI stage 1 changes what the cut posterior means
#'
#' Mean-field ADVI fits a factorized Gaussian in the unconstrained space. It
#' therefore understates marginal variance and drops posterior correlation
#' between baseline parameters. In cut mode that error propagates in a
#' specific direction: the stage-1 components are drawn from too narrow a
#' spread, so the pooled `te` posterior is too narrow too, and its intervals
#' under-cover. Point estimates are usually close; uncertainty is not.
#' @keywords internal
warn_variational_stage1 <- function() {
  cli::cli_warn(c(
    "Cut stage 1 is an ADVI approximation, not a posterior sample.",
    "!" = "Mean-field ADVI understates posterior variance and ignores
           posterior correlation, so the stage-1 components span too narrow a
           range and the pooled treatment-effect intervals will be too tight.",
    i = "Stage-1 R-hat / ESS / divergence gating does not apply; the manifest
         records {.code stage1$converged = NA}, not a pass.",
    i = "Use this to iterate, then re-run with
         {.code cut$stage1_method = \"sample\"} for results you report."
  ))
}

#' Chain-stratified stage-1 draw selection
#'
#' Port of `cut.select_stage1_draws`: quotas split evenly across chains (first
#' chains take the remainder), within-chain picks are seeded, without
#' replacement, and sorted; components are numbered 1..M in ascending
#' (chain, iteration) order for reproducibility.
#'
#' @param n_chains,per_chain Stage-1 chain count and retained draws per chain.
#' @param m Number of draws to select.
#' @param selection_seed RNG seed for the selection stream.
#' @return A data frame with `component`, `chain`, `iteration`, `draw` (the
#'   chain-major global draw index).
#' @keywords internal
select_stage1_draws <- function(n_chains, per_chain, m, selection_seed) {
  total <- n_chains * per_chain
  if (m > total) {
    cli::cli_abort(
      "cut.num_stage1_draws ({m}) exceeds retained stage-1 draws ({total})."
    )
  }
  quotas <- chain_quotas(m, n_chains)
  picks <- withr::with_seed(selection_seed, {
    lapply(seq_len(n_chains), function(c) {
      if (quotas[c] == 0L) {
        integer()
      } else {
        sort(sample.int(per_chain, quotas[c], replace = FALSE))
      }
    })
  })
  refs <- do.call(rbind, lapply(seq_len(n_chains), function(c) {
    if (length(picks[[c]]) == 0) {
      return(NULL)
    }
    data.frame(chain = c, iteration = picks[[c]])
  }))
  refs <- refs[order(refs$chain, refs$iteration), , drop = FALSE]
  refs$component <- seq_len(nrow(refs))
  refs$draw <- as.integer((refs$chain - 1L) * per_chain + refs$iteration)
  rownames(refs) <- NULL
  refs[c("component", "chain", "iteration", "draw")]
}

#' Derived seed for a component's untreated predictive draws
#'
#' The stage-2 sampler seeds occupy `stage2_seed + component`, so the
#' predictive stream is offset by a large multiplier to keep it disjoint.
#' The arithmetic is done in double precision and wrapped into the valid
#' `set.seed()` range: `stage2_seed * 1000` overflows 32-bit integers for any
#' base seed above ~2.1 million (the default `8675309 + 3` among them).
#'
#' @param stage2_seed Base stage-2 seed.
#' @param component Component index (1-based).
#' @return A single integer seed.
#' @keywords internal
predictive_seed <- function(stage2_seed, component) {
  as.integer((as.numeric(stage2_seed) * 1000 + component) %% 2147483647)
}

run_stage2_component <- function(model, sd2, stage2_mcmc, seed_i, quiet = TRUE) {
  ch <- resolve_chains(stage2_mcmc)
  fit <- model$sample(
    data = sd2,
    chains = ch$chains,
    parallel_chains = ch$parallel_chains,
    iter_warmup = stage2_mcmc$iter_warmup,
    iter_sampling = stage2_mcmc$iter_sampling,
    thin = stage2_mcmc$thin,
    adapt_delta = stage2_mcmc$adapt_delta,
    max_treedepth = stage2_mcmc$max_treedepth %||% 10L,
    seed = seed_i,
    refresh = 0,
    show_messages = FALSE,
    show_exceptions = FALSE
  )
  fit
}

#' Fit the countsynth model with two-stage cut inference
#'
#' @param data A `countsynth_data` object.
#' @param rank Factorization rank (defaults to the type's first
#'   `ranks_to_test`).
#' @param config The [countsynth_config()] object; `config$cut` supplies the cut
#'   settings (defaults from [countsynth_cut_opts()] otherwise). Set
#'   `config$cut$stage1_method = "variational"` to fit stage 1 with ADVI
#'   instead of NUTS -- much faster, but approximate and ungated; see
#'   [countsynth_cut_opts()].
#' @param parallel `"none"` (sequential stage-2 loop, default) or `"future"`
#'   (requires \pkg{furrr}; stage-2 components run on the active future plan
#'   with one parallel chain each).
#' @return A `countsynth_cut_fit` object with the pooled tidy `draws` frame
#'   (provenance columns `cut_component`, `stage1_draw`, `stage1_chain`,
#'   `stage1_iteration`), the stage-1 `countsynth_fit`, a `stage1_ppc` draws frame
#'   for the PPC suite, per-component records, and the convergence `manifest`.
#' @export
countsynth_cut_fit <- function(data, rank = NULL, config,
                          parallel = c("none", "future")) {
  checkmate::assert_class(data, "countsynth_data")
  checkmate::assert_class(config, "countsynth_config")
  parallel <- match.arg(parallel)
  if (!identical(config$model$inference_mode, "cut")) {
    cli::cli_warn(
      "Running cut inference although {.field model.inference_mode} is not 'cut'."
    )
  }
  validate_cut_data(data)
  settings <- resolve_cut_settings(config)
  type_spec <- config$model$types[[data$type]]
  rank <- rank %||% type_spec$ranks_to_test[[1]]

  is_nb <- config$model$outcome_distribution == "NB"
  D <- length(data$units)

  # Stage 1: baseline model on control cells, with the counterfactual
  # predictive emitted for the stage-1 PPC frame.
  stage1_method <- settings$stage1_method
  if (stage1_method == "variational") warn_variational_stage1()
  cli::cli_alert_info(
    "Cut stage 1: fitting baseline model (rank {rank}, {fit_method_label(stage1_method)})"
  )
  stage1 <- countsynth_fit(
    data,
    rank = rank, config = config, model_treated = FALSE, gen_ypred = TRUE,
    method = stage1_method, variational = settings$stage1_variational
  )
  stage1_gate <- convergence_gate(stage1)
  ci <- chain_iteration_vectors(stage1$fit)

  mu1 <- extract_var_matrix(stage1$fit, "mu_ctrl")
  disp1 <- if (is_nb && config$model$sample_disp) {
    extract_var_matrix(stage1$fit, "disp")
  } else {
    NULL
  }

  stage1_ppc <- countsynth_draws(stage1)

  refs <- select_stage1_draws(
    ci$n_chains, ci$per_chain, settings$num_stage1_draws,
    settings$selection_seed
  )

  model2 <- countsynth_stan_model("cut_stage2")
  exposed_cell <- as.integer(stage1$stan_data$exp_cell)
  cell_unit <- as.integer(stage1$stan_data$cell_unit)

  # Optional treatment-effect regression: one design shared by every stage-2
  # component (it depends only on the data, not on the stage-1 draw).
  te_design <- if (!is.null(config$model$treatment_effects)) {
    build_te_design(config$model$treatment_effects, data)
  } else {
    NULL
  }

  phi_for_ref <- function(draw) {
    if (!is_nb) {
      return(NULL)
    }
    if (is.null(disp1)) rep(1 / config$model$nb_disp, D) else 1 / disp1[draw, ]
  }

  run_component <- function(i) {
    ref <- refs[i, ]
    phi_unit <- phi_for_ref(ref$draw)
    sd2 <- stan_data_stage2(
      data,
      mu_ctrl_flat = mu1[ref$draw, ],
      phi_unit = phi_unit,
      outcome_distribution = config$model$outcome_distribution,
      adjust_for_missingness = config$model$adjust_for_missingness,
      te_design = te_design
    )
    fit2 <- run_stage2_component(
      model2, sd2, settings$stage2_mcmc,
      seed_i = settings$stage2_seed + ref$component
    )
    diag <- convergence_gate(
      fit2,
      gate_params = config$mcmc$gate_params,
      thresholds = config$mcmc$convergence
    )
    ci2 <- chain_iteration_vectors(fit2)
    te2 <- extract_var_matrix(fit2, "te")

    # Evenly strided per-chain output thinning (deterministic, no RNG).
    quota <- settings$stage2_draws_per_component %||% ci2$per_chain
    keep <- lapply(seq_len(ci2$n_chains), function(c) {
      (c - 1L) * ci2$per_chain + strided_indices(ci2$per_chain, quota)
    })
    keep_idx <- unlist(keep)
    te_out <- te2[keep_idx, , drop = FALSE]
    chain_ids <- rep(seq_len(ci2$n_chains), each = quota)
    iter_ids <- rep(seq_len(quota), times = ci2$n_chains)
    n_out <- nrow(te_out)

    # Regression coefficients ride along on the same thinned draw indices, so
    # a coefficient row and its te row come from the same stage-2 draw.
    te_coefs <- NULL
    if (!is.null(te_design)) {
      m <- te_coef_matrices(fit2, te_design)
      sub <- function(x) if (is.null(x)) NULL else x[keep_idx, , drop = FALSE]
      te_coefs <- te_coef_tidy(
        te_design, sub(m$beta), sub(m$scale), sub(m$z), chain_ids, iter_ids
      )
    }

    # Fresh untreated predictive from the component's frozen baseline.
    withr::with_seed(predictive_seed(settings$stage2_seed, ref$component), {
      mu_grid <- matrix(mu1[ref$draw, ], nrow = n_out, ncol = ncol(mu1), byrow = TRUE)
      phi_row <- if (is_nb) phi_unit[cell_unit] else NULL
      ypred_out <- sample_untreated_predictions(mu_grid, phi_row)
      list(
        ref = ref, diag = diag, te_out = te_out, ypred_out = ypred_out,
        mu_grid = mu_grid, chain_ids = chain_ids, iter_ids = iter_ids,
        n_out = n_out, te_coefs = te_coefs,
        retained_draws = ci2$n_chains * ci2$per_chain
      )
    })
  }

  n_comp <- nrow(refs)
  cli::cli_alert_info(
    "Cut stage 2: {n_comp} conditional fit{?s} ({settings$stage2_mcmc$iter_warmup}+{settings$stage2_mcmc$iter_sampling} iterations each)"
  )
  results <- if (parallel == "future") {
    rlang::check_installed("furrr")
    furrr::future_map(
      seq_len(n_comp), run_component,
      .options = furrr::furrr_options(seed = NULL)
    )
  } else {
    out <- vector("list", n_comp)
    cli::cli_progress_bar("Stage-2 components", total = n_comp)
    for (i in seq_len(n_comp)) {
      out[[i]] <- run_component(i)
      cli::cli_progress_update()
    }
    cli::cli_progress_done()
    out
  }

  # Pool components with provenance; .draw is globally unique, .chain is the
  # real stage-2 chain, .iteration indexes the component's output subsample.
  output_counts <- unique(vapply(results, `[[`, integer(1), "n_out"))
  if (length(output_counts) > 1) {
    cli::cli_abort(
      "Unequal output draw counts across cut components
       ({.val {sort(output_counts)}}); equal counts preserve equal weights."
    )
  }
  draw_offset <- 0L
  component_frames <- vector("list", n_comp)
  component_te_frames <- vector("list", n_comp)
  component_records <- vector("list", n_comp)
  for (i in seq_len(n_comp)) {
    r <- results[[i]]
    df <- draws_frame_core(
      r$mu_grid, r$te_out, r$ypred_out, r$chain_ids, r$iter_ids,
      data, exposed_cell
    )
    df$.draw <- df$.draw + draw_offset
    df$cut_component <- r$ref$component
    df$stage1_draw <- r$ref$draw
    df$stage1_chain <- r$ref$chain
    df$stage1_iteration <- r$ref$iteration
    component_frames[[i]] <- df
    if (!is.null(r$te_coefs)) {
      tc <- r$te_coefs
      tc$.draw <- tc$.draw + draw_offset
      tc$cut_component <- r$ref$component
      tc$stage1_draw <- r$ref$draw
      tc$stage1_chain <- r$ref$chain
      tc$stage1_iteration <- r$ref$iteration
      component_te_frames[[i]] <- tc
    }
    component_records[[i]] <- c(
      list(
        component = r$ref$component, stage1_draw = r$ref$draw,
        stage1_chain = r$ref$chain, stage1_iteration = r$ref$iteration
      ),
      r$diag,
      list(retained_draws = r$retained_draws, output_draws = r$n_out)
    )
    draw_offset <- draw_offset + r$n_out
  }
  draws <- dplyr::bind_rows(component_frames)
  attr(draws, "groups") <- data$groups
  attr(draws, "units") <- data$units
  attr(draws, "times") <- data$times
  attr(draws, "has_denominator") <- "denominator" %in% names(data$df)
  class(draws) <- c("countsynth_draws", class(draws))

  te_draws <- NULL
  if (!is.null(te_design)) {
    te_draws <- dplyr::bind_rows(component_te_frames)
    class(te_draws) <- c("countsynth_te_draws", class(te_draws))
  }

  all_converged <- all(vapply(component_records, function(r) isTRUE(r$converged), logical(1)))
  # An ungated (variational) stage 1 must not silently count as a pass, but
  # it must not fail the run either -- the run-level flag then reports only
  # what was actually gated, and `stage1_gated` says so.
  stage1_gated <- !is.na(stage1_gate$converged)
  manifest <- list(
    inference_mode = "cut",
    stage1_method = stage1_method,
    stage1_gated = stage1_gated,
    converged = all_converged && (!stage1_gated || isTRUE(stage1_gate$converged)),
    stage1 = stage1_gate,
    stage2 = list(
      all_converged = all_converged,
      failed_fits = sum(!vapply(component_records, function(r) isTRUE(r$converged), logical(1))),
      fits = component_records
    )
  )

  new_countsynth_class(
    list(
      draws = draws,
      te_draws = te_draws,
      te_design = te_design,
      stage1 = stage1,
      stage1_ppc = stage1_ppc,
      component_records = component_records,
      manifest = manifest,
      settings = settings,
      config = config,
      data = data,
      rank = as.integer(rank),
      type = data$type,
      inference_mode = "cut",
      stage1_method = stage1_method,
      fit = stage1$fit
    ),
    "countsynth_cut_fit"
  )
}

#' Terminal table of per-component cut diagnostics
#' @param x A `countsynth_cut_fit` object.
#' @export
cut_component_table <- function(x) {
  checkmate::assert_class(x, "countsynth_cut_fit")
  recs <- x$component_records
  fmt <- "%9s %8s %10s %12s %5s %8s  %s"
  cli::cli_h1("Cut stage-2 components")
  cli::cli_verbatim(sprintf(
    fmt, "component", "s1 chain", "max R-hat", "min bulk ESS", "div",
    "div rate", "status"
  ))
  n_pass <- 0L
  for (r in recs) {
    ok <- isTRUE(r$converged)
    n_pass <- n_pass + ok
    line <- sprintf(
      fmt, r$component, r$stage1_chain, sprintf("%.4f", r$rhat_max),
      sprintf("%.0f", r$ess_bulk_min), r$divergences,
      sprintf("%.2f%%", 100 * (r$divergence_fraction %||% 0)),
      if (ok) "PASS" else "FAIL"
    )
    if (ok) cli::cli_verbatim(line) else cli::cli_verbatim(cli::col_red(line))
  }
  cli::cli_alert_info("{n_pass}/{length(recs)} components pass")
  invisible(x)
}
