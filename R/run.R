# Full pipeline: config -> data -> fit (joint or cut) -> draws artifact ->
# convergence JSON -> reports. Mirrors the Python `bpnmf run` command's
# artifact layout:
#   <output_dir>/<type>/{dist}_{outcome}_{type}_{rank}[_cut].csv|.parquet
#   <output_dir>/<type>/{stem}_convergence.json
#   <output_dir>/<type>/df_<type>.csv
#   <output_dir>/<type>/[rank_<rank>/]figs/...

outcome_name <- function(config) {
  if (!is.null(config$outcome)) {
    return(config$outcome)
  }
  p <- config$schema$outcomes_from_prefixes
  if (!is.null(p)) {
    return(sub("_$", "", p$outcome_prefix))
  }
  "outcome"
}

draws_stem <- function(config, type, rank, cut = FALSE) {
  stem <- sprintf(
    "%s_%s_%s_%d",
    config$model$outcome_distribution, outcome_name(config), type, rank
  )
  if (cut) paste0(stem, "_cut") else stem
}

write_draws_file <- function(draws, stem, dir, format) {
  if (format == "parquet") {
    rlang::check_installed("arrow")
    path <- file.path(dir, paste0(stem, ".parquet"))
    arrow::write_parquet(as.data.frame(draws), path)
  } else {
    path <- file.path(dir, paste0(stem, ".csv"))
    utils::write.csv(as.data.frame(draws), path, row.names = FALSE)
  }
  path
}

#' Run the full bpnmf pipeline from a config
#'
#' Fits every requested (type, rank), writes the tidy draws artifact, the
#' convergence gate JSON, and (when `output$figures` selects any) the full
#' figure/table report -- with the same directory layout and filenames as the
#' Python `bpnmf run` command.
#'
#' @param config A [bpnmf_config()] object (or path to a YAML config).
#' @param types Subset of `config$model$types` names to run (`NULL` = all).
#' @param ranks Rank override (`NULL` = each type's `ranks_to_test`).
#' @return Invisible named list of per-(type, rank) results, each holding the
#'   fit object, draws path, and gate.
#' @export
bpnmf_run <- function(config, types = NULL, ranks = NULL) {
  if (is.character(config)) {
    config <- read_bpnmf_config(config)
  }
  checkmate::assert_class(config, "bpnmf_config")
  if (length(config$model$types) == 0) {
    cli::cli_abort("Config must define at least one entry in {.field model.types}.")
  }
  types <- types %||% names(config$model$types)
  checkmate::assert_subset(types, names(config$model$types))
  is_cut <- identical(config$model$inference_mode, "cut")
  out <- config$output

  results <- list()
  for (type in types) {
    type_spec <- config$model$types[[type]]
    type_dir <- file.path(config$output_dir, type)
    if (out$clean && dir.exists(type_dir)) {
      unlink(type_dir, recursive = TRUE)
    }
    dir.create(type_dir, recursive = TRUE, showWarnings = FALSE)

    data <- bpnmf_data(config, type = type)
    utils::write.csv(
      data$df, file.path(type_dir, sprintf("df_%s.csv", type)),
      row.names = FALSE
    )

    rank_list <- ranks %||% type_spec$ranks_to_test
    for (rank in rank_list) {
      cli::cli_h1("{type} @ rank {rank} ({if (is_cut) 'cut' else 'joint'})")
      run_dir <- if (length(rank_list) > 1) {
        file.path(type_dir, sprintf("rank_%d", rank))
      } else {
        type_dir
      }
      dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
      stem <- draws_stem(config, type, rank, cut = is_cut)

      if (is_cut) {
        fit <- bpnmf_cut_fit(data, rank = rank, config = config)
        draws <- fit$draws
        gate <- fit$manifest
        ppc_draws <- fit$stage1_ppc
        write_draws_file(
          ppc_draws, paste0(stem, "_stage1_ppc"), type_dir, "csv"
        )
        cut_component_table(fit)
      } else {
        fit <- bpnmf_fit(data, rank = rank, config = config)
        draws <- bpnmf_draws(fit)
        gate <- convergence_gate(fit)
        ppc_draws <- NULL
      }

      draws_path <- write_draws_file(draws, stem, type_dir, out$draws_format)
      write_convergence_json(
        gate, file.path(type_dir, paste0(stem, "_convergence.json"))
      )
      if (isFALSE(gate$converged)) {
        reasons <- gate_failure_bullets(
          gate, config$mcmc$convergence, if (is_cut) NULL else fit
        )
        advice <- gate_failure_advice(gate, config$mcmc$convergence)
        cli::cli_warn(c(
          "Convergence gate FAILED for {stem}; artifacts still written.",
          stats::setNames(reasons, rep("*", length(reasons))),
          i = "Gated on {.val {config$mcmc$gate_params}}; widen or narrow with
               {.field mcmc.gate_params}.",
          i = "Per-parameter detail: {.code parameter_diagnostics(fit)}.",
          stats::setNames(advice, rep("i", length(advice)))
        ))
      }
      # Non-gating context (divergences on a pass, treedepth either way):
      # meant to replace cmdstanr's own unglossed sampler warnings as the
      # thing a user actually reads.
      notes <- diagnostic_context_notes(
        if (is_cut) gate$stage1 else gate, config$mcmc$convergence
      )
      if (length(notes) > 0) {
        cli::cli_inform(stats::setNames(notes, rep("i", length(notes))))
      }
      if (out$save_traces) {
        saveRDS(
          posterior::as_draws_rvars(fit$fit$draws()),
          file.path(type_dir, paste0(stem, "_draws.rds"))
        )
      }

      if (length(out$figures) > 0) {
        bpnmf_report(
          draws,
          output_dir = run_dir,
          target_unit = out$target_unit,
          groups = out$report_groups,
          figures = out$figures,
          aggregate_units = out$aggregate_units,
          ppc_draws = ppc_draws,
          ppc_units = out$ppc_units,
          ppc_exclude_units = out$ppc_exclude_units,
          ppc_acf_lags = out$ppc_acf_lags,
          ppc_unit_corr_max_time = out$ppc_unit_corr_max_time,
          fit_gap_per_unit = out$fit_gap_per_unit,
          interval_aggregates = out$interval_aggregates,
          print_tables = out$print_tables,
          print_target_table = out$print_target_table,
          html_tables = out$html_tables,
          fit = fit,
          denominator_may_be_affected = out$denominator_may_be_affected
        )
      }

      results[[sprintf("%s_rank%d", type, rank)]] <- list(
        fit = fit, draws_path = draws_path, gate = gate
      )
    }
  }
  invisible(results)
}
