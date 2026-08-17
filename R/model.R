# Stan model compilation with a persistent cache. Models live in inst/stan
# and compile lazily on first use; binaries are cached under
# tools::R_user_dir("bpnmf", "cache") so every fit (and all 25-50 cut stage-2
# runs) reuses one compiled executable.

stan_file_path <- function(name) {
  path <- system.file("stan", paste0(name, ".stan"), package = "bpnmf")
  if (!nzchar(path)) {
    cli::cli_abort("Stan model {.val {name}} not found in the installed package.")
  }
  path
}

#' Compile (or fetch from cache) a bundled Stan model
#'
#' @param name `"joint"` (also the cut stage-1 baseline, via the
#'   `model_treated` data flag) or `"cut_stage2"`.
#' @param quiet Suppress compilation output.
#' @return A `cmdstanr::CmdStanModel`.
#' @export
bpnmf_stan_model <- function(name = c("joint", "cut_stage2"), quiet = TRUE) {
  name <- match.arg(name)
  check_cmdstan()
  cache_dir <- tools::R_user_dir("bpnmf", "cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cmdstanr::cmdstan_model(
    stan_file_path(name),
    include_paths = system.file("stan", "include", package = "bpnmf"),
    dir = cache_dir,
    quiet = quiet
  )
}

#' Verify CmdStan is available, with an actionable install hint
#' @keywords internal
check_cmdstan <- function() {
  ok <- tryCatch(
    nzchar(cmdstanr::cmdstan_path()),
    error = function(e) FALSE
  )
  if (!ok) {
    cli::cli_abort(c(
      "CmdStan was not found.",
      i = "Install it with {.code cmdstanr::install_cmdstan()} (one-time, ~10 min).",
      i = "If CmdStan is already installed, point cmdstanr at it with
           {.code cmdstanr::set_cmdstan_path()}."
    ))
  }
  invisible(TRUE)
}

#' Resolve chain settings from MCMC options
#'
#' Port of the Python auto-parallelism logic for CPU sampling: under
#' `auto_parallelism`, run `min(max_chains, cores)` chains in parallel;
#' otherwise use the manual settings (default 4 chains, sequential).
#' @keywords internal
resolve_chains <- function(mcmc) {
  if (mcmc$auto_parallelism) {
    cores <- max(1L, parallel::detectCores())
    chains <- min(mcmc$max_chains, max(1L, cores))
    list(chains = chains, parallel_chains = chains)
  } else {
    chains <- mcmc$chains %||% 4L
    list(chains = chains, parallel_chains = mcmc$parallel_chains %||% 1L)
  }
}
