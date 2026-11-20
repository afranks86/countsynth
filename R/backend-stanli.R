# The stanli backend: the bundled models run through stanli
# (https://github.com/seantalts/stanli), a Stan interpreter that evaluates
# the same Stan program over precompiled Stan Math kernels with its own
# reverse pass, and samples it with Stan's own NUTS. Nothing is compiled, no
# CSV is written, and the draws never leave memory.
#
# Everything downstream of a fit talks to it through cmdstanr's method names
# (`$draws()`, `$metadata()`, `$sampler_diagnostics()`), which stanli's
# `as_cstanfit()` provides, so the backend is only this model object: a
# `$sample()` taking the cmdstanr arguments countsynth passes.

stanli_stan_model <- function(name) {
  check_stanli()
  file <- stan_file_path(name)
  include <- system.file("stan", "include", package = "countsynth")
  structure(
    list(
      sample = function(...) stanli_sample(file, include, ...),
      # ADVI stays on CmdStan: stanli's variational fit is not a cmdstanr
      # lookalike, and cut stage 1 reads it as one.
      variational = function(...) {
        cli::cli_inform(
          c(i = "ADVI runs on CmdStan; the {.val stanli} backend samples only."),
          .frequency = "regularly", .frequency_id = "countsynth_stanli_advi"
        )
        countsynth_stan_model(name, backend = "cmdstan")$variational(...)
      },
      stan_file = function() file
    ),
    class = "countsynth_stanli_model"
  )
}

check_stanli <- function() {
  rlang::check_installed(
    "stanli",
    reason = "for `mcmc$backend = \"stanli\"`.",
    action = function(...) {
      utils::install.packages(
        "stanli", repos = c("https://seantalts.r-universe.dev",
                            getOption("repos"))
      )
    }
  )
  if (!stanli::stanli_available()) {
    cli::cli_abort(c(
      "The stanli runtime is not installed.",
      i = "Install it with {.code stanli::stanli_install()} (one-time, ~12 MB)."
    ))
  }
  invisible(TRUE)
}

# cmdstanr's `$sample()` for the arguments countsynth passes. Data go through
# cmdstanr's own JSON writer, so both backends read the same bytes.
stanli_sample <- function(file, include, data, chains = 4,
                          parallel_chains = chains, iter_warmup = 1000,
                          iter_sampling = 1000, thin = 1, adapt_delta = 0.8,
                          max_treedepth = 10, seed = NULL, refresh = NULL,
                          init = NULL, show_messages = TRUE,
                          show_exceptions = TRUE, diagnostics = NULL, ...) {
  dots <- list(...)
  if (length(dots) > 0) {
    cli::cli_abort(c(
      "The {.val stanli} backend does not take {.arg {names(dots)}}.",
      i = "Use {.code mcmc$backend = \"cmdstan\"} for CmdStan-only options."
    ))
  }
  seed <- seed %||% sample.int(.Machine$integer.max, 1)
  json <- data
  if (!is.character(data)) {
    json <- tempfile(fileext = ".json")
    on.exit(unlink(json), add = TRUE)
    cmdstanr::write_stan_json(data, json)
  }
  # Building is cheap (well under a second on the bundled panels), so every
  # fit -- each cut stage-2 component included -- builds its own. If part of
  # a program ever has no compiled path, stanli warns, naming it, and runs
  # that part through its much slower interpreter.
  model <- stanli::stanli_model(
    file = file, data = json, include_paths = include, seed = seed
  )
  fit <- stanli::sample_model(
    model,
    chains = chains, parallel_chains = parallel_chains, seed = seed,
    warmup = iter_warmup, samples = iter_sampling, thin = thin,
    delta = adapt_delta, max_depth = max_treedepth,
    init = stanli_init(model, init, chains, seed),
    init_radius = if (is.numeric(init) && length(init) == 1) init else 2,
    refresh = if (isFALSE(show_messages)) 0 else refresh %||% 100
  )
  fit <- cmdstanr_metadata(stanli::as_cstanfit(fit))
  if ("ebfmi" %in% diagnostics) warn_low_ebfmi(fit)
  fit
}

# stanli's metadata lists zero-size variables (the unused te-regression
# block, say), which its own `$draws()` then refuses as unknown; cmdstanr
# reads the list off the CSV header, where they never appear. The gate and
# the trace plots ask for every listed variable, so they are dropped here.
cmdstanr_metadata <- function(fit) {
  metadata <- fit$metadata
  fit$metadata <- function() {
    md <- metadata()
    size <- vapply(md$stan_variable_sizes, prod, numeric(1))
    md$stan_variables <- intersect(md$stan_variables, names(size)[size > 0])
    md$stan_variable_sizes <- md$stan_variable_sizes[md$stan_variables]
    md
  }
  fit
}

# cmdstanr's `init`: NULL or a radius (random starts, handled by the
# sampler), or constrained values -- a list per chain, or a function of
# chain_id -- which stanli takes on the unconstrained scale.
stanli_init <- function(model, init, chains, seed) {
  if (is.null(init) || (is.numeric(init) && length(init) == 1)) {
    return(NULL)
  }
  inits <- if (is.function(init)) {
    lapply(seq_len(chains), function(c) {
      if ("chain_id" %in% names(formals(init))) init(chain_id = c) else init()
    })
  } else {
    init
  }
  if (length(inits) != chains) {
    cli::cli_abort("{.arg init} needs one list per chain ({chains}).")
  }
  # Parameters an init leaves out start at random, as in CmdStan.
  do.call(rbind, lapply(seq_len(chains), function(c) {
    stanli::unconstrain(model, inits[[c]], partial = TRUE,
                        ignore_unknown = TRUE, seed = seed + c)
  }))
}

# cmdstanr's E-BFMI check (`diagnostics = "ebfmi"`), which countsynth keeps
# because it has no gate of its own: same statistic, same message.
warn_low_ebfmi <- function(fit, threshold = 0.3) {
  energy <- posterior::extract_variable_matrix(
    fit$sampler_diagnostics(), "energy__"
  )
  if (nrow(energy) < 3 || anyNA(energy)) {
    return(invisible(NULL))
  }
  ebfmi <- apply(energy, 2, function(x) {
    (sum(diff(x)^2) / length(x)) / stats::var(x)
  })
  low <- sum(ebfmi < threshold)
  if (low > 0) {
    message(
      "Warning: ", low, " of ", length(ebfmi), " chains had an E-BFMI less ",
      "than ", threshold, ".\nSee https://mc-stan.org/misc/warnings for ",
      "details.\n"
    )
  }
  invisible(unname(ebfmi))
}
