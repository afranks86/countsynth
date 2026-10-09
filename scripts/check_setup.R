# Confirm this machine can fit a countsynth model.
#
#   Rscript scripts/check_setup.R
#
# or, without cloning anything (it needs no local copy and no installed
# package -- the checks below name what is missing and how to get it):
#
#   source("https://raw.githubusercontent.com/afranks86/countsynth/main/scripts/check_setup.R")
#
# It modifies nothing. It checks the R version, the C++ toolchain, the Stan
# installation and the countsynth package, then fits a deliberately tiny
# model end to end. If a step fails, copy the whole output and send it on.

cat("countsynth setup check\n")
cat(strrep("=", 40), "\n\n")

results <- list()

check <- function(name, expr, needs_all_prior = FALSE) {
  cat(sprintf("[ ] %s ... ", name))
  # Running the smoke test after a prerequisite failed only produces a second,
  # more confusing error, so skip rather than cascade.
  if (needs_all_prior && !all(unlist(results))) {
    cat("SKIPPED (fix the failure(s) above first)\n")
    return(invisible(NA))
  }
  ok <- tryCatch({
    expr
    TRUE
  }, error = function(e) {
    cat("\n    ERROR:", conditionMessage(e), "\n")
    FALSE
  })
  cat(if (ok) "OK\n" else "FAILED\n")
  results[[name]] <<- ok
  invisible(ok)
}

# 1. R version ----------------------------------------------------------
check("R version >= 4.1", {
  if (getRversion() < "4.1") {
    stop("R ", getRversion(), " found; please install R 4.1 or newer from https://cran.r-project.org")
  }
})

# 2. cmdstanr installed ---------------------------------------------------
check("cmdstanr package installed", {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "cmdstanr is not installed. Run:\n",
      "    install.packages(\"cmdstanr\", repos = c(\"https://stan-dev.r-universe.dev\", getOption(\"repos\")))"
    )
  }
})

# 3. C++ compiler toolchain (checked the way cmdstanr itself needs it) -----
check("C++ compiler toolchain", {
  ok <- tryCatch({
    cmdstanr::check_cmdstan_toolchain(fix = FALSE, quiet = TRUE)
    TRUE
  }, error = function(e) FALSE)
  if (!isTRUE(ok)) {
    stop(
      "No working C++ toolchain found for CmdStan.\n",
      "    Windows: install Rtools matching your R version from https://cran.r-project.org/bin/windows/Rtools/\n",
      "    Mac: run `xcode-select --install` in Terminal\n",
      "    Linux: install g++ via your package manager (e.g. `sudo apt install g++`)\n",
      "    Then try: cmdstanr::check_cmdstan_toolchain(fix = TRUE)"
    )
  }
})

# 4. CmdStan toolchain installed and discoverable --------------------------
check("CmdStan toolchain installed", {
  path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)
  if (is.null(path) || !nzchar(path)) {
    stop("CmdStan not found. Run:\n    cmdstanr::install_cmdstan()   # one-time, ~10 min")
  }
  cat("\n    found at:", path, "\n    ")
})

# 5. countsynth installed ---------------------------------------------------
check("countsynth package installed", {
  if (!requireNamespace("countsynth", quietly = TRUE)) {
    stop(
      "countsynth is not installed. Run:\n",
      "    install.packages(\"remotes\")\n",
      "    remotes::install_github(\"afranks86/countsynth\", dependencies = TRUE)"
    )
  }
})

# 6. End-to-end smoke test ---------------------------------------------------
# Six units, eight yearly periods, rank 2, 50 + 50 iterations: a few seconds
# of sampling, which is all it takes to prove that data prep, Stan and the
# draws frame work on this machine. It is far too little for inference, and
# so reports convergence warnings and a failing gate -- expected here, and
# suppressed, because on a setup check that reads as a broken install.
#
# The model still has to compile, which is the slow part on a first run and
# the part most likely to expose a broken toolchain.
check("End-to-end smoke test (compile + sample)", {
  cat("\n    compiling and sampling ... ")
  elapsed <- system.time({
    suppressWarnings(suppressMessages(invisible(capture.output({
      library(countsynth)
      states <- unique(utils::read.csv(countsynth_example_data())$state)
      keep <- c("Texas", "California", "New York", "Florida", "Ohio",
                "Pennsylvania")
      cfg <- countsynth_example_config(
        model = countsynth_model_opts(
          outcome_distribution = "NB",
          types = list(total = countsynth_type(
            groups = "total", ranks_to_test = 2,
            exclude_units = setdiff(states, keep)
          ))
        ),
        mcmc = countsynth_mcmc_opts(
          iter_warmup = 50, iter_sampling = 50, seed = 1, progress = FALSE
        ),
        time_aggregation = countsynth_time_aggregation(
          enabled = TRUE, period = "yearly"
        )
      )
      dat <- countsynth_data(cfg)
      fit <- countsynth_fit(dat, config = cfg)
      draws <- countsynth_draws(fit)
    }))))
  })
  if (!nrow(draws) > 0) {
    stop("The model ran but produced no posterior draws.")
  }
  cat("\n    ", format(nrow(draws), big.mark = ","), "draw rows in",
      round(elapsed[["elapsed"]]), "seconds\n    ")
}, needs_all_prior = TRUE)

# Summary -------------------------------------------------------------------
cat("\n", strrep("=", 40), "\n", sep = "")
n_fail <- sum(!unlist(results))
if (n_fail == 0) {
  cat("All checks passed. This machine can fit countsynth models.\n")
} else {
  cat(sprintf(
    paste0(
      "%d check(s) failed. Fix the item(s) marked FAILED above, in order,\n",
      "then re-run this script. If you're stuck, email the organizer the\n",
      "full output of this script.\n"
    ),
    n_fail
  ))
}
