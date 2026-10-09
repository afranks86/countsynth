# Run this script before the workshop to confirm your machine is ready.
#
#   Rscript scripts/check_workshop_setup.R
#
# or open it in RStudio and click "Source". It does not modify anything;
# it only checks your R version, compiler toolchain, Stan installation, and
# the countsynth package, then runs a smoke test of the full pipeline. The
# smoke test takes several minutes. If any step fails, copy the full output
# and send it to the workshop organizer.

cat("countsynth workshop setup check\n")
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
# The bundled example runs only 200 warmup + 200 sampling iterations, so it
# WILL report convergence warnings and a failing convergence gate. That is
# expected and irrelevant here: this step only asks whether the pipeline runs
# end to end on this machine, not whether the answers are trustworthy. The
# noise is suppressed so it doesn't look like a setup problem.
check("End-to-end smoke test (compile + sample on bundled example)", {
  cat("\n    this takes several minutes, please wait ... ")
  elapsed <- system.time({
    suppressWarnings(suppressMessages(invisible(capture.output({
      library(countsynth)
      cfg <- countsynth_example_config()
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
  cat("All checks passed. You're ready for the workshop.\n")
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
