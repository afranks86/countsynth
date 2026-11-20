# Does stanli compute the same log density and gradient as CmdStan on
# joint.stan, and what does one gradient cost? CmdStan's side is its
# standalone `log_prob` method (jacobian=1, the density NUTS samples).

source("/Users/afranks/Dropbox/course/countsynth/scripts/bench-stanli/common.R")
library(stanli)

cmdstan_lp_grad <- function(exe, data_json, Q) {
  qfile <- tempfile(fileext = ".json")
  out <- tempfile(fileext = ".csv")
  jsonlite::write_json(list(params_r = Q), qfile, digits = NA)
  status <- system2(exe, c(
    "log_prob", paste0("unconstrained_params=", qfile), "jacobian=1",
    "data", paste0("file=", data_json),
    "output", paste0("file=", out), "sig_figs=18"
  ), stdout = FALSE, stderr = FALSE)
  stopifnot(status == 0)
  m <- as.matrix(utils::read.csv(out, comment.char = "#"))
  list(lp = m[, 1], grad = m[, -1, drop = FALSE])
}

exe <- countsynth_stan_model("joint")$exe_file()

for (name in names(CASES)) {
  for (treated in c(TRUE, FALSE)) {
    json <- case_json(name, treated)
    t0 <- proc.time()[["elapsed"]]
    warns <- character()
    m <- withCallingHandlers(
      stanli_model(file = stan_file(), data = json,
                   include_paths = INCLUDE_DIR),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    t_build <- proc.time()[["elapsed"]] - t0
    n <- m$n_unconstrained

    set.seed(1)
    Q <- matrix(stats::runif(5 * n, -1, 1), nrow = 5)
    cs <- cmdstan_lp_grad(exe, json, Q)
    sl <- lapply(seq_len(nrow(Q)), function(i) log_prob_grad(m, Q[i, ]))
    lp_sl <- vapply(sl, `[[`, numeric(1), "lp")
    g_sl <- do.call(rbind, lapply(sl, `[[`, "grad"))

    # Per-gradient cost, timed at a typical (here: random) point.
    reps <- 300L
    q <- Q[1, ]
    t0 <- proc.time()[["elapsed"]]
    for (i in seq_len(reps)) log_prob_grad(m, q)
    t_grad <- (proc.time()[["elapsed"]] - t0) / reps

    cat(sprintf(
      "\n== %s (%s): n_unc=%d build=%.2fs\n  lp diff (stanli - cmdstan): %s\n  max |grad diff| / max |grad|: %.2e\n  stanli per-gradient: %.3f ms\n",
      name, if (treated) "joint" else "baseline", n, t_build,
      paste(sprintf("%.2e", lp_sl - cs$lp), collapse = " "),
      max(abs(g_sl - cs$grad)) / max(abs(cs$grad)),
      1000 * t_grad
    ))
    if (length(warns)) cat("  WARNINGS:\n", paste("   -", warns, collapse = "\n"), "\n")
  }
}
