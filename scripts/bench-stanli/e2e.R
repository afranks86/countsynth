# End to end through the package API, at the config's own MCMC settings:
# the fit (joint, or both cut stages), the draws frame and the convergence
# gate -- every stage the backend touches. The report is left out; it reads
# the draws frame, which is the same object on either backend.
#   Rscript e2e.R <config> <type> <rank> <cmdstan|stanli>

pkgload::load_all("/Users/afranks/Dropbox/course/countsynth", quiet = TRUE)
OUT_DIR <- "/Users/afranks/Dropbox/course/countsynth/scripts/bench-stanli/out"

args <- commandArgs(trailingOnly = TRUE)
cfg_file <- args[[1]]
type <- args[[2]]
rank <- as.integer(args[[3]])
backend <- args[[4]]

setwd("/Users/afranks/Dropbox/course/test_countsynth")
config <- read_countsynth_config(cfg_file)
config$mcmc$backend <- backend
config$mcmc$progress <- FALSE
data <- countsynth_data(config, type = type)
is_cut <- identical(config$model$inference_mode, "cut")

stage <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, s = proc.time()[["elapsed"]] - t0)
}
if (is_cut) {
  f <- stage(suppressWarnings(
    countsynth_cut_fit(data, rank = rank, config = config)
  ))
  d <- list(value = f$value$draws, s = 0)
  g <- list(value = f$value$manifest, s = 0)
  stage1 <- f$value$stage1$fit
} else {
  f <- stage(countsynth_fit(data, rank = rank, config = config))
  d <- stage(countsynth_draws(f$value))
  g <- stage(convergence_gate(f$value))
  stage1 <- f$value$fit
}
total <- f$s + d$s + g$s
res <- list(
  config = cfg_file, backend = backend, cut = is_cut,
  fit_s = f$s, draws_s = d$s, gate_s = g$s, total_s = total,
  converged = g$value$converged,
  stage1_time = stage1$time()$chains
)
saveRDS(res, file.path(OUT_DIR, sprintf(
  "e2e_%s_%s.rds", tools::file_path_sans_ext(basename(cfg_file)), backend
)))
cat(sprintf(
  "E2E %-28s %-8s fit %7.1fs  draws %5.1fs  gate %5.1fs  total %7.1fs  converged %s\n",
  basename(cfg_file), backend, f$s, d$s, g$s, total, g$value$converged
))
