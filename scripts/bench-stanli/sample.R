# One NUTS run of joint.stan on one engine, timed and summarized.
#   Rscript sample.R <case> <cmdstan|stanli> [warmup] [samples] [seed]
# Same data file, seed, chains, adapt_delta and treedepth on both engines;
# 4 chains in parallel, as countsynth runs them. Writes out/<case>_<engine>.rds.

source("/Users/afranks/Dropbox/course/countsynth/scripts/bench-stanli/common.R")

args <- commandArgs(trailingOnly = TRUE)
name <- args[[1]]
engine <- args[[2]]
warmup <- as.integer(if (length(args) > 2) args[[3]] else 300)
samples <- as.integer(if (length(args) > 3) args[[4]] else 300)
seed_arg <- if (length(args) > 4) as.integer(args[[5]]) else NULL

json <- case_json(name, TRUE)
config <- withr::with_dir(TEST_DIR, read_countsynth_config(CASES[[name]]$config))
mcmc <- config$mcmc
seed <- seed_arg %||% mcmc$seed
delta <- mcmc$adapt_delta
depth <- mcmc$max_treedepth %||% 10L
chains <- 4L

t_start <- proc.time()[["elapsed"]]
if (engine == "cmdstan") {
  model <- countsynth_stan_model("joint")
  t_ready <- proc.time()[["elapsed"]]
  fit <- model$sample(
    data = json, chains = chains, parallel_chains = chains,
    iter_warmup = warmup, iter_sampling = samples, adapt_delta = delta,
    max_treedepth = depth, seed = seed, refresh = 0,
    show_messages = FALSE, show_exceptions = FALSE, diagnostics = NULL
  )
  t_sampled <- proc.time()[["elapsed"]]
  draws <- fit$draws(format = "draws_array")
  t_draws <- proc.time()[["elapsed"]]
  tm <- fit$time()$chains
  chain_warmup <- tm$warmup
  chain_sampling <- tm$sampling
  sdiag <- fit$sampler_diagnostics(format = "draws_array")
} else {
  library(stanli)
  model <- stanli_model(file = stan_file(), data = json,
                        include_paths = INCLUDE_DIR, seed = seed)
  t_ready <- proc.time()[["elapsed"]]
  fit <- sample_model(
    model, chains = chains, parallel_chains = chains, seed = seed,
    warmup = warmup, samples = samples, delta = delta, max_depth = depth,
    refresh = 0
  )
  t_sampled <- proc.time()[["elapsed"]]
  draws <- as_draws_array(fit)
  t_draws <- proc.time()[["elapsed"]]
  chain_warmup <- fit$report$warmup_seconds
  chain_sampling <- fit$report$sampling_seconds
  sdiag <- posterior::as_draws_array(fit$sampler)
}

leap <- posterior::extract_variable_matrix(sdiag, "n_leapfrog__")
div <- posterior::extract_variable_matrix(sdiag, "divergent__")
depth_hit <- posterior::extract_variable_matrix(sdiag, "treedepth__") >= depth

key <- posterior::subset_draws(
  draws, variable = intersect(c("mu_ctrl", "te"),
                              sub("\\[.*", "", posterior::variables(draws))),
  regex = FALSE
)
summ <- posterior::summarise_draws(
  key, "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean"
)

res <- list(
  case = name, engine = engine, warmup = warmup, samples = samples,
  setup_s = t_ready - t_start,
  sample_s = t_sampled - t_ready,
  draws_s = t_draws - t_sampled,
  chain_warmup = chain_warmup, chain_sampling = chain_sampling,
  leap_per_chain = colSums(leap),
  ms_per_grad = 1000 * chain_sampling / colSums(leap),
  divergent = sum(div), treedepth_hits = sum(depth_hit),
  summary = summ
)
res$seed <- seed
saveRDS(res, file.path(OUT_DIR, if (is.null(seed_arg)) {
  sprintf("%s_%s.rds", name, engine)
} else {
  sprintf("%s_%s_seed%d.rds", name, engine, seed)
}))
cat(sprintf(
  "%-18s %-8s seed %d setup %6.1fs  sample %7.1fs  draws %5.1fs  ms/grad %.3f  leapfrog/iter %.0f  div %d  min ESS bulk %.0f  max rhat %.3f\n",
  name, engine, seed, res$setup_s, res$sample_s, res$draws_s,
  mean(res$ms_per_grad), mean(leap), res$divergent,
  min(summ$ess_bulk, na.rm = TRUE), max(summ$rhat, na.rm = TRUE)
))
