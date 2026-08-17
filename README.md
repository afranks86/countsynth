# bpnmf

Estimate causal treatment effects on panel count data (births, low birth
weight, preterm outcomes, …) with a Bayesian hierarchical model built on a
low-rank nonnegative factorization of the log-rate surface. MCMC runs in
**Stan** via [cmdstanr](https://mc-stan.org/cmdstanr/); figures are
**ggplot2**. This is a full R port of the Python package
`bayesian_panel_nmf`, aimed at demographers, social scientists, and political
scientists who work in R.

## What it does

Given a panel of counts indexed by *group × unit × time* (e.g. births by
subgroup, state, and bimonth) with a 0/1 treatment indicator:

- models the untreated log-rate surface as a rank-R nonnegative factorization
  plus unit and time fixed effects (Poisson or Negative Binomial likelihood),
- adds a hierarchical treatment-effect block on exposed cells (staggered
  adoption supported; non-centered parameterization for the weakly-identified
  scales),
- optionally integrates the likelihood over **censored small counts** (values
  1–9 suppressed in the source data),
- supports **joint** inference and two-stage **cut** (modular) inference,
  where the baseline is fit on control data only and treatment effects are
  estimated conditional on frozen baseline draws,
- reports counterfactual expected counts, excess and percent-change effects
  with credible intervals, four posterior-predictive checks, and an R-hat /
  ESS / divergence convergence gate.

## Installation

```r
install.packages("cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()   # one-time, ~10 min

# install bpnmf from source
# install.packages("remotes")
remotes::install_local("path/to/bpnmf-r")
```

## Quick start

```r
library(bpnmf)

cfg  <- bpnmf_example_config()          # bundled US fertility data
dat  <- bpnmf_data(cfg)                 # CSV -> standardized panel + arrays
fit  <- bpnmf_fit(dat, config = cfg)    # NUTS via cmdstanr

summary(fit)                            # convergence gate
parameter_diagnostics(fit)              # per-parameter R-hat / ESS table

draws <- bpnmf_draws(fit)               # tidy long draws frame
bpnmf_summary_table(draws)              # headline observed vs expected
plot(fit, which = "unit_fit", unit = "Texas")
bpnmf_interval_plot(draws, estimand = "ratio")

# or run everything (fits + artifacts + figures) from the config:
bpnmf_run(cfg)
```

Existing Python YAML configs load directly:

```r
cfg <- read_bpnmf_config("configs/fertility_cut_config.yaml")
```

Python's `num_warmup` / `num_samples` / `thinning` / `target_accept` /
`random_seed` translate to cmdstanr's `iter_warmup` / `iter_sampling` /
`thin` / `adapt_delta` / `seed`; everything else keeps its name. Unknown keys
are rejected, like the Python loader.

## Cut (modular) inference

```r
cfg  <- bpnmf_example_config(
  model = bpnmf_model_opts(
    types = list(total = bpnmf_type("total", 3)),
    inference_mode = "cut"
  ),
  cut = bpnmf_cut_opts(num_stage1_draws = 25, stage2_draws_per_component = 100)
)
cfit <- bpnmf_cut_fit(bpnmf_data(cfg), config = cfg)
summary(cfit)               # stage-1 gate + per-component table
draws <- bpnmf_draws(cfit)  # pooled, with cut_component / stage1_* provenance
```

## Parity with the Python implementation

- The Stan models are checked against the NumPyro models by **exact
  log-density comparison** at shared parameter values (fixtures exported from
  Python; `tests/testthat/test-logdensity-parity.R`).
- Cell ordering, chain-stratified cut selection, output thinning, and
  convergence-gate bands are fixture-tested against the Python
  implementation.
- MCMC draws are **not** bitwise-identical across implementations (different
  samplers and RNGs); equivalence is distributional.

## License

MIT.
