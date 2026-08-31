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
remotes::install_github("afranks86/bpnmf")
```

Check the toolchain before a long run — `cmdstanr::cmdstan_path()` should
print a path rather than error:

```r
cmdstanr::cmdstan_path()
```

## Quick start (bundled example)

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

`bpnmf_example_config()` is a smoke-test configuration (200 warmup + 200
sampling iterations): fast enough to verify the pipeline end-to-end, not
suitable for inference. It is also the best starting template — every field
you would set for your own data is set there. Print it (`cfg`) or read its
source (`R/example.R`) to see a complete, working configuration.

---

# Using bpnmf on your own panel data

## 1. Shape your data

One CSV, **one row per (unit, time)**, with a separate column for each
outcome count (and, optionally, each denominator). Outcome columns are the
*groups* the model runs over — the K dimension:

| state   | time       | exposed_births | births_total | births_nhblack | pop_total | pop_nhblack |
|---------|------------|----------------|--------------|----------------|-----------|-------------|
| Alabama | 2016-01-01 | 0              | 9509         | 2837           | 952636    | 292788      |
| Alabama | 2016-03-01 | 0              | 9334         | 2727           | 952636    | 292788      |
| Texas   | 2022-09-01 | 1              | …            | …              | …         | …           |

Requirements and gotchas:

- **Time must be a date**, not a number. A bare `2016` in the time column is
  an error; use `2016-01-01`. With `date_format: auto` the parser tries
  `%Y-%m-%d`, `%m/%d/%y`, `%m/%d/%Y`, `%d-%m-%Y`; otherwise give an explicit
  `strptime` format.
- **Treatment is 0/1 per row**, so *staggered adoption is just the indicator
  turning on at different times for different units* — nothing else to
  declare. The indicator is one column shared by all groups.
- **The panel must be balanced** across (group, unit, time). If some cells
  genuinely do not exist, set `allow_unbalanced_panel: true` and they are
  marked missing (with a warning). Duplicate (group, unit, time) rows are an
  error.
- **Denominators are optional but recommended.** They enter as an exposure
  offset and must be present and strictly positive on every row (a `NA` or
  `0` denominator is an error, not a dropped row). Denominators are divided
  by 10,000 internally, so effects are reported as rates per 10k.
  Without denominators the model works on raw counts.
- **Suppressed small counts should be blank/`NA`** in the outcome column, not
  zero. See `adjust_for_missingness` below.

## 2. Where configs live

There is no configuration directory inside the package — **the config is
yours, and lives in your own project**. It comes in two interchangeable
forms:

| Form | How you make it | When to use it |
|------|-----------------|----------------|
| An R object | `bpnmf_config()` plus the `bpnmf_*_opts()` constructors | interactive work, scripts, programmatic sweeps |
| A YAML file | any path you like, conventionally `configs/myproject.yaml`, loaded with `read_bpnmf_config("configs/myproject.yaml")` | reproducible runs, sharing, configs already written for the Python package |

Both produce the same validated `bpnmf_config` object; every function
(`bpnmf_data()`, `bpnmf_fit()`, `bpnmf_run()`) takes that object, and
`bpnmf_run()` also accepts a path to the YAML directly. Validation happens
the moment you build the config, before any sampling: unknown keys, missing
cross-field requirements, and contradictory options are all errors with the
offending YAML path named.

Relative `input_file` / `output_dir` paths in a YAML config resolve against
your **current working directory**, not the config file's directory — run
from the project root, or use absolute paths.

## 3. A complete annotated config

Copy this to `configs/myproject.yaml` and edit. Every key is shown; all
except `data.input_file`, `data.output_dir`, and `data.schema` are optional
and appear here with their defaults.

```yaml
data:
  input_file: data/my_panel.csv
  output_dir: results/my_panel
  schema:
    unit_col: state              # panel unit
    time_col: time               # date column
    treatment_col: exposed       # 0/1 indicator
    # Option A: discover outcomes by column prefix
    outcomes_from_prefixes:
      outcome_prefix: births_    # births_total, births_nhblack, ...
      denominator_prefix: pop_   # must match label-for-label: pop_total, ...
      include: [total, nhblack]  # omit to take every matching column
    # Option B (mutually exclusive with A): list them explicitly
    # outcomes:
    #   - {outcome_col: births_total, label: total, denominator_col: pop_total}

  date_format: auto              # or an explicit strptime format
  start_date: "2016-01-01"       # inclusive
  end_date: "2024-01-01"         # EXCLUSIVE
  time_aggregation:
    enabled: false
    period: bimonthly            # calendar bins: monthly | bimonthly | quarterly | yearly
    # n_periods: 3               # ...or combine N consecutive periods (any resolution)
  allow_unbalanced_panel: false
  outcome: births                # label used in output filenames

model:
  outcome_distribution: NB       # NB | Poisson
  nb_disp: 0.0001                # fixed dispersion; concentration = 1/nb_disp
  sample_disp: false             # sample per-unit dispersion instead (NB only)
  adjust_for_missingness: true   # integrate over suppressed counts 1-9
  model_treated: true            # false = baseline-only fit, no effects
  inference_mode: joint          # joint | cut
  types:
    total:                       # name of the model type -> output subdirectory
      groups: [total]            # which outcome labels this type models
      ranks_to_test: [3]         # one fit per rank
      # total_from: [nhblack, nhwhite, hisp]   # build a synthetic "total"
      # total_all: false                       # ... or sum every label
      # exclude_units: [Alaska]

mcmc:
  auto_parallelism: true         # pick chains/cores automatically
  max_chains: 4
  # num_chains: 4                # manual (with auto_parallelism: false)
  # chain_method: parallel       # sequential | parallel
  num_warmup: 1000
  num_samples: 2500
  thinning: 10                   # retain num_samples / thinning draws
  target_accept: 0.8             # NUTS adapt_delta
  random_seed: 8675309
  progress_bar: true
  # gate_params: [mu_ctrl, te]   # restrict the convergence gate
  convergence:
    rhat_warn: 1.01
    rhat_fail: 1.05
    ess_min: 400
    ess_fail_fraction: 0.25
    divergence_fail_fraction: 0.01   # share of retained draws; 0 = allow none

output:
  figures: false                 # true/all/none, or a list of figure names:
                                 # unit_fit, unit_gap, raw_rate, interval,
                                 # group_comparison, ppc
  clean: false                   # wipe the type's output dir before writing
  save_traces: false             # also save full draws as .rds
  target_unit: Texas             # highlighted unit
                                 # (default: the aggregate unit, else the
                                 #  treated unit with the most periods)
  report_groups: [total]
  fit_gap_per_unit: false        # fit/gap figure for every treated unit
  print_tables: true
  print_target_table: true
  draws_format: csv              # csv | parquet
  # ppc_units: [Texas]
  # ppc_exclude_units: [Alaska]
  # ppc_acf_lags: [6]
  # ppc_unit_corr_max_time: "2022-01-01"
  interval_aggregates: true      # show aggregate units in interval.png
  # aggregate_units:             # synthetic reporting-only units
  #   - unit: "All treated"
  #     include_treated_units: true   # or include_all_units / include_units

# Only read when model.inference_mode is "cut":
cut:
  num_stage1_draws: 25
  stage2_draws_per_component: 100
  # selection_seed: 8675311      # default: random_seed + 2
  # stage2_seed: 8675312         # default: random_seed + 3
  # stage2_mcmc:                 # overrides merged over the mcmc block
  #   num_warmup: 500
  #   num_samples: 500
  # stage1_method: variational   # sample (default, NUTS) | variational (ADVI)
  # stage1_variational:          # ADVI knobs; only read under "variational"
  #   algorithm: meanfield       # meanfield | fullrank
  #   draws: 1000
  #   iter: 10000
  #   tol_rel_obj: 0.01
```

Then:

```r
cfg <- read_bpnmf_config("configs/myproject.yaml")
cfg                       # prints a summary: input, mode, types, mcmc
dat <- bpnmf_data(cfg)    # prints K groups, D units, N periods, exposed cells
bpnmf_run(cfg)
```

Inspect `dat` before you sample — it is the cheapest check that the schema
did what you meant. Zero exposed cells, an unexpected unit count, or a
surprising number of missing cells all show up here in seconds instead of
after an hour of NUTS.

The same config built in R:

```r
cfg <- bpnmf_config(
  input_file = "data/my_panel.csv",
  output_dir = "results/my_panel",
  schema = bpnmf_schema(
    unit_col = "state", time_col = "time", treatment_col = "exposed",
    outcomes_from_prefixes = bpnmf_prefixes(
      outcome_prefix = "births_", denominator_prefix = "pop_",
      include = c("total", "nhblack")
    )
  ),
  model = bpnmf_model_opts(
    outcome_distribution = "NB",
    types = list(total = bpnmf_type(groups = "total", ranks_to_test = 3))
  ),
  mcmc = bpnmf_mcmc_opts(iter_warmup = 1000, iter_sampling = 2500, thin = 10),
  output = bpnmf_output_opts(figures = TRUE, target_unit = "Texas"),
  start_date = "2016-01-01", end_date = "2024-01-01",
  time_aggregation = bpnmf_time_aggregation(enabled = TRUE, period = "bimonthly")
)
```

Note that the R constructors use cmdstanr's MCMC names while the YAML uses
the Python package's: `num_warmup` → `iter_warmup`, `num_samples` →
`iter_sampling`, `thinning` → `thin`, `target_accept` → `adapt_delta`,
`random_seed` → `seed`, `progress_bar` → `progress`, `num_chains` →
`chains`. Everything else keeps its name in both forms.

## 4. The decisions that actually matter

**Groups and model types.** A *type* is one named bundle of outcome groups
fit together, and gets its own output subdirectory. Splitting outcomes across
types (`total`, `age`, `race`) means separate fits; putting several groups
inside one type means they share the hierarchical treatment block and borrow
strength from each other. A group named `total` that is not itself an outcome
column must be constructed — give the type either `total_from: [a, b, c]` or
`total_all: true`, or you get an error.

**Rank.** `ranks_to_test` is the number of latent factors in the
factorization of the untreated surface. It is the main capacity knob: too low
and the pre-treatment fit is visibly biased, too high and the model can start
absorbing the treatment effect itself. Give it several values
(`ranks_to_test: [2, 3, 4, 5]`) — each is fit separately into
`<type>/rank_<r>/` — and compare pre-period fit and the PPC suite across
them. Rank 3 is a reasonable starting point for panels the size of the
bundled example (≈50 units × ≈50 periods).

**Likelihood.** `NB` is the default and the safe choice for count data with
any overdispersion; `Poisson` is a good deal faster if the mean-variance
relationship really holds. `nb_disp` fixes the dispersion (concentration
`1/nb_disp`); `sample_disp: true` estimates it per unit instead, which costs
sampling time and needs `NB`.

**Missingness.** `adjust_for_missingness: true` is meaningful only for data
that suppresses small counts, i.e. where a blank cell means "somewhere in
1–9" (US vital-statistics conventions, among others). It integrates the
likelihood over exactly those values and additionally tells the model that
every *observed* cell was **not** in 1–9. If your `NA`s mean "genuinely
unknown" or "not collected", set it to `false` — those cells are then simply
dropped from the likelihood.

**Time filtering and aggregation.** `end_date` is exclusive. Aggregation sums
counts, averages denominators, and takes the max of the treatment indicator
within each period — so a period that is partly exposed counts as exposed.
Choose the period to balance signal and length: coarser periods mean less
noise per cell but fewer time points for the factorization to work with.

`time_aggregation` bins the time axis one of two ways, and you set exactly one
of them. `period` bins by the **calendar** — rows are grouped into calendar
months, two-month blocks, quarters, or years, so bins land on calendar
boundaries (Q1 is always Jan–Mar) whatever date the panel starts on. It
assumes monthly-or-finer input, since it bins on year and month.

`n_periods` bins by **position**: every N consecutive time points in the panel
are combined, at whatever resolution the panel actually has. Daily data over a
single month has no calendar bin to fall back on — `period: monthly` would
collapse it to one point — but `n_periods: 7` turns it into weeks. Blocks are
cut from the panel's sorted distinct times rather than per unit, so units stay
aligned even when one is missing a period, and a short trailing block is kept
with its real (shorter) exposure recorded in `start_date`/`end_date`, so
person-year rates stay correct.

The older key name `aggregation` still loads, with a deprecation warning.

**Sampling and the convergence gate.** The defaults (1000 warmup, 2500
sampling, thin 10) are a real run, not a smoke test; expect a substantial
wait on a full panel. Raise `target_accept` toward 0.95 if you see
divergences. `converged` requires a PASS on both R-hat and ESS *and* a
divergent-transition rate at or below `divergence_fail_fraction` — 1% of the
retained draws by default, with 0 demanding none at all. The gate uses the
rate rather than a raw count so that the threshold means the same thing
however long you sample. Note that the factorization is
rotation-non-identifiable by design: the individual factors (`time_fac`,
`unit_weight`) are not expected to mix well, and it is the combined
quantities — the baseline log-rate surface `mu_ctrl` and the treatment
effects `te`, which the draws frame reports as `mu` and `mu_treated` — that
are identified and interpretable. If the gate is dominated by the factor
parameters rather than by anything you report, restrict it with
`gate_params: [mu_ctrl, te]`; divergences still count run-wide.
`parameter_diagnostics(fit)` shows the per-parameter breakdown behind the
gate verdict.

**Target unit.** `target_unit` selects the unit highlighted in tables and
per-unit figures. For a summary across treated units, add an
`aggregate_units` entry with `include_treated_units: true` — these are
reporting-only synthetic units, computed from draws after fitting, so they do
not change the model.

Left unset, `target_unit` is auto-detected: **a configured aggregate unit wins
when there is one**, on the reasoning that if you defined a pooled unit, the
pooled effect is the headline. Failing that, it falls back to the treated unit
with the most treated periods. So a config with `aggregate_units` and no
`target_unit` reports "All treated", not whichever single state happens to
have the longest exposure.

Aggregate units also appear in `interval.png`, in their own band above the
individual units — they pool the same draws as the units they cover, so
ranking them together would read as a peer comparison when it is not. Set
`interval_aggregates: false` to plot only the real units.

## 5. What lands on disk

`bpnmf_run(cfg)` writes, per model type (and per rank when more than one rank
is requested):

```
<output_dir>/<type>/
  df_<type>.csv                              # the standardized long panel
  {NB|Poisson}_{outcome}_{type}_{rank}.csv   # tidy draws (or .parquet)
  ..._convergence.json                       # gate: R-hat, ESS, divergences
  [rank_<rank>/]figs/...                     # figures, if output.figures
```

A failed gate is a warning, not a stop — artifacts are still written so you
can diagnose the run.

## Cut (modular) inference

In cut mode the baseline is fit on control cells only (stage 1), a set of
stage-1 draws is frozen as "cut components", and treatment effects are
estimated conditional on each (stage 2). This stops the treated cells from
informing the baseline surface.

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

Cost scales as `num_stage1_draws` × the stage-2 run, so stage 2 is usually
configured much shorter than stage 1 via `cut.stage2_mcmc`. Retained draw
counts must be equal across components so that pooling weights them equally.
Seeds are derived from `mcmc.random_seed` unless you set `selection_seed` /
`stage2_seed` explicitly; `cut.stage2_mcmc` may not set a seed of its own.

### Fast stage 1 with ADVI

Stage 1 is the expensive half — the full factorization over every control
cell — and full MCMC on it can take hours. `stage1_method = "variational"`
runs Stan's ADVI there instead, turning that into minutes. Everything
downstream is unchanged: ADVI still yields draws of `mu_ctrl`, and the same
seeded, stratified selection promotes some of them to cut components.

```r
cfg <- bpnmf_example_config(
  model = bpnmf_model_opts(
    types = list(total = bpnmf_type("total", 3)),
    inference_mode = "cut"
  ),
  cut = bpnmf_cut_opts(
    num_stage1_draws = 25,
    stage1_method = "variational",
    stage1_variational = list(algorithm = "meanfield", draws = 1000)
  )
)
cfit <- bpnmf_cut_fit(bpnmf_data(cfg), config = cfg)
```

What you give up is real, and it is uncertainty rather than location.
Mean-field ADVI fits a factorized Gaussian in the unconstrained space, so it
understates marginal variance and drops posterior correlation between
baseline parameters. In cut mode that error propagates in one direction: the
stage-1 components span too narrow a range, so the pooled `te` posterior is
too narrow too and its intervals under-cover. Point estimates are usually
close; intervals are not trustworthy.

There is also no convergence gate to lean on — R-hat, ESS and divergences
are all chain-based quantities that ADVI simply does not have. The manifest
records this rather than papering over it:

- `manifest$stage1$converged` is `NA` (“not gated”), never `TRUE`/`FALSE`
- `manifest$stage1_gated` is `FALSE`, and `manifest$converged` then reflects
  only the stage-2 fits, which are still full MCMC and still gated
- `parameter_diagnostics()` and `bpnmf_trace_plot()` error on a variational
  fit instead of returning a meaningless single-chain R-hat

CmdStan’s own run-specific complaints (“the variational approximation may be
poor”, “maximum number of iterations is reached”) are surfaced as R warnings.

Use ADVI to iterate on rank, priors, and data prep; re-run with
`stage1_method = "sample"` for anything you intend to report. The same switch
is available directly on a single fit as `bpnmf_fit(..., method =
"variational")`.

## Parity with the Python implementation

- The Stan models are checked against the NumPyro models by **exact
  log-density comparison** at shared parameter values (fixtures exported from
  Python; `tests/testthat/test-logdensity-parity.R`).
- Cell ordering, chain-stratified cut selection, output thinning, and
  convergence-gate bands are fixture-tested against the Python
  implementation.
- MCMC draws are **not** bitwise-identical across implementations (different
  samplers and RNGs); equivalence is distributional.
- Configs written for the Python package load unchanged via
  `read_bpnmf_config()`; unknown keys are rejected the same way.

## License

MIT.
