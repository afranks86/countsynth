# countsynth

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
  scales), optionally driven by **covariates** through an lme4-style formula
  so the effect can vary with time since treatment or a mediator,
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

# install countsynth from source
# install.packages("remotes")
remotes::install_github("afranks86/countsynth")
```

Check the toolchain before a long run — `cmdstanr::cmdstan_path()` should
print a path rather than error:

```r
cmdstanr::cmdstan_path()
```

## Quick start (bundled example)

```r
library(countsynth)

cfg  <- countsynth_example_config()          # bundled US fertility data
dat  <- countsynth_data(cfg)                 # CSV -> standardized panel + arrays
fit  <- countsynth_fit(dat, config = cfg)    # NUTS via cmdstanr

summary(fit)                            # convergence gate
parameter_diagnostics(fit)              # per-parameter R-hat / ESS table

draws <- countsynth_draws(fit)               # tidy long draws frame
countsynth_summary_table(draws)              # headline observed vs expected
plot(fit, which = "unit_fit", unit = "Texas")
countsynth_interval_plot(draws, estimand = "ratio")

# or run everything (fits + artifacts + figures) from the config:
countsynth_run(cfg)
```

`countsynth_example_config()` is a smoke-test configuration (200 warmup + 200
sampling iterations): fast enough to verify the pipeline end-to-end, not
suitable for inference. It is also the best starting template — every field
you would set for your own data is set there. Print it (`cfg`) or read its
source (`R/example.R`) to see a complete, working configuration.

---

# Using countsynth on your own panel data

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
| An R object | `countsynth_config()` plus the `countsynth_*_opts()` constructors | interactive work, scripts, programmatic sweeps |
| A YAML file | any path you like, conventionally `configs/myproject.yaml`, loaded with `read_countsynth_config("configs/myproject.yaml")` | reproducible runs, sharing, configs already written for the Python package |

Both produce the same validated `countsynth_config` object; every function
(`countsynth_data()`, `countsynth_fit()`, `countsynth_run()`) takes that object, and
`countsynth_run()` also accepts a path to the YAML directly. Validation happens
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
  # factor_variation_pct: 25     # expected swing of the unit-specific temporal
                                 # factors; unset = Gamma(20, 20)
  # time_level_variation_pct: 60 # expected swing of the common time level;
                                 # unset = Gamma(1, 1), which is very diffuse
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
  gate_params: [mu_ctrl, te]     # gate on these prefixes (this is the default;
                                 # "all" gates on every sampled variable)
  convergence:
    rhat_warn: 1.01
    rhat_fail: 1.05
    ess_min: 400
    ess_fail_fraction: 0.25
    divergence_fail_fraction: 0.01   # share of retained draws; 0 = allow none

output:
  figures: false                 # true/all/none, or a list of figure names:
                                 # unit_fit, unit_gap, raw_rate, interval,
                                 # group_comparison, ppc, te_regression
  clean: false                   # wipe the type's output dir before writing
  save_traces: false             # also save full draws as .rds
  target_unit: Texas             # highlighted unit
                                 # (default: the aggregate unit, else the
                                 #  treated unit with the most periods)
  report_groups: [total]
  fit_gap_per_unit: false        # fit/gap figure for every treated unit
  print_tables: true             # by-unit summary table in the terminal
  print_target_table: false      # also print the target unit's own table
  html_tables: true              # write gt HTML tables (needs the gt package)
  draws_format: csv              # csv | parquet
  # ppc_units: [Texas]
  # ppc_exclude_units: [Alaska]
  # ppc_acf_lags: [6]
  # ppc_unit_corr_max_time: "2022-01-01"
  interval_aggregates: true      # show aggregate units in interval.png
  # aggregate_units:             # synthetic reporting-only units
  #   - unit: "All treated"
  #     include_treated_units: true   # or include_all_units / include_units

  # Optional: covariate model for the treatment effect (see the
  # "Treatment-effect covariates" vignette). Omit for the default
  # group / unit / group:unit hierarchy.
  # treatment_effects:
  #   formula: ~ 1 + event_time + dist_change + (1 + event_time | unit)
  #   covariates_file: clinic_distance.csv
  #   standardize: true

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
cfg <- read_countsynth_config("configs/myproject.yaml")
cfg                       # prints a summary: input, mode, types, mcmc
dat <- countsynth_data(cfg)    # prints K groups, D units, N periods, exposed cells
countsynth_run(cfg)
```

Inspect `dat` before you sample — it is the cheapest check that the schema
did what you meant. Zero exposed cells, an unexpected unit count, or a
surprising number of missing cells all show up here in seconds instead of
after an hour of NUTS.

The same config built in R:

```r
cfg <- countsynth_config(
  input_file = "data/my_panel.csv",
  output_dir = "results/my_panel",
  schema = countsynth_schema(
    unit_col = "state", time_col = "time", treatment_col = "exposed",
    outcomes_from_prefixes = countsynth_prefixes(
      outcome_prefix = "births_", denominator_prefix = "pop_",
      include = c("total", "nhblack")
    )
  ),
  model = countsynth_model_opts(
    outcome_distribution = "NB",
    types = list(total = countsynth_type(groups = "total", ranks_to_test = 3))
  ),
  mcmc = countsynth_mcmc_opts(iter_warmup = 1000, iter_sampling = 2500, thin = 10),
  output = countsynth_output_opts(figures = TRUE, target_unit = "Texas"),
  start_date = "2016-01-01", end_date = "2024-01-01",
  time_aggregation = countsynth_time_aggregation(enabled = TRUE, period = "bimonthly")
)
```

Each option group has its own constructor so that ~80 options do not collapse
into one unreadable signature. You do not have to call them, though: anywhere
a `countsynth_*` object is expected you can pass a **plain named list** of the same
arguments, which is then run through that very constructor — same validation,
same defaults, same result. So the config above can be written as one call
with the same shape as the YAML:

```r
cfg <- countsynth_config(
  input_file = "data/my_panel.csv",
  output_dir = "results/my_panel",
  schema = list(
    unit_col = "state", time_col = "time", treatment_col = "exposed",
    outcomes_from_prefixes = list(
      outcome_prefix = "births_", denominator_prefix = "pop_",
      include = c("total", "nhblack")
    )
  ),
  model = list(
    outcome_distribution = "NB",
    types = list(total = list(groups = "total", ranks_to_test = 3))
  ),
  mcmc = list(iter_warmup = 1000, iter_sampling = 2500, thin = 10),
  output = list(figures = TRUE, target_unit = "Texas"),
  start_date = "2016-01-01", end_date = "2024-01-01",
  time_aggregation = list(enabled = TRUE, period = "bimonthly")
)
```

Mix the two freely. A misspelled name in a list is an error naming the field
and listing the valid arguments, so the list form is not a way to smuggle a
typo past validation. Reach for the constructors when you want argument
completion and `?countsynth_model_opts` at your fingertips; reach for lists when
you want one call that mirrors the YAML.

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

**Latent factor scale.** `time_fac`, the rank-`R` temporal basis, has a
`Gamma(shape, shape)` prior: mean 1, entering the log-rate as
`log(time_fac)`, so `sd(log time_fac) ≈ 1/√shape`. It sets how far the
low-rank temporal structure may swing multiplicatively away from a unit's own
level — the *unit × time interaction*, since unit levels are carried by
`unit_fe` and the common trend by `time_fe`.

`factor_variation_pct` sets it in interpretable units: an expected swing of
`p` percent means `shape = 1 / log(1 + p/100)²`.

| `factor_variation_pct` | shape | central 95% multiplier |
|---|---|---|
| 10 | 110 | ×0.83 – ×1.20 |
| 25 | 20 | ×0.61 – ×1.48 |
| 50 | 6.1 | ×0.40 – ×2.05 |
| 100 | 2.1 | ×0.13 – ×2.75 |

Left unset it is `shape = 20`, carried over from the fertility application —
about a ±25% swing. Whether that is loose or tight depends entirely on the
grouping: on the bundled fertility panel the observed unit × time interaction
is ~3% for the race groups (`shape` ≈ 1000, so the prior is ~7× looser than
the data needs) but ~21% for the marital-status groups (`shape` ≈ 22, which
the default matches almost exactly). That order-of-magnitude spread within one
dataset is why the default is a fixed weakly-informative value rather than
something estimated: a fixed multiple of the empirical scale would be
reasonable for one grouping and degenerate for another.

**The common time level.** `time_fe` carries the group-wide temporal level,
shared by every unit, and has the same `Gamma(shape, shape)` form set by
`time_level_variation_pct`. Unset it is `Gamma(1, 1)`, which is far more
diffuse than it looks: `sd(log time_fe) = 1.28` — a central 95% range of ×0.03
to ×3.7 — making it the loosest prior in the model.

That matters beyond prior belief. `unit_fe_mu` has a flat prior and the
likelihood sees only `log(time_fe) + unit_fe_mu`, so the two trade off exactly
and this prior is the *only* thing separating them. On the bundled example
their posterior correlation is −0.99, and they are the worst-mixing parameters
in the model. Tightening `time_level_variation_pct` narrows that ridge: at a
50% swing (shape 6.6) the correlation falls to about −0.91 — and keeps
falling as you tighten further, to about −0.82 at 25% — while `unit_fe_mu`'s
worst R-hat drops from ~1.27 to ~1.02, and `time_fe`'s from ~1.24 to ~1.10
with bulk ESS improving by something between 1.5× and 3×. Those figures come
from short two-chain runs and move around a fair bit between replicates — the
direction is reliable, the magnitudes are not.

What it does *not* do is change the answer. `mu_ctrl`'s worst R-hat does not
improve materially at any setting — across replicates it wanders between
about 1.10 and 1.19 with no clear relation to the prior — because only the sum
was ever identified — tightening reallocates
it between two parameters rather than learning anything new. So this knob
buys cleaner diagnostics on the nuisance parameters, not a better posterior,
which is why the loose default is left alone. It is worth reaching for if you
gate on all parameters rather than the `mu_ctrl` / `te` default, where the
split otherwise dominates the verdict.

Both knobs are read as an expected multiplicative swing and inverted exactly
(`trigamma(shape) = log(1 + p/100)²`), not through the usual `1/sd²`
approximation — which is fine above shape 10 but 28% off at shape 1, exactly
where `time_fe` sits.

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
however long you sample.

**The gate covers `mu_ctrl` and `te` by default**, not every sampled
variable. The factorization is rotation-non-identifiable by design: the
individual factors (`time_fac`, `unit_weight`) are not expected to mix well,
because chains that settle on different factor labelings give them an
enormous R-hat while nothing you report has moved. The identified,
interpretable quantities are the baseline log-rate surface `mu_ctrl` and the
treatment effects `te` — which the draws frame reports as `mu` and
`mu_treated` — so those are what the verdict is built from. Widen it by
listing more prefixes, or set `gate_params: all` to gate on every sampled
variable. Divergences always count run-wide, whatever the gate covers, and
`parameter_diagnostics(fit)` reports every variable with a `gated` column
showing which ones counted.

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

**Tables.** Every run writes `summary_table_by_unit.csv` (display-formatted:
counts, rates per 1,000 person-years, pre-formatted CIs, `*` for a two-sided
posterior p < 0.05) and `post_treatment_summary.csv` (the same estimands as
plain numeric columns, for joining and plotting). With the `gt` package
installed you also get `summary_table.html` and `summary_table_by_unit.html`,
which are the display tables rendered for publication — counts and rates under
their own spanners, one row group per unit.

Build one yourself from a draws frame:

```r
draws <- countsynth_draws(fit)
countsynth_gt_table(draws)                  # headline unit, one row per group
countsynth_gt_table(draws, by_unit = TRUE)  # every treated unit, grouped
gt::gtsave(countsynth_gt_table(draws, by_unit = TRUE), "summary.html")
```

The terminal prints the by-unit table only. `print_target_table: true` adds
the target unit's own table above it — its rows are already in the by-unit
table, so it is off by default.

## 5. What lands on disk

`countsynth_run(cfg)` writes, per model type (and per rank when more than one rank
is requested):

```
<output_dir>/<type>/
  df_<type>.csv                              # the standardized long panel
  {NB|Poisson}_{outcome}_{type}_{rank}.csv   # tidy draws (or .parquet)
  ..._convergence.json                       # gate: R-hat, ESS, divergences
  [rank_<rank>/]figs/
    summary_table_by_unit.csv                # display-formatted, per unit
    summary_table*.html                      # the same, via gt
    post_treatment_summary.csv               # numeric estimands + CIs
    expected_vs_observed.csv                 # per (unit, time, group) detail
    *.png                                    # figures, if output.figures
```

A failed gate is a warning, not a stop — artifacts are still written so you
can diagnose the run.

## Cut (modular) inference

In cut mode the baseline is fit on control cells only (stage 1), a set of
stage-1 draws is frozen as "cut components", and treatment effects are
estimated conditional on each (stage 2). This stops the treated cells from
informing the baseline surface.

```r
cfg  <- countsynth_example_config(
  model = countsynth_model_opts(
    types = list(total = countsynth_type("total", 3)),
    inference_mode = "cut"
  ),
  cut = countsynth_cut_opts(num_stage1_draws = 25, stage2_draws_per_component = 100)
)
cfit <- countsynth_cut_fit(countsynth_data(cfg), config = cfg)
summary(cfit)               # stage-1 gate + per-component table
draws <- countsynth_draws(cfit)  # pooled, with cut_component / stage1_* provenance
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
cfg <- countsynth_example_config(
  model = countsynth_model_opts(
    types = list(total = countsynth_type("total", 3)),
    inference_mode = "cut"
  ),
  cut = countsynth_cut_opts(
    num_stage1_draws = 25,
    stage1_method = "variational",
    stage1_variational = list(algorithm = "meanfield", draws = 1000)
  )
)
cfit <- countsynth_cut_fit(countsynth_data(cfg), config = cfg)
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
- `parameter_diagnostics()` and `countsynth_trace_plot()` error on a variational
  fit instead of returning a meaningless single-chain R-hat

CmdStan’s own run-specific complaints (“the variational approximation may be
poor”, “maximum number of iterations is reached”) are surfaced as R warnings.

Use ADVI to iterate on rank, priors, and data prep; re-run with
`stage1_method = "sample"` for anything you intend to report. The same switch
is available directly on a single fit as `countsynth_fit(..., method =
"variational")`.

## Treatment-effect covariates

By default the treatment effect is a group / unit / group:unit hierarchy,
which says how large the effect is but not what moves it. A
`treatment_effects` formula replaces that hierarchy with a regression
surface, in lme4-style syntax, shared by joint and cut inference:

```r
cfg <- countsynth_example_config(
  model = countsynth_model_opts(
    types = list(total = countsynth_type("total", 3)),
    treatment_effects = countsynth_te_opts(
      # Effect trending with time since treatment, plus a mediator, with
      # per-unit intercepts and slopes shrunk toward the surface.
      formula = ~ 1 + event_time + dist_change + (1 + event_time | unit),
      covariates = dist   # data frame or CSV, joined by unit / time / group
    )
  )
)
fit <- countsynth_fit(countsynth_data(cfg), config = cfg)

countsynth_te_regression_plot(fit, predictor = "event_time")            # surface
countsynth_te_regression_plot(fit, predictor = "event_time", by = "unit")  # by level
countsynth_te_coef_plot(fit, terms = "all")                             # forest
countsynth_te_coef_table(fit)                                           # summary
```

`event_time` (periods since the unit's first treated period), `time_idx`,
`group`, and `unit` are always available; covariates are joined onto the
exposed cells. All existing tables and figures are unchanged. Leaving
`treatment_effects` unset reproduces earlier results exactly, down to the
random-number stream.

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
  `read_countsynth_config()`; unknown keys are rejected the same way.

## License

MIT.
