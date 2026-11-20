# stanli vs CmdStan

Benchmarks for `mcmc.backend: stanli` (see `R/backend-stanli.R`), which runs the
bundled Stan programs through [stanli](https://github.com/seantalts/stanli), a
Stan interpreter: same program and data, Stan's own NUTS, its own reverse-mode
pass over precompiled Stan Math kernels, and no C++ compilation or CSV
round-trip.

Data come from the panels in `~/Dropbox/course/test_countsynth`, built exactly
as `countsynth_fit()` builds them (`common.R`), written once to JSON, and read
by both engines.

| script | what it does |
|---|---|
| `parity.R` | log density + gradient at random unconstrained points, stanli vs CmdStan's standalone `log_prob`; stanli's per-gradient cost |
| `sample.R` | one NUTS run of `joint.stan` on one engine; timings, leapfrogs, ESS, posterior summaries to `out/` |
| `run_all.sh` | `sample.R` over all cases, engines interleaved (300+300, 4 parallel chains) |
| `run_seeds.sh` | the same at three more seeds (posterior agreement against a same-engine null) |
| `e2e.R`, `run_e2e.sh` | through the package API at each config's own settings: fit, draws frame, gate; cut runs both stages |
| `analyze.R` | summary table from `out/*.rds` |

Run long jobs on AC power with the lid open: on battery with the lid closed
the Mac enters Maintenance Sleep and every run spanning it is garbage
(`caffeinate` does not prevent it). Check `pmset -g log | grep -E " (Sleep|Wake)  "`.

## Results (2026-11-20, M-series Mac, 8 performance cores, CmdStan 2.39.0 at stanc O1, stanli 0.30.1)

**Parity.** On all four panels, joint and baseline: log density within 3e-8
(float rounding on |lp| ~ 1e4), gradient within 3e-16 relative. No part of
either model falls back to stanli's interpreter. Model build 0.06-0.4 s.

**Sampling** (300 warmup + 300 sampling, 4 parallel chains, config's
`adapt_delta`):

| case | model | CmdStan | stanli | wall | per gradient | leapfrog ratio | median ESS (CS / sl) |
|---|---|---|---|---|---|---|---|
| fertility_full | NB, rank 3, Dirichlet weights | 19.3 s | 9.4 s | 2.06x | 2.14x | 1.03 | 1021 / 1004 |
| fertility_full (seed 101) | | 20.9 s | 8.9 s | 2.36x | 2.13x | 0.96 | 1054 / 1075 |
| mortality_con | Poisson, rank 3 | 24.7 s | 10.8 s | 2.28x | 2.30x | 1.00 | 1254 / 1294 |
| mortality_full | Poisson, rank 5, rank shrinkage, heavy censoring | 63.8 s | 29.7 s | 2.15x | 2.13x | 1.01 | 1379 / 1439 |
| fertility_te_race | NB, rank 4, rank shrinkage, te regression | 395.9 s | 193.4 s | 2.05x | 2.03x | 1.00 | 1216 / 1214 |

The whole gain is per-gradient cost: leapfrog counts and ESS match, so ESS per
second improves by the same 2.0-2.4x. Minimum ESS varied run to run on both
engines (fertility_full's stanli minimum was 60 at one seed and 215 at
another).

**Agreement.** Share of `mu_ctrl`/`te` cells whose posterior means differ by
more than 3 combined MCSEs: stanli vs CmdStan 0-2.1% across cases; on
fertility_full, CmdStan vs CmdStan at another seed 0.9%, stanli vs stanli
2.6%, stanli vs CmdStan 0.7%. No sign the engines target different
posteriors.

Not yet run cleanly: `run_seeds.sh` beyond seed 101 on fertility_full, and
`run_e2e.sh` (package-level, including the cut stage-2 loop). A sleep
contaminated the first attempt; those results are in `out/contaminated/`.
