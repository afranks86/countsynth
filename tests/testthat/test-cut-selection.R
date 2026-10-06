test_that("selection is chain-stratified, sorted, without replacement", {
  refs <- select_stage1_draws(
    n_chains = 4, per_chain = 100, m = 10, selection_seed = 7
  )
  expect_equal(nrow(refs), 10)
  expect_identical(refs$component, 1:10)
  # quotas 3,3,2,2
  expect_equal(as.integer(table(refs$chain)), c(3L, 3L, 2L, 2L))
  # sorted within chain, no duplicates
  for (c in 1:4) {
    iters <- refs$iteration[refs$chain == c]
    expect_identical(iters, sort(unique(iters)))
  }
  # global draw index is chain-major
  expect_identical(refs$draw, (refs$chain - 1L) * 100L + refs$iteration)
  # deterministic under the same seed
  refs2 <- select_stage1_draws(4, 100, 10, 7)
  expect_identical(refs, refs2)
  # different under a different seed
  refs3 <- select_stage1_draws(4, 100, 10, 8)
  expect_false(identical(refs$iteration, refs3$iteration))
})

test_that("selection errors when m exceeds retained draws", {
  expect_error(select_stage1_draws(2, 3, 10, 1), "exceeds")
})

test_that("cut settings derive seeds and merge stage2 overlay", {
  cfg <- countsynth_config(
    input_file = "x.csv", output_dir = tempdir(),
    schema = countsynth_schema(
      "u", "t", "tr",
      outcomes = list(countsynth_outcome("o", "g"))
    ),
    model = countsynth_model_opts(
      types = list(g = countsynth_type("g", 2)), inference_mode = "cut"
    ),
    mcmc = countsynth_mcmc_opts(seed = 100, iter_warmup = 1000, iter_sampling = 2500),
    cut = countsynth_cut_opts(
      num_stage1_draws = 5,
      stage2_mcmc = list(num_warmup = 200, num_samples = 400, thinning = 2)
    )
  )
  s <- resolve_cut_settings(cfg)
  expect_equal(s$selection_seed, 102L)
  expect_equal(s$stage2_seed, 103L)
  expect_equal(s$stage2_mcmc$iter_warmup, 200)
  expect_equal(s$stage2_mcmc$iter_sampling, 400)
  expect_equal(s$stage2_mcmc$thin, 2)
  # non-overlaid settings inherit from the top-level mcmc config
  expect_equal(s$stage2_mcmc$adapt_delta, 0.8)
  expect_equal(s$stage2_mcmc$seed, 100L)
})

test_that("seed collisions warn", {
  cfg <- countsynth_config(
    input_file = "x.csv", output_dir = tempdir(),
    schema = countsynth_schema(
      "u", "t", "tr",
      outcomes = list(countsynth_outcome("o", "g"))
    ),
    model = countsynth_model_opts(
      types = list(g = countsynth_type("g", 2)), inference_mode = "cut"
    ),
    mcmc = countsynth_mcmc_opts(seed = 100),
    cut = countsynth_cut_opts(selection_seed = 100)
  )
  expect_warning(resolve_cut_settings(cfg), "collides")
})

test_that("predictive seed stays in integer range and is component-distinct", {
  # the default stage-2 seed: 8675309 * 1000 overflows a 32-bit integer
  s <- predictive_seed(8675312L, 1L)
  expect_true(is.integer(s))
  expect_false(is.na(s))
  expect_no_error(withr::with_seed(s, runif(1)))
  expect_false(predictive_seed(8675312L, 1L) == predictive_seed(8675312L, 2L))
  # small seeds keep the plain arithmetic
  expect_identical(predictive_seed(7L, 3L), 7003L)
})
