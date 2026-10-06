# The claim rank shrinkage makes is a prior claim: the amplitude of the
# low-rank term, and the number of components the units actually load on,
# should stop depending on the truncation R. These tests check it against the
# real Stan program by sampling it with an empty likelihood -- every
# observation removed, so the posterior *is* the prior -- at a small and a
# large rank.
#
# Testing the Stan model rather than an R re-derivation of the same algebra
# is the point: the invariance lives in how the stick-breaking profile and
# the per-unit softmax deviation interact, which is exactly what a
# re-implementation would get wrong in the same way twice.

# The joint model's data list with every likelihood term dropped.
prior_only_data <- function(rank, rank_shrinkage) {
  data <- make_test_data() |> (\(df) {
    arrays <- build_model_arrays(df, c("g1", "g2"))
    arrays$df <- df
    arrays$type <- "both"
    structure(arrays, class = c("countsynth_data", "list"))
  })()
  sd <- stan_data_joint(
    data,
    rank = rank, model_treated = FALSE, adjust_for_missingness = FALSE,
    gen_ypred = FALSE, rank_shrinkage = rank_shrinkage
  )
  empty_int <- as.array(integer(0))
  sd$n_obs <- 0L
  sd$obs_cell <- empty_int
  sd$y <- empty_int
  sd$obs_unit <- empty_int
  sd$n_cens <- 0L
  sd$cens_cell <- empty_int
  sd$cens_unit <- empty_int
  sd$n_notcens <- 0L
  sd$notcens_cell <- empty_int
  sd$notcens_unit <- empty_int
  sd
}

# E[sum_r unit_weight[k, d, r]^2] over draws, groups and units. This is the
# factor that multiplies the temporal prior's variance in the low-rank term,
# so it is what `factor_variation_pct` is implicitly scaled by; its inverse
# is the effective number of components a unit loads on.
prior_weight_concentration <- function(rank, rank_shrinkage, seed = 4242) {
  model <- countsynth_stan_model("joint")
  fit <- model$sample(
    data = prior_only_data(rank, rank_shrinkage),
    chains = 1, iter_warmup = 400, iter_sampling = 600, seed = seed,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE
  )
  # Under shrinkage the sampled simplex is empty and the realized weights
  # come from generated quantities instead; whichever is non-empty is the one
  # the likelihood used.
  var <- if (is.null(rank_shrinkage)) "unit_weight" else "unit_weight_fitted"
  w <- posterior::draws_of(
    posterior::as_draws_rvars(fit$draws(variables = var))[[var]]
  )
  # [draw, group, unit, component] -> mean over everything but the component
  # axis, which is summed.
  mean(apply(w^2, c(1, 2, 3), sum))
}

test_that("the flat weight prior makes the low-rank amplitude depend on R", {
  skip_on_cran()
  skip_if_no_cmdstan()
  # Under Dirichlet(1, ..., 1), E[sum_r w^2] = 2 / (R + 1) exactly: offering
  # more components silently tightens the prior on the structure they are
  # there to express. This is the behaviour rank shrinkage replaces, and it
  # is asserted here so that the comparison below has a baseline that is
  # measured rather than assumed.
  for (R in c(4, 12)) {
    expect_equal(
      prior_weight_concentration(R, NULL), 2 / (R + 1),
      tolerance = 0.08
    )
  }
})

test_that("rank shrinkage makes the low-rank amplitude R-invariant", {
  skip_on_cran()
  skip_if_no_cmdstan()
  opts <- countsynth_rank_shrinkage_opts()
  small <- prior_weight_concentration(6, opts)
  large <- prior_weight_concentration(18, opts)
  # Tripling the truncation moves this by a few percent, where the flat prior
  # would cut it by two thirds (2/19 vs 2/7). The tolerance is loose because
  # 600 draws of a heavy-tailed functional is a noisy estimate; the effect
  # being tested is a factor of three.
  expect_equal(small, large, tolerance = 0.2)
  expect_gt(large, 2 * 2 / (18 + 1))

  # And the sparsity is real: the effective number of components in the
  # shared profile stays well inside the truncation at both ranks.
  expect_lt(1 / large, 6)
})

test_that("rank_shrink = 0 reproduces the pre-shrinkage model draw for draw", {
  skip_on_cran()
  skip_if_no_cmdstan()
  # The parameters rank shrinkage adds are all zero-size when it is off, so
  # the unconstrained vector is unchanged and a given seed must retrace the
  # same path. test-logdensity-parity.R is the sharper check of the same
  # property but needs `compile_model_methods`, whose Rcpp build is broken on
  # some toolchains; this one only needs a plain fit, so it always runs.
  data <- make_test_data() |> (\(df) {
    arrays <- build_model_arrays(df, c("g1", "g2"))
    arrays$df <- df
    arrays$type <- "both"
    structure(arrays, class = c("countsynth_data", "list"))
  })()
  model <- countsynth_stan_model("joint")
  sd <- stan_data_joint(data, rank = 3, gen_ypred = FALSE)
  draws_at <- function(seed) {
    fit <- model$sample(
      data = sd, chains = 1, iter_warmup = 120, iter_sampling = 60,
      seed = seed, refresh = 0, show_messages = FALSE,
      show_exceptions = FALSE
    )
    posterior::as_draws_matrix(fit$draws(variables = c("mu_ctrl", "te")))
  }
  # Same seed twice through the same code path: exact equality is the
  # baseline this test's sensitivity rests on. If Stan were nondeterministic
  # here, the cross-version claim would be untestable this way.
  expect_identical(draws_at(4242), draws_at(4242))
  # And the shrinkage fields really are inert when off: perturbing every one
  # of them must not move a single draw.
  sd_perturbed <- sd
  sd_perturbed$group_mass_shape <- 7
  sd_perturbed$group_mass_rate <- 0.3
  sd_perturbed$unit_sd_scale <- 5
  fit_perturbed <- model$sample(
    data = sd_perturbed, chains = 1, iter_warmup = 120, iter_sampling = 60,
    seed = 4242, refresh = 0, show_messages = FALSE, show_exceptions = FALSE
  )
  expect_identical(
    posterior::as_draws_matrix(
      fit_perturbed$draws(variables = c("mu_ctrl", "te"))
    ),
    draws_at(4242)
  )
})
