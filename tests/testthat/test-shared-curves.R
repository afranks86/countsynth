# Shared temporal curves across groups, the crossed weight effects that come
# with them under rank shrinkage, and the rank-1 model without curves. The
# load-bearing checks are reconstructions: Stan's mu_ctrl must equal the
# surface rebuilt in R from the sampled pieces, which is the only thing tying
# "shared" and "dropped" to what the likelihood actually used.

sc_data <- function(groups = c("g1", "g2")) {
  df <- make_test_data()
  df <- df[df$group %in% groups, ]
  arrays <- build_model_arrays(df, groups)
  arrays$df <- df
  arrays$type <- "both"
  structure(arrays, class = c("bpnmf_data", "list"))
}

# Draws come back through CmdStan's CSV, so rebuilt values agree only to its
# output precision (~1e-7 here). Any structural error -- a wrong curve set,
# weights the likelihood didn't use -- shows up at 1e-2 or worse.
SC_TOL <- 1e-5

sc_fit <- function(sd) {
  bpnmf_stan_model("joint")$sample(
    data = sd, chains = 1, iter_warmup = 120, iter_sampling = 40, seed = 7,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE
  )
}

sc_draws <- function(fit, v) {
  posterior::draws_of(posterior::as_draws_rvars(fit$draws(variables = v))[[v]])
}

# max |mu_ctrl - rebuilt| over draws and cells. `weights` names the variable
# holding the weights the likelihood used; NULL means no curves (rank 1).
sc_reconstruction_error <- function(fit, sd, weights = NULL, shared = FALSE) {
  K <- sd$K
  D <- sd$D
  N <- sd$N
  lt <- log(sc_draws(fit, "time_fe"))
  fm <- sc_draws(fit, "unit_fe_mu")
  fs <- sc_draws(fit, "unit_fe_sigma")
  fz <- sc_draws(fit, "unit_fe_z")
  mu <- unclass(posterior::as_draws_matrix(fit$draws(variables = "mu_ctrl")))
  if (!is.null(weights)) {
    tf <- sc_draws(fit, "time_fac")
    w <- sc_draws(fit, weights)
  }
  err <- 0
  for (i in seq_len(nrow(mu))) {
    for (k in seq_len(K)) {
      for (d in seq_len(D)) {
        lf <- if (is.null(weights)) {
          rep(0, N)
        } else {
          log(matrix(tf[i, if (shared) 1 else k, , ], N) %*% w[i, k, d, ])
        }
        cells <- ((k - 1) * D + (d - 1)) * N + seq_len(N)
        rebuilt <- lf + lt[i, , k] + fm[i, k] + fs[i, k] * fz[i, d, k] +
          sd$log_denom[cells]
        err <- max(err, abs(rebuilt - mu[i, cells]))
      }
    }
  }
  err
}

test_that("share_fac is set only when there is something to share", {
  two <- sc_data()
  expect_equal(stan_data_joint(two, rank = 3)$share_fac, 0L)
  expect_equal(stan_data_joint(two, rank = 3, shared_curves = TRUE)$share_fac, 1L)
  # One group: the shared and per-group models are the same model.
  expect_equal(
    stan_data_joint(sc_data("g1"), rank = 3, shared_curves = TRUE)$share_fac, 0L
  )
  # Rank 1: there are no curves to share.
  expect_equal(stan_data_joint(two, rank = 1, shared_curves = TRUE)$share_fac, 0L)

  sd <- stan_data_joint(
    two, rank = 3, shared_curves = TRUE,
    rank_shrinkage = bpnmf_rank_shrinkage_opts(
      group_sd_prior = 0.3, shared_unit_sd_prior = 0.7
    )
  )
  expect_equal(sd$group_sd_scale, 0.3)
  expect_equal(sd$shared_unit_sd_scale, 0.7)
})

test_that("shared_curves and the crossed-effect priors validate", {
  expect_false(bpnmf_model_opts()$shared_curves)
  expect_true(bpnmf_model_opts(shared_curves = TRUE)$shared_curves)
  expect_error(bpnmf_model_opts(shared_curves = "yes"), "flag")
  o <- bpnmf_rank_shrinkage_opts()
  expect_equal(o$group_sd_prior, 0.5)
  expect_equal(o$shared_unit_sd_prior, 1)
  expect_error(bpnmf_rank_shrinkage_opts(group_sd_prior = 0), "not >= ")
  # The YAML round trip lives in test-config.R, beside the helpers it needs.
})

test_that("rank 1 has no curves: mu is the time level, unit level and offset", {
  skip_on_cran()
  skip_if_no_cmdstan()
  sd <- stan_data_joint(sc_data(), rank = 1, gen_ypred = FALSE)
  fit <- sc_fit(sd)
  # At rank 1 the lone curve would be exactly absorbed by time_fe, so the
  # model no longer carries it (a zero-size variable has no sizes entry).
  expect_null(fit$metadata()$stan_variable_sizes$time_fac)
  expect_lt(sc_reconstruction_error(fit, sd), SC_TOL)
})

test_that("shared curves: one curve set, and the likelihood uses it", {
  skip_on_cran()
  skip_if_no_cmdstan()
  sd <- stan_data_joint(sc_data(), rank = 3, gen_ypred = FALSE, shared_curves = TRUE)
  fit <- sc_fit(sd)
  expect_equal(fit$metadata()$stan_variable_sizes$time_fac, c(1, sd$N, 3))
  expect_lt(
    sc_reconstruction_error(fit, sd, weights = "unit_weight", shared = TRUE),
    SC_TOL
  )
})

test_that("shared curves with rank shrinkage: crossed effects around one profile", {
  skip_on_cran()
  skip_if_no_cmdstan()
  sd <- stan_data_joint(
    sc_data(), rank = 4, gen_ypred = FALSE, shared_curves = TRUE,
    rank_shrinkage = bpnmf_rank_shrinkage_opts()
  )
  fit <- sc_fit(sd)
  sizes <- fit$metadata()$stan_variable_sizes
  expect_equal(sizes$time_fac, c(1, sd$N, 4))
  expect_equal(sizes$stick, c(1, 3))             # one global stick set
  expect_equal(sizes$group_profile_z, c(4, sd$K))
  expect_equal(sizes$unit_shared_z, c(4, sd$D))

  w <- sc_draws(fit, "unit_weight_fitted")
  expect_lt(max(abs(apply(w, c(1, 2, 3), sum) - 1)), SC_TOL)
  expect_lt(
    sc_reconstruction_error(fit, sd, weights = "unit_weight_fitted", shared = TRUE),
    SC_TOL
  )

  # Each group's profile is the global one perturbed on the log-ratio scale.
  lg <- sc_draws(fit, "log_global_weight")
  lk <- sc_draws(fit, "log_group_weight")
  gsd <- sc_draws(fit, "group_profile_sd")
  gz <- sc_draws(fit, "group_profile_z")
  err <- 0
  for (i in seq_len(dim(lg)[1])) {
    for (k in seq_len(sd$K)) {
      a <- lg[i, 1, ] + gsd[i, 1] * gz[i, , k]
      err <- max(err, abs(a - log(sum(exp(a))) - lk[i, k, ]))
    }
  }
  expect_lt(err, SC_TOL)

  # The report gains a global row ahead of the groups.
  obj <- new_bpnmf_class(
    list(fit = fit, stan_data = sd, data = sc_data(), rank = 4L), "bpnmf_fit"
  )
  eff <- bpnmf_eff_rank_summary(obj)
  expect_equal(eff$group, c("(global)", "g1", "g2"))
  expect_true(all(eff$median >= 1 & eff$median <= 4))
})
