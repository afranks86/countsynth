# The strongest cross-implementation check: evaluate the Stan models'
# log-density at constrained parameter values exported from Python (with
# numpyro's log_density at the same values) and require the two to agree up
# to an additive constant per configuration. Pairwise differencing cancels
# normalizing constants dropped by Stan's `~` statements and the Dirichlet
# uniform constant; anything data- or parameter-dependent must match exactly.
# The fixture keeps the Python side's original `state_*` field names; the
# left-hand names here are this package's Stan variables, which were renamed
# to `unit_*`. Both sides of `unit_fe_mu = p$state_fe_mu` are correct.

py_to_stan_params <- function(p, treated, sample_disp) {
  N <- dim(p$time_fac)[1] %||% length(p$time_fac)
  out <- list(
    # numpyro (N, R, K) -> Stan array[K] matrix[N, R]
    time_fac = aperm(p$time_fac, c(3, 1, 2)),
    unit_fe_mu = p$state_fe_mu,
    unit_fe_sigma = p$state_fe_sigma,
    unit_fe_z = p$state_fe_z,
    # The Python model has the uncentered Gamma time level, so parity is
    # checked against time_level = "uncentered", whose parameter this is.
    time_fe_free = p$time_fe,
    # numpyro (D, K, R) -> Stan array[K, D] simplex[R]
    unit_weight = aperm(p$unit_weight, c(2, 1, 3))
  )
  if (treated) {
    out <- c(out, list(
      treatment_it_scale = as.array(p$treatment_it_scale),
      treatment_unit_scale = as.array(p$treatment_state_scale),
      treatment_category_scale = as.array(p$treatment_category_scale),
      unit_category_scale = as.array(p$state_category_scale),
      treatment_kt_z = p$treatment_kt_z,
      unit_treatment_effect_z = p$state_treatment_effect_z,
      unit_category_te_z = p$state_category_te_z,
      category_treatment_effect = p$category_treatment_effect
    ))
  }
  if (sample_disp) {
    out$disp <- p$disp
  }
  out
}

# tol is absolute spread across cases. Log-densities here are O(1e3-1e4);
# the lgamma-heavy NB terms accumulate ~1e-8 relative float noise between
# the two implementations, so 1e-4 absolute is tight enough to catch any
# real mis-port (which shows up at O(0.1) or more) without false alarms.
expect_constant_offset <- function(lp_stan, lp_py, label, tol = 1e-4) {
  offsets <- lp_stan - lp_py
  spread <- max(offsets) - min(offsets)
  expect_lt(spread, tol, label = sprintf("%s: offset spread %.2e", label, spread))
}

test_that("Stan joint model log-density matches numpyro up to a constant", {
  skip_on_cran()
  skip_if_no_cmdstan()
  fx <- jsonlite::fromJSON(fixture_path("logdensity_parity.json"),
    simplifyVector = FALSE
  )
  d <- fx$data
  K <- d$K
  D <- d$D
  N <- d$N
  as_arr <- function(x, dims) {
    array(unlist(x), dim = rev(dims)) |> aperm(rev(seq_along(dims)))
  }
  data <- structure(
    list(
      Y = as_arr(d$Y, c(K, D, N)),
      denominators = as_arr(d$denominators, c(K, D, N)),
      control_idx_array = as_arr(d$control, c(K, D, N)) == 1,
      missing_idx_array = as_arr(d$missing, c(K, D, N)) == 1,
      groups = paste0("g", seq_len(K)),
      units = LETTERS[seq_len(D)],
      times = seq(as.Date("2020-01-01"), by = "month", length.out = N),
      df = NULL, type = "test"
    ),
    class = c("countsynth_data", "list")
  )

  # Compile into a per-session dir: a pre-existing cached executable cannot
  # be reused with model methods (log_prob / unconstrain_variables).
  model <- cmdstanr::cmdstan_model(
    stan_file_path("joint"),
    include_paths = system.file("stan", "include", package = "countsynth"),
    dir = withr::local_tempdir(),
    compile_model_methods = TRUE, quiet = TRUE
  )

  cases <- fx$cases
  joint_cases <- Filter(function(cs) cs$model == "joint", cases)
  configs <- unique(vapply(joint_cases, `[[`, character(1), "config"))

  for (cfg_name in configs) {
    group <- Filter(function(cs) cs$config == cfg_name, joint_cases)
    first <- group[[1]]
    sd <- stan_data_joint(
      data,
      rank = d$rank,
      model_treated = first$model_treated,
      outcome_distribution = first$outcome_dist,
      nb_disp = 1e-4,
      time_level = "uncentered",
      sample_disp = first$sample_disp,
      adjust_for_missingness = first$adjust,
      gen_ypred = FALSE
    )
    fit <- model$sample(
      data = sd, chains = 1, iter_warmup = 30, iter_sampling = 5,
      refresh = 0, show_messages = FALSE, show_exceptions = FALSE, seed = 1
    )
    fit$init_model_methods(verbose = FALSE)

    lp_stan <- numeric(length(group))
    lp_py <- numeric(length(group))
    for (i in seq_along(group)) {
      cs <- group[[i]]
      p <- cs$params
      pj <- list(
        time_fac = as_arr(p$time_fac, c(N, d$rank, K)),
        unit_fe_mu = unlist(p$state_fe_mu),
        unit_fe_sigma = unlist(p$state_fe_sigma),
        unit_fe_z = as_arr(p$state_fe_z, c(D, K)),
        time_fe = as_arr(p$time_fe, c(N, K)),
        unit_weight = as_arr(p$unit_weight, c(D, K, d$rank))
      )
      if (cs$model_treated) {
        pj <- c(pj, list(
          treatment_it_scale = p$treatment_it_scale,
          treatment_unit_scale = p$treatment_state_scale,
          treatment_category_scale = p$treatment_category_scale,
          unit_category_scale = p$state_category_scale,
          treatment_kt_z = unlist(p$treatment_kt_z),
          unit_treatment_effect_z = unlist(p$state_treatment_effect_z),
          unit_category_te_z = as_arr(p$state_category_te_z, c(K, D)),
          category_treatment_effect = unlist(p$category_treatment_effect)
        ))
      }
      if (cs$sample_disp) {
        pj$disp <- unlist(p$disp)
      }
      stan_pars <- py_to_stan_params(pj, cs$model_treated, cs$sample_disp)
      upars <- fit$unconstrain_variables(stan_pars)
      lp_stan[i] <- fit$log_prob(upars, jacobian = FALSE)
      lp_py[i] <- cs$log_density
    }
    expect_constant_offset(lp_stan, lp_py, cfg_name)
  }
})

test_that("Stan stage-2 model log-density matches numpyro up to a constant", {
  skip_on_cran()
  skip_if_no_cmdstan()
  fx <- jsonlite::fromJSON(fixture_path("logdensity_parity.json"),
    simplifyVector = FALSE
  )
  d <- fx$data
  K <- d$K
  D <- d$D
  N <- d$N
  as_arr <- function(x, dims) {
    array(unlist(x), dim = rev(dims)) |> aperm(rev(seq_along(dims)))
  }
  data <- structure(
    list(
      Y = as_arr(d$Y, c(K, D, N)),
      denominators = as_arr(d$denominators, c(K, D, N)),
      control_idx_array = as_arr(d$control, c(K, D, N)) == 1,
      missing_idx_array = as_arr(d$missing, c(K, D, N)) == 1,
      groups = paste0("g", seq_len(K)), units = LETTERS[seq_len(D)],
      times = seq(as.Date("2020-01-01"), by = "month", length.out = N),
      df = NULL, type = "test"
    ),
    class = c("countsynth_data", "list")
  )
  mu_flat <- flatten_kdn(as_arr(d$mu_ctrl_fixed, c(K, D, N)))
  sd2 <- stan_data_stage2(
    data, mu_flat,
    phi_unit = unlist(d$phi_fixed),
    outcome_distribution = "NB", adjust_for_missingness = TRUE
  )
  model <- cmdstanr::cmdstan_model(
    stan_file_path("cut_stage2"),
    include_paths = system.file("stan", "include", package = "countsynth"),
    dir = withr::local_tempdir(),
    compile_model_methods = TRUE, quiet = TRUE
  )
  fit <- model$sample(
    data = sd2, chains = 1, iter_warmup = 30, iter_sampling = 5,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE, seed = 1
  )
  fit$init_model_methods(verbose = FALSE)

  s2 <- Filter(function(cs) cs$model == "stage2", fx$cases)
  lp_stan <- numeric(length(s2))
  lp_py <- numeric(length(s2))
  for (i in seq_along(s2)) {
    p <- s2[[i]]$params
    stan_pars <- list(
      treatment_it_scale = p$treatment_it_scale,
      treatment_unit_scale = p$treatment_state_scale,
      treatment_category_scale = p$treatment_category_scale,
      unit_category_scale = p$state_category_scale,
      treatment_kt_z = unlist(p$treatment_kt_z),
      unit_treatment_effect_z = unlist(p$state_treatment_effect_z),
      unit_category_te_z = as_arr(p$state_category_te_z, c(K, D)),
      category_treatment_effect = unlist(p$category_treatment_effect)
    )
    upars <- fit$unconstrain_variables(stan_pars)
    lp_stan[i] <- fit$log_prob(upars, jacobian = FALSE)
    lp_py[i] <- s2[[i]]$log_density
  }
  expect_constant_offset(lp_stan, lp_py, "stage2_nb_adjust")
})
