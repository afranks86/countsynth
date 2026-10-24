# The strongest cross-implementation check: evaluate the Stan models'
# log-density at constrained parameter values exported from Python (with
# numpyro's log_density at the same values) and require the two to agree up
# to an additive constant per configuration. Pairwise differencing cancels
# normalizing constants dropped by Stan's `~` statements and the Dirichlet
# uniform constant; anything data- or parameter-dependent must match exactly.
# The fixture keeps the Python side's original `state_*` field names, and so
# does the per-case list handed to py_to_stan_params(); its left-hand names
# are this package's Stan variables, which were renamed to `unit_*`. Both
# sides of `unit_fe_mu = p$state_fe_mu` are correct.
#
# One term is deliberately not shared. With the censoring adjustment on, the
# Python model also adds log P(y not in 1..9) for every observed cell, which
# counts the fact that the cell was not suppressed twice (P(y) already
# implies it); the Stan models dropped it on 2026-10-23. The comparison adds
# that term back on the Stan side, computed in R at the same parameters, so
# everything else still has to match the Python density exactly.
#
# Densities come from CmdStan's standalone `log_prob` method and the cell
# means from a one-draw `fixed_param` run, both on the plainly compiled
# package model. cmdstanr's model methods (`compile_model_methods = TRUE`)
# need an Rcpp build that some toolchains cannot do, and they are not needed.

parity_fixture <- function() {
  fx <- jsonlite::fromJSON(fixture_path("logdensity_parity.json"),
    simplifyVector = FALSE
  )
  d <- fx$data
  dims <- c(d$K, d$D, d$N)
  data <- structure(
    list(
      Y = parity_arr(d$Y, dims),
      denominators = parity_arr(d$denominators, dims),
      control_idx_array = parity_arr(d$control, dims) == 1,
      missing_idx_array = parity_arr(d$missing, dims) == 1,
      groups = paste0("g", seq_len(d$K)), units = LETTERS[seq_len(d$D)],
      times = seq(as.Date("2020-01-01"), by = "month", length.out = d$N),
      df = NULL, type = "test"
    ),
    class = c("countsynth_data", "list")
  )
  list(fx = fx, d = d, data = data)
}

parity_arr <- function(x, dims) {
  array(unlist(x), dim = rev(dims)) |> aperm(rev(seq_along(dims)))
}

py_to_stan_params <- function(p, treated, sample_disp) {
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

# Stan's log density (no Jacobian) at constrained `params`.
stan_log_prob <- function(model, data, params) {
  files <- withr::local_tempfile(fileext = c(".json", ".json", ".csv"))
  cmdstanr::write_stan_json(data, files[1])
  cmdstanr::write_stan_json(params, files[2])
  status <- system2(
    model$exe_file(),
    c(
      "log_prob", paste0("constrained_params=", files[2]), "jacobian=0",
      "data", paste0("file=", files[1]),
      "output", paste0("file=", files[3]), "sig_figs=15"
    ),
    stdout = FALSE, stderr = FALSE
  )
  if (status != 0) {
    stop("CmdStan log_prob failed")
  }
  utils::read.csv(files[3], comment.char = "#")$lp__
}

# Model quantities at constrained `params`: a fixed_param draw initialized
# there reports the parameters' own transformed parameters.
stan_quantities_at <- function(model, data, params, vars) {
  # Zero-size parameters (the unused te-regression block) need no init.
  withr::local_options(cmdstanr_warn_inits = FALSE)
  fit <- model$sample(
    data = data, init = list(params), fixed_param = TRUE, chains = 1,
    iter_sampling = 1, refresh = 0, show_messages = FALSE, sig_figs = 15
  )
  stats::setNames(lapply(vars, function(v) {
    as.numeric(posterior::as_draws_matrix(fit$draws(variables = v)))
  }), vars)
}

# The term the Stan models no longer have: log P(y not in 1..9) summed over
# the observed cells at log-means `eta` (NB2 with concentration `phi`, or
# Poisson when `phi` is NULL).
not_suppressed_term <- function(eta, phi = NULL) {
  mass <- vapply(seq_along(eta), function(i) {
    mu <- exp(eta[i])
    sum(if (is.null(phi)) {
      stats::dpois(1:9, mu)
    } else {
      stats::dnbinom(1:9, size = phi[i], mu = mu)
    })
  }, numeric(1))
  sum(log1p(-mass))
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
  pf <- parity_fixture()
  d <- pf$d
  K <- d$K
  D <- d$D
  N <- d$N
  model <- countsynth_stan_model("joint")

  joint_cases <- Filter(function(cs) cs$model == "joint", pf$fx$cases)
  configs <- unique(vapply(joint_cases, `[[`, character(1), "config"))

  for (cfg_name in configs) {
    group <- Filter(function(cs) cs$config == cfg_name, joint_cases)
    first <- group[[1]]
    sd <- stan_data_joint(
      pf$data,
      rank = d$rank,
      model_treated = first$model_treated,
      outcome_distribution = first$outcome_dist,
      nb_disp = 1e-4,
      time_level = "uncentered",
      sample_disp = first$sample_disp,
      adjust_for_missingness = first$adjust,
      gen_ypred = FALSE
    )

    lp_stan <- numeric(length(group))
    lp_py <- numeric(length(group))
    for (i in seq_along(group)) {
      cs <- group[[i]]
      p <- cs$params
      # Still in the fixture's (Python) names and shapes; py_to_stan_params()
      # renames and reorders them.
      pj <- list(
        time_fac = parity_arr(p$time_fac, c(N, d$rank, K)),
        state_fe_mu = unlist(p$state_fe_mu),
        state_fe_sigma = unlist(p$state_fe_sigma),
        state_fe_z = parity_arr(p$state_fe_z, c(D, K)),
        time_fe = parity_arr(p$time_fe, c(N, K)),
        unit_weight = parity_arr(p$unit_weight, c(D, K, d$rank))
      )
      if (cs$model_treated) {
        pj <- c(pj, list(
          treatment_it_scale = p$treatment_it_scale,
          treatment_state_scale = p$treatment_state_scale,
          treatment_category_scale = p$treatment_category_scale,
          state_category_scale = p$state_category_scale,
          treatment_kt_z = unlist(p$treatment_kt_z),
          state_treatment_effect_z = unlist(p$state_treatment_effect_z),
          state_category_te_z = parity_arr(p$state_category_te_z, c(K, D)),
          category_treatment_effect = unlist(p$category_treatment_effect)
        ))
      }
      if (cs$sample_disp) {
        pj$disp <- unlist(p$disp)
      }
      stan_pars <- py_to_stan_params(pj, cs$model_treated, cs$sample_disp)
      lp_stan[i] <- stan_log_prob(model, sd, stan_pars)
      if (cs$adjust) {
        q <- stan_quantities_at(
          model, sd, stan_pars, c("mu_ctrl", if (cs$model_treated) "te")
        )
        eta <- q$mu_ctrl
        if (cs$model_treated) {
          eta[sd$exp_cell] <- eta[sd$exp_cell] + q$te
        }
        units <- sd$cell_unit[sd$obs_cell]
        phi <- if (cs$outcome_dist != "NB") {
          NULL
        } else if (cs$sample_disp) {
          1 / pj$disp[units]
        } else {
          rep(1 / 1e-4, length(units))
        }
        lp_stan[i] <- lp_stan[i] + not_suppressed_term(eta[sd$obs_cell], phi)
      }
      lp_py[i] <- cs$log_density
    }
    expect_constant_offset(lp_stan, lp_py, cfg_name)
  }
})

test_that("Stan stage-2 model log-density matches numpyro up to a constant", {
  skip_on_cran()
  skip_if_no_cmdstan()
  pf <- parity_fixture()
  d <- pf$d
  K <- d$K
  D <- d$D
  N <- d$N
  mu_flat <- flatten_kdn(parity_arr(d$mu_ctrl_fixed, c(K, D, N)))
  phi_fixed <- unlist(d$phi_fixed)
  sd2 <- stan_data_stage2(
    pf$data, mu_flat,
    phi_unit = phi_fixed,
    outcome_distribution = "NB", adjust_for_missingness = TRUE
  )
  model <- countsynth_stan_model("cut_stage2")

  s2 <- Filter(function(cs) cs$model == "stage2", pf$fx$cases)
  lp_stan <- numeric(length(s2))
  lp_py <- numeric(length(s2))
  for (i in seq_along(s2)) {
    p <- s2[[i]]$params
    stan_pars <- list(
      treatment_it_scale = p$treatment_it_scale,
      treatment_unit_scale = as.array(p$treatment_state_scale),
      treatment_category_scale = as.array(p$treatment_category_scale),
      unit_category_scale = as.array(p$state_category_scale),
      treatment_kt_z = unlist(p$treatment_kt_z),
      unit_treatment_effect_z = unlist(p$state_treatment_effect_z),
      unit_category_te_z = parity_arr(p$state_category_te_z, c(K, D)),
      category_treatment_effect = unlist(p$category_treatment_effect)
    )
    lp_stan[i] <- stan_log_prob(model, sd2, stan_pars)
    eta <- stan_quantities_at(model, sd2, stan_pars, "mu_exposed")$mu_exposed
    lp_stan[i] <- lp_stan[i] + not_suppressed_term(
      eta[sd2$obs_e], phi_fixed[sd2$exp_d[sd2$obs_e]]
    )
    lp_py[i] <- s2[[i]]$log_density
  }
  expect_constant_offset(lp_stan, lp_py, "stage2_nb_adjust")
})
