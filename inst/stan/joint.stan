// Joint Bayesian panel NMF model (port of models/joint.py). With
// model_treated=0 this is also the cut stage-1 baseline model
// (models/cut_baseline.py is scientifically identical to that branch; here a
// data flag replaces the Python file duplication).
//
// Cells of the (K groups, D units, N times) panel are addressed by 1-based
// row-major flat indices, cell = ((k-1)*D + (d-1))*N + n, matching numpy's
// reshape(-1); all index vectors are precomputed in R (stan-data.R).
functions {
  #include censoring.stanfunctions
  // Log weights of a truncated stick-breaking draw, the final stick taking
  // the whole remainder so the weights sum to one. Built in logs because a
  // component nothing uses can underflow to exactly zero, and log(0) would
  // hand autodiff an infinite gradient.
  vector log_stick_breaking(vector nu) {
    int R = rows(nu) + 1;
    vector[R] out;
    real log_rem = 0;
    for (r in 1 : (R - 1)) {
      out[r] = log_rem + log(nu[r]);
      log_rem += log1m(nu[r]);
    }
    out[R] = log_rem;
    return out;
  }
}
data {
  int<lower=1> K;                        // groups
  int<lower=1> D;                        // units
  int<lower=1> N;                        // time periods
  int<lower=1> R;                        // factorization rank
  int<lower=1> KDN;                      // K * D * N
  vector[KDN] log_denom;                 // log scaled exposure denominators
  array[KDN] int<lower=1, upper=D> cell_unit; // unit index of each flat cell

  // Direct-likelihood cells: !missing & (model_treated | control).
  int<lower=0> n_obs;
  array[n_obs] int<lower=1, upper=KDN> obs_cell;
  array[n_obs] int<lower=0> y;
  array[n_obs] int<lower=1, upper=D> obs_unit;

  // Exposed (treated) cells in canonical flat order; defines treatment_kt_z.
  int<lower=0> n_exposed;
  array[n_exposed] int<lower=1, upper=K> exp_k;
  array[n_exposed] int<lower=1, upper=D> exp_d;
  array[n_exposed] int<lower=1, upper=K * D> exp_kd; // (d-1)*K + k gather index
  array[n_exposed] int<lower=1, upper=KDN> exp_cell;

  // Censoring adjustment sets (all cells if model_treated, else control).
  int<lower=0> n_cens;                   // suppressed: add log P(y in 1..9)
  array[n_cens] int<lower=1, upper=KDN> cens_cell;
  array[n_cens] int<lower=1, upper=D> cens_unit;
  int<lower=0> n_notcens;                // observed: add log P(y not in 1..9)
  array[n_notcens] int<lower=1, upper=KDN> notcens_cell;
  array[n_notcens] int<lower=1, upper=D> notcens_unit;

  int<lower=0, upper=1> model_treated;
  int<lower=0, upper=1> is_nb;           // 1 = NB2, 0 = Poisson
  int<lower=0, upper=1> sample_disp;
  int<lower=0, upper=1> adjust_missing;
  real<lower=0> nb_disp;                 // fixed dispersion (phi = 1/nb_disp)
  int<lower=0, upper=1> gen_ypred;       // emit counterfactual predictive
  // Gamma(shape, shape) on the temporal basis: mean 1, and sd(log) ~
  // 1/sqrt(shape), so this sets how far the low-rank temporal structure may
  // swing multiplicatively away from a unit's own level. 20 is the historical
  // default (about a +/-25% swing).
  real<lower=0> time_fac_shape;
  // Gamma(shape, shape) on the common time level. 1 is the historical
  // default (Gamma(1, 1)) and is very diffuse: sd(log) = 1.28. Since
  // unit_fe_mu is flat and the likelihood sees only log(time_fe) +
  // unit_fe_mu, this prior alone separates the two.
  real<lower=0> time_fe_shape;

  // Optional rank shrinkage (rank_shrink == 1): the per-unit component
  // weights stop being iid uniform on the simplex and become multiplicative
  // perturbations of a shared per-group component-popularity profile, itself
  // a stick-breaking draw truncated at R sticks -- a two-level HDP over the
  // R candidate temporal curves, with a logistic-normal lower level (see the
  // parameters block for why not the textbook Dirichlet). Every parameter it adds is zero-size when
  // rank_shrink == 0, so the unconstrained vector -- and therefore any draw
  // from a given seed -- is bit-identical to the pre-shrinkage model.
  //
  // Why it matters beyond rank selection: under Dirichlet(1,...,1),
  // E[sum_r unit_weight^2] = 2 / (R + 1), so the prior variance of the whole
  // low-rank term, v * sum_r w_r^2, falls off like 1/R. `time_fac_shape`
  // therefore does not mean what it says unless the rank is held fixed. A
  // stick-breaking profile has E[sum_r profile^2] = 1 / (1 + mass) with no R
  // in it at all, and the per-unit weights inherit that.
  int<lower=0, upper=1> rank_shrink;
  real<lower=0> group_mass_shape;        // Gamma prior on group_weight_mass
  real<lower=0> group_mass_rate;
  real<lower=0> unit_sd_scale;           // HalfNormal prior on unit_weight_sd

  // Shared temporal curves (share_fac == 1; R side sets it only for K >= 2
  // and R >= 2): one set of R curves for every group instead of one set per
  // group, so component r is the same temporal pattern in each group and its
  // loadings are comparable across groups. With rank shrinkage the weights
  // then get crossed effects on the log-ratio scale around one global
  // profile: group, unit (shared across groups), and group x unit -- the
  // (1 | group) + (1 | unit) + (1 | group:unit) structure the legacy
  // treatment effect uses, applied to the loadings.
  int<lower=0, upper=1> share_fac;
  real<lower=0> group_sd_scale;          // HalfNormal prior on group_profile_sd
  real<lower=0> shared_unit_sd_scale;    // HalfNormal prior on unit_shared_sd

  // Optional treatment-effect regression (te_reg == 1, requires
  // model_treated == 1): the legacy group/unit/group:unit hierarchy is
  // replaced by a linear fixed-effect surface X * te_beta plus J ragged
  // random-effect terms whose per-level coefficients shrink toward that
  // surface. The iid exposed-cell effect (treatment_kt_z) is kept in both
  // parameterizations. All sizes below are zero when te_reg == 0.
  int<lower=0, upper=1> te_reg;
  int<lower=0> P;                        // fixed-effect design columns
  matrix[n_exposed, P] X;
  int<lower=0> J;                        // random-effect terms
  int<lower=0> Qtot;                     // total term predictors, sum(Q)
  int<lower=0> Utot;                     // total coefficients, sum(L .* Q)
  array[J] int<lower=1> L;               // levels per term
  array[J] int<lower=1> Q;               // predictors per term
  matrix[n_exposed, Qtot] Z;             // column-concatenated term designs
  array[J, n_exposed] int<lower=1> re_level; // level of each exposed cell
  vector<lower=0>[P] te_beta_prior_scale;
  vector<lower=0>[Qtot] te_re_prior_scale;
}
transformed data {
  int DN = D * N;
  // Legacy treatment hierarchy is active only without the regression design.
  int te_legacy = (model_treated == 1 && te_reg == 0) ? 1 : 0;
  vector[9] sup;                         // censored-count support 1..9
  vector[9] lgamma_sup1;
  for (v in 1 : 9) {
    sup[v] = v;
    lgamma_sup1[v] = lgamma(v + 1);
  }
  vector[9] neg_lgamma_sup1 = -lgamma_sup1;
  // Fixed-dispersion fast path: when sample_disp == 0 the concentration is
  // data, so the likelihood's phi and the censoring coefficients carry no
  // autodiff cost at all.
  real phi_fixed = inv(nb_disp);
  real log_phi_fixed = log(phi_fixed);
  vector[9] b_fixed = nb_censor_coeff(phi_fixed, sup, lgamma_sup1);

  // Curves and weights exist only for R >= 2. At R = 1 every weight is 1 and
  // the lone curve multiplies every unit's rate alike, so the likelihood sees
  // only its product with time_fe and it is exactly absorbed: the model is
  // the rank-1 surface time_fe[n, k] * exp(unit_fe[d, k]) either way, and
  // dropping the curve removes a factor nothing but the priors could split.
  // More generally R counts the rank of the multiplicative rate surface:
  // time_fe supplies one direction and R curves on a simplex add R - 1.
  int R_fac = R >= 2 ? R : 0;
  int K_fac = share_fac == 1 ? 1 : K;    // curve sets
  int K_stick = rank_shrink == 1 ? K_fac : 0;   // one global or one per group
  int crossed = (rank_shrink == 1 && share_fac == 1) ? 1 : 0;
}
parameters {
  // One curve set per group, or one shared set; none at R = 1.
  array[K_fac] matrix<lower=0>[N, R_fac] time_fac; // Gamma(shape, shape)
  vector[K] unit_fe_mu;                     // improper flat (no statement)
  vector<lower=0>[K] unit_fe_sigma;         // HalfNormal(0.5)
  matrix[D, K] unit_fe_z;                   // std normal (non-centered)
  matrix<lower=0>[N, K] time_fe;             // Gamma(shape, shape), logged
  // Flat weights, Dirichlet(1,...,1) per (k, d). Rank shrinkage replaces
  // them wholesale (zero-size here, weights built in transformed parameters
  // from the blocks below), so with rank_shrink == 0 and R >= 2 this
  // declaration -- name, size and position -- is exactly the pre-shrinkage
  // model's. At R = 1 a simplex[1] has no free coordinates anyway.
  array[(rank_shrink == 1 || R_fac == 0) ? 0 : K,
        (rank_shrink == 1 || R_fac == 0) ? 0 : D] simplex[R] unit_weight;

  // Rank shrinkage, all zero-size when rank_shrink == 0.
  //
  // Stick fractions rather than a sparse symmetric Dirichlet on the simplex.
  // The two are close in distribution (both approximate a DP truncated at R
  // atoms) and far apart in geometry, so they are not interchangeable: on
  // the simplex a near-zero profile weight demands a coordinate pinned to
  // the boundary, whose transformed tail decays only like exp(alpha * z) and
  // so runs out to z ~ -1/alpha, and NUTS diverges getting there. Here the
  // same near-zero weight is a *product* of moderate stick fractions -- an
  // unused eighth component needs no extreme coordinate anywhere.
  array[K_stick] vector<lower=0, upper=1>[R - 1] stick;
  array[rank_shrink] real<lower=0> group_weight_mass; // DP mass, Beta(1, mass)
  // Per-unit deviation from the group profile: logistic-normal rather than
  // Dirichlet, and non-centered. unit_weight_sd is the multiplicative spread
  // of a unit's loadings about its group's profile (an sd on the log scale),
  // the same knob a Dirichlet concentration c would set at about
  // sqrt(2 / c).
  //
  // Dirichlet(unit_conc * profile) is the textbook HDP form and was measured
  // at 5.4x the wall time of this on the 51-unit fertility panel (1193s vs
  // 223s for 300 iterations at rank 8) with the same leapfrog count and no
  // divergences either way -- so the price is per gradient, not geometry: a
  // Dirichlet whose concentration is itself sampled needs an lgamma and a
  // digamma per component per unit, 408 of each per gradient there, where
  // the softmax below needs an exp. Non-centering is the second reason and
  // the one that should scale: here the deviation's prior is std_normal
  // whatever the profile does, so an unused component's z is merely
  // unidentified -- flat likelihood, proper prior, no boundary -- while
  // under the Dirichlet its width is set by profile[r], which is sampled.
  //
  // It also drops a fudge factor the Dirichlet form needs. There, unused
  // components sit at concentration ~0, where the simplex transform funnels,
  // so a floor is required; but a per-component floor puts back weight that
  // grows with R (breaking the invariance this exists to buy) and a floor
  // spread as total/R shrinks back toward the boundary. No Dirichlet, no
  // floor, no tension.
  array[rank_shrink] real<lower=0> unit_weight_sd;
  // Columns are (k, d) pairs at (k - 1) * D + d, the same group-major order
  // the flat cell index uses; R rows so each unit's deviation is one
  // contiguous column.
  matrix[rank_shrink == 1 ? R : 0, rank_shrink == 1 ? K * D : 0] unit_weight_z;
  // Crossed effects, zero-size unless curves are shared under shrinkage.
  // Group k's profile departs from the global one by group_profile_sd *
  // group_profile_z[, k]; unit d departs the same way in every group by
  // unit_shared_sd * unit_shared_z[, d]; unit_weight_sd * unit_weight_z is
  // then the group x unit remainder. Each has the softmax's flat shift
  // direction, proper under std_normal and harmless for the same reason.
  array[crossed] real<lower=0> group_profile_sd;
  matrix[crossed == 1 ? R : 0, crossed == 1 ? K : 0] group_profile_z;
  array[crossed] real<lower=0> unit_shared_sd;
  matrix[crossed == 1 ? R : 0, crossed == 1 ? D : 0] unit_shared_z;

  // Treatment block; zero-size when model_treated == 0. Declaration order
  // matches the pre-regression model, and the legacy hierarchy is zero-size
  // exactly when the regression replaces it, so the legacy unconstrained
  // vector is unchanged -- same seed, same draws.
  array[model_treated] real<lower=0> treatment_it_scale;      // HalfNormal(0.1)
  array[te_legacy] real<lower=0> treatment_unit_scale;       // HalfNormal(1)
  array[te_legacy] real<lower=0> treatment_category_scale;    // HalfNormal(1)
  array[te_legacy] real<lower=0> unit_category_scale;        // HalfNormal(1)
  vector[model_treated == 1 ? n_exposed : 0] treatment_kt_z;
  vector[te_legacy == 1 ? D : 0] unit_treatment_effect_z;
  matrix[te_legacy == 1 ? K : 0, te_legacy == 1 ? D : 0] unit_category_te_z;
  vector[te_legacy == 1 ? K : 0] category_treatment_effect;   // centered

  // Per-unit dispersion; Uniform(0,1) via constraint + factor prior below.
  vector<lower=0, upper=1>[sample_disp == 1 ? D : 0] disp;

  // Regression surface; all sizes are zero when te_reg == 0.
  vector[P] te_beta;
  vector<lower=0>[Qtot] te_re_scale;     // HalfNormal(te_re_prior_scale)
  vector[Utot] te_re_z;                  // non-centered level coefficients
}
transformed parameters {
  vector[KDN] mu_ctrl;                   // untreated log-count surface
  // Component-popularity profile per group, in logs (generated quantities
  // reports the natural scale). Without shared curves each group has its own
  // stick-breaking draw; with them there is one global draw, and each group's
  // profile is a log-ratio deviation from it -- so a component the global
  // profile has zeroed out stays near zero in every group.
  array[rank_shrink == 1 ? K : 0] vector[R] log_group_weight;
  array[crossed] vector[R] log_global_weight;
  // The realized per-unit weights under rank shrinkage, where the sampled
  // simplex is empty: each unit's group profile perturbed multiplicatively
  // and renormalized, so every column still sums to one and the level stays
  // with unit_fe exactly as in the flat model. Computed once here, for the
  // likelihood and the output alike.
  array[rank_shrink == 1 ? K : 0, rank_shrink == 1 ? D : 0]
    vector[R] unit_weight_fitted;
  vector[model_treated == 1 ? n_exposed : 0] te;
  vector[is_nb == 1 ? D : 0] phi_unit;   // NB2 concentration per unit

  if (crossed == 1) {
    log_global_weight[1] = log_stick_breaking(stick[1]);
    for (k in 1 : K) {
      log_group_weight[k] = log_softmax(
        log_global_weight[1] + group_profile_sd[1] * col(group_profile_z, k)
      );
    }
  } else {
    for (k in 1 : K_stick) {
      log_group_weight[k] = log_stick_breaking(stick[k]);
    }
  }
  for (k in 1 : (rank_shrink == 1 ? K : 0)) {
    for (d in 1 : D) {
      vector[R] dev = unit_weight_sd[1] * col(unit_weight_z, (k - 1) * D + d);
      if (crossed == 1) {
        dev += unit_shared_sd[1] * col(unit_shared_z, d);
      }
      unit_weight_fitted[k, d] = softmax(log_group_weight[k] + dev);
    }
  }

  for (k in 1 : K) {
    row_vector[D] fe_k = unit_fe_mu[k]
                         + unit_fe_sigma[k] * to_row_vector(unit_fe_z[, k]);
    int base = (k - 1) * DN;
    // to_vector is column-major: column d contributes N contiguous entries,
    // exactly the row-major (k, d, n) flat layout of this group's cells.
    if (R_fac > 0) {
      // Natural-scale low-rank product: log(TF_k %*% W_k) equals the Python
      // log-sum-exp assembly of log time_fac + log weights; both factors are
      // positive and O(1), so the direct product is stable and faster.
      matrix[R, D] W;
      if (rank_shrink == 1) {
        for (d in 1 : D) {
          W[, d] = unit_weight_fitted[k, d];
        }
      } else {
        for (d in 1 : D) {
          W[, d] = unit_weight[k, d];
        }
      }
      matrix[N, D] factor_kd = time_fac[share_fac == 1 ? 1 : k] * W;
      // Summation order kept exactly as before shared curves existed, so a
      // flat-weight configuration reproduces its draws bit for bit. (A
      // rank_shrinkage one reproduces the same log density but not the same
      // bits: restructuring how its weights are built changed the autodiff
      // graph, and a last-ulp gradient difference eventually flips a
      // trajectory decision -- verified identical lp__ for 12 iterations.)
      mu_ctrl[(base + 1) : (base + DN)]
        = to_vector(log(factor_kd) + rep_matrix(log(time_fe[, k]), D)
                    + rep_matrix(fe_k, N))
          + log_denom[(base + 1) : (base + DN)];
    } else {
      mu_ctrl[(base + 1) : (base + DN)]
        = to_vector(rep_matrix(log(time_fe[, k]), D) + rep_matrix(fe_k, N))
          + log_denom[(base + 1) : (base + DN)];
    }
  }

  if (model_treated == 1) {
    te = treatment_kt_z * treatment_it_scale[1];
    if (te_reg == 0) {
      te += unit_treatment_effect_z[exp_d] * treatment_unit_scale[1]
            + category_treatment_effect[exp_k]
            + to_vector(unit_category_te_z)[exp_kd] * unit_category_scale[1];
    } else {
      if (P > 0) {
        te += X * te_beta;
      }
      // Ragged gather: term j predictor q holds L[j] level coefficients laid
      // out term-major then predictor-major in te_re_z.
      int uo = 0;
      int qo = 0;
      for (j in 1 : J) {
        for (q in 1 : Q[j]) {
          te += col(Z, qo + q)
                .* (segment(te_re_z, uo + 1, L[j]) * te_re_scale[qo + q])[re_level[j]];
          uo += L[j];
        }
        qo += Q[j];
      }
    }
  }

  if (is_nb == 1) {
    phi_unit = sample_disp == 1 ? inv(disp) : rep_vector(inv(nb_disp), D);
  }
}
model {
  // mu with the treatment effect scattered onto exposed cells (local: the
  // reporting layer reconstructs it from mu_ctrl + te).
  vector[KDN] mu = mu_ctrl;
  if (model_treated == 1) {
    mu[exp_cell] = mu_ctrl[exp_cell] + te;
  }

  // Baseline priors (joint.py:32-63).
  for (k in 1 : (R_fac > 0 ? K_fac : 0)) {
    to_vector(time_fac[k]) ~ gamma(time_fac_shape, time_fac_shape);
  }
  unit_fe_sigma ~ normal(0, 0.5);           // half-normal via <lower=0>
  to_vector(unit_fe_z) ~ std_normal();
  to_vector(time_fe) ~ gamma(time_fe_shape, time_fe_shape);
  // Component weights. Without rank shrinkage, unit_weight ~
  // Dirichlet(1,...,1) is uniform on the simplex and the declaration alone
  // supplies the prior (a dirichlet statement would add only a constant).
  // With it, the two-level HDP: a sparse shared profile per group, and each
  // unit's weights Dirichlet-distributed around that profile.
  //
  // The sticks are the whole trick. stick ~ Beta(1, mass) makes the profile
  // a DP(mass) draw truncated at R atoms, for which
  // E[sum_r profile^2] = 1 / (1 + mass) exactly -- the number of components
  // a group occupies, and with it the prior amplitude of the low-rank term,
  // is set by the mass and not by the truncation. Components past that count
  // get group_weight[k, r] near zero, and every unit in the group inherits
  // that near-zero weight instead of independently finding a use for the
  // spare column. That is what makes the fit robust to an over-specified
  // rank.
  //
  // unit_weight_sd sets how far a unit may depart from its group's profile:
  // near zero pins every unit to the shared profile, large lets each unit
  // load on its own few components. Both it and the mass are sampled, not
  // fixed, because how much cross-unit sharing the panel supports is exactly
  // what one does not know in advance -- and their posteriors say more about
  // that than a rank sweep does.
  //
  // Both are shared across groups: two scalars estimated from all K * D
  // weight vectors, where per-group versions would be 2K scalars each
  // informed by D. The profile itself stays per-group.
  if (rank_shrink == 1) {
    group_weight_mass[1] ~ gamma(group_mass_shape, group_mass_rate);
    unit_weight_sd[1] ~ normal(0, unit_sd_scale); // half-normal via <lower=0>
    to_vector(unit_weight_z) ~ std_normal();
    for (k in 1 : K_stick) {
      stick[k] ~ beta(1, group_weight_mass[1]);
    }
    // The group level is informed by K profiles, so group_sd_scale should
    // usually be tighter than unit_sd_scale: with K = 4 its prior matters.
    if (crossed == 1) {
      group_profile_sd[1] ~ normal(0, group_sd_scale);
      to_vector(group_profile_z) ~ std_normal();
      unit_shared_sd[1] ~ normal(0, shared_unit_sd_scale);
      to_vector(unit_shared_z) ~ std_normal();
    }
  }

  // Treatment priors (joint.py:164-207).
  if (model_treated == 1) {
    treatment_it_scale[1] ~ normal(0, 0.1);
    treatment_kt_z ~ std_normal();
    if (te_reg == 0) {
      treatment_unit_scale[1] ~ normal(0, 1);
      treatment_category_scale[1] ~ normal(0, 1);
      unit_category_scale[1] ~ normal(0, 1);
      unit_treatment_effect_z ~ std_normal();
      to_vector(unit_category_te_z) ~ std_normal();
      category_treatment_effect ~ normal(0, treatment_category_scale[1]);
    } else {
      te_beta ~ normal(0, te_beta_prior_scale);
      te_re_scale ~ normal(0, te_re_prior_scale);
      te_re_z ~ std_normal();
    }
  }

  // Dispersion factor prior (joint.py:224-232): disp ~ Uniform(0,1) via the
  // constraint, plus factor -1/2 log(disp) - 100 sqrt(disp) per unit.
  if (sample_disp == 1) {
    target += sum(-0.5 * log(disp) - 100 * sqrt(disp));
  }

  // Direct likelihood on observed cells. With fixed dispersion phi is a
  // data scalar, letting the vectorized lpmf hoist every phi-only term.
  if (is_nb == 1) {
    if (sample_disp == 1) {
      y ~ neg_binomial_2_log(mu[obs_cell], phi_unit[obs_unit]);
    } else {
      y ~ neg_binomial_2_log(mu[obs_cell], phi_fixed);
    }
  } else {
    y ~ poisson_log(mu[obs_cell]);
  }

  // Censored small-count marginalization (likelihood.py): suppressed cells
  // contribute log P(y in 1..9); observed cells log P(y not in 1..9).
  // Concentration-only terms are precomputed once per unit (or in
  // transformed data when dispersion is fixed); each cell then costs a
  // single log_sum_exp over the nine-point support.
  if (adjust_missing == 1) {
    vector[n_cens] mass_c;
    vector[n_notcens] mass_nc;
    if (is_nb == 1 && sample_disp == 1) {
      vector[D] log_phi = log(phi_unit);
      array[D] vector[9] b;
      for (d in 1 : D) {
        b[d] = nb_censor_coeff(phi_unit[d], sup, lgamma_sup1);
      }
      for (i in 1 : n_cens) {
        int d = cens_unit[i];
        mass_c[i] = suppressed_mass_nb(mu[cens_cell[i]], phi_unit[d],
                                       log_phi[d], b[d], sup);
      }
      for (i in 1 : n_notcens) {
        int d = notcens_unit[i];
        mass_nc[i] = suppressed_mass_nb(mu[notcens_cell[i]], phi_unit[d],
                                        log_phi[d], b[d], sup);
      }
    } else if (is_nb == 1) {
      for (i in 1 : n_cens) {
        mass_c[i] = suppressed_mass_nb(mu[cens_cell[i]], phi_fixed,
                                       log_phi_fixed, b_fixed, sup);
      }
      for (i in 1 : n_notcens) {
        mass_nc[i] = suppressed_mass_nb(mu[notcens_cell[i]], phi_fixed,
                                        log_phi_fixed, b_fixed, sup);
      }
    } else {
      for (i in 1 : n_cens) {
        mass_c[i] = suppressed_mass_pois(mu[cens_cell[i]],
                                         neg_lgamma_sup1, sup);
      }
      for (i in 1 : n_notcens) {
        mass_nc[i] = suppressed_mass_pois(mu[notcens_cell[i]],
                                          neg_lgamma_sup1, sup);
      }
    }
    target += sum(mass_c);
    target += sum(log1m_exp(mass_nc));
  }
}
generated quantities {
  // The shared component-popularity profile on the natural scale (the model
  // itself only needs its log), and the effective number of components each
  // group is using: the inverse Simpson index 1 / sum_r group_weight[k, r]^2.
  // eff_rank reads m when m components split the mass evenly and the other
  // R - m are unused, so it is the rank actually in use -- with no cutoff to
  // pick, and invariant to the component relabeling that makes group_weight
  // useless to monitor coordinate by coordinate across chains (stick-breaking
  // orders the profile only in expectation, not within a draw). Its prior
  // mean is about 1 + mass. If it sits well below R the truncation was
  // generous enough; if it presses against R, raise the rank and refit.
  //
  // Because time_fe carries the common direction, eff_rank near 1 means the
  // units share one curve that time_fe absorbs -- no deviation beyond the
  // trend -- and eff_rank tracks the rank of the rate surface, trend
  // included. With shared curves, global_eff_rank is the same summary of the
  // global profile: the rank the panel uses overall, of which each group's
  // eff_rank uses a subset.
  array[rank_shrink == 1 ? K : 0] vector[R] group_weight;
  vector[rank_shrink == 1 ? K : 0] eff_rank;
  array[crossed] vector[R] global_weight;
  vector[crossed] global_eff_rank;
  for (k in 1 : (rank_shrink == 1 ? K : 0)) {
    group_weight[k] = exp(log_group_weight[k]);
    eff_rank[k] = inv(dot_self(group_weight[k]));
  }
  if (crossed == 1) {
    global_weight[1] = exp(log_global_weight[1]);
    global_eff_rank[1] = inv(dot_self(global_weight[1]));
  }

  // Counterfactual untreated posterior predictive: rate exp(mu_ctrl), never
  // exp(mu) -- matches Python's Predictive(model_treated=False) recomputing
  // from latents. Log-mean clamped so the count RNG cannot overflow (2^30).
  vector[gen_ypred == 1 ? KDN : 0] ypred;
  if (gen_ypred == 1) {
    for (c in 1 : KDN) {
      real lm = fmin(mu_ctrl[c], 20.79);
      ypred[c] = is_nb == 1
                 ? neg_binomial_2_log_rng(lm, phi_unit[cell_unit[c]])
                 : poisson_log_rng(lm);
    }
  }
}
