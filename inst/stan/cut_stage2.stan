// Cut stage-2 model (port of models/cut_treatment.py): the treatment block
// conditioned on one frozen stage-1 draw of the baseline surface. mu_ctrl and
// the matched NB concentration enter as data, so exposed outcomes can never
// feed back into the baseline -- the cut boundary.
functions {
  #include censoring.stanfunctions
}
data {
  int<lower=1> K;
  int<lower=1> D;
  int<lower=1> KDN;
  vector[KDN] mu_ctrl;                   // one fixed stage-1 draw

  int<lower=1> n_exposed;
  array[n_exposed] int<lower=1, upper=K> exp_k;
  array[n_exposed] int<lower=1, upper=D> exp_d;
  array[n_exposed] int<lower=1, upper=K * D> exp_kd; // (d-1)*K + k gather index
  array[n_exposed] int<lower=1, upper=KDN> exp_cell;

  // Likelihood subsets are positions within the exposed-cell list.
  int<lower=0> n_obs;                    // exposed & !missing
  array[n_obs] int<lower=1, upper=n_exposed> obs_e;
  array[n_obs] int<lower=0> y;
  int<lower=0> n_cens;                   // exposed & missing (suppressed)
  array[n_cens] int<lower=1, upper=n_exposed> cens_e;

  int<lower=0, upper=1> is_nb;
  int<lower=0, upper=1> adjust_missing;
  vector<lower=0>[is_nb == 1 ? D : 0] phi_unit; // matched stage-1 concentration

  // Optional treatment-effect regression (te_reg == 1): the legacy
  // group/unit/group:unit hierarchy is replaced by a linear fixed-effect
  // surface X * te_beta plus J ragged random-effect terms whose per-level
  // coefficients shrink toward that surface. The iid exposed-cell effect
  // (treatment_kt_z) is kept in both parameterizations.
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
  vector[9] sup;                         // censored-count support 1..9
  vector[9] lgamma_sup1;
  for (v in 1 : 9) {
    sup[v] = v;
    lgamma_sup1[v] = lgamma(v + 1);
  }
  vector[9] neg_lgamma_sup1 = -lgamma_sup1;
  // phi is data here, so the per-unit censoring coefficients are free.
  vector[is_nb == 1 ? D : 0] log_phi = log(phi_unit);
  array[is_nb == 1 ? D : 0] vector[9] b_cens;
  for (d in 1 : (is_nb == 1 ? D : 0)) {
    b_cens[d] = nb_censor_coeff(phi_unit[d], sup, lgamma_sup1);
  }
}
parameters {
  // Declaration order matches the pre-regression model, and the blocks below
  // are zero-size in the parameterization they do not belong to, so the
  // legacy unconstrained vector is unchanged -- same seed, same draws.
  real<lower=0> treatment_it_scale;      // HalfNormal(0.1)
  array[te_reg == 1 ? 0 : 1] real<lower=0> treatment_unit_scale; // HalfNormal(1)
  array[te_reg == 1 ? 0 : 1] real<lower=0> treatment_category_scale; // HalfNormal(1)
  array[te_reg == 1 ? 0 : 1] real<lower=0> unit_category_scale; // HalfNormal(1)
  vector[n_exposed] treatment_kt_z;
  vector[te_reg == 1 ? 0 : D] unit_treatment_effect_z;
  matrix[te_reg == 1 ? 0 : K, te_reg == 1 ? 0 : D] unit_category_te_z;
  vector[te_reg == 1 ? 0 : K] category_treatment_effect; // centered
  // Regression surface; all sizes are zero when te_reg == 0.
  vector[P] te_beta;
  vector<lower=0>[Qtot] te_re_scale;     // HalfNormal(te_re_prior_scale)
  vector[Utot] te_re_z;                  // non-centered level coefficients
}
transformed parameters {
  vector[n_exposed] te = treatment_kt_z * treatment_it_scale;
  vector[n_exposed] mu_exposed;
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
  mu_exposed = mu_ctrl[exp_cell] + te;
}
model {
  treatment_it_scale ~ normal(0, 0.1);
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

  if (is_nb == 1) {
    y ~ neg_binomial_2_log(mu_exposed[obs_e], phi_unit[exp_d[obs_e]]);
  } else {
    y ~ poisson_log(mu_exposed[obs_e]);
  }

  if (adjust_missing == 1) {
    if (is_nb == 1) {
      for (i in 1 : n_cens) {
        int d = exp_d[cens_e[i]];
        target += suppressed_mass_nb(mu_exposed[cens_e[i]], phi_unit[d],
                                     log_phi[d], b_cens[d], sup);
      }
      if (n_obs > 0) {
        vector[n_obs] mass;
        for (i in 1 : n_obs) {
          int d = exp_d[obs_e[i]];
          mass[i] = suppressed_mass_nb(mu_exposed[obs_e[i]], phi_unit[d],
                                       log_phi[d], b_cens[d], sup);
        }
        target += sum(log1m_exp(mass));
      }
    } else {
      for (i in 1 : n_cens) {
        target += suppressed_mass_pois(mu_exposed[cens_e[i]],
                                       neg_lgamma_sup1, sup);
      }
      if (n_obs > 0) {
        vector[n_obs] mass;
        for (i in 1 : n_obs) {
          mass[i] = suppressed_mass_pois(mu_exposed[obs_e[i]],
                                         neg_lgamma_sup1, sup);
        }
        target += sum(log1m_exp(mass));
      }
    }
  }
}
