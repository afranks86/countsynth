// Cut stage-2 model (port of models/cut_treatment.py): the treatment block
// conditioned on one frozen stage-1 draw of the baseline surface. mu_ctrl and
// the matched NB concentration enter as data, so exposed outcomes can never
// feed back into the baseline — the cut boundary.
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
}
parameters {
  real<lower=0> treatment_it_scale;      // HalfNormal(0.1)
  real<lower=0> treatment_state_scale;   // HalfNormal(1)
  real<lower=0> treatment_category_scale; // HalfNormal(1)
  real<lower=0> state_category_scale;    // HalfNormal(1)
  vector[n_exposed] treatment_kt_z;
  vector[D] state_treatment_effect_z;
  matrix[K, D] state_category_te_z;
  vector[K] category_treatment_effect;   // centered
}
transformed parameters {
  vector[n_exposed] te;
  vector[n_exposed] mu_exposed;
  for (e in 1:n_exposed) {
    te[e] = treatment_kt_z[e] * treatment_it_scale
            + state_treatment_effect_z[exp_d[e]] * treatment_state_scale
            + category_treatment_effect[exp_k[e]]
            + state_category_te_z[exp_k[e], exp_d[e]] * state_category_scale;
    mu_exposed[e] = mu_ctrl[exp_cell[e]] + te[e];
  }
}
model {
  treatment_it_scale ~ normal(0, 0.1);
  treatment_state_scale ~ normal(0, 1);
  treatment_category_scale ~ normal(0, 1);
  state_category_scale ~ normal(0, 1);
  treatment_kt_z ~ std_normal();
  state_treatment_effect_z ~ std_normal();
  to_vector(state_category_te_z) ~ std_normal();
  category_treatment_effect ~ normal(0, treatment_category_scale);

  if (is_nb == 1) {
    y ~ neg_binomial_2_log(mu_exposed[obs_e], phi_unit[exp_d[obs_e]]);
  } else {
    y ~ poisson_log(mu_exposed[obs_e]);
  }

  if (adjust_missing == 1) {
    for (i in 1:n_cens) {
      target += suppressed_mass(mu_exposed[cens_e[i]],
                                is_nb == 1 ? phi_unit[exp_d[cens_e[i]]] : 1.0,
                                is_nb);
    }
    for (i in 1:n_obs) {
      target += log1m_exp(suppressed_mass(mu_exposed[obs_e[i]],
                                          is_nb == 1 ? phi_unit[exp_d[obs_e[i]]] : 1.0,
                                          is_nb));
    }
  }
}
