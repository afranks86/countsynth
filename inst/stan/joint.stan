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
}
parameters {
  array[K] matrix<lower=0>[N, R] time_fac;   // Gamma(20, 20), logged in use
  vector[K] state_fe_mu;                     // improper flat (no statement)
  vector<lower=0>[K] state_fe_sigma;         // HalfNormal(0.5)
  matrix[D, K] state_fe_z;                   // std normal (non-centered)
  matrix<lower=0>[N, K] time_fe;             // Gamma(1, 1), logged in use
  array[K, D] simplex[R] unit_weight;        // Dirichlet(1,...,1) per (k, d)

  // Treatment block; zero-size when model_treated == 0.
  array[model_treated] real<lower=0> treatment_it_scale;      // HalfNormal(0.1)
  array[model_treated] real<lower=0> treatment_state_scale;   // HalfNormal(1)
  array[model_treated] real<lower=0> treatment_category_scale; // HalfNormal(1)
  array[model_treated] real<lower=0> state_category_scale;    // HalfNormal(1)
  vector[model_treated == 1 ? n_exposed : 0] treatment_kt_z;
  vector[model_treated == 1 ? D : 0] state_treatment_effect_z;
  matrix[model_treated == 1 ? K : 0, model_treated == 1 ? D : 0] state_category_te_z;
  vector[model_treated == 1 ? K : 0] category_treatment_effect; // centered

  // Per-unit dispersion; Uniform(0,1) via constraint + factor prior below.
  vector<lower=0, upper=1>[sample_disp == 1 ? D : 0] disp;
}
transformed parameters {
  vector[KDN] mu_ctrl;                   // untreated log-count surface
  vector[model_treated == 1 ? n_exposed : 0] te;
  vector[is_nb == 1 ? D : 0] phi_unit;   // NB2 concentration per unit

  for (k in 1:K) {
    // Natural-scale low-rank product: log(TF_k %*% W_k) equals the Python
    // log-sum-exp assembly of log time_fac + log unit_weight; both factors
    // are positive and O(1), so the direct product is stable and faster.
    matrix[R, D] W;
    for (d in 1:D) {
      W[, d] = unit_weight[k, d];
    }
    matrix[N, D] factor_kd = time_fac[k] * W;
    for (d in 1:D) {
      real fe_kd = state_fe_mu[k] + state_fe_sigma[k] * state_fe_z[d, k];
      int base = ((k - 1) * D + (d - 1)) * N;
      for (n in 1:N) {
        mu_ctrl[base + n] = log(factor_kd[n, d]) + fe_kd + log(time_fe[n, k])
                            + log_denom[base + n];
      }
    }
  }

  if (model_treated == 1) {
    for (e in 1:n_exposed) {
      te[e] = treatment_kt_z[e] * treatment_it_scale[1]
              + state_treatment_effect_z[exp_d[e]] * treatment_state_scale[1]
              + category_treatment_effect[exp_k[e]]
              + state_category_te_z[exp_k[e], exp_d[e]] * state_category_scale[1];
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
    for (e in 1:n_exposed) {
      mu[exp_cell[e]] += te[e];
    }
  }

  // Baseline priors (joint.py:32-63).
  for (k in 1:K) {
    to_vector(time_fac[k]) ~ gamma(20, 20);
  }
  state_fe_sigma ~ normal(0, 0.5);           // half-normal via <lower=0>
  to_vector(state_fe_z) ~ std_normal();
  to_vector(time_fe) ~ gamma(1, 1);
  // unit_weight ~ Dirichlet(1,...,1) is uniform on the simplex: the
  // declaration alone supplies the prior (a dirichlet statement would add
  // only a constant).

  // Treatment priors (joint.py:164-207).
  if (model_treated == 1) {
    treatment_it_scale[1] ~ normal(0, 0.1);
    treatment_state_scale[1] ~ normal(0, 1);
    treatment_category_scale[1] ~ normal(0, 1);
    state_category_scale[1] ~ normal(0, 1);
    treatment_kt_z ~ std_normal();
    state_treatment_effect_z ~ std_normal();
    to_vector(state_category_te_z) ~ std_normal();
    category_treatment_effect ~ normal(0, treatment_category_scale[1]);
  }

  // Dispersion factor prior (joint.py:224-232): disp ~ Uniform(0,1) via the
  // constraint, plus factor -1/2 log(disp) - 100 sqrt(disp) per unit.
  if (sample_disp == 1) {
    target += sum(-0.5 * log(disp) - 100 * sqrt(disp));
  }

  // Direct likelihood on observed cells.
  if (is_nb == 1) {
    y ~ neg_binomial_2_log(mu[obs_cell], phi_unit[obs_unit]);
  } else {
    y ~ poisson_log(mu[obs_cell]);
  }

  // Censored small-count marginalization (likelihood.py): suppressed cells
  // contribute log P(y in 1..9); observed cells log P(y not in 1..9).
  if (adjust_missing == 1) {
    for (i in 1:n_cens) {
      target += suppressed_mass(mu[cens_cell[i]],
                                is_nb == 1 ? phi_unit[cens_unit[i]] : 1.0,
                                is_nb);
    }
    for (i in 1:n_notcens) {
      target += log1m_exp(suppressed_mass(mu[notcens_cell[i]],
                                          is_nb == 1 ? phi_unit[notcens_unit[i]] : 1.0,
                                          is_nb));
    }
  }
}
generated quantities {
  // Counterfactual untreated posterior predictive: rate exp(mu_ctrl), never
  // exp(mu) — matches Python's Predictive(model_treated=False) recomputing
  // from latents. Log-mean clamped so the count RNG cannot overflow (2^30).
  vector[gen_ypred == 1 ? KDN : 0] ypred;
  if (gen_ypred == 1) {
    for (c in 1:KDN) {
      real lm = fmin(mu_ctrl[c], 20.79);
      ypred[c] = is_nb == 1
                 ? neg_binomial_2_log_rng(lm, phi_unit[cell_unit[c]])
                 : poisson_log_rng(lm);
    }
  }
}
