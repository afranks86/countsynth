"""Export cross-language parity fixtures from the Python implementation.

Run from the Python repo so its venv resolves:
    cd ~/hierarchical-bayesian-NMF-refactor
    uv run python ~/countsynth/tools/export_fixtures.py

Writes JSON fixtures into countsynth/tests/testthat/fixtures/:
  - flatten_parity.json      numpy reshape(-1) ordering of arrays/masks
  - logdensity_parity.json   constrained parameter draws + numpyro log
                             densities for several model configurations
  - cut_mechanics.json       chain quotas + linspace-round subsample indices
"""

import json
from pathlib import Path

import numpy as np

FIXTURE_DIR = Path(__file__).resolve().parent.parent / "tests" / "testthat" / "fixtures"
FIXTURE_DIR.mkdir(parents=True, exist_ok=True)

K, D, N, RANK = 2, 3, 5, 2
rng = np.random.default_rng(20260817)

# --- synthetic panel -------------------------------------------------------
denominators = rng.uniform(0.5, 3.0, size=(K, D, N))
Y = rng.poisson(50.0 * denominators).astype(float)
control = np.ones((K, D, N), dtype=bool)
control[:, 2, 3:] = False  # unit 3 treated from period 4 on
missing = np.zeros((K, D, N), dtype=bool)
missing[0, 0, 1] = True
missing[1, 1, 4] = True
Y[missing] = 0.0

flatten = {
    "K": K,
    "D": D,
    "N": N,
    "rank": RANK,
    "Y_flat": Y.reshape(-1).tolist(),
    "denominators_flat": denominators.reshape(-1).tolist(),
    "control_flat": control.reshape(-1).astype(int).tolist(),
    "missing_flat": missing.reshape(-1).astype(int).tolist(),
    "exposed_cells_1based": (np.where(~control.reshape(-1))[0] + 1).tolist(),
}
(FIXTURE_DIR / "flatten_parity.json").write_text(json.dumps(flatten))
print("wrote flatten_parity.json")

# --- cut mechanics ---------------------------------------------------------
def chain_quotas(m, n_chains):
    base, extra = divmod(m, n_chains)
    return [base + (1 if c < extra else 0) for c in range(n_chains)]


def strided(n, quota):
    return (np.round(np.linspace(0, n - 1, quota)).astype(int) + 1).tolist()


cut_mech = {
    "quotas": [
        {"m": m, "chains": c, "result": chain_quotas(m, c)}
        for m, c in [(10, 4), (25, 4), (7, 3), (4, 4), (5, 2), (1, 4), (3, 8)]
    ],
    "strided": [
        {"n": n, "quota": q, "result_1based": strided(n, q)}
        for n, q in [(10, 4), (250, 100), (100, 100), (7, 3), (9, 4), (5, 1), (6, 4)]
    ],
}
(FIXTURE_DIR / "cut_mechanics.json").write_text(json.dumps(cut_mech))
print("wrote cut_mechanics.json")

# --- log-density parity ----------------------------------------------------
from jax import numpy as jnp  # noqa: E402
from numpyro.infer.util import log_density  # noqa: E402

from bayesian_panel_nmf.models.joint import model as joint_model  # noqa: E402
from bayesian_panel_nmf.models.cut_treatment import stage2_model  # noqa: E402

M = 8  # parameter sets per configuration
n_exposed = int((~control).sum())


def sample_params(with_treatment, with_disp):
    p = {
        "time_fac": rng.gamma(20.0, 1 / 20.0, size=(N, RANK, K)),
        "state_fe_mu": rng.normal(0, 1, size=(K,)),
        "state_fe_sigma": np.abs(rng.normal(0, 0.5, size=(K,))) + 0.05,
        "state_fe_z": rng.normal(0, 1, size=(D, K)),
        "time_fe": rng.gamma(1.0, 1.0, size=(N, K)) + 0.05,
        "unit_weight": rng.dirichlet(np.ones(RANK), size=(D, K)),
    }
    if with_treatment:
        p.update(
            {
                "treatment_it_scale": np.abs(rng.normal(0, 0.1)) + 0.01,
                "treatment_state_scale": np.abs(rng.normal(0, 1)) + 0.05,
                "treatment_category_scale": np.abs(rng.normal(0, 1)) + 0.05,
                "state_category_scale": np.abs(rng.normal(0, 1)) + 0.05,
                "treatment_kt_z": rng.normal(0, 1, size=(n_exposed,)),
                "state_treatment_effect_z": rng.normal(0, 1, size=(D,)),
                "state_category_te_z": rng.normal(0, 1, size=(K, D)),
                "category_treatment_effect": rng.normal(0, 0.5, size=(K,)),
            }
        )
    if with_disp:
        p["disp"] = rng.uniform(0.05, 0.95, size=(D,))
    return p


def jaxify(params):
    return {k: jnp.asarray(v) for k, v in params.items()}


configs = [
    # name, outcome_dist, model_treated, adjust, sample_disp
    ("nb_treated_adjust", "NB", True, True, False),
    ("poisson_treated_adjust", "Poisson", True, True, False),
    ("nb_baseline_adjust", "NB", False, True, False),
    ("nb_treated_noadjust", "NB", True, False, False),
    ("nb_treated_sampledisp", "NB", True, True, True),
]

cases = []
for name, dist_name, treated, adjust, sample_disp in configs:
    for i in range(M):
        params = sample_params(with_treatment=treated, with_disp=sample_disp)
        kwargs = dict(
            y=jnp.asarray(Y),
            rank=RANK,
            outcome_dist=dist_name,
            adjust_for_missingness=adjust,
            nb_disp=1e-4,
            sample_disp=sample_disp,
            model_treated=treated,
        )
        ld, _ = log_density(
            joint_model,
            (jnp.asarray(denominators), control, missing),
            kwargs,
            jaxify(params),
        )
        cases.append(
            {
                "config": name,
                "model": "joint",
                "outcome_dist": dist_name,
                "model_treated": treated,
                "adjust": adjust,
                "sample_disp": sample_disp,
                "params": {k: np.asarray(v).tolist() for k, v in params.items()},
                "log_density": float(ld),
            }
        )
    print(f"  joint config {name}: done")

# Stage-2: fixed mu_ctrl, treatment block only.
mu_ctrl_fixed = rng.normal(np.log(50.0 * denominators), 0.1)
phi_fixed = np.full(D, 1e4)
for i in range(M):
    params = sample_params(with_treatment=True, with_disp=False)
    params = {
        k: v
        for k, v in params.items()
        if k.startswith(("treatment", "state_category", "category", "state_treatment"))
    }
    ld, _ = log_density(
        stage2_model,
        (jnp.asarray(mu_ctrl_fixed), control),
        dict(
            missing_idx_array=missing,
            y=jnp.asarray(Y),
            outcome_dist="NB",
            nb_concentration=jnp.asarray(phi_fixed),
            adjust_for_missingness=True,
        ),
        jaxify(params),
    )
    cases.append(
        {
            "config": "stage2_nb_adjust",
            "model": "stage2",
            "outcome_dist": "NB",
            "model_treated": True,
            "adjust": True,
            "sample_disp": False,
            "params": {k: np.asarray(v).tolist() for k, v in params.items()},
            "log_density": float(ld),
        }
    )
print("  stage2 config: done")

payload = {
    "data": {
        "K": K,
        "D": D,
        "N": N,
        "rank": RANK,
        "Y": Y.tolist(),
        "denominators": denominators.tolist(),
        "control": control.astype(int).tolist(),
        "missing": missing.astype(int).tolist(),
        "mu_ctrl_fixed": np.asarray(mu_ctrl_fixed).tolist(),
        "phi_fixed": phi_fixed.tolist(),
    },
    "cases": cases,
}
(FIXTURE_DIR / "logdensity_parity.json").write_text(json.dumps(payload))
print("wrote logdensity_parity.json")
