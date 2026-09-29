"""
Numerical verification of Theorems 1-3 and the entropic corollaries.

Run:  python verify_bounds.py
Every check prints PASS/FAIL; the script exits with status 1 on any failure.
"""

from __future__ import annotations

import sys

import numpy as np

import majorization as mj

rng = np.random.default_rng(2026)
STATES = mj.state_grid(n_sphere=30000)
RANDOM = mj.random_ball(20000, rng)
failures = 0


def check(name: str, ok: bool, detail: str = "") -> None:
    global failures
    failures += not ok
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))


def random_povm(rng: np.random.Generator, sharp: bool = False) -> mj.TwoOutcomePOVM:
    axis = rng.normal(size=3)
    if sharp:
        return mj.projective(axis)
    eta = rng.random()
    mu = (1 - eta) * rng.uniform(-1, 1)
    return mj.TwoOutcomePOVM(axis=axis, eta=eta, mu=mu)


# --- 0. Bloch parametrisation agrees with the Born rule ---------------------
m = random_povm(rng)
r = mj.random_ball(1, rng)[0]
check("Bloch formula == Born rule", np.allclose(m.probs(r), m.probs_born(r)))

# --- 1. Theorem 1: projective, direct sum -----------------------------------
worst = 0.0
for _ in range(40):
    M, N = random_povm(rng, sharp=True), random_povm(rng, sharp=True)
    c = mj.overlap_cmax(M, N)
    w_an = mj.omega_direct_sum_projective(c)
    w_num = mj.numerical_omega(M, N, "direct", STATES)
    worst = max(worst, np.abs(w_an - w_num).max())
    vecs = mj.direct_sum(M.probs(RANDOM), N.probs(RANDOM))
    assert all(mj.is_majorized(v, w_an) for v in vecs[:2000])
check("Thm 1: omega_plus = (1, sqrt c, 1-sqrt c, 0) matches lattice supremum",
      worst < 2e-3, f"max deviation {worst:.2e}")

# --- 2. Theorem 2: projective, tensor product -------------------------------
worst = 0.0
for _ in range(40):
    M, N = random_povm(rng, sharp=True), random_povm(rng, sharp=True)
    c = mj.overlap_cmax(M, N)
    w_an = mj.omega_tensor_projective(c)
    w_num = mj.numerical_omega(M, N, "tensor", STATES)
    worst = max(worst, np.abs(w_an - w_num).max())
    vecs = mj.tensor(M.probs(RANDOM), N.probs(RANDOM))
    assert all(mj.is_majorized(v, w_an) for v in vecs[:2000])
check("Thm 2: omega_times = (w, 1-w, 0, 0) matches lattice supremum",
      worst < 2e-3, f"max deviation {worst:.2e}")

# --- 3. Theorem 3: general two-outcome POVMs, direct sum --------------------
worst = 0.0
for _ in range(60):
    M, N = random_povm(rng), random_povm(rng)
    w_an = mj.omega_direct_sum(M, N)
    w_num = mj.numerical_omega(M, N, "direct", STATES)
    worst = max(worst, np.abs(w_an - w_num).max())
    assert np.all(np.diff(w_an) <= 1e-12), "omega must be non-increasing"
check("Thm 3: POVM formula matches lattice supremum (and is already ordered)",
      worst < 2e-3, f"max deviation {worst:.2e}")

# --- 3b. Theorem 3 again, via the exact eigenvalue characterisation ---------
worst = 0.0
for _ in range(200):
    M, N = random_povm(rng), random_povm(rng)
    worst = max(worst, np.abs(mj.omega_direct_sum(M, N)
                              - mj.omega_direct_sum_operators(M.effects, N.effects)).max())
check("Thm 3 == eigenvalue formula S_k = max lambda_max(sum of k effects)",
      worst < 1e-10, f"max deviation {worst:.1e}")

# --- 3c. Theorem 4: two-outcome projective measurements in dimension d ------
worst = 0.0
for _ in range(200):
    dim = int(rng.integers(2, 7))
    P = mj.random_projector(dim, int(rng.integers(1, dim)), rng)
    Q = mj.random_projector(dim, int(rng.integers(1, dim)), rng)
    I = np.eye(dim)
    c = mj.max_overlap_projectors(P, Q)
    w = mj.omega_direct_sum_operators([P, I - P], [Q, I - Q])
    worst = max(worst, np.abs(w - mj.omega_direct_sum_projective(c)).max())
check("Thm 4: omega_plus = (1, sqrt c, 1-sqrt c, 0) in dimension 2..6, c = max||P_i Q_j||^2",
      worst < 1e-8, f"max deviation {worst:.1e}")

# --- 4. Non-attainability (Proposition 1) -----------------------------------
M, N = mj.projective(mj.axis_in_xz_plane(0)), mj.projective(mj.axis_in_xz_plane(np.pi / 3))
c = mj.overlap_cmax(M, N)
w = mj.omega_direct_sum_projective(c)
vecs = mj.direct_sum(M.probs(STATES), N.probs(STATES))
gap = np.min(np.max(np.abs(mj.sort_desc(vecs) - w), axis=1))
check("Prop 1: no state attains omega_plus for non-commuting M, N", gap > 1e-2,
      f"min sup-distance {gap:.3f}")

# --- 5. Entropic corollaries -------------------------------------------------
ok = True
for theta in np.linspace(0.01, np.pi / 2, 25):
    M, N = mj.projective(mj.axis_in_xz_plane(0)), mj.projective(mj.axis_in_xz_plane(theta))
    c = mj.overlap_cmax(M, N)
    exact, _ = mj.min_entropy_sum(M, N)
    b_plus = mj.shannon(mj.omega_direct_sum_projective(c))
    b_times = mj.shannon(mj.omega_tensor_projective(c))
    ok &= exact >= b_plus - 1e-9 and exact >= b_times - 1e-9 and exact >= mj.maassen_uffink(c) - 1e-9
    # min-entropy: tight Landau-Pollak/Deutsch bound
    hinf, _ = mj.min_entropy_sum(M, N, alpha=np.inf)
    ok &= abs(hinf - (-2 * np.log2((1 + np.sqrt(c)) / 2))) < 1e-6
check("Shannon bounds are valid; H_inf bound -2 log2((1+sqrt c)/2) is tight", ok)

ok = True
for _ in range(200):
    M, N = random_povm(rng), random_povm(rng)
    exact, _ = mj.min_entropy_sum(M, N, num=801)
    ok &= exact >= mj.shannon(mj.omega_direct_sum(M, N)) - 1e-8
check("Cor. 3: H(M)+H(N) >= H(omega_plus) for random unsharp/biased POVMs", ok)

t1 = mj.crossover_sqrt_c(lambda c: mj.shannon(mj.omega_direct_sum_projective(c)))
t2 = mj.crossover_sqrt_c(lambda c: mj.shannon(mj.omega_tensor_projective(c)))
print(f"\nCrossover with Maassen-Uffink: sqrt(c) = {t1:.4f} (direct sum), {t2:.4f} (tensor)")
print(f"   i.e. Bloch angle theta = {np.degrees(2*np.arccos(t1)):.2f} deg, "
      f"{np.degrees(2*np.arccos(t2)):.2f} deg")

sys.exit(1 if failures else 0)
