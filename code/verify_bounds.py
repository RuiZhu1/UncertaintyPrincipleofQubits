"""
Numerical verification of every theorem, proposition and corollary of the paper.

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


def pair(theta: float, eta: float = 1.0):
    """Unbiased measurements along z and at Bloch angle theta from z, sharpness eta."""
    return (mj.TwoOutcomePOVM(mj.axis_in_xz_plane(0.0), eta=eta),
            mj.TwoOutcomePOVM(mj.axis_in_xz_plane(theta), eta=eta))


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

# --- 3. Theorem 4: general two-outcome POVMs, direct sum --------------------
worst = 0.0
for _ in range(60):
    M, N = random_povm(rng), random_povm(rng)
    w_an = mj.omega_direct_sum(M, N)
    w_num = mj.numerical_omega(M, N, "direct", STATES)
    worst = max(worst, np.abs(w_an - w_num).max())
    assert np.all(np.diff(w_an) <= 1e-12), "omega must be non-increasing"
check("Thm 4: POVM formula matches lattice supremum (and is already ordered)",
      worst < 2e-3, f"max deviation {worst:.2e}")

# --- 3b. Theorem 4 again, via the exact eigenvalue characterisation ---------
worst = 0.0
for _ in range(200):
    M, N = random_povm(rng), random_povm(rng)
    worst = max(worst, np.abs(mj.omega_direct_sum(M, N)
                              - mj.omega_direct_sum_operators(M.effects, N.effects)).max())
check("Thm 4 == eigenvalue formula S_k = max lambda_max(sum of k effects)",
      worst < 1e-10, f"max deviation {worst:.1e}")

# --- 3c. Theorem 5: binary projective measurements in dimension d ------------
# A common eigenvector forces c = 1.  It always exists for odd d or unequal
# ranks, so non-trivial instances (c < 1) need even d and rank(P) = rank(Q) = d/2.
worst, cs = 0.0, []
for _ in range(150):
    dim = int(rng.choice([2, 4, 6]))
    P = mj.random_projector(dim, dim // 2, rng)
    Q = mj.random_projector(dim, dim // 2, rng)
    I = np.eye(dim)
    c = mj.max_overlap_projectors(P, Q)
    cs.append(c)
    w = mj.omega_direct_sum_operators([P, I - P], [Q, I - Q])
    worst = max(worst, np.abs(w - mj.omega_direct_sum_projective(c)).max())
check("Thm 5: omega_plus in d = 2, 4, 6 with c = max||P_i Q_j||^2 (all 150 instances have c < 1)",
      worst < 1e-8 and max(cs) < 1 - 1e-6,
      f"max deviation {worst:.1e}, c in [{min(cs):.3f}, {max(cs):.3f}]")

worst = 0.0
for _ in range(50):
    dim = int(rng.integers(3, 8))
    rank_p, rank_q = int(rng.integers(1, dim)), int(rng.integers(1, dim))
    if 2 * rank_p == 2 * rank_q == dim:
        continue  # generic c < 1 case, covered above
    P = mj.random_projector(dim, rank_p, rng)
    Q = mj.random_projector(dim, rank_q, rng)
    I = np.eye(dim)
    w = mj.omega_direct_sum_operators([P, I - P], [Q, I - Q])
    worst = max(worst, abs(mj.max_overlap_projectors(P, Q) - 1), np.abs(w - [1, 1, 0, 0]).max())
check("Thm 5, degenerate case: common eigenvector => c = 1 and omega_plus = (1, 1, 0, 0)",
      worst < 1e-8, f"max deviation {worst:.1e}")

# --- 4. Non-attainability (Proposition 1) -----------------------------------
M, N = pair(np.pi / 3)
c = mj.overlap_cmax(M, N)
w = mj.omega_direct_sum_projective(c)
vecs = mj.direct_sum(M.probs(STATES), N.probs(STATES))
gap = np.min(np.max(np.abs(mj.sort_desc(vecs) - w), axis=1))
check("Prop 1: no state attains omega_plus for non-commuting M, N", gap > 1e-2,
      f"min sup-distance {gap:.3f}")

# --- 5. Shannon: validity, exact minimum, critical angle --------------------
ok = True
for theta in np.linspace(0.01, np.pi / 2, 25):
    M, N = pair(theta)
    c = mj.overlap_cmax(M, N)
    exact, _ = mj.min_entropy_sum(M, N)
    b_plus = mj.shannon(mj.omega_direct_sum_projective(c))
    b_times = mj.shannon(mj.omega_tensor_projective(c))
    ok &= exact >= b_plus - 1e-9 and exact >= b_times - 1e-9 and exact >= mj.maassen_uffink(c) - 1e-9
    # min-entropy: tight Landau-Pollak/Deutsch bound
    hinf, _ = mj.min_entropy_sum(M, N, alpha=np.inf)
    ok &= abs(hinf - (-2 * np.log2((1 + np.sqrt(c)) / 2))) < 1e-6
check("Cor. 1-2: Shannon bounds are valid; H_inf bound -2 log2((1+sqrt c)/2) is tight", ok)

c_bar = mj.critical_c()
theta_bar = 2 * np.arccos(np.sqrt(c_bar))
ok = True
for deg in (5, 20, 40, 60, 65, 66.5):
    M, N = pair(np.radians(deg))
    c = mj.overlap_cmax(M, N)
    ok &= abs(mj.min_entropy_sum(M, N, num=20001)[0] - mj.bisector_shannon(c)) < 1e-9
for deg in (68, 70, 80, 90):
    M, N = pair(np.radians(deg))
    c = mj.overlap_cmax(M, N)
    ok &= mj.min_entropy_sum(M, N, num=20001)[0] < mj.bisector_shannon(c) - 1e-6
check("Shannon minimum = 2h((1+sqrt c)/2) exactly iff theta <= theta_bar", ok,
      f"sqrt(c_bar) artanh sqrt(c_bar) = 1: c_bar = {c_bar:.4f}, "
      f"theta_bar = {np.degrees(theta_bar):.2f} deg")

# --- 6. Renyi and Tsallis families (Corollaries 1-2, Proposition 2) ----------
ok_valid = ok_order = True
for deg in (10, 30, 45, 60, 75, 85, 90):
    M, N = pair(np.radians(deg))
    c = mj.overlap_cmax(M, N)
    for a in (0.2, 0.5, 0.8, 1.0):
        rd = mj.renyi_direct_sum_bound(c, a)
        ok_valid &= mj.min_entropy_sum(M, N, alpha=a, num=2001)[0] >= rd - 1e-9
        ok_order &= rd >= mj.renyi_tensor_bound(c, a) - 1e-12
    for a in (1.5, 2.0, 5.0, 20.0, np.inf):
        ok_valid &= (mj.min_entropy_sum(M, N, alpha=a, num=2001)[0]
                     >= mj.renyi_tensor_bound(c, a) - 1e-9)
    for q in (0.3, 0.5, 2.0, 3.0):
        ok_valid &= (mj.min_tsallis_sum(M, N, q, num=2001)[0]
                     >= mj.tsallis_direct_sum_bound(c, q) - 1e-9)
check("Cor. 1-2: direct-sum Renyi (alpha <= 1), Tsallis (q > 0) and tensor Renyi bounds hold",
      ok_valid)
check("Prop 2(i): direct-sum Renyi bound >= tensor bound for alpha <= 1", ok_order)

ok = True
for deg in (30, 60, 85):
    M, N = pair(np.radians(deg))
    c = mj.overlap_cmax(M, N)
    # violation 2 log2[(1 + sqrt c) / (2 c^(1/4))]: 4.2e-4 at 30 deg, larger at wider angles
    ok &= mj.min_entropy_sum(M, N, alpha=np.inf)[0] < mj.renyi_direct_sum_bound(c, np.inf) - 1e-4
check("Prop 2(ii): the direct-sum form fails at alpha = inf", ok)

# --- 7. Theorem 3: strict for 0 < alpha < inf, tight at alpha = 0 and inf ---
gaps, tight = [], 0.0
for deg in (30, 60, 85):
    M, N = pair(np.radians(deg))
    c = mj.overlap_cmax(M, N)
    for a in (0.5, 1.0, 2.0):
        exact = mj.min_entropy_sum(M, N, alpha=a, num=2001)[0]
        gaps.append(exact - mj.renyi_tensor_bound(c, a))
        if a <= 1:
            gaps.append(exact - mj.renyi_direct_sum_bound(c, a))
    for q in (0.5, 2.0):
        gaps.append(mj.min_tsallis_sum(M, N, q, num=2001)[0] - mj.tsallis_direct_sum_bound(c, q))
    for a in (0.0, np.inf):
        tight = max(tight, abs(mj.min_entropy_sum(M, N, alpha=a, num=2001)[0]
                               - mj.renyi_tensor_bound(c, a)))
check("Thm 3: all bounds strict for 0 < alpha < inf and q > 0; tight at alpha = 0 and inf",
      min(gaps) > 1e-4 and tight < 1e-9, f"smallest gap {min(gaps):.1e}, "
      f"largest deviation at alpha in {{0, inf}}: {tight:.1e}")

# --- 8. POVMs: entropic corollary and the joint-measurability example -------
ok = True
for _ in range(200):
    M, N = random_povm(rng), random_povm(rng)
    exact, _ = mj.min_entropy_sum(M, N, num=801)
    ok &= exact >= mj.shannon(mj.omega_direct_sum(M, N)) - 1e-8
check("Thm 4: H(M)+H(N) >= H(omega_plus) for random unsharp/biased POVMs", ok)

eta = 0.7  # orthogonal axes, eta = zeta <= 1/sqrt(2): jointly measurable
M, N = pair(np.pi / 2, eta)
busch = (np.linalg.norm(eta * M.axis + eta * N.axis)
         + np.linalg.norm(eta * M.axis - eta * N.axis))
w = mj.omega_direct_sum(M, N)
vecs = mj.direct_sum(M.probs(STATES), N.probs(STATES))
gap = np.min(np.max(np.abs(mj.sort_desc(vecs) - w), axis=1))
s2 = w[0] + w[1]
check("Sec. VI: jointly measurable pair still has S_2 < 1 + eta and an unattained bound",
      busch <= 2 and abs(s2 - (1 + eta / np.sqrt(2))) < 1e-12 and gap > 1e-3,
      f"Busch sum {busch:.3f} <= 2, S_2 = {s2:.3f} < {1 + eta:.1f}, min sup-distance {gap:.3f}")

t1 = mj.crossover_sqrt_c(lambda c: mj.shannon(mj.omega_direct_sum_projective(c)))
t2 = mj.crossover_sqrt_c(lambda c: mj.shannon(mj.omega_tensor_projective(c)))
print(f"\nCrossover with Maassen-Uffink: sqrt(c) = {t1:.4f} (direct sum), {t2:.4f} (tensor)")
print(f"   i.e. Bloch angle theta = {np.degrees(2*np.arccos(t1)):.2f} deg, "
      f"{np.degrees(2*np.arccos(t2)):.2f} deg")

sys.exit(1 if failures else 0)
