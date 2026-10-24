"""
Numerical verification of every theorem, proposition and corollary of the paper,
and of the numerical statements made in the text.

Run:  python verify_bounds.py
Every check prints PASS/FAIL; the script exits with status 1 on any failure.
"""

from __future__ import annotations

import sys

import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq

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

# c from the anticommutator of the reflections A = 2P - 1 and B = 2Q - 1
rng_ac = np.random.default_rng(5)  # separate stream: the numbers printed below stay unchanged
worst = 0.0
for _ in range(300):
    dim = int(rng_ac.integers(2, 9))
    P = mj.random_projector(dim, int(rng_ac.integers(1, dim)), rng_ac)
    Q = mj.random_projector(dim, int(rng_ac.integers(1, dim)), rng_ac)
    A, B = 2 * P - np.eye(dim), 2 * Q - np.eye(dim)
    worst = max(worst, abs(mj.max_overlap_projectors(P, Q)
                           - (0.5 + 0.25 * np.linalg.norm(A @ B + B @ A, 2))))
check("Thm 5 remark: c = 1/2 + ||{A, B}||/4 with A = 2P - 1, B = 2Q - 1 (d = 2..8, all ranks)",
      worst < 1e-10, f"max deviation {worst:.1e}")

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

# Beyond theta_bar the minimiser moves continuously from the bisector to an eigenstate
fracs = []
for deg in (68, 70, 75, 80, 85, 89, 90):
    theta = np.radians(deg)
    M, N = pair(theta)
    r_star = mj.min_entropy_sum(M, N, num=20001)[1]
    psi = np.arctan2(r_star[0], r_star[2]) % np.pi  # polar angle in the x-z plane (r ~ -r)
    dist = min(abs(psi - theta / 2), np.pi - abs(psi - theta / 2))
    fracs.append(dist / (theta / 2))  # 0 at the bisector, 1 at an eigenstate
check("Sec. V C: beyond theta_bar the Shannon minimiser moves from the bisector to an eigenstate",
      fracs[0] > 0.05 and bool(np.all(np.diff(fracs) > 0)) and abs(fracs[-1] - 1) < 1e-4,
      "fraction of the way at 68, 70, 75, 80, 85, 89, 90 deg: "
      + ", ".join(f"{f:.2f}" for f in fracs))


# Renyi orders: range of angles in which the bisector is optimal
# (known: Zozor, Bosyk and Portesi, J. Phys. A 46, 465301 (2013), Cor. 1-2)
def bisector_renyi(c: float, a: float) -> float:
    t = np.sqrt(c)
    return 2 * float(mj.renyi(np.array([(1 + t) / 2, (1 - t) / 2]), a))


def bisector_optimal(deg: float, a: float) -> bool:
    M, N = pair(np.radians(deg))
    return (mj.min_entropy_sum(M, N, alpha=a)[0]
            > bisector_renyi(mj.overlap_cmax(M, N), a) - 1e-10)


alpha_dag = brentq(lambda a: bisector_renyi(0.5, a) - 1, 1.1, 2.0)  # MUB: bisector = eigenstate
ok = (bisector_optimal(42, 0.75) and not bisector_optimal(44, 0.75)
      and bisector_optimal(79.5, 1.2) and not bisector_optimal(81, 1.2)
      and not bisector_optimal(90, alpha_dag - 1e-3)
      and all(bisector_optimal(deg, a) for a in (alpha_dag + 1e-3, 2.0, 5.0, np.inf)
              for deg in (10, 30, 50, 70, 85, 89, 90)))
th = np.radians(np.linspace(0.01, 90, 9000))
cs_grid = np.cos(th / 2) ** 2
eig_dist = np.stack([cs_grid, 1 - cs_grid], axis=-1)  # N-statistics of an eigenstate of M
bis_dist = np.stack([(1 + np.sqrt(cs_grid)) / 2, (1 - np.sqrt(cs_grid)) / 2], axis=-1)
for a in (0.1, 0.3, 0.5):  # the eigenstates beat the bisector at every angle
    ok &= bool(np.all(mj.renyi(eig_dist, a) < 2 * mj.renyi(bis_dist, a)))
ok &= (np.allclose(mj.renyi(eig_dist, 0.5), np.log2(1 + np.sin(th)))
       and np.allclose(2 * mj.renyi(bis_dist, 0.5), 2 * np.log2(1 + np.sin(th / 2))))
check("Sec. V C: bisector optimal up to ~43 deg (alpha = 3/4) and ~80 deg (alpha = 1.2), "
      "at all angles for alpha >= alpha_dagger only, at no angle for alpha <= 1/2", ok,
      f"alpha_dagger = {alpha_dag:.4f}")

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
check("Cor. 1-2: direct-sum Renyi (alpha <= 1), Tsallis (gamma > 0) and tensor Renyi bounds hold",
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
check("Thm 3: all bounds strict for 0 < alpha < inf and gamma > 0; tight at alpha = 0 and inf",
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

# Photon counting: detector efficiencies set the bias and the sharpness
H_proj, V_proj = np.diag([1.0, 0.0]), np.diag([0.0, 1.0])  # |H> has Bloch vector +z
ok = True
for e_h, e_v in [(0.9, 0.9), (0.9, 0.7), (0.6, 0.3), (0.2, 0.95)]:
    E_H = e_h * H_proj + 0.5 * ((1 - e_h) * H_proj + (1 - e_v) * V_proj)  # no click: coin toss
    model = mj.TwoOutcomePOVM([0, 0, 1], eta=(e_h + e_v) / 2, mu=(e_h - e_v) / 2)
    ok &= np.allclose(E_H, model.effects[0])
e_d = 0.8
single = mj.TwoOutcomePOVM([0, 0, 1], eta=e_d, mu=e_d - 1)  # no click counted as outcome 2
ok &= np.allclose(e_d * H_proj, single.effects[0]) and abs(single.kappa - 1) < 1e-12
check("Sec. VI: efficiencies (eta_H, eta_V) give mu = (eta_H - eta_V)/2, eta = (eta_H + eta_V)/2; "
      "a single detector gives |mu| + eta = 1", ok)

# --- 9. Covariance form of Lemma 2 and the quantum-optics remarks (Sec. VIII) --
rng_qo = np.random.default_rng(9)
theta = np.radians(60)
M, N = pair(theta)
A, B = mj.sigma_dot(M.axis), mj.sigma_dot(N.axis)
e2 = mj.unit(N.axis - np.cos(theta) * M.axis)
robertson = max(abs(np.trace(mj.density_matrix(np.cos(f) * M.axis + np.sin(f) * e2)
                             @ (A @ B - B @ A))) for f in np.linspace(0, 2 * np.pi, 361))
x, y = rng_qo.uniform(-1.2, 1.2, size=(2, 20000))
inside = x**2 - 2 * x * y * np.cos(theta) + y**2 <= np.sin(theta) ** 2
gamma = np.stack([np.stack([1 - x * x, np.cos(theta) - x * y], -1),
                  np.stack([np.cos(theta) - x * y, 1 - y * y], -1)], -2)
psd = np.linalg.eigvalsh(gamma)[:, 0] >= 0
check("Sec. III: (x, y) in E_theta iff the covariance matrix is PSD; "
      "Sec. VIII: the Robertson bound vanishes in the plane of the axes",
      bool(np.array_equal(inside, psd)) and robertson < 1e-12,
      f"max |<[A, B]>| over in-plane states {robertson:.1e}")

NF = 600  # Fock-space truncation: photon numbers 0..599
a_op = np.diag(np.sqrt(np.arange(1, NF)), 1).astype(complex)
vac = np.zeros(NF, complex)
vac[0] = 1
parity = np.diag((-1.0) ** np.arange(NF))


def displacement(beta: complex) -> np.ndarray:
    return expm(beta * a_op.conj().T - np.conj(beta) * a_op)


def squeezed_vacuum(r: float) -> np.ndarray:
    """S(r)|0> with S(z) = exp[(z* a^2 - z a^dag^2)/2]: x squeezed, p anti-squeezed."""
    return expm(0.5 * r * (a_op @ a_op - a_op.conj().T @ a_op.conj().T)) @ vac


D = displacement(1j)  # |beta| = 1 along the anti-squeezed axis
ok, vals = True, []
for r in (0.0, 1.0, 2.0):
    psi = squeezed_vacuum(r)
    dpar = np.vdot(psi, D @ parity @ D.conj().T @ psi).real
    ok &= (abs(np.vdot(psi, parity @ psi).real - 1) < 1e-8
           and abs(dpar - np.exp(-2 * np.exp(-2 * r))) < 1e-8)
    vals.append(dpar)
check("Sec. VIII (v): squeezed vacua approach simultaneous eigenstates of parity and displaced parity",
      ok, "<D Pi D^+> = " + ", ".join(f"{v:.3f}" for v in vals) + " for r = 0, 1, 2")

ok = True
for al in (0.5, 1.0, 2.0):
    k1, k2 = displacement(al) @ vac, displacement(-al) @ vac
    lam, vec = np.linalg.eigh(np.outer(k1, k1.conj()) + np.outer(k2, k2.conj()))
    cat = (k1 + k2) / np.linalg.norm(k1 + k2)
    ok &= (abs(lam[-1] - (1 + np.exp(-2 * al**2))) < 1e-10
           and abs(abs(np.vdot(cat, vec[:, -1])) - 1) < 1e-10)
check("Sec. VIII (v): pi[Q(a) + Q(-a)] <= 1 + exp(-2|a|^2), with equality for the even cat state", ok)

s_z = -parity  # +1 on odd, -1 on even photon numbers
s_m = np.zeros((NF, NF))
s_m[np.arange(0, NF, 2), np.arange(1, NF, 2)] = 1  # sum_n |2n><2n+1|
s_x, s_y = s_m + s_m.T, -1j * (s_m.T - s_m)
psi = squeezed_vacuum(0.8)
ok = (np.allclose(s_x @ s_x, np.eye(NF)) and np.allclose(s_x @ s_z + s_z @ s_x, 0)
      and np.allclose(s_x @ s_y - s_y @ s_x, 2j * s_z)
      and abs(np.vdot(psi, s_z @ psi).real + 1) < 1e-8)
check("Sec. VIII (v): pseudospin operators obey the Pauli algebra; the squeezed vacuum is an "
      "s_z eigenstate", ok)

t1 = mj.crossover_sqrt_c(lambda c: mj.shannon(mj.omega_direct_sum_projective(c)))
t2 = mj.crossover_sqrt_c(lambda c: mj.shannon(mj.omega_tensor_projective(c)))
print(f"\nCrossover with Maassen-Uffink: sqrt(c) = {t1:.4f} (direct sum), {t2:.4f} (tensor)")
print(f"   i.e. Bloch angle theta = {np.degrees(2*np.arccos(t1)):.2f} deg, "
      f"{np.degrees(2*np.arccos(t2)):.2f} deg")

sys.exit(1 if failures else 0)
