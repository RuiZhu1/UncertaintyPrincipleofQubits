"""
Core routines for majorization uncertainty relations of two-outcome qubit
measurements.

Conventions
-----------
* A qubit state is described by its Bloch vector r (|r| <= 1),
  rho = (I + r . sigma) / 2.
* A general two-outcome qubit POVM {E0, E1 = I - E0} is written as
      E0 = [ (1 + mu) I + eta * (a . sigma) ] / 2,
  with a unit vector a, eta >= 0 and |mu| + eta <= 1.
  Projective (sharp, rank-one) measurements have mu = 0, eta = 1.
* Outcome distribution: p(r) = ( (1 + mu + eta a.r)/2 , (1 - mu - eta a.r)/2 ).
* Majorization x < y (x is majorized by y) for vectors with equal total means
  every partial sum of the decreasingly ordered x is <= that of y.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq, minimize

# ---------------------------------------------------------------------------
# Pauli matrices and states
# ---------------------------------------------------------------------------

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = (SX, SY, SZ)


def sigma_dot(v: np.ndarray) -> np.ndarray:
    """Return v . sigma for a real 3-vector v."""
    return v[0] * SX + v[1] * SY + v[2] * SZ


def density_matrix(r: np.ndarray) -> np.ndarray:
    """Qubit density matrix with Bloch vector r."""
    return 0.5 * (I2 + sigma_dot(np.asarray(r, dtype=float)))


def unit(v) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def axis_in_xz_plane(angle: float) -> np.ndarray:
    """Unit Bloch vector at polar angle `angle` from +z in the x-z plane."""
    return np.array([np.sin(angle), 0.0, np.cos(angle)])


# ---------------------------------------------------------------------------
# Two-outcome POVMs
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TwoOutcomePOVM:
    """E0 = [(1+mu) I + eta a.sigma]/2,  E1 = I - E0."""

    axis: np.ndarray
    eta: float = 1.0
    mu: float = 0.0

    def __post_init__(self):
        object.__setattr__(self, "axis", unit(self.axis))
        if self.eta < 0 or abs(self.mu) + self.eta > 1 + 1e-12:
            raise ValueError("need eta >= 0 and |mu| + eta <= 1")

    @property
    def effects(self) -> tuple[np.ndarray, np.ndarray]:
        e0 = 0.5 * ((1 + self.mu) * I2 + self.eta * sigma_dot(self.axis))
        return e0, I2 - e0

    @property
    def kappa(self) -> float:
        """Largest achievable outcome probability is (1 + kappa)/2."""
        return abs(self.mu) + self.eta

    def probs(self, r: np.ndarray) -> np.ndarray:
        """Outcome distribution(s) for Bloch vector(s) r of shape (..., 3)."""
        x = np.asarray(r) @ self.axis
        p0 = 0.5 * (1 + self.mu + self.eta * x)
        return np.stack([p0, 1 - p0], axis=-1)

    def probs_born(self, r: np.ndarray) -> np.ndarray:
        """Same as `probs`, but via the Born rule tr(rho E_i) (for checks)."""
        rho = density_matrix(r)
        return np.array([np.trace(rho @ e).real for e in self.effects])


def projective(axis) -> TwoOutcomePOVM:
    return TwoOutcomePOVM(axis=axis, eta=1.0, mu=0.0)


def overlap_cmax(m: TwoOutcomePOVM, n: TwoOutcomePOVM) -> float:
    """c = max_{i,j} |<m_i|n_j>|^2 = (1 + |m.n|)/2 for projective qubit measurements."""
    return 0.5 * (1 + abs(float(m.axis @ n.axis)))


# ---------------------------------------------------------------------------
# Majorization utilities
# ---------------------------------------------------------------------------


def sort_desc(v: np.ndarray) -> np.ndarray:
    return -np.sort(-np.asarray(v, dtype=float), axis=-1)


def lorenz(v: np.ndarray) -> np.ndarray:
    """Partial sums of the decreasingly ordered vector (without the leading 0)."""
    return np.cumsum(sort_desc(v), axis=-1)


def is_majorized(x: np.ndarray, y: np.ndarray, tol: float = 1e-10) -> bool:
    """True iff x is majorized by y (x < y)."""
    lx, ly = lorenz(x), lorenz(y)
    return bool(abs(lx[-1] - ly[-1]) < 1e-8 and np.all(lx <= ly + tol))


def least_concave_majorant(partial_sums: np.ndarray) -> np.ndarray:
    """
    Least concave majorant of the points (k, S_k), k = 0..d, with S_0 = 0.

    Given the pointwise maxima S_k of the Lorenz curves of a family of vectors,
    the supremum of the family in the majorization lattice is the vector whose
    Lorenz curve is the least concave majorant of (k, S_k)
    [Cicalese & Vaccaro, IEEE Trans. Inf. Theory 48, 933 (2002)].
    """
    s = np.concatenate([[0.0], np.asarray(partial_sums, dtype=float)])
    d = len(s) - 1
    # Upper convex hull (monotone chain).
    hull: list[int] = []
    for k in range(d + 1):
        while len(hull) >= 2:
            i, j = hull[-2], hull[-1]
            # remove j if it lies on or below the segment i -> k
            if (s[j] - s[i]) * (k - i) <= (s[k] - s[i]) * (j - i) + 1e-15:
                hull.pop()
            else:
                break
        hull.append(k)
    return np.interp(np.arange(d + 1), hull, s[hull])[1:]


def vector_from_lorenz(partial_sums: np.ndarray) -> np.ndarray:
    return np.diff(np.concatenate([[0.0], partial_sums]))


def lattice_supremum(vectors: np.ndarray) -> np.ndarray:
    """Supremum (least upper bound) of a family of vectors in the majorization lattice."""
    s = lorenz(vectors).max(axis=0)
    return vector_from_lorenz(least_concave_majorant(s))


# ---------------------------------------------------------------------------
# Entropies
# ---------------------------------------------------------------------------


def shannon(v: np.ndarray, axis: int = -1) -> np.ndarray:
    """Shannon 'entropy' sum -v log2 v (also used for vectors with total 2)."""
    v = np.clip(np.asarray(v, dtype=float), 0, None)
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = np.where(v > 0, -v * np.log2(v), 0.0)
    return terms.sum(axis=axis)


def renyi(p: np.ndarray, alpha: float, axis: int = -1) -> np.ndarray:
    """Renyi entropy (bits) of normalized distribution(s) p."""
    p = np.clip(np.asarray(p, dtype=float), 0, None)
    if np.isclose(alpha, 1.0):
        return shannon(p, axis=axis)
    if np.isinf(alpha):
        return -np.log2(p.max(axis=axis))
    return np.log2((p**alpha).sum(axis=axis)) / (1 - alpha)


def binary_entropy(t: np.ndarray) -> np.ndarray:
    t = np.asarray(t, dtype=float)
    return shannon(np.stack([t, 1 - t], axis=-1))


# ---------------------------------------------------------------------------
# Analytic bounds (the results of the paper)
# ---------------------------------------------------------------------------


def omega_direct_sum(m: TwoOutcomePOVM, n: TwoOutcomePOVM) -> np.ndarray:
    """
    Theorem 3: optimal direct-sum bound for two general two-outcome qubit POVMs.

        kappa  = max(|mu|+eta, |nu|+zeta)
        Lambda = max_{s,t = +-1} [ s mu + t nu + | s eta a + t zeta b | ]
        omega  = ( (1+kappa)/2, (1+Lambda-kappa)/2, (1+kappa-Lambda)/2, (1-kappa)/2 )
    """
    kappa = max(m.kappa, n.kappa)
    lam = max(
        s * m.mu + t * n.mu + np.linalg.norm(s * m.eta * m.axis + t * n.eta * n.axis)
        for s in (1, -1)
        for t in (1, -1)
    )
    return 0.5 * np.array([1 + kappa, 1 + lam - kappa, 1 + kappa - lam, 1 - kappa])


def omega_direct_sum_projective(c: float) -> np.ndarray:
    """Theorem 1: omega_plus = (1, sqrt c, 1 - sqrt c, 0)."""
    t = np.sqrt(c)
    return np.array([1.0, t, 1 - t, 0.0])


def omega_tensor_projective(c: float) -> np.ndarray:
    """Theorem 2: omega_times = (w, 1 - w, 0, 0) with w = (1 + sqrt c)^2 / 4."""
    w = 0.25 * (1 + np.sqrt(c)) ** 2
    return np.array([w, 1 - w, 0.0, 0.0])


def maassen_uffink(c: float) -> float:
    """Maassen-Uffink bound H(M) + H(N) >= -log2 c."""
    return -np.log2(c)


# ---------------------------------------------------------------------------
# Numerics over the state space
# ---------------------------------------------------------------------------


def fibonacci_sphere(n: int) -> np.ndarray:
    """Quasi-uniform points on the unit sphere."""
    k = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * k / n)
    theta = np.pi * (1 + 5**0.5) * k
    return np.stack(
        [np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], axis=-1
    )


def random_ball(n: int, rng: np.random.Generator) -> np.ndarray:
    """Bloch vectors uniformly distributed in the unit ball."""
    v = rng.normal(size=(n, 3))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    return v * rng.random((n, 1)) ** (1 / 3)


def state_grid(n_sphere: int = 20000, radii=(1.0, 0.9, 0.7, 0.5, 0.25, 0.0)) -> np.ndarray:
    pts = fibonacci_sphere(n_sphere)
    return np.concatenate([r * pts for r in radii])


def direct_sum(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    return np.concatenate([p, q], axis=-1)


def tensor(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    return (p[..., :, None] * q[..., None, :]).reshape(*p.shape[:-1], -1)


def numerical_omega(m: TwoOutcomePOVM, n: TwoOutcomePOVM, kind: str = "direct",
                    states: np.ndarray | None = None) -> np.ndarray:
    """Supremum of {p(rho) (+|x) q(rho)} over a grid of states (lower estimate of omega)."""
    if states is None:
        states = state_grid()
    p, q = m.probs(states), n.probs(states)
    vecs = direct_sum(p, q) if kind == "direct" else tensor(p, q)
    return lattice_supremum(vecs)


def great_circle(m: TwoOutcomePOVM, n: TwoOutcomePOVM, num: int = 4001) -> np.ndarray:
    """Pure states in the plane spanned by the two measurement axes."""
    e1 = m.axis
    e2 = n.axis - (n.axis @ e1) * e1
    if np.linalg.norm(e2) < 1e-12:  # commuting case: any orthogonal direction
        e2 = np.cross(e1, [1.0, 0, 0]) if abs(e1[0]) < 0.9 else np.cross(e1, [0, 1.0, 0])
    e2 = unit(e2)
    phi = np.linspace(0, 2 * np.pi, num, endpoint=False)
    return np.cos(phi)[:, None] * e1 + np.sin(phi)[:, None] * e2


def min_entropy_sum(m: TwoOutcomePOVM, n: TwoOutcomePOVM, alpha: float = 1.0,
                    num: int = 4001) -> tuple[float, np.ndarray]:
    """
    Exact minimum of H_alpha(M) + H_alpha(N) over all qubit states.

    The outcome statistics depend on r only through (x, y) = (a.r, b.r), whose
    range is an ellipse E.  (i) Shannon entropy is concave in (x, y), so its
    minimum sits on the extreme points of E, i.e. on the great circle spanned by
    the two axes (valid for arbitrary two-outcome POVMs).  (ii) For unbiased
    measurements (mu = nu = 0) scaling r outwards makes both distributions more
    ordered, so the same holds for every Schur-concave functional, e.g. any
    Renyi entropy (Lemma 2 of the paper).  A dense 1D grid is refined with a
    local optimizer.
    """
    e1 = m.axis
    e2 = n.axis - (n.axis @ e1) * e1
    e2 = unit(e2) if np.linalg.norm(e2) > 1e-12 else great_circle(m, n, 5)[1]

    def f(phi):
        r = np.cos(phi) * e1 + np.sin(phi) * e2
        return float(renyi(m.probs(r), alpha) + renyi(n.probs(r), alpha))

    phis = np.linspace(0, 2 * np.pi, num, endpoint=False)
    vals = np.array([f(ph) for ph in phis])
    k = int(np.argmin(vals))
    res = minimize(lambda z: f(z[0]), x0=[phis[k]], method="Nelder-Mead",
                   options={"xatol": 1e-12, "fatol": 1e-14})
    best = min(res.fun, vals[k])
    phi_star = res.x[0] if res.fun <= vals[k] else phis[k]
    return best, np.cos(phi_star) * e1 + np.sin(phi_star) * e2


def crossover_sqrt_c(bound) -> float:
    """sqrt(c) at which a majorization bound equals the Maassen-Uffink bound."""
    g = lambda t: bound(t**2) - maassen_uffink(t**2)
    return brentq(g, 1 / np.sqrt(2) + 1e-9, 1 - 1e-9)


# ---------------------------------------------------------------------------
# Exact direct-sum bound for arbitrary POVMs in any dimension
# ---------------------------------------------------------------------------


def omega_direct_sum_operators(*povms: list[np.ndarray]) -> np.ndarray:
    """
    Optimal direct-sum majorization bound for any finite family of POVMs.

    The sum of any fixed subset of outcome probabilities equals tr(rho F) with F
    the sum of the corresponding effects, whose maximum over states is
    lambda_max(F).  Hence S_k = max_{|J| = k} lambda_max(sum_{j in J} F_j), and
    omega is read off the least concave majorant of (k, S_k).  Exponential in
    the total number of outcomes, which is harmless for the cases studied here.
    """
    from itertools import combinations

    effects = [e for povm in povms for e in povm]
    d = len(effects)
    s = np.empty(d)
    for k in range(1, d + 1):
        s[k - 1] = max(np.linalg.eigvalsh(sum(effects[j] for j in idx))[-1]
                       for idx in combinations(range(d), k))
    return vector_from_lorenz(least_concave_majorant(s))


def random_projector(dim: int, rank: int, rng: np.random.Generator) -> np.ndarray:
    z = rng.normal(size=(dim, rank)) + 1j * rng.normal(size=(dim, rank))
    q, _ = np.linalg.qr(z)
    return q @ q.conj().T


def max_overlap_projectors(P: np.ndarray, Q: np.ndarray) -> float:
    """c = max_{i,j} ||P_i Q_j||_inf^2 for the two-outcome measurements {P, I-P}, {Q, I-Q}."""
    I = np.eye(P.shape[0])
    return max(np.linalg.norm(a @ b, 2) ** 2 for a in (P, I - P) for b in (Q, I - Q))
