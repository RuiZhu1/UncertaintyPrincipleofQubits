"""
Generate all figures of the paper into ../paper/figures/.

Run:  python make_figures.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")  # headless: no GUI backend needed
import matplotlib.pyplot as plt
import numpy as np

import majorization as mj

OUT = Path(__file__).resolve().parent.parent / "paper" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# Categorical slots in fixed order (validated reference palette, light mode).
C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

# Journal style: all text typeset by LaTeX (Computer Modern, as in the paper),
# closed axes boxes, and ticks pointing inwards on all four sides of every panel.
mpl.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{amsmath}\usepackage{amssymb}",
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 9,
    "axes.labelsize": 9,
    "legend.fontsize": 7.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "text.color": INK,
    "axes.edgecolor": INK,
    "axes.labelcolor": INK,
    "axes.linewidth": 0.8,
    "axes.grid": False,
    "axes.spines.top": True,
    "axes.spines.right": True,
    "xtick.color": INK,
    "ytick.color": INK,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "xtick.minor.size": 2.0,
    "ytick.minor.size": 2.0,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "lines.linewidth": 1.5,
    "legend.frameon": False,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})
COL = 3.4  # single-column width (inches) for Phys. Rev. A


def pair(theta: float):
    return mj.projective(mj.axis_in_xz_plane(0.0)), mj.projective(mj.axis_in_xz_plane(theta))


# ---------------------------------------------------------------------------
# Fig. 1 - the achievable set in the (x, y) = (<m.sigma>, <n.sigma>) plane
# ---------------------------------------------------------------------------
def fig_achievable_set():
    fig, ax = plt.subplots(figsize=(COL, COL))
    phi = np.linspace(0, 2 * np.pi, 600)
    for theta, col, ls in [(np.pi / 6, C1, "-"), (np.pi / 3, C2, "--"), (np.pi / 2, C3, "-.")]:
        x, y = np.cos(phi), np.cos(phi - theta)
        ax.fill(x, y, color=col, alpha=0.08, lw=0)
        ax.plot(x, y, color=col, ls=ls, label=rf"$\theta={np.degrees(theta):.0f}^\circ$")

    # Saturating states for theta = 60 deg
    theta = np.pi / 3
    t = np.cos(theta / 2)
    ax.plot([t, -t], [t, -t], "o", ms=5, color=C2, mec="white", mew=1.2, zorder=5,
            label=r"$\hat r=\pm(\hat m+\hat n)/|\hat m+\hat n|$")
    ax.plot([1, -1, np.cos(theta), -np.cos(theta)], [np.cos(theta), -np.cos(theta), 1, -1],
            "s", ms=4.5, color=C2, mec="white", mew=1.2, zorder=5,
            label=r"eigenstates of $M$, $N$")

    ax.set_xlabel(r"$x=\hat m\cdot\vec r$")
    ax.set_ylabel(r"$y=\hat n\cdot\vec r$")
    ax.set_xlim(-1.15, 1.15)
    ax.set_ylim(-1.15, 1.15)
    ax.set_aspect("equal")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, handlelength=2.2,
              fontsize=7, columnspacing=1.2)
    fig.savefig(OUT / "achievable_set.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Fig. 2 - Shannon bounds vs Bloch angle
# ---------------------------------------------------------------------------
def fig_entropic_bounds():
    thetas = np.linspace(1e-3, np.pi / 2, 181)
    exact, mu, bplus, btimes = [], [], [], []
    for th in thetas:
        M, N = pair(th)
        c = mj.overlap_cmax(M, N)
        exact.append(mj.min_entropy_sum(M, N, num=1201)[0])
        mu.append(mj.maassen_uffink(c))
        bplus.append(mj.shannon(mj.omega_direct_sum_projective(c)))
        btimes.append(mj.shannon(mj.omega_tensor_projective(c)))
    deg = np.degrees(thetas)

    fig, ax = plt.subplots(figsize=(COL, 2.5))
    ax.plot(deg, exact, color=INK, lw=1.8, label="exact minimum")
    ax.plot(deg, bplus, color=C1, label=r"$H(\omega_\oplus)$ (direct sum)")
    ax.plot(deg, btimes, color=C3, ls=(0, (1.2, 1.2)), label=r"$H(\omega_\otimes)$ (tensor)")
    ax.plot(deg, mu, color=C2, ls="--", label=r"Maassen--Uffink $-\log_2 c$")

    t1 = mj.crossover_sqrt_c(lambda c: mj.shannon(mj.omega_direct_sum_projective(c)))
    th1 = np.degrees(2 * np.arccos(t1))
    ax.axvline(th1, color=INK2, lw=0.6, ls=":")
    ax.text(th1 - 1.5, 0.08, rf"${th1:.1f}^\circ$", ha="right", fontsize=7, color=INK2)

    ax.set_xlabel(r"Bloch angle $\theta$ between $\hat m$ and $\hat n$ (deg)")
    ax.set_ylabel(r"lower bound on $H(M)+H(N)$ (bits)")
    ax.set_xlim(0, 90)
    ax.set_ylim(0, 1.05)
    ax.set_xticks([0, 15, 30, 45, 60, 75, 90])
    ax.legend(loc="upper left")
    fig.savefig(OUT / "entropic_bounds.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Fig. 3 - Lorenz curves at theta = 60 deg
# ---------------------------------------------------------------------------
def fig_lorenz():
    theta = np.pi / 3
    M, N = pair(theta)
    c = mj.overlap_cmax(M, N)
    w = mj.omega_direct_sum_projective(c)
    k = np.arange(5)
    L = lambda v: np.concatenate([[0], mj.lorenz(v)])

    rng = np.random.default_rng(1)
    states = mj.random_ball(150, rng)
    vecs = mj.direct_sum(M.probs(states), N.probs(states))

    fig, ax = plt.subplots(figsize=(COL, 2.5))
    for v in vecs:
        ax.plot(k, L(v), color=GRID, lw=0.6, zorder=1)
    t = np.cos(theta / 2)
    eig = mj.direct_sum(M.probs(M.axis), N.probs(M.axis))
    bis = mj.direct_sum(M.probs(mj.unit(M.axis + N.axis)), N.probs(mj.unit(M.axis + N.axis)))
    ax.plot(k, L(eig), color=C2, ls="--", marker="s", ms=4, label=r"eigenstate $\hat r=\hat m$")
    ax.plot(k, L(bis), color=C3, ls="-.", marker="^", ms=4.5, label=r"bisector $\hat r\propto\hat m+\hat n$")
    ax.plot(k, L(w), color=C1, lw=1.8, marker="o", ms=4.5, label=r"optimal bound $\omega_\oplus$", zorder=4)
    ax.plot([], [], color=GRID, lw=2, label="random states")

    ax.set_xlabel(r"$k$")
    ax.set_ylabel(r"$\sum_{i\leq k}(p\oplus q)^\downarrow_i$")
    ax.set_xticks(k)
    ax.set_ylim(0, 2.08)
    ax.legend(loc="lower right")
    ax.text(0.03, 0.95, rf"$\theta=60^\circ,\ \sqrt{{c}}={t:.3f}$", transform=ax.transAxes,
            fontsize=7.5, color=INK2, va="top")
    fig.savefig(OUT / "lorenz.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Fig. 4 - the whole Renyi family at fixed Bloch angle
# ---------------------------------------------------------------------------
def fig_renyi_family():
    alphas = np.logspace(-1, np.log10(50), 141)
    valid = alphas <= 1
    fig, axes = plt.subplots(2, 1, figsize=(COL, 4.3), sharex=True)
    for ax, deg, tag in zip(axes, (60, 85), ("(a)", "(b)")):
        M, N = pair(np.radians(deg))
        c = mj.overlap_cmax(M, N)
        exact = np.array([mj.min_entropy_sum(M, N, alpha=a, num=1201)[0] for a in alphas])
        dsum = np.array([mj.renyi_direct_sum_bound(c, a) for a in alphas])
        tens = np.array([mj.renyi_tensor_bound(c, a) for a in alphas])
        ax.plot(alphas, exact, color=INK, lw=1.8, label="exact minimum")
        ax.plot(alphas[valid], dsum[valid], color=C1,
                label=r"direct sum $H_\alpha(\sqrt{c},1-\sqrt{c})$, $\alpha\leq1$")
        ax.plot(alphas[~valid], dsum[~valid], color=C1, lw=0.9, ls=(0, (4, 2)),
                label=r"same expression, $\alpha>1$ (not guaranteed)")
        ax.plot(alphas, tens, color=C3, ls=(0, (1.2, 1.2)), label=r"tensor product $H_\alpha(\omega_\otimes)$")
        ax.plot(alphas[valid], np.full(valid.sum(), mj.maassen_uffink(c)), color=C2, ls="--",
                label=r"Maassen--Uffink, $\alpha\leq1$")
        ax.axvline(1, color=INK2, lw=0.6, ls=":")
        ax.text(0.97, 0.93, rf"{tag}\quad$\theta={deg}^\circ$", transform=ax.transAxes, ha="right", va="top",
                fontsize=8, color=INK2)
        ax.set_ylabel(r"bound on $H_\alpha(M)+H_\alpha(N)$", fontsize=8)
        ax.set_ylim(0, 1.08)
    axes[-1].set_xscale("log")
    axes[-1].set_xlim(alphas[0], alphas[-1])
    axes[-1].set_xticks([0.1, 0.5, 1, 2, 5, 10, 50])
    axes[-1].set_xticklabels(["0.1", "0.5", "1", "2", "5", "10", "50"])
    axes[-1].set_xlabel(r"R\'enyi order $\alpha$")
    axes[1].legend(loc="lower left", fontsize=6.8, handlelength=2.4)
    fig.savefig(OUT / "renyi_family.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Fig. 5 - unsharp measurements: bound vs exact as a function of sharpness
# ---------------------------------------------------------------------------
def fig_povm():
    etas = np.linspace(0.02, 1.0, 50)
    fig, ax = plt.subplots(figsize=(COL, 2.5))
    for theta, col, lab in [(np.pi / 2, C1, r"90^\circ"), (np.pi / 4, C2, r"45^\circ")]:
        exact, bound = [], []
        for eta in etas:
            M = mj.TwoOutcomePOVM(mj.axis_in_xz_plane(0), eta=eta)
            N = mj.TwoOutcomePOVM(mj.axis_in_xz_plane(theta), eta=eta)
            exact.append(mj.min_entropy_sum(M, N, num=1201)[0])
            bound.append(mj.shannon(mj.omega_direct_sum(M, N)))
        ax.plot(etas, exact, color=col, label=rf"exact, $\theta={lab}$")
        ax.plot(etas, bound, color=col, ls="--", label=rf"$H(\omega_\oplus)$, $\theta={lab}$")
    ax.set_xlabel(r"sharpness $\eta=\zeta$ (unbiased, $\mu=\nu=0$)")
    ax.set_ylabel(r"$H(M)+H(N)$ (bits)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0.3, 2.1)
    ax.legend(loc="lower left", ncol=1)
    fig.savefig(OUT / "povm.pdf")
    plt.close(fig)


if __name__ == "__main__":
    fig_achievable_set()
    fig_entropic_bounds()
    fig_lorenz()
    fig_renyi_family()
    fig_povm()
    print("figures written to", OUT)
