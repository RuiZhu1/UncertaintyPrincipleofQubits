# Material for future quantum-optics work

This document preserves material removed from the focused Discussion of the current manuscript. It is not part of the manuscript's claims or numerical evidence. Citation keys refer to `references.bib`.

## Relation to quantum optics: removed discussion material

These results contrast instructively with the quadratures of a light mode. There, coherent and squeezed states are minimum-uncertainty states: for their principal quadratures they saturate both the Heisenberg relation and the entropic relation of Bia{\l}ynicki-Birula and Mycielski~\cite{BBM1975}; moreover, all pure Gaussian states are equivalent under continuous majorization in phase space and are conjectured to majorize every Wigner-positive state~\cite{VanHerstraeten2023}. For noncommuting two-outcome qubit measurements, and for binary projective measurements with $c<1$, no such state exists (Propositions 1 and 3), and in these cases the optimal state can depend on the uncertainty measure, in line with the finding of Korzekwa and Lostaglio~\cite{KorzekwaLostaglio2016} that no qubit state minimizes all information-theoretic uncertainty measures at once, and with the measure dependence of qubit uncertainty discussed in Ref.~\cite{LuisBosykPortesi2016}.

The boundary of $\mathcal{E}_\theta$, traced by the pure states in the plane of the two axes, plays the role of the family of squeezed states, trading predictability of $N$ for that of $M$. The balanced bisector state, the analogue of a coherent state, minimizes the Shannon sum up to $\bar\theta\approx67^\circ$; for larger angles the minimizer moves towards an eigenstate of either measurement, the analogue of an infinitely squeezed state, and the range in which the balanced state is optimal shrinks for smaller R\'enyi orders and grows for larger ones (Sec. V and Ref.~\cite{ZozorBosykPortesi2013}).

The ellipse itself resembles the covariance ellipse of a squeezed vacuum only loosely: it describes the pair of measurements rather than a single state, and although it is a disk for mutually unbiased measurements and flattens as they become compatible, its area $\pi\sin\theta$ shrinks instead of being conserved. A single qubit cannot be squeezed in the usual sense~\cite{KitagawaUeda1993}, as all its pure states are spin coherent states; the essential difference is algebraic. For quadratures the commutator is a $c$-number, which makes the Heisenberg bound state independent. For binary observables it is the anticommutator that is, blockwise, a $c$-number, and it fixes $c=\frac12+\frac14\|\{A,B\}\|$.

The Robertson bound $\Delta A\,\Delta B\ge|\langle(\hat m\times\hat n)\cdot\vec\sigma\rangle|$, by contrast, vanishes for every state whose Bloch vector lies in the plane of the two axes, which contains all optimal states. Commutator-based relations are therefore uninformative here, and state-independent descriptions are needed: tight variance relations~\cite{Abbott2016} or, for all Schur-concave measures at once, majorization. The information missing from the Robertson bound is carried by the anticommutator term of the Robertson--Schr\"odinger relation,

\[
\det\Gamma\ge\langle(\hat m\times\hat n)\cdot\vec\sigma\rangle^2,
\]

which every pure qubit state saturates for these Pauli observables; in the plane of the axes it reduces to $\det\Gamma=0$, the boundary of $\mathcal{E}_\theta$.

## Extended photonic context

The results apply directly to photonic qubits. Restricted to single photons, the Stokes operators of polarization are the Pauli matrices, and two linear polarization analyzers at a relative angle $\vartheta$ correspond to $\theta=2\vartheta$, so that $c=\cos^2\vartheta$ is Malus's law; the crossover with MU and the angle $\bar\theta$ correspond to $\vartheta\approx40.7^\circ$ and $\vartheta\approx33.5^\circ$, respectively. Coherent-light experiments have tested the corresponding polarization-ensemble statistics, with the state reconstructed from measured Stokes parameters~\cite{Wang2020}, rather than the single-photon biased-POVM or steering tests introduced in the current manuscript.

For many photons the Stokes operators form a collective spin that can be squeezed (polarization squeezing~\cite{Korolkova2002,Bowen2002}); for a single photon, which cannot be squeezed, only the trade-off along the boundary of $\mathcal{E}_\theta$ remains. Imperfect detectors lead to the unsharp and biased measurements of Theorem 4, with the bias set by the mismatch of detection efficiencies; by Proposition 3(iii), such biases can create a trade-off even between two measurements in the same polarization basis, and by Proposition 4 they can invalidate the ideal-detector bounds.

## Open question removed from the manuscript: binary measurements on light modes

Which experimentally accessible binary measurements on a light mode have $c<1$, and what do the finite-dimensional binary bounds imply for them? In infinite dimension $c=\sup_{i,j}\|P_iQ_j\|^2$ can equal $1$ without a common eigenvector.

For the photon-number parity operator $\Pi$ and its suitably displaced counterpart, whose expectation values sample the Wigner function~\cite{Royer1977}, the squeezed vacua $\ket{\psi_r}$ contain only even photon numbers. For a fixed displacement $D(\alpha)$ with $\alpha=i\alpha_p$ along the anti-squeezed quadrature,

\[
\begin{aligned}
  \bra{\psi_r}\Pi\ket{\psi_r}&=1,\\
  \bra{\psi_r}D(\alpha)\Pi D^\dagger(\alpha)\ket{\psi_r}
  &=e^{-2\alpha_p^2e^{-2r}}\longrightarrow1.
\end{aligned}
\]

Thus these states form an approximate common-eigenstate sequence, and no non-trivial state-independent relation exists.

Tests for coherent states, $\{\ket{\nu}\!\bra{\nu},\mathds{1}-\ket{\nu}\!\bra{\nu}\}$, realized by a displacement followed by on/off detection with no-click probability $\pi Q(\nu)=\bra{\nu}\rho\ket{\nu}$, with $Q$ the Husimi function, also have $c=1$, because all states orthogonal to both coherent states are common eigenvectors.

Individual partial sums remain informative: the identity

\[
\lambda_{\max}(\ket{\nu_1}\!\bra{\nu_1}+\ket{\nu_2}\!\bra{\nu_2})
=1+|\langle\nu_1|\nu_2\rangle|
\]

gives

\[
\pi[Q(\nu_1)+Q(\nu_2)]\le1+e^{-|\nu_1-\nu_2|^2/2},
\]

with equality for the phase-aligned optical cat state $\ket{\nu_1}+e^{-i\varphi}\ket{\nu_2}$, $\varphi=\arg\langle\nu_1|\nu_2\rangle$~\cite{Ourjoumtsev2007,Lewenstein2021} (the even cat state for $\nu_2=-\nu_1$), which thus plays the role of the bisector state. Non-trivial relations require $\|\{A,B\}\|<2$, as for the pseudospin operators of Ref.~\cite{Chen2002}, or for polarization measurements on single photons.

## Numerical-methods notes removed from Appendix A

The following statement was removed from the current manuscript:

> The statements about light modes in Sec. VIII were checked in a Fock space of dimension 600.

The associated Fock-space checks remain in `code/verify_bounds.py` as material for future work. They cover squeezed-vacuum parity and displaced parity, the cat-state partial-sum identity, and the pseudospin algebra. They are not presented as numerical evidence for the focused current manuscript.
