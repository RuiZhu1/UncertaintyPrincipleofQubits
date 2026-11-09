# Revision of October 26: what changed in `main.tex` and why

All changes are in `paper/main.tex`. The structure, format and title are unchanged: revtex4-2 PRA style, Secs. I–VIII plus Appendix A, proofs after each result. The changes answer the referee-style review of the same day.

## New content

**Proposition 4, Corollary 3 and Fig. 6** — new paragraph *"Imperfect detectors and uncertainty-based tests"* at the end of Sec. VI.

- (i) Unbiased detector noise is a doubly stochastic map on the statistics, so every bound derived for ideal detectors stays valid.
- (ii) Biased measurements obey the ideal relation p⊕q ≺ (1, √c, 1−√c, 0) for all states **iff Λ ≤ 2√c**. When the condition fails, some state violates the ideal Landau–Pollak sum.
- (iii) If the same detectors read both settings, the condition becomes |μ| ≤ √c(1−η). Click/no-click detection violates it for every η_d < 1.
- Shannon thresholds (orthogonal axes, random no-click assignment):

  | η_H | ideal relation fails for | 1-bit Shannon bound fails for |
  |---|---|---|
  | 0.99 | η_V < 0.942 | η_V < 0.578 |
  | 0.95 | η_V < 0.708 | η_V < 0.211 |
  | 0.9 | η_V < 0.417 | never |

- Click/no-click detection at 90°: the Shannon minimum is h(η_d/2) < 1.
- Corollary 3 (steering): in a local-hidden-state model, H(M|A₁)+H(N|A₂) ≥ min[H(M)+H(N)] ≥ H(ω⊕(M,N)). With biased detectors, the ideal bound can certify steering for a product state.

Other new content:

- **Sec. V A, after Proposition 2:** for qubits the direct-sum Rényi form is a valid bound for every α ≤ 2. This is shown numerically at all angles and analytically at α = 2, where it reduces to (√c−1)³(√c+3) ≤ 0. It fails only above an angle-dependent order between 2.00 and 2.58.
- **Sec. VI, after Theorem 4:** fixed-spectrum (|r| = r) bound ½(1+r, 1−r+2r√c, 1+r−2r√c, 1−r), not attained for r > 0.
- **Sec. VIII:** a paragraph on noise versus bias, the testable non-attainment distance min_ρ max_k[S_k − L_k], and a proposed experiment with mismatched detectors.

## Review comments → changes

| Comment | Change |
|---|---|
| M1 Partovi's priority | Introduction credits Partovi and the 2026 Conlon et al. review. Theorem 2 header cites Partovi (the vector for mutually unbiased bases). After Proposition 1 the text says the proof is short once the saturating states are known. |
| M2 attributions | KTW credited for c = ½ + ¼‖{A,B}‖ in any dimension (Sec. VII and Introduction). Yuan et al. credited for the general optimality of the flattened direct-sum vector. Abbott et al. cited at Theorem 4. Sec. VII separates the assembled parts from the new ones. |
| M3(a) Rényi beyond α = 1 | See "New content" above. Text and caption of Fig. 4 updated. |
| M3(b) eigenstate range | Eigenstates optimal for α ≤ α₀(c), with α₀ > ½ (0.52, 0.56, 0.66 at 20°, 60°, 85°). |
| M4 operational meaning | Proposition 4, Corollary 3, Fig. 6, and the Sec. VIII paragraph (including Wang et al.'s "touching" Lorenz curves). |
| M5 fixed spectrum | Sec. VI. |
| m1 "22%" in Sec. VIII | "Never below 78%, exceeds by up to 28%". |
| m2 bound (19) wording | "Strictly stronger than h(√c), at least as strong as MU" (Sec. V C and Fig. 3 caption). |
| m3 tensor ratio | h(w) ≥ 0.761 × the exact minimum, worst at 70.6°. |
| m4 Proposition 1(ii) proof | Two-vanishing-entries argument. |
| m5 Proposition 3(i) | Sign condition completed. |
| m6 Corollary 4 (binary POVMs, dimension d) | max{1, ·} removed, with justification in the proof. |
| m8 Appendix A | "Unless both ranks equal d/2". |
| m9 pseudospin | Claim that the squeezed vacuum is an eigenstate removed. |
| m10 citations | Luis–Bosyk–Portesi 2016, Puchała et al. 2018, Conlon et al. 2026, Yang–Qiao 2026, Walborn 2011, Schneeloch 2013, Uola 2020 (all in `references.bib`). |

Numbering: the new Proposition 4 follows Proposition 3. The steering corollary is Corollary 3, so the former Corollary 3 (binary POVMs in dimension d) is now Corollary 4.

## Code

- `majorization.py`: `big_lambda`, `ideal_relation_holds`, `detector_pair`, `click_detector`, `shannon_min_in_plane`.
- `verify_bounds.py`: 8 new checks; **46/46 PASS**.
- `make_figures.py`: `fig_detector_bias()` writes `figures/detector_bias.pdf`. Without cm-super the script falls back to Matplotlib's built-in Computer Modern. **Re-run `python make_figures.py` on the workstation** so Fig. 6 is LaTeX-typeset like the others.

## Still open (unchanged from before)

- Gilad's affiliation and funding (TODOs in `main.tex`).
- Marcelo's approval of the abstract, which now has one more sentence on detector bias.
- Journal choice.

# Revision of November 1: peer-review corrections and additions

## Scientific corrections

- Rewrote the abstract and the introduction's contribution statement to distinguish known constructions from the paper's actual contributions.
- Corrected the higher-dimensional tensor-product attribution: the vector is the binary case of the general FGG construction. The new elements are its explicit evaluation, the observation that it is already ordered, and the attainment analysis.
- Replaced Proposition 3(ii) with a necessary and sufficient attainment criterion for every commuting two-outcome qubit POVM. The same endpoint state must saturate both $S_1$ and $S_2$; attainment of $S_3$ then follows automatically.
- Limited the direct-sum Rényi statement through $\alpha=2$ to numerical evidence, except for the analytical endpoint proof at $\alpha=2$.
- Restricted the Sánchez-Ruiz–KTW strict comparison to noncommuting angles, the Robertson–Schrödinger equality to pure qubit states with these Pauli observables, and the squeezed-parity argument to displacement along the anti-squeezed quadrature.

## Steering addition

Added a closed-form Werner-state example with
$\mathcal S(w,\eta_d)=h(\eta_d(1+w)/2)+h(\eta_d(1-w)/2)$ and corrected bound
$B(\eta_d)=h(1-\eta_d+\eta_d/\sqrt2)$. The Bell-state violation threshold is
$\eta_*=(4+\sqrt2)/7\approx0.773459$. The text explicitly states that nonviolation
of this sufficient witness does not prove an LHS model. Added references on
loss-tolerant and measurement-setting-error-tolerant steering.

## Code and verification

- Added `commuting_attainment`, `werner_conditional_entropy`, and `click_mub_bound` to `majorization.py`.
- Added a 2000-pair random test of the commuting attainment criterion and Born-rule checks of the Werner formula and efficiency threshold to `verify_bounds.py`.
- Added the squeezed-direction parity check to document why the displacement qualification is necessary.
- Updated theorem numbering and numerical-minimum terminology in code and documentation.
- The complete suite now has **48/48 PASS**. `main.pdf` was rebuilt successfully and checked visually. One non-fatal REVTeX float-placement warning remains, with all six figures present.

## Submission target and remaining approvals

The target venue is Communications Physics. Remaining author-side items include approval of the revised abstract, corresponding-author and code-availability decisions, and the previously noted affiliation and funding information.

# Follow-up minor revision after the second peer review

- Clarified in the introduction, tightness discussion, and conclusion that Rényi tightness at orders zero and infinity refers specifically to the tensor-product bounds.
- Changed the appendix's remaining “exact minima” wording to “numerical minima”.
- Corrected the README's Sánchez-Ruiz–KTW comparison and the binary-POVM roadmap label to Corollary 4, and added `detector_bias.pdf` to the figure list.
- Retained the original Data Availability Statement: the data are available from the corresponding author upon reasonable request. No repository link is included in the manuscript.

# Clarifications following review of the merged revision

- Unified the introduction's Shannon-gap statement to the ratio of the direct-sum bound to the exact minimum (about 78%), and explicitly identified the bound itself as the denominator when the Discussion quotes the approximately 28% excess.
- Clarified that $\eta_*=(4+\sqrt2)/7$ is a limitation of the present steering witness, not a fundamental steering threshold for Werner states.
- Added the geometric interpretation of the commuting attainment criterion: one extremal state must simultaneously maximize the largest single outcome probability and the largest cross-pair sum, whereas biased measurements can attain these maxima at opposite endpoints.

These changes affect presentation only; no analytical result or numerical code is changed.

# Scope corrections following the external review

- Corrected the Discussion's non-attainment statement: noncommutation is not sufficient in arbitrary finite dimension. The text now covers noncommuting projective qubit measurements and, more generally, binary projective measurements with $c<1$.
- Limited the direct-sum Rényi validity statement to sampled angles and sampled orders.
- Changed “The attaining state” to “An attaining state” in Proposition 3(ii) to allow degenerate alternatives.
- Described the closest-state Lorenz distance as experimentally accessible through optimized state preparation and outcome-statistics estimation, rather than as a single directly measurable observable, and noted finite sampling and precision.
- Stated the steering model's calibration assumptions: efficiencies, biases, and no-click assignment rules must be independently characterized and fixed across Alice's settings; dark counts, drift, and setting-dependent losses require an extended model.

No theorem formula or numerical code is changed.

# Counterexample-driven corrections and regression coverage

- Corrected the remaining quantum-optics statement: nonattainment is now claimed for noncommuting two-outcome qubit measurements and for binary projective measurements with $c<1$; the optimal state can, rather than must, depend on the uncertainty measure.
- Corrected the Fig. 6(a) discussion. The quoted thresholds refer to the lower-$\eta_V$ branch. At $(\eta_H,\eta_V)=(0.9,1)$, the state $-\hat a$ gives $H(M)+H(N)=h(0.475)=0.998196<1$, so the previous whole-square “never violated” claim was false.
- Expanded the Fig. 6 verification to an $81\times81$ grid over the full efficiency square, including detector-label symmetry, validity of the corrected bound, the explicit counterexample, and the high-efficiency branch. For $\eta_H=0.9$, the high branch begins at $\eta_V\approx0.999671$ in the current numerical optimization.
- Added a noncommuting rank-two projector regression in $d=3$ with $c=1$ and both optimal vectors attained.
- Replaced “all bounds of Sec. V” with the explicit direct-sum, tensor-product, strictness, and tightness results proved in that dimension.
- Reworded the Rényi statement so that validity through $\alpha=2$ and the earliest failure beyond it are explicitly numerical observations over sampled angles and orders.
- Changed Proposition 3(ii) to “An attaining state can be chosen as ...” and the introduction to “spuriously certify steering”.

The Data Availability Statement remains unchanged.

# Steering optimality and scope refinements

- Added the exact Shannon bound $E(\eta_d)$ from the first inequality of the steering corollary and compared it with the majorization bound. The Bell-state threshold is $\eta_{\rm S}\approx0.680568$ for the exact Shannon witness, whereas $\eta_*=(4+\sqrt2)/7\approx0.773459$ is specifically the majorization-witness threshold. For $\eta_{\rm S}<\eta_d\le\eta_*$, the exact Shannon witness certifies the Bell state while the majorization witness does not.
- Added numerical verification of the exact threshold and the two thresholds' relative positions to `verify_bounds.py`; the total remains 48 checks.
- Limited the $\Lambda\le2\sqrt c$ necessity claim to the ideal direct-sum majorization relation. It is sufficient for bounds inherited from that relation but generally not necessary for an individual entropic bound.
- Distinguished state-dependent preparation trade-offs from randomness intrinsic to unsharp POVMs or classical no-click post-processing.
- Specified probability $1/2$ for each no-click assignment and excluded the jointly zero projector pair from $\lambda_{\max}(P+Q)=1+\|PQ\|$, which does not affect the maxima.
- Updated the abstract, introduction, Discussion, README, and numerical verification accordingly. The Data Availability Statement remains unchanged.

# Tsallis application and optics clarifications

- Added a concise non-Shannon application: for every Tsallis order $\gamma>0$, the calibrated POVM bound gives $T_\gamma(p)+T_\gamma(q)\ge(2-\sum_k\omega_k^\gamma)/(\gamma-1)$. Random biased-POVM checks were added to the existing Theorem 4 verification; the suite remains at 48 checks.
- Replaced the qualitative squeezed-parity statement with the explicit approximate common-eigenstate sequence: for fixed $\alpha=i\alpha_p$ along the anti-squeezed quadrature, $\langle\Pi\rangle=1$ and $\langle D(\alpha)\Pi D^\dagger(\alpha)\rangle=e^{-2\alpha_p^2e^{-2r}}\to1$.
- Clarified that the coherent-light reference tests polarization-ensemble statistics via Stokes reconstruction, not the single-photon biased-POVM or steering tests introduced here.
- Compressed the cat-state explanation into the analytic identity $\lambda_{\max}(|\nu_1\rangle\langle\nu_1|+|\nu_2\rangle\langle\nu_2|)=1+|\langle\nu_1|\nu_2\rangle|$, from which the quoted partial-sum bound follows. The existing coherent-state common-eigenvector argument is retained.
- Added an explicit theorem-by-theorem literature-map explanation to README, distinguishing known formulas and optimality from new proofs, attainment criteria, closed-form evaluations, and physical applications.

No data-availability text or repository disclosure is changed. No dark-count, efficiency-drift, or setting-dependent-loss model is added, and the existing 78%/28% presentation is retained.

# Focused Discussion and Fig. 6 layout

- Moved the Fig. 6(b) legend below the axes in a two-column layout so that it no longer overlaps the curves, and regenerated `figures/detector_bias.pdf`.
- Removed the in-plot “ideal detectors” text annotations after visual review. The y-axis now marks $60^\circ$ and $90^\circ$ at the corresponding horizontal dashed-line heights, with matching colors.
- Removed the italic “Relation to quantum optics” heading and substantially compressed the material. The Discussion now contains one photonic-qubit paragraph covering Malus's law, coherent-light polarization-ensemble tests, and detector-efficiency mismatch.
- Deleted the former open question (v) on binary measurements on light modes.
- Removed the Robertson–Schrödinger and Fock-space statements from the numerical appendix, including the statement about a Fock space of dimension 600. Updated the Werner appendix sentence to mention both $\eta_*$ and $\eta_{\rm S}$.
- Preserved the removed quantum-optics discussion, extended photonic context, former open question (v), and numerical-methods notes in `paper/quantum_optics_future.md` for a future article. The associated Fock-space checks remain in the code for that future work, but are not presented as evidence for the focused manuscript.
- Synchronized README's Discussion summary, repository layout, and roadmap. No data-availability wording is changed.

# Five-point mathematical review (November 9)

## 1. Rényi minimizers and numerical stability

- Replaced the former angle-dependent eigenstate threshold $\alpha_0(c)$ and its numerical values by the sharp threshold $\alpha=1/2$ for every non-MUB angle. For $\alpha>1/2$, the small-rotation expansion $F_\alpha(\varphi)-F_\alpha(0)=-f'_\alpha(\theta)\varphi+o(\varphi)$ proves that an eigenstate cannot remain a minimizer. The MUB exception remains separate.
- Replaced the old grid-based threshold check with stable $\sin^2(\varphi/2)$ probabilities and logarithmically spaced small rotations. The alpha-star (bisector threshold) analysis is unchanged.
- Retained the numerical, rather than certified-global, status of the quoted direct-sum ratios and failure orders.

## 2. Literature and novelty positioning

- Added Yang–Qiao, *Geometric Optimization over Quantum State Spaces: Tight Uncertainty Relations and Resource Certification*, Phys. Rev. Lett. **137**, 160201 (2026), and contrasted its certified support-function outer approximation with the present closed-form binary vectors and attainment criteria.
- Added Tavakoli (2024) and Svegborn et al. (2025) to the discussion of steering with imperfect trusted measurements. The theorem-by-theorem literature map now includes the new results.

## 3. Direct-sum Rényi bounds for all orders above one

- Added the induced bound $H_\alpha(p)+H_\alpha(q)\ge[2/(1-\alpha)]\log_2[(\sum_k\omega_{\oplus,k}^\alpha)/2]$ for $\alpha>1$, with the AM–GM argument and attribution to Baek–Nha–Son (2019), Cor. 1.
- Stated its applicability to arbitrary biased two-outcome qubit POVMs through Theorem 4. Replaced claims that only the tensor product remains at large orders: it is the particular expression $H_\alpha(\sqrt c,1-\sqrt c)$ that fails, not all direct-sum Rényi consequences.
- Added projective and random biased-POVM validity checks to the existing corollary and Theorem 4 checks.

## 4. Calibrated steering benchmark

- Retained the exact-Shannon and majorization thresholds $\eta_{\rm S}\approx0.680568$ and $\eta_*=(4+\sqrt2)/7$, and added the calibrated linear two-setting witness for click/no-click detection.
- For the Werner family, the calibrated linear expression is $2w$ independently of efficiency, so exact-statistics detection gives $w>1/\sqrt2$ for every $\eta_d>0$. This is explicitly a comparison of witnesses, not a fundamental steering boundary.
- Added explicit Born-rule matrix checks of the linear expression at several Werner parameters and efficiencies.

## 5. Results replacing parts of the former open questions

- Added Proposition 5 for unbiased equal-sharpness qubit POVMs: $S_1=(1+\eta s)^2/4$, $S_2=(1+\eta)/2$, $S_3=1-m_\eta(s)$, and $S_4=1$, where $m_\eta(s)=(1-\eta s)^2/4$ for $\eta\le s$ and $(1-s^2)(1-\eta^2)/4$ for $\eta\ge s$. The proof uses AM–GM, the exact $S_2$ ordering, and a quadratic minimization on the ellipse. The ordered increments require no flattening.
- Added Proposition 6: for binary POVMs in any dimension, optimal direct-sum attainment is equivalent to nonzero intersection of the top eigenspaces of a maximizing effect and a maximizing cross-sum.
- Narrowed the open questions to biased/unequal-sharpness tensor products, more explicit higher-dimensional spectral characterizations and tensor products, and the remaining multi-measurement and mixed-order problems.
- Added 70 grid/random comparisons for Proposition 5, branch-continuity checks, the explicit $S_3=0.955$ example, and trivial/MUB regressions for Proposition 6.

## Author information and scope

- Added Gilad Gour's e-mail and his Technion Faculty of Mathematics / Helen Diller Quantum Center and GTIIT Mathematics with Computer Science Program affiliations. Only his funding acknowledgement remains open.
- No dark-count, drift, or setting-dependent-loss model is added. The Data Availability Statement remains unchanged and contains no repository link.

Final verification for this revision: all 48 checks pass. `main.pdf` was rebuilt and visually checked (17 pages), with no unresolved references or overfull boxes; the existing REVTeX float-placement warnings remain.
