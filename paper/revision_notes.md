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
