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
