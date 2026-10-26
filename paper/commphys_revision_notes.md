# Communications Physics revision: what changed and why

Branch `commphys-revision` (local, not pushed). The PRA-style manuscript `paper/main.tex` is unchanged; the new version is `paper/commphys.tex` (compiles with `pdflatex commphys && bibtex commphys && pdflatex commphys && pdflatex commphys`).

## 1. Format (Communications Physics Article)

- Structure: unheaded abstract (143 words, no references), Introduction / Results (subsections) / Discussion / Methods. All proofs moved to Methods, which does not count towards the ~5,000-word main-text guideline.
- Nature-style superscript citations (`natbib` + `unsrtnat`; swap in `naturemag.bst` if available), line numbers and 1.5 spacing for review.
- Six figures (≤10 display items). Every legend starts with a bold title sentence.
- New end statements: Data availability, Code availability, Author contributions, Competing interests (the last three contain TODOs for the authors).
- New title: *Majorization uncertainty bounds for incompatible binary measurements cannot be saturated* (the old title is a fine fallback).

## 2. New scientific content (the review's main request: operational significance)

**Proposition 4, Corollary 1 and Fig. 6 (new subsection "Imperfect detectors and uncertainty-based tests").**

- (i) Unbiased noise is a doubly stochastic map on the statistics, so every bound derived for ideal detectors remains valid.
- (ii) Biased measurements obey the ideal relation p⊕q ≺ (1, √c, 1−√c, 0) for all states **iff Λ ≤ 2√c**. Otherwise a state violates the ideal Landau–Pollak sum. This follows from Theorems 1 and 4 together with the lattice-supremum property.
- (iii) When both settings are read out by the same detectors, the condition becomes |μ| ≤ √c(1−η). A click/no-click detector violates it for every η_d < 1.
- Shannon level (numerics, Fig. 6):
  - Random no-click assignment, η_H = 0.99: the ideal relation fails for η_V < 0.942, and the 1-bit Shannon bound fails for η_V < 0.578.
  - η_H = 0.95: the thresholds are 0.708 and 0.211.
  - Click/no-click detection at 90°: the minimum is h(η_d/2) < 1 (analytic).
- Corollary 1 (steering): for a local-hidden-state model, H(M|A₁)+H(N|A₂) ≥ min[H(M)+H(N)] ≥ H(ω⊕(M,N)). With biased detectors, the ideal bound can certify steering for a *product* state.
- Also new: the spectrum-dependent bound for states with |r| = r follows from Theorem 4 with η = ζ = r, and Proposition 3 shows it is not attained. This answers review comment 5.

## 3. Review comments → changes

| Comment | Change in `commphys.tex` |
|---|---|
| M1 Partovi's priority on non-attainment | Introduction credits Partovi (general remark, MUB tensor vector) and the 2026 Conlon et al. review. Theorem 2 cites Partovi. After Proposition 1 the text says the proof is short once the saturating states are known. |
| M2 attributions in Secs. VI–VII | KTW credited for c = ½ + ¼‖{A,B}‖ in any dimension. Yuan et al. credited for the general optimality of the flattened direct-sum vector. Abbott et al. cited at Theorem 4. The dimension-d text separates assembled parts from new ones. |
| M3(a) Rényi direct sum beyond α = 1 | New paragraph: valid for all α ≤ 2 at every angle (numerical). Analytic at α = 2: (√c−1)³(√c+3) ≤ 0, proof in Methods. Fails above an angle-dependent order between 2.00 and 2.58. Fig. 4 legend updated. |
| M3(b) eigenstate range | Eigenstates are optimal for α ≤ α₀(c) with α₀ > ½ (0.52, 0.56, 0.66 at 20°, 60°, 85°). |
| M4 operational meaning | Section 2 above. The Discussion also notes that Wang et al.'s "touching" Lorenz curves can only touch individual partial sums, and proposes min_ρ max_k[S_k − L_k] as a measurable quantity. |
| M5 fixed spectrum | Done (section 2). |
| m1 "22%" in Discussion | Now "never below 78%" / "exceeds by up to 28%". The 22% wording appears only where it is correct. |
| m2 bound (19) wording | "Strictly stronger than h(√c), at least as strong as MU, equal at 90°" (text and Fig. 3 legend). |
| m3 tensor ratio | Added: h(w) ≥ 0.761 × the exact minimum, worst at 70.6°. |
| m4 Proposition 1(ii) proof | Replaced by the two-vanishing-entries argument. |
| m5 sign condition in Proposition 3(i) | Completed. |
| m6 max{1, ·} in Corollary 3 | Removed, with justification. |
| m8 common eigenvectors | Now "unless both ranks equal d/2". |
| m9 pseudospin / squeezed vacuum | The squeezed-vacuum claim was dropped. Only the requirement ‖{A,B}‖ < 2 remains. |
| m10 citations | Luis–Bosyk–Portesi 2016, Puchała et al. 2018, Conlon et al. 2026, Yang–Qiao 2026, Walborn 2011, Schneeloch 2013, Uola 2020 added to `references.bib`. |
| m11 format | Section 1. |
| m12 code availability | Statement added; repository URL or DOI left as a TODO. |

## 4. Code

- `majorization.py`: `big_lambda`, `ideal_relation_holds`, `detector_pair`, `click_detector`, `shannon_min_in_plane`.
- `verify_bounds.py`: 8 new checks covering Proposition 4 (i)–(iii), Fig. 6, the detector thresholds, the α ≤ 2 validity, α₀(c), the tensor ratio and the fixed-spectrum bound. **46/46 PASS.**
- `make_figures.py`: `fig_detector_bias()` writes `figures/detector_bias.pdf`. Without `type1ec.sty` (cm-super), the script falls back to Matplotlib's built-in Computer Modern. The copy of `detector_bias.pdf` in the repository was made that way, so **re-run `python make_figures.py` on the group workstation** to get LaTeX-typeset text matching the other figures.

## 5. For the authors to decide (TODOs in the PDF are red)

1. Gilad's affiliation and funding; corresponding author.
2. Author contributions and the competing-interests statement.
3. Code availability: public repository or Zenodo DOI vs "on request". The review recommended a public archive; Marcelo's earlier choice was "on request".
4. Whether to keep the new title.
5. Whether the detector-bias result (Proposition 4, Corollary 1) is the operational message you want to lead with. If yes, consider a short experimental proposal with a photonics collaborator. Reviewers at this journal will still weigh significance, and an experiment would strengthen the case most.
