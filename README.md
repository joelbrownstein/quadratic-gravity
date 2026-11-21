# Perfect-square sector of quadratic gravity

This repository is an Overleaf-ready JCAP manuscript scaffold for a numerical study of the perfect-square (PS) conformal sector of quadratic gravity.

## Start on Overleaf

1. In Overleaf, choose **New Project → Upload Project**.
2. Upload a ZIP archive of this repository, or connect the GitHub repository and select `main.tex` as the root document.
3. Compile with pdfLaTeX. The repository includes the JCAP style file (`jcappub.sty`) and bibliography style (`JHEP.bst`), so it does not depend on a nonstandard Overleaf package.

## Local build

```sh
latexmk -pdf main.tex
```

The first LaTeX run requires BibTeX; `latexmk` handles the required sequence.  Build products are ignored by Git.

## Manuscript map

Start with [LITERATURE_MAP.md](LITERATURE_MAP.md) for the claim-to-source evidence map and retrieval limits. Read [MANUSCRIPT_DESIGN.md](MANUSCRIPT_DESIGN.md) for the structural audit and staged narrative design.

1. `sections/introduction.tex` — motivation, author perspective, and guide to the argument.
2. `sections/classical_gravity.tex` — classical curvature-squared operators, modes and Newtonian response.
3. `sections/parent_theory.tex` — QQG renormalization, independent couplings, running-coupling debate and physical states.
4. `sections/model.tex` — controlled PS scalar model, geometric restriction, perturbative anchors and scope.
5. `sections/numerics.tex` — proposed sampling, validation and evidence protocol.
6. `sections/predictions.tex` — predictions, decision criteria and first numerical deliverable.
7. `sections/black_holes.tex` — Euclidean geometry, solution branches, stability and singularity-resolution criteria.
8. `sections/cosmology.tex` — no-boundary probabilities, arrows of time, the nonrotating interior/cosmology bridge, infrared matching and dated observational tests.
9. `sections/ai_exploration.tex` — literature-backed geometry–growth frontier, corrected learned sampling, calibration and matched-cost comparisons.
10. `sections/conclusions.tex` — present scope and staged outlook.

Appendices define the candidate lattice action and finite-lattice bound, observables and step scaling, and cosmological identities. `references.bib` contains the literature anchors. `BLACK_HOLE_READING_LIST.md` provides an annotated specialist reading order. `HARTLE_COSMOLOGY_READING_LIST.md` records the no-boundary and arrow-of-time follow-up, its evidence levels, and Bond queries.

## Deterministic regulator checks

```sh
python3 scripts/validate_scalar_regulator.py
```

This standard-library-only script checks the proposed action derivative against finite differences, shift symmetry, discrete integration by parts, the free bilaplacian eigenvalue, and finite-lattice inequalities. It does not run HMC or produce interacting ensembles. `regulator-validation.json` records the development check; rerun the script after changing the stencil.

## Submission guardrails

- Keep the abstract and conclusions consistent with the research-design status until interacting measurements exist. No cosmological fits or production ensembles have been run.
- Do not state that the PS model resolves the spin-two ghost problem; it does not contain the Weyl/tensor sector.
- Treat the PS scalar calculation as the first independently useful computational vertical slice; infrared/cosmological matching and the Weyl/tensor sector are later stages with separate resource needs.
- Do not equate a smooth Euclidean horizon, bounded Ricci scalar, or horizon avoidance with resolution of a Lorentzian singularity. The PS sector excludes Schwarzschild’s nonzero Weyl curvature.
- Keep the cosmological discussion prospective: the PS ensembles do not yet predict dark energy, and an f(R) control model is not a derived QQG prediction.
- Record AI assistance and verify citations and numerical claims independently; no company participation or computational allocation is implied by the methods discussion.
- Keep source data, figure-generation scripts, and analysis settings in version control.
- Before submission, build a clean PDF, run a spell check, verify every citation, and include the generated `.bbl` file with any JCAP submission archive.

The manuscript uses the JCAP `jcappub` package and `JHEP.bst`, in line with JCAP's author guidance.  The package is included locally so the project imports cleanly into Overleaf.

`AI_RESEARCH_MAP.md` records the focused Bond search, primary-source checks and AI experiment design.
