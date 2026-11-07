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

- `main.tex` — JCAP title page, abstract, section order, bibliography setup.
- `sections/model.tex` — the PS action and its relation to conformally flat quadratic gravity.
- `sections/numerics.tex` — numerical strategy and reproducibility commitments.
- `sections/results.tex` — placeholders for measured results only.
- `sections/lattice_action.tex` and `sections/observables.tex` — reproducibility appendices.
- `figures/` — source-controlled PDF, PNG, or JPG figures.
- `references.bib` — bibliography entries.

## Submission guardrails

- Replace all bracketed placeholders and draft-only text before circulation.
- Do not state that the PS model resolves the spin-two ghost problem; it does not contain the Weyl/tensor sector.
- Keep source data, figure-generation scripts, and analysis settings in version control.
- Before submission, build a clean PDF, run a spell check, verify every citation, and include the generated `.bbl` file with any JCAP submission archive.

The manuscript uses the JCAP `jcappub` package and `JHEP.bst`, in line with JCAP's author guidance.  The package is included locally so the project imports cleanly into Overleaf.
