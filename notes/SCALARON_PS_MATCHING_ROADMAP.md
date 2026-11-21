# Scalaron–PS matching roadmap

Standalone research note; excluded from the TeX manuscript.

Packet B update (21 November 2026): [SCALARON_MODE_DECOMPOSITION.md](SCALARON_MODE_DECOMPOSITION.md) derives the invariant quadratic operator, its spin-projector inverse, and the 2+5+1 spectrum. Curvature contraction reproduces Packet A. A direct local-equation test shows that conformal restriction before variation admits a spurious massless branch at finite Planck mass. Independent review and all nonlinear/quantum matching gates remain pending.

Packet A update (21 November 2026): [SCALARON_BENCHMARK_DERIVATION.md](SCALARON_BENCHMARK_DERIVATION.md) contains Codex's classical auxiliary-field derivation, Einstein-frame potential, independent trace-equation mass check, curvature source response and leading PS comparison. Internal analytical checks agree; separate expert review and the remaining quantum gates are pending.

Gate 0 update (21 November 2026): see [SCALARON_CONVENTIONS.md](SCALARON_CONVENTIONS.md). Working conventions and a ready-to-run packet A are complete for the classical finite-Planck-mass scalaron benchmark. Full quantum measure/contour choices remain open. Use the conventions sheet in preference to schematic signs and source-functional notation below; in particular distinguish Lorentzian alpha_L from Euclidean alpha_E and connected W[J] from Legendre-transformed Gamma[r].

**Status:** research plan, 21 November 2026.  This document is not part of the manuscript and makes no claim that the perfect-square (PS) scalar is the physical scalaron of quadratic quantum gravity (QQG).

## The question to make precise

Start from a declared four-dimensional parent action, for example

\[
S[g]=\int d^4x\sqrt{g}\left[-\frac{M_{\rm Pl}^{2}}{2}R+
\alpha R^{2}+\beta C_{\mu\nu\rho\sigma}^{2}+\theta E+\cdots\right].
\]

The target is not a vague identification of fields.  It is a conditional statement about specified scalar-curvature observables in a specified energy window:

\[
\Gamma^{(0)}_{\rm full}[J_R]
=\Gamma_{\rm candidate}[J_R]+\text{controlled remainder}.
\]

Here `full` includes a declared gauge fixing, measure, Euclidean contour and renormalization scheme; `candidate` may be the existing PS action or an improved scalar action.  The remainder must name both its expansion parameter and the observables for which it is bounded.  Without these choices, there is no well-posed proof to ask of a person or a model.

The central warning is structural.  Conformal flatness sets the Weyl tensor to zero.  Therefore the PS sector cannot be exactly equivalent to generic QQG with nonzero `beta` unless the omitted tensor sector is independently shown to decouple for the stated observables and scale range.  Integrating out that sector may generate nonlocal terms or additional local operators.

## Fixed terminology

* **Scalaron:** the physical spin-zero curvature degree of freedom in a specified `R+R^2` or `f(R)` theory.
* **Conformal factor:** the field `sigma` in `g_mu_nu = exp(2 lambda sigma) gbar_mu_nu`; it is a parametrization, not automatically the scalaron.
* **PS scalar:** the conformally flat restricted theory with
  \[
  S_{\rm PS}=\frac12\int[\Box\sigma+\lambda(\partial\sigma)^2]^2.
  \]
* **Matching:** equality of selected renormalized observables up to a declared, tested error.  A shared classical expression is not quantum matching.

## Decision gates

Each wave must produce a retained derivation and a list of assumptions.  Do not begin the next wave merely because a language model supplied fluent text.

| Wave | Core question | Required output | Gate to continue |
|---|---|---|---|
| 0 | What is the parent theory? | One-page conventions sheet | Human chooses action, signature, background, gauge, contour, regulator, and match observables. |
| 1 | What is the scalaron when `beta=0`? | Auxiliary-field derivation and source checks | Independent hand check of signs, dimensions, and background solution. |
| 2 | Which modes exist in the parent theory? | Linearized York-decomposition spectrum | Gauge-invariant scalar pole and tensor poles identified in one convention. |
| 3 | What exactly does PS retain? | Conformal reduction of every parent operator | The Einstein-induced two-derivative scalar term and all omitted Weyl information are explicit. |
| 4 | Can the tensor sector be removed? | EFT power counting and an integration-out calculation | A mass hierarchy, contour prescription, and locality range are stated; otherwise record failure. |
| 5 | Does a candidate scalar action match? | Correlator and running comparison | Agreement or disagreement is demonstrated for a frozen observable set. |
| 6 | Does the lattice theory support it? | Regulator, counterterm, and finite-volume program | The proposed lattice calculation has an observable tied to the matching claim. |

## Wave 0 — conventions before algebra

Prepare `SCALARON_CONVENTIONS.md` before asking Claude for long derivations.  It must fix:

1. Lorentzian or Euclidean starting action, curvature-sign convention, and metric signature.
2. Whether the target is the Wilsonian action, 1PI effective action, or connected curvature correlators.
3. Background geometry: initially flat Euclidean space is the simplest comparison point; a constant-curvature background is a separate task.
4. Gauge fixing, ghosts, treatment of the constant conformal mode, and integration contour.
5. Operator normalization for `R^2`, `C^2`, and the Euler density.
6. The observable set: start with the scalar curvature two-point function, scalar pole and residue, then selected three- and four-point functions.

**Human decision:** the desired physical regime.  A low-energy Einstein-plus-scalaron EFT, a high-momentum four-derivative regime, and a proposed UV fixed point are different matching questions.

## Wave 1 — the clean scalaron benchmark

First set `beta=0`.  Replace the curvature-square term by an auxiliary field `chi`, schematically

\[
\alpha R^2\;\longleftrightarrow\;2\alpha\chi R-\alpha\chi^2,
\]

with every coefficient derived in the chosen convention.  Then perform the metric rescaling that produces an Einstein-frame scalar–tensor action, if its assumptions are satisfied.

**Claude packet A — derivation audit**

> Work from the attached convention sheet only.  Derive the auxiliary-field representation of the specified `R+R^2` action, vary with respect to the auxiliary field, and substitute it back.  Derive the Einstein-frame kinetic term and potential.  Track boundary terms, signs, dimensions, domain restrictions, and every assumption of the Weyl rescaling.  Return a line-numbered LaTex derivation, a table of conventions, and a list of claims that need primary-source verification.  Do not call the conformal factor the scalaron without demonstrating the field relation.

**Deliverable:** a short derivation which can be checked by hand and against two primary or standard reference sources.  This is the correct first place to use formal symbolic algebra or proof-assistant checks for identities.

## Wave 2 — spectrum and constraints

Linearize the chosen parent action about the declared background and use a York decomposition.  Separate gauge variables and constraints from propagating fields.  Determine the scalar curvature pole and residue in the `beta=0` benchmark; only then include the `C^2` sector and its spin-two pole.

**Claude packet B — adversarial spectrum check**

> Starting from the attached action and conventions, independently derive the quadratic fluctuation operator and count physical degrees of freedom.  Give the scalar and spin-two propagating poles in the declared normalization.  Search for sign, gauge, and boundary-condition ambiguities.  Attempt to disprove any claimed scalar-only reduction by identifying couplings to the Weyl/tensor sector.  Cite sources for each nontrivial statement and label direct derivation separately from literature summary.

**Gate:** a human can point to the gauge-invariant scalar observable that will later be matched.  The model output is not evidence until its algebra is independently checked.

## Wave 3 — compare the PS restriction honestly

Set

\[
g_{\mu\nu}=e^{2\lambda\sigma}\delta_{\mu\nu}
\]

and reduce the complete parent action, not only `R^2`.  The `R^2` term yields the existing PS square after normalization.  The Einstein–Hilbert term yields a two-derivative contribution schematically proportional to

\[
M_{\rm Pl}^2\int e^{2\lambda\sigma}(\partial\sigma)^2,
\]

while the Weyl term vanishes on this restricted ansatz.  The latter fact is a loss of dynamical information, not a proof of decoupling.

**Claude packet C — operator inventory**

> Reduce each operator of the specified parent action under the conformally flat ansatz.  Produce a table with: parent operator; exact reduced expression; derivative order; symmetry; coefficient; whether it is retained by the PS action; and the physical information lost by the restriction.  Use integration by parts only with stated boundary conditions.  Search specifically for relevant two-derivative and cosmological-constant operators that the PS action omits.

The output should motivate one of three candidate actions: the original PS theory, PS plus the relevant operators demanded by matching, or the auxiliary scalaron action.  Do not select one by aesthetic preference.

## Wave 4 — conditional EFT elimination of the tensor sector

An integration-out calculation is justified only if there is a stated hierarchy, such as `E << M_2`, and a defined contour/physical prescription for the spin-two sector.  Its generic output is an effective action containing higher operators and potentially nonlocal kernels.  A local scalar action can be adequate only after an expansion and error estimate are established.

**Claude packet D — matching-or-obstruction search**

> Assume the supplied hierarchy and action.  Derive the tree-level effective scalar action obtained after eliminating the designated heavy field, keeping the first requested order in `E/M_2`.  Identify all nonlocal terms and the conditions under which they can be expanded locally.  If the hierarchy or contour is insufficient, give a precise obstruction rather than inventing a local action.  Compare the resulting operator basis with PS and list missing operators.

## Wave 5 — observable matching and renormalization

Freeze a minimal comparison set before calculating:

\[
\langle R R\rangle_{\rm ren},\qquad
\text{scalar pole and residue},\qquad
\langle R R R\rangle_{\rm ren},\qquad
\text{one declared running coupling}.
\]

Match at a scale `mu` in one scheme.  Compare not just the action form but the values and running of renormalized quantities.  The desired statement is conditional, for example “agreement through order `p^4/M_2^4` for these correlators,” rather than “PS proves QQG.”

**Claude packet E — falsification-first comparison**

> Given the two actions, conventions, scheme, and observable set, derive the lowest-order predictions on both sides.  Seek the first observable or operator at which equivalence fails.  Report the calculation as a table: observable; PS result; candidate-scalar result; full-theory result or EFT expansion; matching coefficient; remainder; required verification.  Do not infer equality outside the stated order and regime.

## Wave 6 — consequence for the lattice project

Only after Waves 1–5 should the lattice action be changed.  The current PS regulator remains a useful independent nonperturbative study.  A matching-driven extension may require an additional two-derivative operator, a cosmological term, and counterterms, each with a clear continuum observable.  The immediate numerical prerequisite remains the existing HMC validation at fixed finite volume.

For any candidate action, preserve:

* exact definition of the zero-mode treatment and measure;
* free-theory covariance test where available;
* force, reversibility, and acceptance checks for HMC;
* independent seeds, autocorrelation analysis, and retained failed runs;
* a direct link from every measured lattice observable to a matching target.

## Evidence and provenance protocol

Claude can search, propose algebra, generate code and challenge derivations.  It cannot be cited as scientific authority.  For every task retain the prompt, model/version, references, symbolic notebooks or scripts, assumptions, intermediate expressions, and an independent check.  Use a proof assistant only for a fully specified mathematical identity; it certifies the encoded statement, not whether the encoded action, contour, or physical interpretation is correct.

The useful success outcomes are deliberately asymmetric:

* **Positive:** a controlled scalar effective action and an explicit matching range are found.
* **Negative:** PS fails a stated correlation-function or operator-basis test; this precisely defines its scope.
* **Inconclusive:** no justified hierarchy or contour is available; record the missing physics rather than replacing it with a speculative action.

## Immediate next task

Do not ask Claude to solve generic QQG in one prompt.  Complete Wave 0, then run packet A on the `beta=0` scalaron benchmark.  Its deliverable is a compact, checkable derivation and a proposed notation map between the auxiliary scalar, the scalar curvature, and the PS conformal field.
