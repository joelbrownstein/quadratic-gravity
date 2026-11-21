# Gate 0: scalaron–PS conventions

21 November 2026. Working choices for the first analytical benchmark, adopted to make packet A executable. This is a standalone note, excluded from TeX. Gate 0 is complete for a **classical, tree-level comparison**; a nonperturbative QQG measure and ghost prescription remain open research choices.

## 1. Target and order of work

First derive the scalaron of Lorentzian R+R² gravity at finite positive Planck mass. Then compare its scalar-curvature response with the restricted Euclidean PS model. Do not assume a common quantum measure or interchange removal of the Einstein term, conformal restriction, and variation.

The immediate target is classical auxiliary-field equivalence and the linearized scalar response around Minkowski space. No integration of tensor quantum fluctuations is required for packet A. The eventual target is matching specified correlators; it is not equality of entire theories.

## 2. Lorentzian benchmark

Use units hbar=c=1, signature (-,+,+,+), and

\[
R^\rho{}_{\sigma\mu\nu}=\partial_\mu\Gamma^\rho_{\nu\sigma}
-\partial_\nu\Gamma^\rho_{\mu\sigma}
+\Gamma^\rho_{\mu\lambda}\Gamma^\lambda_{\nu\sigma}
-\Gamma^\rho_{\nu\lambda}\Gamma^\lambda_{\mu\sigma},
\qquad R_{\sigma\nu}=R^\rho{}_{\sigma\rho\nu}.
\]

The parent bulk action is

\[
S_L=\int d^4x\sqrt{-g}\left[\frac{M^2}{2}R+\alpha_L R^2+\beta_L C^2\right],
\qquad M^2=(8\pi G)^{-1}>0.
\]

Set cosmological constant and matter to zero; background g=eta, R=0. For packet A set beta_L=0 and alpha_L>0. This is a scalaron benchmark, not a renormalizably closed replacement for generic QQG. Retain finite M throughout. For a later tensor comparison choose beta_L<0 and translate to the manuscript's alpha_L=1/(6 f0²), beta_L=-1/(2 f2²). No hierarchy between the resulting masses is assumed.

Use perturbations with boundary terms vanishing in the local calculation. Finite-region gravitational boundary actions require a separate specification; do not claim a boundary-value equivalence from this bulk exercise. Constant Euler terms are omitted only for this fixed-topology bulk calculation.

Dimensions: [x]=-1, [R]=[chi]=2, [M]=[canonical scalar]=1; alpha_L, beta_L, lambda and the PS field sigma are dimensionless.

## 3. Auxiliary variable and admissible frame transformation

Define the classical auxiliary action by replacing alpha_L R² with

\[
2\alpha_L\chi R-\alpha_L\chi^2.
\]

Packet A must establish chi=R by variation and derive the frame transformation, rather than borrow a potential with different conventions. Define F(chi)=M²+4 alpha_L chi and restrict the Einstein-frame construction to F>0; the proposed transformed metric is gE=(F/M²)g. Require the map, kinetic term, potential, vacuum and scalar mass to be independently checked.

This replacement is a classical stationary-action identity. A quantum auxiliary-field Gaussian identity needs a normalization and integration contour and is not asserted here.

## 4. Gauge, states and quantum limits

No gauge fixing is needed to establish the covariant classical auxiliary identity. For a subsequent flat-background propagator calculation use g=eta+h and the linear de Donder condition partial^mu h_mu_nu - (1/2) partial_nu h=0. If inverting a gauge-fixed operator, record its gauge-fixing coefficient explicitly and check that the linearized curvature response is independent of it. R^(1) is invariant under linearized diffeomorphisms about R=0; nonlinear R(x) at a coordinate point is not automatically a diffeomorphism-invariant observable.

Use the perturbative Minkowski vacuum and Feynman boundary conditions for the beta_L=0 benchmark. No contour or probability prescription is selected for the massive spin-two ghost. No nonperturbative parent metric measure is selected. These omissions block later quantum-equivalence claims, but do not block packet A.

## 5. Euclidean PS comparison model

Retain the manuscript's independently defined model:

\[
S_{PS}=\frac12\int d^4x[\Delta\sigma+\lambda(\partial\sigma)^2]^2,
\qquad g_{\mu\nu}=e^{2\lambda\sigma}\delta_{\mu\nu},
\qquad \Delta=\sum_\mu\partial_\mu^2.
\]

Thus Delta acting on exp(ipx) gives -p². On a periodic lattice use the existing real, zero-mean scalar subspace and its flat finite-dimensional measure; weight exp(-S_lat). This zero-mode convention belongs to PS and must not be imposed on the massive scalaron without justification.

The manuscript's exact restricted identity is

\[
R=-6\lambda e^{-2\lambda\sigma}[\Delta\sigma+\lambda(\partial\sigma)^2],
\qquad \alpha_E\int\sqrt g R^2=72\alpha_E\lambda^2 S_{PS}.
\]

Normalization requires alpha_E=1/(72 lambda²) for nonzero lambda. Keep alpha_E distinct from alpha_L. A positive Euclidean PS action does not establish analytic continuation from the chosen Lorentzian branch. Packet A must not silently set alpha_E=alpha_L or combine this identity with the Lorentzian f0 relation. The lambda=0 free PS theory exists independently; the displayed geometric normalization is singular there.

## 6. Observable and matching contract

Start with the linearized curvature R^(1), its sourced response and propagator pole. Introduce a scalar source through integral sqrt(-g) J_R R, and specify source transformations under any field redefinition. Derive pole position and residue with that fixed normalization; separate polynomial contact terms from poles. Compare with the PS curvature expression, not merely the sigma propagator 1/p⁴.

If a quantum calculation is later undertaken, define Z[J], connected W[J] and the Legendre-transformed Gamma[r] separately, where r=delta W/delta J. The earlier roadmap's notation Gamma[J_R] is schematic and must be replaced by W[J_R] for a source-functional comparison. Composite-operator renormalization and source contact counterterms must then be included. At tree level there is no subtraction scheme to fix; a future loop packet must specify one before reporting running.

Regime: finite-M linearized scalaron physics first. A UV PS limit and a low-energy elimination of a massive spin-two mode are distinct proposed limits and require separate demonstrations. In particular, removing the Weyl term leaves the ordinary graviton; scalar–tensor equivalence does not eliminate Einstein gravity.

## 7. Packet A: ready to give Claude

> Read this conventions sheet and sections/classical_gravity.tex and sections/model.tex. Derive the classical auxiliary-field and Einstein-frame formulations of the beta_L=0 Lorentzian benchmark step by step. Verify the auxiliary equation, the F>0 domain, kinetic normalization, potential, Minkowski vacuum and scalar mass. Independently obtain the linearized curvature equation from the original metric equations and compare its pole with the Einstein-frame result. Track boundary terms and source normalization. Explain precisely which steps would need a contour or Jacobian at the quantum level. Do not identify alpha_E with alpha_L or the PS conformal factor with the scalaron. Return a compact derivation, a notation map, independent checks and unresolved questions in notes/SCALARON_BENCHMARK_DERIVATION.md. Cite exact source sections for imported claims. Do not alter the manuscript or run simulations.

Acceptance checks: recovering the original bulk action after eliminating chi; dimensionally correct Einstein-frame quantities; vanishing potential and first derivative at the chosen vacuum; consistent scalar pole from two derivations; explicit domain of invertibility. Failure of a check must remain visible.

## Sources and evidence

Local anchors: sections/classical_gravity.tex (Lorentzian action), sections/parent_theory.tex (f0/f2 normalization), sections/model.tex (independent Euclidean PS normalization), and sections/lattice_action.tex (regulated measure). These were inspected for this gate.

Background references already in the bibliography: [Salvio, Quadratic Gravity](https://arxiv.org/abs/1804.09944) and [De Felice and Tsujikawa, f(R) theories](https://arxiv.org/abs/1002.4928). Their bibliographic records were rechecked for this gate; exact sections and calculations must be verified in packet A. The choices above define the task and are not a newly proved matching result.
