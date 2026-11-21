# Weyl-sector ringdown: spectrum, sourced response, then rotation

21 November 2026. Research roadmap; excluded from the TeX manuscript. W0 and W1 are now completed as working calculations, including a GR axial benchmark solver. No massive-mode spectrum, sourced response, waveform fit or coupling bound has been completed. Later gates remain proposed tasks. This work is separate from the scalaron–PS packets.

## Objective and primary sources

Build a reproducible classical Einstein–Weyl response benchmark before attempting rotating remnants or an observational bound on the parent QQG coupling. The required chain is **action → constrained perturbations → spectrum → excitation and propagation → detector likelihood → coupling matching**.

- **Spectrum:** Antoniou, Gualtieri and Pani, *Gravitational quasinormal modes of black holes in quadratic gravity*, PRD **111**, 064059 (2025), [arXiv:2412.15037v2](https://arxiv.org/html/2412.15037v2), [DOI](https://doi.org/10.1103/PhysRevD.111.064059). Read II, IV.1 and appendices B–C for the Schwarzschild axial benchmark. The non-Schwarzschild branch is a later, separate target.
- **Response:** the same authors, *Nonperturbative suppression of beyond-general-relativity effects in quadratic gravity*, PRD **113**, 124068 (2026), [arXiv:2604.01631v2](https://arxiv.org/html/2604.01631v2), [DOI](https://doi.org/10.1103/8r69-19zy). Read II.4–II.8, III and appendices A–E for polar radial-infall perturbations, numerical methods and fluxes. Distinguish this classical suppression result from quantum nonperturbative dynamics.
- **Stability:** Held and Zhang, [arXiv:2209.01867](https://arxiv.org/abs/2209.01867), PRD **107**, 064060 (2023). Stability must be checked separately from finding damped nonradial modes.

Local anchors: `SCALARON_CONVENTIONS.md`, `SCALARON_MODE_DECOMPOSITION.md`, and manuscript sections `classical_gravity`, `parent_theory`, `black_holes`. Record exact source versions and equation/figure identifiers in each delivered calculation. Independently derive suspect formulae rather than copying rendered equations without checking the PDF and conventions.

## Gate W0 — conventions and retained sectors

Deliver `WEYL_RINGDOWN_CONVENTIONS.md` before implementing a solver.

Use the existing signature and curvature convention, and write

\[
S=\int\sqrt{-g}\left[\frac{M_{\rm Pl}^2}{2}R+\alpha_LR^2+\beta_LC^2\right]+S_m,
\quad M_{\rm Pl}^2=(8\pi G)^{-1},\quad
\beta_L=-\frac{1}{2f_2^2}.
\]

The initial benchmark sets **alpha_L = 0 explicitly**. It retains the massless graviton and the massive spin-two sector, with classical retarded boundary conditions. It does not select a quantum ghost prescription. Do not infer absence of scalar perturbations merely from background R=0. Reintroducing finite alpha_L and a traceful source needs its own scalar equation and matching check.

Translate the papers' dimensionful coefficient in `(M_Pl²/2)[R − a_W C²]`:

\[
a_W=-2\beta_L/M_{\rm Pl}^2=1/(f_2^2M_{\rm Pl}^2),
\quad m_2^2=1/(2a_W),\quad
r_g=GM_f,\quad r_h=2r_g,\quad x_2=r_gm_2.
\]

Reserve `m2_pole` for the pole mass, `M_BH` for physical remnant mass, and `r_g` for geometric mass in code. The papers' spin-two mass symbol must not be confused with the renormalization scale. Adopt time dependence exp(−i omega t); a decaying mode has Im(omega)<0. State the tortoise coordinate and frequency normalization before comparing tables. For a future detector model, account for redshift between source and observer frequencies.

Acceptance: recover the flat-space pole mass from the original action and auxiliary formulation; verify all units and factors of two; distinguish a_W from the manuscript's R² coefficient; list retained constraints and the reason for any omitted field.

## Gate W1 — recover GR and derive the axial system

Deliver a derivation note and a small, reproducible solver with documented input/output units.

1. Start on Schwarzschild with the auxiliary background field set to its branch value. Derive the perturbation equations and constraints before eliminating variables.
2. Separate axial harmonics, treating exceptional low multipoles separately. Identify the GR massless sector and massive spin-two degrees of freedom without identifying helicity-zero spin-two content with the scalaron.
3. Recover the Regge–Wheeler benchmark. For the l=2 fundamental, verify `r_g omega ≈ 0.37367 − 0.08896 i`; a finer comparison should use a cited high-precision reference.
4. Derive horizon and far-zone series. At infinity the massive channel involves `k = sqrt(omega² − m2_pole²)`; document the analytic branch used for outgoing conditions. Quasinormal resonances and decaying-at-infinity quasibound states are different boundary problems.

Acceptance: constraint residuals converge toward zero; the GR result is insensitive to numerical boundaries and expansion order; the solver explicitly identifies which boundary-value problem each root solves.

## Gate W2 — reproduce the Schwarzschild massive spectrum

Use the 2025 axial calculation as the first literature comparison. Choose named published points and mode families before running; store a table containing the source location, x2, complex frequency, computed value and numerical error. Do not use agreement with a plotted curve as the only evidence.

Use two independently formulated methods where practical, such as continued fractions and boundary matching. Track roots continuously in x2 so that changes of mode ordering are not mistaken for physical discontinuities. Check precision, truncation, radial matching and asymptotic-order dependence. Proposed acceptance targets are relative error below 10^-4 in resolved benchmark frequencies and agreement within documented convergence error; near zero damping or instability use absolute errors and an explicitly revised tolerance.

Identify the stable parameter range for the chosen use. Include the radial stability question in the ledger even though an axial solver cannot settle it. A growing mode must not be discarded because a damped-ringdown template is desired. If a transient calculation uses an unstable background, justify its duration relative to the growth time.

## Gate W3 — derive and validate the polar system

Radial infall does not supply the axial source needed to excite the W2 modes. Derive the polar equations independently before reproducing the 2026 response. Resolve l=0, l=1 and l>=2 constraints separately. Compare appropriate vacuum polar frequencies with the 2026 appendix B; do not impose GR axial–polar isospectrality on the extra sector without proof.

For the source, specify particle energy, mass-ratio limit, trajectory, harmonic normalization and Fourier convention. Begin with the paper's ultrarelativistic radial-infall benchmark, stating which quantity is held fixed in the limit. Check conservation of the point-particle stress tensor on the background. Before moving to finite-energy massive particles in generic R²+C² gravity, explicitly revisit trace sourcing of the scalaron.

Acceptance: source-free equations reproduce W0 masses and constraints; reconstructed source obeys conservation; exceptional multipoles do not rely on divisions that vanish for their l values. Any inconsistent reduction remains a blocking issue for the response stage.

## Gate W4 — reproduce a sourced waveform and suppression

Construct a retarded matrix Green function or independently validated time-domain evolution. Retain source normalization and both the metric and auxiliary-field response. Reproduce at least one identified spectrum and one waveform from the 2026 paper before surveying new parameters.

Vary frequency resolution, integration range, precision and extraction radius. Explain how each asymptotic component propagates: a massive channel is not automatically an ordinary null-infinity strain observable. Reconstruct the metric quantity relevant to a detector, and separate it from an auxiliary-field amplitude. Check the signed energy flux appropriate to the classical ghost sector rather than silently replacing every contribution by a positive GR flux.

To test the reported suppression, scan a justified x2 interval and compare convergence-controlled amplitudes or energies with candidate asymptotic behavior. State which observable is fitted. Demonstrate that an apparent exponential decline lies above numerical noise and is not caused by a frequency cutoff or lost roots. Quote a fitted range and uncertainty; do not extrapolate a numerical fit into an unrestricted theorem.

Acceptance: GR limit for the resolved massless response; source normalization and units verified; matched benchmarks within a predeclared numerical tolerance; independent convergence checks; a visible noise floor. This delivers a sourced Schwarzschild benchmark, not an astrophysical merger prediction.

## Gate W5 — extend to rotation

Proceed only after W1–W4 pass. Start with a controlled slow-rotation expansion on the Kerr branch, deriving angular/parity mixing and constraints. Recover both the Schwarzschild limit and the GR rotating limit. Check the rotation expansion against a higher order or an independent method before applying it to a remnant with spin near 0.68; such spin is not automatically within a first-order approximation's domain.

Then assess whether a two-dimensional frequency-domain solver or time-domain evolution is required. Derive the horizon condition with omega−m Omega_H, investigate massive-mode trapping and possible instabilities, and calculate excitation for a specified plunge. Comparable-mass merger amplitudes require an additional nonlinear calculation or a justified calibration. A rotating 2-2-hole is a separate background project; do not substitute a reflective boundary into a Kerr solver and label the result a QQG prediction.

## Gate W6 — observational inference and quantum interpretation

Build strain templates only after the spectrum and transfer to the detector are known. Fit mass, spin, x2, excitation amplitudes, phases, orientation and ringdown start time jointly. Use injection–recovery checks to expose degeneracies, model mismatch and coverage problems. A free extra-mode amplitude that can vanish generally prevents nondetection alone from excluding its coupling. Any informative excitation prior needs physical justification.

Translate a conditional pole-mass posterior with `f2 = sqrt(2) m2_pole/M_Pl` only at the stated classical/tree-level accuracy. Do not presume a one-sided bound: inaccessible frequencies, suppressed amplitudes and branch structure can make sensitivity nonmonotonic. A bound on running f2 needs a scheme, scale, radiative matching and a specified physical-state prescription. Existing phenomenological LVK deviation posteriors cannot be relabeled as this coupling posterior.

## Deliverables and stopping rules

Each gate should produce a concise derivation, machine-readable benchmark data where relevant, numerical error evidence, source locations and an unresolved-issues list. Keep working calculations in notes and a separately named analysis directory; do not add exploratory results to the manuscript as validated findings. No computing allocation, detection, stability proof or coupling bound is implied by this plan.

W0 completed on 21 November 2026: see [conventions and retained-sector derivation](WEYL_RINGDOWN_CONVENTIONS.md). W1 completed the same day: see [axial derivation and GR benchmark](WEYL_RINGDOWN_AXIAL_DERIVATION.md), the solver in `scripts/weyl_ringdown/`, and results in `notes/validation/`. These remain working calculations for scientific review. Next executable task is **W2**, reproducing the massive Schwarzschild spectrum with independent methods and explicit stability checks. No sourced response has been computed. Stop dependent work if a required boundary condition, retained sector or source normalization is unresolved. The existing PS program remains independently useful but cannot determine these Weyl-mode spectra, excitation amplitudes or f2 bounds.
