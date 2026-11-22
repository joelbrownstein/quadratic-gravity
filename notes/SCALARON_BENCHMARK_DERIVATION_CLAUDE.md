# Packet A: independent derivation and review of the scalaron benchmark

21 November 2026 (America/Denver). Derived by Claude from the public `main` branch of `joelbrownstein/quadratic-gravity`, commit `2b204f1`, for Joel's review. This is an analytical working note, excluded from TeX. It supplements and does not replace Codex's [`SCALARON_BENCHMARK_DERIVATION.md`](SCALARON_BENCHMARK_DERIVATION.md); equation tags (A·) below refer to that note. Every algebraic step was checked symbolically with SymPy 1.14 in [`check_scalaron_benchmark.py`](../scripts/check_scalaron_benchmark.py) (31 checks, all pass). That is a consistency check, not peer review or proof-assistant certification. No simulation was run and no manuscript file was altered.

## 0. Summary

The independent derivation **agrees with every result in Codex's note**: the auxiliary identity, the F>0 chart, the canonical normalization, the potential, the Minkowski vacuum, m0² = M²/(12α_L), and the curvature response with its contact term and pole coupling c_R² = m0²/(2α_L).

It adds five things:

1. **Off-shell Einstein-frame check.** The full linear curvature response, including the contact term, follows from the Einstein frame. Codex checked only the on-shell residue there (§4).
2. **Third route.** Restricting to the conformal mode *before* variation reproduces the same tree-level curvature response, but not the conformal-mode propagator (§4.3).
3. **Euclidean sign map.** Under the standard Wick rotation, the Euclidean R² coefficient is −α_L. The sheet's positive α_E therefore corresponds to α_L < 0, the tachyonic flat-vacuum branch of this benchmark, unless another contour is chosen (§5.1).
4. **Gaussian consistency with PS.** With that sign map, the ultraviolet (p ≫ m0) or M→0 limit of the benchmark's Euclidean curvature response equals the PS Gaussian curvature correlator 36λ² when α_E = 1/(72λ²). This is a statement about a contact term, not an identification of fields (§5.2).
5. **One qualification of Codex's note.** The Euclidean auxiliary-field Gaussian *does* converge for real χ on the healthy branch α_L > 0. It needs a rotated contour on the PS sign (§6).

## 1. Setup imported from the repo

| Item | Source |
|---|---|
| Units ħ=c=1, signature (−,+,+,+), Riemann and Ricci convention, M² = (8πG)⁻¹ > 0 | `notes/SCALARON_CONVENTIONS.md` §2 |
| Lorentzian bulk action S = ∫√−g [M²R/2 + α_L R² + β_L C²]; packet A takes β_L = 0, α_L > 0, Λ = 0, no matter | Conventions §2; `sections/classical_gravity.tex`, eq. `qg-cosmological-action` (there with M_Pl, Λ, α, β) |
| Auxiliary replacement α_L R² → 2α_L χR − α_L χ², F = M² + 4α_L χ, g_E = (F/M²)g on F > 0 | Conventions §3 |
| Translation α_L = 1/(6f0²), β_L = −1/(2f2²), on the branch α > 0, β < 0 | `sections/parent_theory.tex`, eq. `parent-couplings` |
| Gauge condition for later propagator work; R^(1) invariant under linearized diffeomorphisms | Conventions §4 |
| PS action ½∫[Δσ + λ(∂σ)²]², g = e^{2λσ}δ, identity √g R² = 36λ²[…]², α_E = 1/(72λ²) | Conventions §5; `sections/model.tex`, eqs. `ps-action`, `conformal-metric`, `conformal-r-squared` (§`cosmological-action`) |
| Source coupling ∫√−g J_R R; contact terms separated from poles | Conventions §6 |

**Notation clash to fix in the manuscript.** `model.tex` writes the Euclidean Laplacian as `\Box`, while `classical_gravity.tex` and this note use □ for the Lorentzian d'Alembertian. The sheet's Δ = Σ_μ ∂_μ² is unambiguous; adopting it in `model.tex` would prevent sign errors (Δ e^{ipx} = −p² e^{ipx}).

The curvature sign was checked on a sphere: the sheet's convention gives R = +2 on a unit two-sphere (check 0).

## 2. Classical auxiliary-field form

Write α = α_L in §§2–4. Start from

\[
\mathcal L_J=\tfrac{M^2}{2}R+\alpha R^2,\qquad
\mathcal L_{\rm aux}=\tfrac{M^2}{2}R+2\alpha\chi R-\alpha\chi^2
=\tfrac12F(\chi)R-\alpha\chi^2 .
\]

1. **Completing the square (checks 2a, 2c):** 𝓛_aux = 𝓛_J − α(χ−R)².
2. **Auxiliary equation (check 2b):** ∂𝓛_aux/∂χ = 2α(R−χ) = 0, so χ = R for α ≠ 0. Substituting returns 𝓛_J. Because χ enters algebraically, no boundary term arises in this step.
3. **Status.** This is a classical stationary-point identity at fixed metric. Its quantum counterpart is in §6.

## 3. Einstein frame, step by step

**3.1 Frame map.** On F > 0 define Ω² = F/M², u = ln Ω and g_E = Ω²g. In four dimensions

\[
\sqrt{-g}=\Omega^{-4}\sqrt{-g_E},\qquad
R[g]=\Omega^{2}\big[R_E+6\Box_Eu-6(\nabla_Eu)^2\big].
\]

The curvature formula was checked by computing both sides for g_E = diag(−1, a², a², a²) with arbitrary a(t) and u(t, x) (check 1).

**3.2 Curvature term.**

\[
\sqrt{-g}\,\tfrac12FR=\sqrt{-g_E}\Big[\tfrac{M^2}{2}R_E+3M^2\Box_Eu-3M^2(\nabla_Eu)^2\Big].
\]

The 3M²□_E u term is the divergence ∇_E·(3M²∇_E u). It integrates to a flux that is dropped only under the boundary assumption of §7.

**3.3 Kinetic normalization (check 3b).** With u = ½ ln(F/M²) and

\[
\phi=\sqrt{3/2}\,M\ln(F/M^2),\qquad b=\sqrt{2/3}/M,\qquad F=M^2e^{b\phi},
\]

one finds −3M²(∇u)² = −½(∇φ)². Check 3a confirms the inversion; it requires M > 0, as the sheet does.

**3.4 Potential (check 3c).** −√−g αχ² becomes −√−g_E M⁴αχ²/F². With χ = M²(e^{bφ}−1)/(4α),

\[
V(\phi)=\frac{M^4}{16\alpha}\big(1-e^{-b\phi}\big)^2 .
\]

**3.5 Domain.** dφ/dχ = 4α√(3/2) M/F is finite and nonzero on F > 0 (check 3h). The chart χ ∈ (−M²/4α, ∞) maps bijectively onto φ ∈ ℝ. The boundary F → 0⁺ lies at φ → −∞, where V grows like e^{−2bφ}, so it is at infinite field distance. The region F ≤ 0 is not covered by this chart.

**3.6 Vacuum and mass (checks 3d–3g).**
- V(0) = V′(0) = 0, with χ = R = φ = 0, F = M², and g_E = g = η.
- V″(0) = M²/(12α) ≡ m0², which equals f0²M²/2 under α = 1/(6f0²).
- Linearizing gives (□_η − m0²)φ = 0. For α > 0 the scalar has the healthy kinetic sign and is not tachyonic.
- For α < 0 the same chart gives V ≤ 0 with a maximum at φ = 0, so m0² < 0: the tachyonic branch.

## 4. Independent pole derivations and source normalization

**4.1 From the original metric equations.** Varying ∫√−g 𝓛(R) with 𝓛 = M²R/2 + αR² + JR, and using δR = R_μν δg^μν + (g_μν□ − ∇_μ∇_ν)δg^μν in this convention, gives

\[
\mathcal L_RR_{\mu\nu}-\tfrac12\mathcal Lg_{\mu\nu}+(g_{\mu\nu}\Box-\nabla_\mu\nabla_\nu)\mathcal L_R=0 .
\]

The trace is 𝓛_R R − 2𝓛 + 3□𝓛_R = 6α□R − M²R/2 − JR + 3□J = 0. The R² pieces cancel identically (check 4a), as in Codex's (A11) and (A13). At first order in J around Minkowski space (checks 4b, 4c):

\[
R^{(1)}=-\frac{1}{2\alpha}\frac{\Box}{\Box-m_0^2}J
=-\frac{J}{2\alpha}-\frac{c_R^2}{\Box-m_0^2}J,\qquad
c_R^2=\frac{m_0^2}{2\alpha}=\frac{M^2}{24\alpha^2}.
\]

The general f(R) variation step itself was not machine-checked. Its sign structure is confirmed indirectly by route 4.3, which does not use it.

**4.2 From the Einstein frame, off shell (new; check 5).** Use the source-independent map of §3, under which the source becomes Codex's (A16). At first order in J only two terms survive:
- the scalar equation, (□ − m0²)φ = −3b□J, from the J·3b□_Eφ coupling;
- the Einstein-frame trace with the J R_E coupling, −M²R_E/2 + 3□J = 0, so R_E^(1) = 6□J/M².

The terms J e^{−bφ}R_E − J and J(∇φ)² are second order. Mapping back with R = R_E + 3b□φ at linear order (§3.1 with u = bφ/2), and using 9b² = 6/M² and 6m0²/M² = 1/(2α), gives

\[
R^{(1)}=\frac{6}{M^2}\Box J-\frac{9b^2\,\Box^2}{\Box-m_0^2}J
=-\frac{1}{2\alpha}\frac{\Box}{\Box-m_0^2}J .
\]

This reproduces the contact term as well as the pole. Codex's frame check (3b m0² = c_R, confirmed as check 6c) covered the residue only. The pole position, residue and contact coefficient now agree between the two frames.

**4.3 Restricting before variation (new; checks 7a–7c).** For g = e^{2ψ}η (Lorentzian):
- R = −6e^{−2ψ}[□ψ + (∂ψ)²] (check 7a).
- (M²/2)√−g R = 3M²e^{2ψ}(∂ψ)² plus a total derivative (check 7c).
- To quadratic order with the source, 𝓛₂ = 3M²(∂ψ)² + 36α(□ψ)² − 6J□ψ.

Its equation, □[(72α□ − 6M²)ψ − 6J] = 0, gives R^(1) = −6□ψ equal to the result of §4.1, contact term included (check 7b).

What this does and does not show:
- The ψ propagator ∝ 1/[□(72α□ − 6M²)] keeps a massless conformal pole whose kinetic sign is wrong (the 3M²(∂ψ)² term).
- In the full theory that pole is removed by the constraints. Only the curvature correlator agrees between the two treatments.
- This agreement is specific to the linear response to a scalar source at tree level about Minkowski space. It does not show that restriction commutes with variation in general, and it says nothing about nonlinear or tensor dynamics.
- It supports the sheet's rule (§6) to compare curvature correlators, never σ or ψ propagators.

**4.4 Source normalization.** J has mass dimension 2 and couples to the Jordan-frame R.
- The Einstein-frame source is (A16), not J R_E.
- Absorbing J into F instead makes the frame map source-dependent, and derivatives of W[J] then pick up additional terms. §4.2 avoids this by keeping the map fixed.
- Local terms ∫J² can be added as source counterterms. They shift only the contact coefficient, never m0² or c_R².

## 5. Euclidean continuation and the PS comparison

This section does **not** identify α_E with α_L, or σ with the scalaron. It records what a stated continuation implies.

**5.1 Sign map (derivation; assumption stated).** Rotate t = −iτ in the bulk Lagrangian and continue the metric so that √−g → √g with R unchanged. This is the rotation that gives the Gibbons–Hawking sign −(M²/2)∫√g R. Then iS_L = −S_E with

\[
S_E=-\int d^4x\sqrt g\,\Big[\tfrac{M^2}{2}R+\alpha_LR^2\Big]-\int\sqrt g\,JR ,
\]

so the Euclidean R² coefficient is −α_L. On conformally flat metrics g = e^{2ψ}δ:
- S_E = −3M²∫e^{2ψ}(∂ψ)² − 36α_L∫[Δψ + (∂ψ)²]².
- For the healthy branch α_L > 0, both terms are negative: the Einstein-term conformal-factor problem plus a negative R² term.
- A real, positive PS action (α_E > 0) is reached under this rotation only from α_L = −α_E < 0, which is the tachyonic flat-vacuum branch of §3.6.
- Rotating the conformal factor instead (ψ → iψ, à la Gibbons–Hawking–Perry) makes both quadratic terms positive for α_L > 0. It also makes the cubic and quartic PS terms complex: [iΔψ − (∂ψ)²]² is not real.

Neither of the two prescriptions examined here maps the healthy Lorentzian branch onto the real, positive PS action. This is my inference and needs a literature check. It is consistent with the branch distinction that `parent_theory.tex` (§`renormalization-couplings`) attributes to Avramidi's thesis.

**5.2 Gaussian curvature correlators (checks 8, 9, 10).**
- The `model.tex` identities hold exactly (checks 8a, 8b), and α_E = 1/(72λ²) is the unique normalization (check 8c).
- At Gaussian order and p ≠ 0, PS gives ⟨R^(1)R^(1)⟩ = 36λ² p⁴ · p⁻⁴ = 36λ² = 1/(2α_E), a pure contact term (check 9). This agrees with Codex's (A17).
- The benchmark's Euclidean tree-level response, from §4.1 with □ → −p², is R/J = −(1/(2α_L)) · p²/(p² + m0²) (check 10a).

So:
- **For p ≫ m0, or M → 0 at fixed p,** R/J → −1/(2α_L) (checks 10b, 10c). The M → 0 limit of this Jordan-frame linear response is smooth, although the Einstein-frame variables are singular there (b ∝ 1/M).
- **With α_E = −α_L from §5.1,** this limit equals the PS correlator 1/(2α_E) = 36λ² exactly.
- **At finite M,** the scalaron pole, with weight ∝ m0²/(p² + m0²), has no PS counterpart. This agrees with Codex's §6.
- **Consistency with Barker–Glavan.** The PS curvature correlator is pure contact, with no propagating curvature mode. That fits the empty linearized pure-R² spectrum cited in `model.tex` (subsection "Perturbative anchors and the scope of the reduction", ref. `BarkerGlavan:2026`), but it is not a proof of it.

**Interpretation.** The Gaussian PS curvature correlator is the ultraviolet contact term of the benchmark's curvature response, *on the α_L < 0 sign under naive Wick rotation*. It is a tree-level, Gaussian-order consistency of normalizations. It does not identify fields: the scalaron is φ ≈ c_R⁻¹R^(1) ∝ □ψ, not the conformal factor ψ = λσ. Nor is it a quantum equivalence.

## 6. Steps that need a contour, Jacobian or regulator at the quantum level

| Step | Classical status | Quantum requirement |
|---|---|---|
| Auxiliary χ (§2) | Exact algebraic elimination | Lorentzian: an oscillatory Fresnel integral needing an iε prescription and a metric-dependent ultralocal normalization ∏ₓ(α√−g)^{−1/2}. Euclidean with naive Wick: the exponent is −α_L∫√g(χ−R)² + α_L∫√g R², so the real-χ Gaussian **converges for α_L > 0** (qualifying Codex's §1 remark) and diverges for α_L < 0, the PS sign, which needs χ → iχ. |
| Chart restriction F > 0 (§3.5) | Branch choice | Real-χ integration runs over F ≤ 0 too. Restricting the domain or choosing a contour that avoids F = 0 changes the functional integral and must be specified. |
| χ → φ (§3.3) | Bijective on F > 0 | Jacobian ∏ₓ dχ/dφ ∝ ∏ₓ e^{bφ}, an ultralocal δ⁴(0) term: zero in dimensional regularization, nonzero on a lattice. |
| g → g_E (§3.1) | Field redefinition | Jacobian of the metric measure (DeWitt or another choice); gauge fixing and Faddeev–Popov determinants in each frame. Perturbative S-matrix equivalence of frames does not make off-shell effective actions or composite correlators equal. |
| Wick rotation (§5.1) | Not used classically | Contour for the conformal factor; the sign map α_E = −α_L under naive rotation; complex interactions under conformal rotation. |
| Composite source J R (§4.4) | Linear response | Renormalization of the composite operator and source contact counterterms; source transformation under each field redefinition. |
| Dropped divergences (§3.2, §7) | Boundary assumption | Boundary action and boundary data in the functional integral. |
| Ghost sector | Absent (β_L = 0) | No prescription selected for the massive spin-2 state of generic QQG (conventions §4). |

## 7. Boundary terms

- **Divergences dropped.** The 3M²□_E u term of §3.2 and the integrations by parts behind the trace equation and §4.3 are dropped. This needs fields that vanish, or are periodic, at infinity, the sheet's local-calculation assumption.
- **Finite regions.** The original R + αR² action needs extra boundary data: R, equivalently χ, on the boundary. The auxiliary form makes this explicit through a boundary term ∮√|h| F K. The frame change maps it to M²∮√|h_E| K_E plus normal-derivative terms in u, which cancel the flux of 3M²□_E u.
- **Scope.** None of this establishes a boundary-value equivalence. The sheet (§2) already excludes that claim.
- **Euler term.** Constant Euler terms are irrelevant at fixed topology and are omitted (conventions §2; `classical_gravity.tex` eq. `quadratic-basis` and the text after it).

## 8. Notation map

| Symbol | Meaning | Relation, domain or caution |
|---|---|---|
| M | Reduced Planck mass | M² = (8πG)⁻¹ > 0; `M_Pl` in `classical_gravity.tex` |
| α_L, β_L | Lorentzian R², C² coefficients | α_L = α of eq. `qg-cosmological-action`; α_L = 1/(6f0²) on α > 0 (eq. `parent-couplings`) |
| α_E | Coefficient of the Euclidean ∫√g R² functional in `model.tex` | 1/(72λ²). Equals −α_L **only** under naive Wick rotation; never set equal to α_L |
| χ | Auxiliary curvature variable | χ = R classically |
| F | M² + 4α_L χ | Chart F > 0; F/M² = f_R in De Felice–Tsujikawa's normalization |
| u | ln Ω = ½ ln(F/M²) | Weyl-rescaling exponent |
| φ | Canonical scalaron | √(3/2) M ln(F/M²); dimension 1 |
| b | √(2/3)/M | Singular as M → 0 |
| m0 | Scalaron mass | m0² = M²/(12α_L) = f0²M²/2; tachyonic if α_L < 0 |
| c_R | Curvature–scalaron coupling | R^(1) = c_R φ on shell; c_R² = m0²/(2α_L) |
| J | Source for Jordan-frame R | Dimension 2; in the Einstein frame, eq. (A16) |
| ψ | Jordan-frame conformal factor, g = e^{2ψ}η | Used only in §4.3 and §5.1 |
| σ, λ | PS field and coupling | ψ = λσ only on the restricted Euclidean configuration space. σ is **not** the scalaron: φ ∝ □ψ at linear order |
| Δ vs □ | Euclidean Laplacian vs Lorentzian d'Alembertian | `model.tex` currently writes □ for Δ |

## 9. Independent checks

The script `check_scalaron_benchmark.py` contains 31 symbolic checks, all of which pass.

| Group | Check numbers | Content |
|---|---|---|
| Conventions | 0 | Sign of R on a sphere |
| Frame map | 1 | Conformal transformation on a non-flat g_E |
| Auxiliary | 2a–2c | Identity, χ = R, form of F |
| Einstein frame | 3a–3h | Inversion, kinetic −½, V, V(0) = V′(0) = 0, m0², f0 form, chart Jacobian |
| Trace equation | 4a–4c | Algebraic cancellation, response, contact/pole split |
| Off-shell frame | 5 | Einstein-frame response equals the Jordan-frame response |
| On-shell coupling | 6a–6c | c_R, c_R² = m0²/(2α), 3b m0² = c_R |
| Restricted mode | 7a–7c | Exact conformal R, response matching, Einstein-term reduction |
| PS identity | 8a–8c | `model.tex` identities and α_E normalization |
| PS correlator | 9 | ⟨RR⟩ = 36λ² = 1/(2α_E) |
| Euclidean response | 10a–10c | Continuation, ultraviolet limit, M → 0 limit |
| Euclidean auxiliary | 11 | Split of the exponent and convergence condition |

**Acceptance criteria from conventions §7:**
- Original bulk action recovered: **pass** (2a).
- Dimensions of Einstein-frame quantities: **pass** (by inspection: [V] = 4, [φ] = 1, [m0²] = 2).
- V and V′ vanish at the vacuum: **pass** (3d, 3e).
- Scalar pole consistent from two derivations: **pass**, and from three routes including the contact term (4b, 5, 7b).
- Explicit invertibility domain: **pass** (§3.5).

No check failed. One step was not machine-verified: the general f(R) metric variation of §4.1.

## 10. Unresolved questions

1. **Branch and continuation.** Is there a contour that maps the healthy Lorentzian branch (α_L > 0) onto a real, positive Euclidean PS action? If not, the manuscript should say plainly that the PS sign corresponds, under naive rotation, to the tachyonic flat-vacuum branch. This must be checked against Avramidi's thesis and Anderson et al.'s definition of CFQG, in full text.
2. **Limits.** In the benchmark the M → 0 linear response is smooth and matches PS at Gaussian order. Does this survive interactions? The cubic PS vertex has no counterpart at finite M in the conformal restriction without the e^{2ψ} Einstein-term factors.
3. **Constraints.** §4.3 agreed because only the trace sector responds to a scalar source at linear order. Packet B's constrained mode decomposition should test where restriction and variation first disagree, at nonlinear order or with tensor sources.
4. **Source scheme.** At loop level the composite J R needs a stated renormalization and contact-term scheme before W[J] can be compared with PS lattice correlators. The lattice zero-mode subtraction (conventions §5) belongs to PS only.
5. **Measure.** No nonperturbative metric measure, frame Jacobian or gauge-fixed determinant is specified. This blocks any quantum-equivalence claim (conventions §4).
6. **Literature.** De Felice and Tsujikawa (arXiv:1002.4928), §§2.1–2.3 and eq. (3.16), are cited by Codex for the standard f(R) results. I did not read their full text; the results above were rederived, not imported.

## Sources

- Repo files at commit `2b204f1`: `notes/SCALARON_CONVENTIONS.md` §§2–7; `sections/classical_gravity.tex` (eqs. `qg-cosmological-action`, `quadratic-basis`); `sections/parent_theory.tex` (eq. `parent-couplings`, §`renormalization-couplings`); `sections/model.tex` (eqs. `ps-action`, `conformal-metric`, `conformal-r-squared`; §`cosmological-action`); `notes/SCALARON_BENCHMARK_DERIVATION.md` (A1–A17).
- External results were not consulted in full text for this note. The f(R) structure and the Gibbons–Hawking Euclidean sign are standard and rederived here. The attribution to Avramidi is the manuscript's own (`parent_theory.tex`), at abstract level in our earlier Bond pass.
