# Packet A — classical scalaron benchmark

21 November 2026. Derived by Codex for Joel's review, following [the conventions sheet](SCALARON_CONVENTIONS.md). This is an analytical working note, not a Claude-generated result or an independently certified proof. It is excluded from the TeX manuscript.

## Result and scope

For finite M>0 and alpha_L>0, the Lorentzian R+R² bulk action is classically equivalent, on the branch F>0, to Einstein gravity plus a canonically normalized scalar with

\[
V(\phi)=\frac{M^4}{16\alpha_L}(1-e^{-\sqrt{2/3}\phi/M})^2,
\qquad m_0^2=\frac{M^2}{12\alpha_L}.
\tag{A1}
\]

The mass follows independently from the original metric trace equation. This identifies the physical scalaron of the benchmark. It does not establish a quantum equivalence with PS. In fact, the finite-mass curvature pole is absent from the leading Gaussian curvature correlator of the pure PS model with its stated observable dictionary.

All calculations below use signature (-,+,+,+), the curvature convention in the conventions sheet, beta_L=0, zero cosmological constant and no matter. Write alpha=alpha_L only within this note. Boundary variations and total divergences are assumed to vanish; a finite-boundary variational principle needs its own boundary action.

## 1. Auxiliary field: an algebraic identity before a frame change

Start with

\[
S_J=\int\sqrt{-g}\,[M^2R/2+\alpha R^2].
\tag{A2}
\]

Replace it by

\[
S_{aux}=\int\sqrt{-g}\,[F(\chi)R/2-U(\chi)],
\quad F=M^2+4\alpha\chi,\quad U=\alpha\chi^2.
\tag{A3}
\]

Variation at fixed g gives delta S_aux/delta chi=2 alpha sqrt(-g)(R-chi). Therefore chi=R for alpha nonzero, and substitution recovers (A2). More explicitly,

\[
2\alpha\chi R-\alpha\chi^2
=\alpha R^2-\alpha(\chi-R)^2.
\tag{A4}
\]

This proves classical elimination. It is not a convergent Euclidean Gaussian identity over a real auxiliary field.

## 2. Einstein-frame transformation with the derivative terms retained

Define Omega²=F/M²>0, u=ln Omega and gE=Omega² g. In four dimensions,

\[
\sqrt{-g}=\Omega^{-4}\sqrt{-g_E},\qquad
R[g]=\Omega^2[R_E+6\Box_Eu-6(\nabla_Eu)^2].
\tag{A5}
\]

Consequently the curvature term becomes

\[
\sqrt{-g}\,FR/2=\sqrt{-g_E}
[M^2R_E/2+3M^2\Box_Eu-3M^2(\nabla_Eu)^2].
\tag{A6}
\]

The integral of Box_E u is the surface flux of its gradient; drop it only under the stated boundary assumption. Since u=(1/2)ln(F/M²), the remaining kinetic coefficient is -3M²/(4F²). Define

\[
\phi=\sqrt{3/2}\,M\ln(F/M^2),\qquad
b=\sqrt{2/3}/M,\qquad
F=M^2e^{b\phi},\quad
\chi=\frac{M^2}{4\alpha}(e^{b\phi}-1).
\tag{A7}
\]

Then -3M²(grad F)²/(4F²)=-(grad phi)²/2, and

\[
S_E=\int\sqrt{-g_E}\,[M^2R_E/2-(\nabla_E\phi)^2/2-V(\phi)],
\quad V=\frac{M^4U}{F^2}
=\frac{M^4}{16\alpha}(1-e^{-b\phi})^2.
\tag{A8}
\]

The Einstein metric remains dynamical. This transformation does not integrate out the ordinary graviton.

## 3. Domain, vacuum and first mass calculation

F>0 means chi>-M²/(4 alpha). This open branch maps bijectively to real phi: dphi/dchi=4 alpha sqrt(3/2) M/F is nonzero. F=0 is a singular boundary, not another point in this chart. At the vacuum chi=R=phi=0, F=M² and both frames have the same Minkowski metric.

With A=M⁴/(16 alpha),

\[
V'=2Ab(e^{-b\phi}-e^{-2b\phi}),\qquad
V''=2Ab^2(-e^{-b\phi}+2e^{-2b\phi}).
\tag{A9}
\]

Thus V(0)=V'(0)=0 and V''(0)=2Ab²=M²/(12 alpha)>0. The scalar equation is Box_E phi-V'=0, so its linearization is (Box_eta-m0²)phi=0. Around this vacuum the scalar kinetic energy has the healthy sign. Dimensions agree: V has mass dimension four, phi dimension one, V'' dimension two.

## 4. Independent calculation from the metric equations

Let L(R)=M²R/2+alpha R² and L_R=M²/2+2 alpha R. Metric variation gives

\[
L_R R_{\mu\nu}-\frac12Lg_{\mu\nu}
+(g_{\mu\nu}\Box-\nabla_\mu\nabla_\nu)L_R=0.
\tag{A10}
\]

Tracing in four dimensions,

\[
L_RR-2L+3\Box L_R
=-M^2R/2+6\alpha\Box R=0.
\tag{A11}
\]

The algebraic R² contributions cancel exactly. Equation (A11) is an exact vacuum trace equation, although its Box still depends on the metric. Linearizing at Minkowski gives

\[
(\Box_\eta-m_0^2)R^{(1)}=0,
\qquad m_0^2=M^2/(12\alpha),
\quad R^{(1)}=\partial_\mu\partial_\nu h^{\mu\nu}-\Box_\eta h.
\tag{A12}
\]

The mass agrees with (A9). Under a linearized diffeomorphism delta h_mu_nu=partial_mu xi_nu+partial_nu xi_mu, delta R^(1)=0 by commuting flat derivatives. No gauge-fixing coefficient entered this calculation. The trace equation alone is not the entire metric system or a mode-counting proof.

In the manuscript's Lorentzian normalization alpha_L=1/(6 f0²), this is m0²=f0² M²/2. No analogous identification with lambda or alpha_E has been made.

## 5. Curvature source, pole normalization and contact term

Introduce a fixed external scalar source through S[J]=S_J+integral sqrt(-g) J R, with [J]=2. In the original frame the metric equation uses L_J=L+JR and L_J,R=L_R+J. Its trace is

\[
6\alpha\Box R-M^2R/2-JR+3\Box J=0.
\tag{A13}
\]

At linear order in both the source and response, the JR term is second order. With homogeneous response set to zero and the Feynman inverse chosen,

\[
R^{(1)}[J]=-\frac{1}{2\alpha}
\frac{\Box_\eta}{\Box_\eta-m_0^2}J
=-\frac{J}{2\alpha}
-\frac{m_0^2}{2\alpha}(\Box_\eta-m_0^2)^{-1}J.
\tag{A14}
\]

The first term is a local contact response. The second contains a pole with squared curvature-to-canonical-scalar coupling

\[
c_R^2=\frac{m_0^2}{2\alpha}=\frac{M^2}{24\alpha^2},
\qquad c_R=\frac{M}{2\sqrt6\alpha}.
\tag{A15}
\]

This also follows from chi=R and (A7) on source-free solutions: R^(1)=c_R phi. It is an on-shell relation; replacing the geometric curvature operator by c_R phi off shell loses contact information.

For clarity about signs, put s=omega²-|k|², the eigenvalue of Box_eta on exp(-i omega t+i k.x). The canonical time-ordered propagator is i/(s-m0²+i0); the curvature pole part is i c_R²/(s-m0²+i0). The classical susceptibility in (A14) is instead -1/(2 alpha)-c_R²/(s-m0²+i0). If Z[J]=exp(iW[J]), W''=i times the connected time-ordered correlator, explaining this difference of i and sign. The physical scalar spectral weight c_R² is positive. Local source counterterms can change contact pieces without moving this pole.

Using the source-independent transformation (A7), the exact transformed source is

\[
S_{source}=\int\sqrt{-g_E}\,J e^{-b\phi}
\left[R_E+3b\Box_E\phi-\frac32b^2(\nabla_E\phi)^2\right].
\tag{A16}
\]

One may instead absorb J into F, but then the frame change depends on the source and its derivatives must be included when differentiating W. Equation (A16) avoids that ambiguity. On a free scalaron solution R_E^(1)=0 and Box phi=m0² phi, so 3b m0²=c_R, a second check on the pole normalization. Off shell, retain R_E and the derivative coupling.

## 6. What this says about the PS comparison

Use the existing Euclidean PS Gaussian expansion at nonzero momentum, with canonical free covariance <sigma(p)sigma(-p)>=1/p⁴ and geometric linearization R_PS^(1)=-6 lambda Delta sigma. At this leading order,

\[
\langle R_{PS}^{(1)}(p)R_{PS}^{(1)}(-p)\rangle_0
=36\lambda^2,\qquad p\ne0.
\tag{A17}
\]

This is momentum independent, hence contact-like (with a zero-mode subtraction in finite volume). It has no finite scalaron pole. This elementary comparison already rules out equating the two bare Gaussian observable dictionaries at finite M without additional ingredients. It does not rule out an interacting limit, a modified dictionary, or a matched effective action.

Equation (A17) is an expansion around sigma=0; it does not set lambda=0 in the geometric definition. The positive Euclidean PS coefficient alpha_E remains independent of alpha_L. Sending M to zero in (A7) is singular and cannot be used as a proof of equivalence with pure R² or PS. The scalaron benchmark also retains metric constraints that a restriction before variation can remove.

## 7. Quantum questions left open

The identity (A4), when inserted into a Lorentzian functional integral, needs a regulated oscillatory auxiliary Gaussian and its normalization. A Wick rotation must specify its contour. The metric-dependent Gaussian normalization cannot simply be discarded before specifying the measure.

The frame and scalar redefinitions require their functional Jacobians, consistent gauge fixing and ghost determinants, regulator and counterterms. Quantum composite curvature insertions require renormalization and source contact terms. Nonlinear local curvature at a fixed coordinate is not automatically a gauge-invariant observable. Nothing here selects a prescription for the additional spin-two sector of generic QQG.

## 8. Notation and verification record

| Symbol | Meaning | Relation or domain |
|---|---|---|
| M | Reduced Planck mass | M²=(8 pi G)^(-1)>0 |
| alpha=alpha_L | Lorentzian R² coefficient | Positive, dimensionless |
| chi | Auxiliary curvature variable | chi=R on classical auxiliary solutions |
| F | Einstein coefficient times two | M²+4 alpha chi>0 |
| phi | Canonical scalaron | sqrt(3/2) M ln(F/M²) |
| m0 | Scalaron mass | M/sqrt(12 alpha) |
| sigma, lambda | Euclidean PS field and coupling | Dimensionless; not identified with phi, f0 |
| J | Source for Jordan-frame curvature | Mass dimension two |

Analytical checks completed: auxiliary elimination (A4); canonical normalization (A6–A8); vacuum and dimensions (A9); independent mass from metric variation (A12); curvature pole normalization by source response and on-shell frame map (A14–A16). This is internal algebraic cross-checking, not independent peer review or proof-assistant certification. No simulations were run.

## Literature cross-check and next gate

[De Felice and Tsujikawa, f(R) theories](https://arxiv.org/html/1002.4928v2), sections 2.1–2.3, supplies the metric, auxiliary-field and conformal frameworks. Their action is (M²/2) integral sqrt(-g) f(R); their dimensionless derivative f_R corresponds to this note's F/M² on shell. Thus use f(R)=R+(2 alpha/M²)R² when translating. Section 3.2, equation (3.16), gives the Starobinsky potential in a different mass notation. These are established results, rederived here in the manuscript's normalization; the source-response and PS comparison are explicit calculations in this note.

Next: have a separate reviewer check the derivation and source contact terms, then undertake packet B's constrained mode decomposition before trying to remove tensor modes or alter the PS action.
