# Packet B — flat-background modes and the PS restriction

21 November 2026. Codex analytical working derivation for review; follows [Gate 0](SCALARON_CONVENTIONS.md) and [Packet A](SCALARON_BENCHMARK_DERIVATION.md). Excluded from TeX. Classical quadratic order only; no assertion about a nonperturbative quantum contour or ghost resolution.

## 1. Action and conventions

Use g=eta+h (h is dimensionless), signature (-,+,+,+), and the curvature convention of Gate 0:

\[
S=\int\sqrt{-g}\,[M^2R/2+\alpha R^2+\beta C^2],
\quad M>0,\quad\alpha=\alpha_L>0,\quad\beta=\beta_L<0.
\tag{B1}
\]

The cosmological constant vanishes. Surface terms are dropped under the stated vanishing-boundary assumption. Define Box=eta^{mu nu} partial_mu partial_nu; its Fourier eigenvalue is s=omega²-|k|². Lorentzian alpha and beta are not the Euclidean PS coefficient alpha_E.

## 2. York decomposition and gauge variables

For nonzero, non-null momenta, decompose

\[
h_{\mu\nu}=h^{TT}_{\mu\nu}+\partial_\mu v_\nu+\partial_\nu v_\mu
+(\partial_\mu\partial_\nu-\eta_{\mu\nu}\Box/4)\psi
+\eta_{\mu\nu}h/4,
\quad\partial^\mu h^{TT}_{\mu\nu}=0,\quad h^{TT\mu}{}_{\mu}=0,
\quad\partial^\mu v_\mu=0.
\tag{B2}
\]

Under xi_mu=xiT_mu+partial_mu epsilon, delta v=xiT, delta psi=2 epsilon and delta h=2 Box epsilon. Thus t=h-Box psi is invariant. Reorganizing the scalar terms,

\[
h_{\mu\nu}=h^{TT}_{\mu\nu}+\frac13\theta_{\mu\nu}t
+\partial_\mu v_\nu+\partial_\nu v_\mu
+\partial_\mu\partial_\nu(\psi+t/(3\Box)),
\quad\theta_{\mu\nu}=\eta_{\mu\nu}-\frac{\partial_\mu\partial_\nu}{\Box}.
\tag{B3}
\]

The last two terms are pure gauge in this non-null decomposition. Linear curvature is

\[
R^{(1)}=\partial_\mu\partial_\nu h^{\mu\nu}-\Box h=-\Box t.
\tag{B4}
\]

Inverse Box needs boundary conditions and is singular on null and zero modes. Therefore the off-shell decomposition cannot by itself count massless physical states. The massive spin-two field has five polarizations in the 4D TT block, including its helicity-zero and helicity-one components. They are not the gauge vector v or the invariant spin-zero t. This distinction prevents confusing a spatial scalar perturbation with the scalaron.

## 3. Quadratic action derived on the invariant blocks

The Einstein quadratic bulk Lagrangian can be written

\[
\mathcal L_{EH}^{(2)}=-\frac{M^2}{8}
[\partial_\rho h_{\mu\nu}\partial^\rho h^{\mu\nu}
-2\partial_\mu h^{\mu\nu}\partial^\rho h_{\rho\nu}
+2\partial_\mu h^{\mu\nu}\partial_\nu h-(\partial h)^2].
\tag{B5}
\]

On the transverse representative h=TT+theta t/3 its integrated form is

\[
S_{EH}^{(2)}=\int[M^2h^{TT}_{\mu\nu}\Box h^{TT\mu\nu}/8
-M^2t\Box t/12].
\tag{B6}
\]

At flat background R² contributes alpha(Box t)². For Weyl curvature use the four-dimensional identity C²=E4+2 Ricci²-(2/3)R². The Euler term is a surface term here. For TT, Ricci^(1)=-Box hTT/2 and R^(1)=0, giving beta(Box hTT)²/2. The theta t/3 scalar representative is conformal plus pure gauge, so its linear Weyl curvature vanishes. Consequently

\[
S^{(2)}=\frac12\int h^{TT}_{\mu\nu}
[M^2\Box/4+\beta\Box^2]h^{TT\mu\nu}
+\frac12\int t[-M^2\Box/6+2\alpha\Box^2]t.
\tag{B7}
\]

There is no scalar–TT cross term at this order. This is diagonalization of the free operator; it does not establish nonlinear or quantum decoupling.

## 4. Projector inverse and physical poles

On symmetric tensors define P0=theta_mu_nu theta_rho_sigma/3 and

\[
P_2=\tfrac12(\theta_{\mu\rho}\theta_{\nu\sigma}+\theta_{\mu\sigma}\theta_{\nu\rho})-P_0.
\]

With S²=(1/2) integral h K h the invariant kernel is

\[
K=(M^2s/4+\beta s^2)P_2+(-M^2s/2+6\alpha s^2)P_0.
\tag{B8}
\]

The factor of three in the scalar block follows from (theta t/3)²=t²/3. Define

\[
m_2^2=-M^2/(4\beta)>0,\qquad m_0^2=M^2/(12\alpha)>0.
\tag{B9}
\]

The invariant part of the inverse, before choosing pole contours, is

\[
K^{-1}=\frac4{M^2}\left[
\left(\frac1s-\frac1{s-m_2^2}\right)P_2
+\left(-\frac1{2s}+\frac1{2(s-m_0^2)}\right)P_0\right].
\tag{B10}
\]

This is an inverse on the gauge quotient, not a full gauge-fixed propagator; no gauge-fixing coefficient is implicit. For conserved external tensor sources, omitted longitudinal terms do not contribute. A full metric inverse can be constructed with the Gate 0 de Donder condition and an explicit gauge-fixing parameter in a separate calculation. The algebraic ghost residue is fixed here without choosing a nonperturbative ghost contour.

For a conserved symmetric source T, contractions give T P2 T=T_mu_nu T^mu_nu-T²/3 and T P0 T=T²/3. Thus the massless numerator is T_mu_nu T^mu_nu-T²/2, the standard Einstein combination. Its negative P0 coefficient is not an extra propagating scalar ghost. Gauge constraints and residual gauge freedom leave two massless graviton helicities. At the timelike massive spin-two pole there are five independent TT polarizations and an opposite residue relative to the graviton. The scalar pole has positive residue and one polarization.

| Physical mode | Polarizations | Mass squared | Conventional perturbative residue |
|---|---:|---|---|
| Einstein graviton | 2 | 0 | Healthy |
| Massive spin two | 5 | -M²/(4 beta) | Opposite sign (ghost) |
| Scalaron | 1 | M²/(12 alpha) | Healthy |

Total: eight for generic finite nonzero alpha,beta on this branch. At beta=0 the massive spin-two pole is absent and the benchmark has three physical modes. At alpha=0 the scalaron is absent. Neither limiting statement licenses taking M=0 in formulas using the above partial fractions and gauge argument. With alpha=1/(6f0²), beta=-1/(2f2²), the two masses are m0²=f0²M²/2 and m2²=f2²M²/2. Independence of couplings does not itself supply a hierarchy of masses.

## 5. Curvature contraction checks Packet A

Because R^(1)=-s theta^{mu nu}h_mu_nu and theta P0 theta=3, contraction of (B10) yields

\[
D_{RR}(s)=i\frac{6m_0^2}{M^2}\frac{s}{s-m_0^2+i0}
=\frac{i}{2\alpha}+\frac{iM^2}{24\alpha^2}\frac1{s-m_0^2+i0}.
\tag{B11}
\]

Here the healthy scalar Feynman prescription is used; the contact term is understood distributionally. The Weyl block drops out exactly at quadratic order. The curvature pole mass, squared coupling cR²=M²/(24 alpha²), and local contact term agree with Packet A. In its Z=exp(iW) convention, the response kernel is i D_RR, reproducing -s/[2 alpha(s-m0²)]. This is an independent check using the metric Hessian rather than the auxiliary potential or trace equation.

The absence of beta from this free curvature response does not imply absence of beta from interacting curvature correlators. Expanding the full action includes scalar–metric interactions; loops and exchange can involve the tensor sector.

## 6. Restricting the metric before variation loses equations

Set g=e^{2 omega}eta, omega=lambda sigma, in the Lorentzian action only for this diagnostic. At linear order h=2 omega eta, t=6 omega, R^(1)=-6 Box omega, and linear Weyl curvature is zero. The restricted quadratic action is

\[
S_{conf}^{(2)}=\int[-3M^2\omega\Box\omega+36\alpha(\Box\omega)^2].
\tag{B12}
\]

Its Euler equation is Box(Box-m0²)omega=0. It appears to allow both a massless conformal wave and the massive scalaron branch. To test whether either branch solves the parent equations, retain their full tensor form. Weyl variation vanishes on a conformally flat metric; at linear order the remaining equation is

\[
\mathcal E^{(1)}_{\mu\nu}
=(\eta_{\mu\nu}\Box-\partial_\mu\partial_\nu)
[M^2\omega-12\alpha\Box\omega]=0.
\tag{B13}
\]

Its trace is precisely the restricted scalar equation up to a nonzero factor. For Box omega=m0² omega, the bracket vanishes: this branch is a valid linear scalaron representative. For a nonconstant null plane wave with Box omega=0, however, E_mu_nu=-M² partial_mu partial_nu omega is generally nonzero. The spurious massless branch of the restricted scalar equation fails the parent tensor equations. Special constant/affine configurations require the stated boundary treatment and do not establish a massless particle.

This gives a concrete obstruction: imposing conformal flatness before variation retains the trace equation while discarding traceless equations. It is not automatically a consistent truncation of the physical solution space. The non-null projector construction was not used at the null pole in this test; (B13) follows directly from local metric equations.

At the nonlinear level the Einstein-frame scalar stress tensor is quadratic and generally sources the metric, so maintaining a fixed flat Einstein metric is also not a generic nonlinear scalaron solution.

## 7. Consequences for PS and next calculations

The existing Euclidean PS action retains the conformal R² square and omits the Einstein-induced two-derivative term. It has its own flat zero-mean scalar measure. Equations (B12–B13) explain why adding the missing lower-derivative term alone still would not restore all parent constraints.

The positive Euclidean PS measure must not be assigned the Lorentzian pole signs of (B12) by inspection. A Wick-rotation and measure prescription is required. The leading Euclidean PS curvature correlator remains 36 lambda² at nonzero momentum, as computed in Packet A, whereas (B11) contains a finite-mass scalaron pole. Matching needs a stated regime, observable map and constraint treatment.

If one wants a physical scalaron benchmark now, the Einstein-frame canonical scalar plus its dynamical Einstein metric is the established classical description. If one wants to test the PS model's own interacting continuum limit, its existing lattice question remains well defined. Establishing a relation between these tasks requires additional work.

Next priorities: independently review the normalization and null-branch test; derive the nonlinear conformal restriction including the Einstein term (packet C); then test whether any tensor-elimination hierarchy and prescribed quantum contour justify a scalar effective action (packet D). A heavy massive spin-two field can sometimes admit a derivative expansion, but the massless Einstein graviton cannot generically be eliminated into a local scalar-only action.

## 8. Evidence and checks

This note explicitly derives the invariant Hessian, its inverse, the curvature contraction, and the conformal branch test. Checks are analytical cross-checks, not peer review or formal verification. No simulations, loop calculations, or quantum decoupling proof are reported.

Literature anchor: [Salvio, Quadratic Gravity, sections 2.1 and 2.3.1–2.3.3](https://arxiv.org/html/1804.09944v2) identifies the curvature basis, masses and the helicity decomposition of the massive spin-two ghost. Its metric conventions differ; the signs here are derived from (B1), not copied from its Einstein-term sign. The reported spectrum agrees after translating the couplings. Packet A supplies the independent beta=0 scalaron benchmark and its f(R) references.
