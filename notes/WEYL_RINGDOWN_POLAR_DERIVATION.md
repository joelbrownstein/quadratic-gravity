# Gate W3 — polar Schwarzschild equations and radial-infall sources

21 November 2026. Working derivation and reproducible classical checks, excluded from TeX. This follows [W0](WEYL_RINGDOWN_CONVENTIONS.md), [W1](WEYL_RINGDOWN_AXIAL_DERIVATION.md), and [W2](WEYL_RINGDOWN_MASSIVE_SPECTRUM.md). The deliverable is a constrained polar ODE, reconstruction code, conserved source convention, and selected spectral validation. It is not yet a sourced waveform or a bound on f2.

## 1. Retained theory and sourced auxiliary equations

Use the Einstein–Weyl benchmark: alpha_L=0, beta_L=−1/(2 f2²), signature (−+++), W0 curvature conventions, kappa=M_Pl^−2=8 pi G, and

\[
\mu^2=m_2^2=-\frac{M_{\rm Pl}^2}{4\beta_L}
=\frac{f_2^2M_{\rm Pl}^2}{2}.
\]

Dimensionless code variables set r_g=GM_BH=1; hence code mu=r_g m2, omega=r_g omega_physical, horizon r=2. All Fourier modes have exp(−i omega t). The Schwarzschild auxiliary background vanishes. Write the physical metric perturbation as h=H+2q. The linear auxiliary equations are

\[
G^{(1)}_{\mu\nu}(H)=\kappa T_{\mu\nu},\qquad
G^{(1)}_{\mu\nu}(q)+\frac{\mu^2}{2}(q_{\mu\nu}-g_{\mu\nu}q)
=-\frac{\kappa}{2}T_{\mu\nu}.
\tag{W3.1}
\]

Background stress conservation and the linear Bianchi identity imply divergence(q)=gradient(q). On a Ricci-flat background this makes R^(1)[q]=0. Taking the trace of the second equation therefore gives

\[
q=\frac{\kappa T}{3\mu^2},\qquad
\nabla^\mu q_{\mu\nu}=\frac{\kappa}{3\mu^2}\nabla_\nu T.
\tag{W3.2}
\]

Substitution into the linearized Einstein tensor gives the sourced Lichnerowicz equation

\[
\mathcal E_{\mu\nu}\equiv(\Box-\mu^2)q_{\mu\nu}
+2R_{\mu\rho\nu\sigma}q^{\rho\sigma}
=J_{\mu\nu},\qquad
J_{\mu\nu}=\kappa\left(T_{\mu\nu}-\frac13g_{\mu\nu}T
+\frac{\nabla_\mu\nabla_\nu T}{3\mu^2}\right).
\tag{W3.3}
\]

Vacuum and exactly null infall have T=0 and the familiar transverse/traceless constraints. A massive particle has T≠0: imposing those vacuum constraints would give an inconsistent source. The asymptotically flat principal part is (Box−mu²)q, recovering W0's pole mass without redefining f2.

For generic alpha_L≠0, the additional equation is (12 alpha_L Box−M_Pl²)R^(1)=T. A finite-energy massive source then also excites the scalaron. The present finite-energy implementation is explicitly Einstein–Weyl; it is not a complete R²+C² response. The polar helicity-zero component of massive spin two is not the scalaron, and neither is established by the PS scalar toy action.

## 2. Harmonics and exceptional multipoles

Let F=1−2/r, L=ell(ell+1), and normalize complex scalar harmonics by integral |Y_lm|² dOmega=1. Define Y_AB=D_A D_B Y+(L/2)Omega_AB Y. The massive tensor is

\[
q_{tt}=FaY,\quad q_{tr}=bY,\quad q_{rr}=cY/F,
\quad q_{tA}=dD_AY,\quad q_{rA}=eD_AY,
\quad q_{AB}=r^2(K\Omega_{AB}Y+2G Y_{AB}).
\tag{W3.4}
\]

These are coefficients of q, not gauge choices for h. Because the background q vanishes, its linear perturbation is invariant under infinitesimal coordinate gauge transformations. A Regge–Wheeler gauge choice for H does not remove q's angular fields.

Writing Q=kappa T/(3 mu²), the projected constraints are

\[
\begin{aligned}
-a+c+2K&=Q,\\
i\omega a+Fb'+(F'+2F/r)b-Ld/r^2&=-i\omega Q,\\
i\omega b/F+c'+F'(a+c)/(2F)+2(c-K)/r-Le/r^2&=Q',\\
i\omega d/F+Fe'+(F'+2F/r)e+K-(L-2)G&=Q.
\end{aligned}\tag{W3.5}
\]

The last equation exists only for ell≥1.

| Multipole | Fields present | Independent massive radial fields | GR radiative fields |
|---|---|---|---|
| ell=0 | a,b,c,K | K | None; mass/constraint perturbations remain |
| ell=1 | a,b,c,d,e,K | K,e | None; dipole constraint/gauge pieces remain |
| ell≥2 | a,b,c,d,e,K,G | K,e,G | One polar Zerilli field |

The ell=0 and ell=1 coordinate projections are performed directly. No division by L or L−2 is used to obtain them. For the massless sector, the usual Zerilli potential is retained; the suite checks its exact Darboux relation to Regge–Wheeler, connecting it to W1's GR benchmark. This relation does not imply massive axial–polar isospectrality.

## 3. Explicit polar equations

Define B=omega²/F−mu²−L/r². The following left-hand sides use the projections tt/Y, tr/Y, rr/Y, tA/D_A Y, rA/D_A Y, angular trace/(2r²Y), and angular trace-free coefficient/(2r²Y_AB), respectively:

\[
\begin{aligned}
\mathcal E_{tt}={}&F^2a''+\frac{2F(r-1)}{r^2}a'
+(FB-2/r^4)a-\frac{4i\omega}{r^2}b
-\frac{2(2r-3)}{r^4}c+\frac{4F}{r^3}K,\\
\mathcal E_{tr}={}&Fb''+\frac{2(r-1)}{r^2}b'
+\left(B-\frac2{r^2}-\frac4{r^4F}\right)b
-\frac{2i\omega(a+c)}{r^2F}+\frac{2L}{r^3}d,\\
\mathcal E_{rr}={}&c''+\frac{2(r-1)}{r^2F}c'
+\left(\frac BF-\frac4{r^2}-\frac2{r^4F^2}\right)c
-\frac{2(2r-3)}{r^4F^2}a-\frac{4i\omega}{r^2F^2}b
+\frac{4(r-3)}{r^3F}K+\frac{4L}{r^3}e,\\
\mathcal E_{tA}={}&Fd''+(B+4/r^3)d+2Fb/r-2i\omega e/r^2,\\
\mathcal E_{rA}={}&Fe''+4e'/r^2+(B-4F/r^2)e
+2(c-K)/r+2(L-2)G/r-2i\omega d/(r^2F^2),\\
\mathcal E_K={}&FK''+2(r-1)K'/r^2+[B-2(r-4)/r^3]K
+2a/r^3+2(r-3)c/r^3-2LFe/r^3,\\
\mathcal E_G={}&FG''+2(r-1)G'/r^2+(B+2/r^2)G+2Fe/r^3.
\end{aligned}\tag{W3.6}
\]

Set each equal to the corresponding projection of J in (W3.3). Drop the nonexistent vector and tensor equations/harmonics at low ell. The generator constructs the connection, covariant tensor Laplacian and curvature contraction directly, rather than importing master-equation coefficients. It checks angular cancellation exactly at ell=0,1,2,3,4. The angular Laplacian gives affine L dependence; the general operator extracted from ell=2,3 is checked against the independently constructed ell=4 projection. The seven displayed formulas are checked symbolically against that operator.

The generated JSON stores all wave equations, the trace and its first two derivatives, and divergences and their first derivatives, as a matrix acting on (fields, fields', fields''). The reducer keeps (K,K') for ell=0, (K,e,K',e') for ell=1, and (K,e,G,K',e',G') otherwise. Solving the remaining algebraic linear system supplies accelerations and reconstructs every eliminated field and derivative. It checks all residual rows and rejects a singular/inconsistent reduction. Local derivative compatibility is tested separately: satisfying an algebraic matrix alone is insufficient to establish an ODE reduction.

The sampled eliminated ranks are 10,14,15 for ell=0,1,≥2. They establish the expected 1,2,3 massive second-order channels on those patches. This is not a proof that a chosen variable basis remains nonsingular at every complex frequency. Horizon, zero-mass, static and reduction-denominator limits require a regular basis or separate limiting analysis before a global response calculation.

## 4. Independent monopole equation and boundary problem

For comparison with Antoniou–Gualtieri–Pani (AGP2026), Eqs. (32)–(33), put

\[
\begin{gathered}
D=r(2\mu^2F+\mu^4r^2+4\omega^2),\qquad K''+P K'+V K=0,\\
P=\frac{2[4\mu^2rF^2+\mu^4r^2(r-1)+4\omega^2(2r-3)]}{rFD},\\
V=\frac{-F(\mu^6r^6+6\mu^2r^4\omega^2)
+4\omega^2(-12+6r+r^4\omega^2)
+\mu^4r^2[r^4\omega^2-4(8-6r+r^2)]}{DF^2r^3}.
\end{gathered}\tag{W3.7}
\]

Coordinate-derived elimination agrees with these independently transcribed coefficients. The higher-multipole G row also agrees with AGP2026 Appendix A, including its coupling to e.

The numerical spectral check selects the unstable monopole branch discussed in Appendix B.2, at mu=0.3 and 0.4. It does not use axial frequencies as polar reference values. For omega=i Omega, Omega>0, the correct boundaries are an ingoing/regular future-horizon solution and decay at infinity:

\[
K\sim(r-2)^{2\Omega}[1+b_1(r-2)+\cdots],\qquad
K\sim e^{-\sqrt{\mu^2+\Omega^2}\,r}
r^{-1-(\mu^2+2\Omega^2)/\sqrt{\mu^2+\Omega^2}}.
\tag{W3.8}
\]

The first horizon correction follows by substituting this series into (W3.7); it is implemented explicitly. Two-sided integration matches logarithmic derivatives at r=4 or 5, varying the inner offset from 10^−5 to 10^−7, outer radius from 60 to 100, and integration tolerance. An independent finite-difference eigenproblem uses the different master potential in Brito–Cardoso–Pani (BCP2013), Eq. (30), with grid refinement and Richardson extrapolation. These are growing bound states, not damped QNMs.

The results are approximately GM omega=i 0.03427946 and i 0.01122045. Numerical records retain all refinements and comparisons; the 2026 plots do not supply exact author values, so no such values are invented. Appendix B.1's damped scalar-, vector-, and tensor-led polar QNM branches have not all been reproduced. Their branch-resolved convergence is an explicit preflight requirement for any W4 waveform that relies on those resonances.

## 5. Conserved radial-infall source and normalization

Use a test particle with conserved energy E=m_p gamma and E/M_BH≪1, initially infalling from infinity with gamma≥1. Let v=sqrt(gamma²−F). Its background four-velocity and trajectory satisfy

\[
u^t=\gamma/F,\quad u^r=-v,\quad u_t=-\gamma,\quad u_r=-v/F,
\qquad t_p'=-\gamma/(Fv).
\tag{W3.9}
\]

For a fixed angular position, the covariant stress tensor is

\[
T_{\mu\nu}(t,r,\Omega)=\frac{m_p u_\mu u_\nu}{u^t r^2}
\delta[r-r_p(t)]\delta^{(2)}(\Omega,\Omega_p).
\tag{W3.10}
\]

Define hat T=integral dt exp(+i omega t) T and inverse integral d omega/(2 pi). The delta-function transformation gives **1/|dr_p/dt|**, hence a positive Jacobian. With A_lm=Y_lm*(Omega_p) exp[i omega t_p(r)], the raw covariant harmonic amplitudes are

\[
\widehat T_{tt}=\frac{E\gamma}{r^2v}A_{lm},\quad
\widehat T_{tr}=\frac{E}{r^2F}A_{lm},\quad
\widehat T_{rr}=\frac{Ev}{\gamma r^2F^2}A_{lm},\quad
\widehat T=-\frac{E}{\gamma r^2v}A_{lm}.
\tag{W3.11}
\]

Both lowered velocity components are negative, so the covariant tr amplitude is positive. This is not a statement about the sign of a mixed component or a differently normalized tensor harmonic. Choosing the north pole makes only m=0 nonzero, with Y_l0*(0)=sqrt[(2ell+1)/(4pi)]. Radial motion supplies no axial source.

The raw amplitudes satisfy stress conservation exactly, including differentiation of the phase via t_p'. The library also computes the trace's first and second derivatives for (W3.3). It projects the scalar Hessian explicitly; consequently angular components of J may be nonzero for finite-mass infall even though the corresponding raw T components vanish. The trace/divergence right-hand sides must be supplied together with the wave source.

The ultrarelativistic benchmark is gamma→infinity **at fixed E**, with m_p→0. Then t_p'=−1/F, T=0, and (T_tt,T_tr,T_rr)=E A_lm(1,1/F,1/F²)/r². This limit preserves the test-particle energy scale. Letting gamma diverge at fixed particle mass would not.

AGP2026 uses a symmetric Fourier normalization (inverse 1/sqrt(2pi)); our raw Fourier amplitudes are sqrt(2pi) times symmetric-transform amplitudes of the same components. Tensor-harmonic normalization must additionally be converted before comparing paper source coefficients or emitted energies. An arbitrary trajectory time origin changes phases; the caller supplies t_p(r). The code does not integrate a retarded waveform or choose a physical excitation prescription.

## 6. Validation record and next gate

See [operator artifact](validation/w3_polar_operators.json), [validation results](validation/w3_validation.json), and [commands/code guide](../scripts/weyl_ringdown/README.md). Checks include exact stress conservation and null-limit identities; all seven compact equations; angular projections; massless GR intertwining; local reconstruction at ell=0,1,2,3,4,7 and r=2.3,4,11 for vacuum, null and massive sources; independent monopole spectral calculations; and numerical refinements. The sampled source and homogeneous reconstruction derivatives agree to about 10^−11 in the stated scaled norms. This is a numerical check, not a uniform error theorem.

The saved run passes all 54 local cases. The largest equation residual is 5.93×10^−13; the largest homogeneous jet-derivative residual is 2.90×10^−11. Monopole coefficient agreement is better than 9×10^−15 at the sampled points. The spectral acceptance threshold is an absolute frequency difference of 2×10^−8, with the same threshold on the final boundary refinement:

| mu | Omega from K matching | Independent extrapolated Omega | Absolute difference |
|---|---:|---:|---:|
| 0.3 | 0.034279461746 | 0.034279461738 | 7.70×10^−12 |
| 0.4 | 0.011220446774 | 0.011220446746 | 2.79×10^−11 |

These extra digits record numerical agreement, not physical precision or rigorous error bounds. The finite-difference reference uses spacings 0.1, 0.05, 0.025 on tortoise-coordinate bounds [−900,150], with an additional inner-boundary check at −700. The JSON records versions, the generated-operator hash, residuals and all individual settings.

W3 is complete at this documented scope. W4 must still establish regular boundary bases and damped polar branch benchmarks for the chosen response, construct the retarded Green function, reconstruct h=H+2q including the GR source, and validate a spectrum/waveform and flux normalization. No source amplitude, energy suppression, detector strain, rotation correction, quantum ghost prescription, or f2 posterior is inferred here. The PS scalar sector cannot supply these spin-two excitation coefficients.

## References

1. G. Antoniou, L. Gualtieri and P. Pani, *Nonperturbative suppression of beyond-General-Relativity effects in quadratic gravity*, [arXiv:2604.01631v2](https://arxiv.org/abs/2604.01631v2), 18 August 2026. Sections II.4–II.6, Eqs. (32)–(33), Appendix A, Appendix B.2, and Appendix C. This is the benchmark, not a substitute for the independent coordinate derivation.
2. R. Brito, V. Cardoso and P. Pani, *Massive spin-2 fields on black hole spacetimes: Instability of the Schwarzschild and Kerr solutions and bounds on the graviton mass*, [arXiv:1304.6725v3](https://arxiv.org/abs/1304.6725v3), Phys. Rev. D **88**, 023514 (2013), especially Eq. (30).
3. A. Held and J. Zhang, *Instability of spherically symmetric black holes in quadratic gravity*, [arXiv:2209.01867](https://arxiv.org/abs/2209.01867). Alternative radial-master formulation and instability context.
4. G. Antoniou, L. Gualtieri and P. Pani, *Gravitational quasinormal modes of black holes in quadratic gravity*, [arXiv:2412.15037v3](https://arxiv.org/abs/2412.15037v3), Phys. Rev. D **111**, 064059 (2025). Axial comparison and auxiliary-field conventions; axial results alone do not validate a radial-infall response.
