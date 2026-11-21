# Gate W1 — axial Schwarzschild equations and the GR benchmark

21 November 2026. Classical working derivation and reproducible numerical benchmark, following [W0](WEYL_RINGDOWN_CONVENTIONS.md). Excluded from TeX. W1 recovers a GR mode and derives the massive axial system; it does **not** compute a massive-mode spectrum, a sourced response, a rotating waveform or a bound on f2.

## 1. Covariant starting point and field choice

Use W0's Einstein–Weyl action with alpha_L=0, bar q=0, Ricci-flat Schwarzschild and finite positive spin-two mass mu=m2_pole. Set r_g=GM_BH=1 in all equations below; restore dimensions with r→r/r_g, omega→r_g omega, mu→r_g m2_pole. Thus the horizon is r=2, not r=1. Define F=1−2/r and d/dr_*=F d/dr. All fields have time dependence exp(−i omega t).

The linear equations and constraints are

\[
\delta G_{\mu\nu}(h)+\mu^2q_{\mu\nu}=0,\qquad
(\bar\Box-\mu^2)q_{\mu\nu}
+2\bar R_{\mu\rho\nu\sigma}q^{\rho\sigma}=0,
\quad q^\mu{}_{\mu}=0,\quad\bar\nabla^\mu q_{\mu\nu}=0.
\tag{W1.1}
\]

Using G^(1)(q)=−mu² q/2 on these solutions, H=h−2q obeys G^(1)(H)=0. We choose Regge–Wheeler gauge for **H**, not simultaneously for H and the physical metric h. The latter is reconstructed as h=H+2q. This keeps the GR sector homogeneous while preserving the massive contribution to the physical metric. A gauge-fixed equation for the original h, such as the one in the literature, can instead have explicit q sources; the two descriptions must not be mixed coefficient by coefficient without transforming their gauge choices.

## 2. Axial harmonics and the constraint before elimination

Let D_A be the unit-sphere derivative, L=l(l+1), and use

\[
S_A=(\csc\theta\,\partial_\phi Y,-\sin\theta\,\partial_\theta Y),
\quad S_{AB}=D_{(A}S_{B)},\quad
D^AS_A=0,\quad D^BS_{AB}=-\frac{L-2}{2}S_A.
\tag{W1.2}
\]

Suppress the l,m labels and set

\[
q_{tA}=F_0(r)S_A,\quad q_{rA}=F_1(r)S_A,
\quad q_{AB}=-2F_2(r)S_{AB},\quad q_{ab}=0.
\tag{W1.3}
\]

The trace vanishes identically. The time and radial components of the divergence vanish by the angular identities. Its angular component gives

\[
F_1'+\left(\frac{F'}F+\frac2r\right)F_1
+\frac{i\omega}{F^2}F_0+\frac{L-2}{r^2F}F_2=0.
\tag{W1.4}
\]

Define Q=F F1 and Z=F2/r. For omega≠0 the constraint fixes

\[
F_0=\frac{iF}{\omega}\left[Q'+\frac{2Q}{r}+\frac{L-2}{r}Z\right].
\tag{W1.5}
\]

Static perturbations must be treated directly with (W1.4), not by dividing by omega. In particular, the stationary angular-momentum perturbation is not covered by the reconstruction below.

## 3. Projected massive equations

Project the tensor wave operator in (W1.1) on S_A and −2S_AB, substitute (W1.5), and collect Q,Z. This gives

\[
\frac{d^2}{dr_*^2}\begin{pmatrix}Q\\Z\end{pmatrix}
+\left[\omega^2\mathbf1-F\mathbf W(r)\right]
\begin{pmatrix}Q\\Z\end{pmatrix}=0,
\tag{W1.6}
\]
\[
\mathbf W=\mu^2\mathbf1+\frac{\mathbf A}{r^2}+\frac{\mathbf B}{r^3},\quad
\mathbf A=\begin{pmatrix}L+4&2(L-2)\\2&L-2\end{pmatrix},\quad
\mathbf B=\begin{pmatrix}-16&-6(L-2)\\0&2\end{pmatrix}.
\tag{W1.7}
\]

In particular, the off-diagonal terms in the differential equation are −2(L−2)F(r−3)Z/r³ and −2FQ/r². They are not a scalaron coupling. The matrix need not be symmetric in these noncanonically normalized radial variables.

For a check that retains the original equations, call the two left-hand sides in (W1.6) E_Q,E_Z and call the covariant operator in (W1.1) E_mu_nu. Their projected components are

\[
E_{rA}/S_A=E_Q/F^2,\qquad
E_{AB}/(-2S_{AB})=rE_Z/F.
\tag{W1.8}
\]

The tA equation vanishes after using E_Q=E_Z=0, their required radial derivative, and the constraint. The saved symbolic verifier constructs the Schwarzschild connection and curvature directly, differentiates (W1.3), and checks these identities for m=0 at l=1,2,3. It does not import a radial potential from a QNM package. The angular identities in (W1.2) supply the general-l reduction; the coordinate checks are independent finite-l checks, not an exhaustive symbolic proof for every harmonic.

Exceptional multipoles:

- l=0: S_A and S_AB vanish; there is no axial monopole. The known monopolar stability question is polar and remains outside W1.
- l=1: S_AB vanishes. Set F2=Z=0 **before** forming a two-field system. A single massive equation remains with W=mu²+6/r²−16/r³. There is no second equation E_Z=0 to impose. The GR axial dipole has no radiative degree of freedom; stationary angular momentum and gauge perturbations require separate handling.
- l>=2: two massive axial radial fields and one massless GR axial field remain. These are the axial subset of W0's seven polarizations. No independent scalaron has been restored.

The covariant setup and literature comparison are Antoniou–Gualtieri–Pani (2025), section IV.1 and appendix B, and Brito–Cardoso–Pani (2013). The potentials here were checked from the covariant operator rather than accepted from rendered formulae alone.

## 4. GR reduction and metric reconstruction

In Regge–Wheeler gauge take H_tA=h0 S_A, H_rA=h1 S_A, H_AB=0. The tensor-angular Einstein equation gives

\[
i\omega h_0+F(Fh_1)'=0.
\tag{W1.9}
\]

Define psi=Fh1/r, which yields

\[
h_1=\frac rF\psi,\qquad
h_0=\frac{iF}{\omega}(r\psi)',\qquad
\frac{d^2\psi}{dr_*^2}+
\left[\omega^2-F\left(\frac L{r^2}-\frac6{r^3}\right)\right]\psi=0.
\tag{W1.10}
\]

The last equation is the remaining dynamical Einstein equation after substitution. Its potential has dimension inverse length squared before rescaling. Setting q=0 selects this GR branch; it is not the same operation as taking mu→0 in the massive system. The solver below searches only the l=2 fundamental of (W1.10).

## 5. Horizon and infinity series

These boundary expansions apply to (W1.6) and can also be specialized to the scalar GR equation by using W=L/r²−6/r³. They are derived by substitution; the helper `axial_system.py` implements the massive recurrences for l>=2.

**Horizon.** Put u=r−2 and U=exp(−i omega r_*) S(u), S=sum a_n u^n. The factored equation is

\[
F S''+(F'-2i\omega)S'-\mathbf W S=0.
\tag{W1.11}
\]

Write F=sum_{j>=1} f_j u^j, f_j=(−1)^(j−1)/2^j, and W=sum W_j u^j. Matching u^(n−1) gives

\[
n\left(\frac n2-2i\omega\right)a_n
=\sum_{j=0}^{n-1}\mathbf W_j a_{n-1-j}
-n\sum_{m=1}^{n-1}m f_{n+1-m}a_m.
\tag{W1.12}
\]

There are two independent ingoing leading vectors a0. In particular a1=W(2)a0/(1/2−2i omega). Exceptional frequencies at which the denominator vanishes require a separate Frobenius treatment. F1=Q/F has an extra inverse-u factor; F2=rZ does not. Regularity must be judged after the ingoing time coordinate is used, not from an individual Schwarzschild-coordinate component.

**Infinity.** Let k²=omega²−mu² and

\[
U=e^{ikr}r^\chi\sum_{n\geq0}\frac{b_n}{r^n},\qquad
\chi=\frac{i(2\omega^2-\mu^2)}{k}.
\tag{W1.13}
\]

The constant term fixes k; cancellation of the 1/r term fixes chi. The first correction is

\[
b_1=\frac{[(\chi^2-\chi-8ik\chi-4k^2+2ik)\mathbf1-\mathbf A]}{2ik}b_0.
\tag{W1.14}
\]

For higher orders set S=sum b_n/r^n, ell_p=ik+chi/r and substitute

\[
F^2S''+(2F^2\ell_p+FF')S'
+[F^2(\ell_p^2+\ell_p')+FF'\ell_p+\omega^2\mathbf1-F\mathbf W]S=0.
\tag{W1.15}
\]

At order r^(−n−1), the new coefficient is −2ik n b_n; lower coefficients determine the rest. The code implements this finite Laurent-polynomial recurrence. It is a local asymptotic expansion, not guaranteed convergent at arbitrary finite r or increasing order. The massless limit of the outgoing factor is exp(i omega r) r^(2i omega), providing a phase check. The threshold k=0 is excluded.

The retarded branch is fixed as in W0 by continuation from the upper omega half-plane. The numerical helper deliberately implements only continuation from real omega>mu into the positive-frequency above-threshold region; it rejects other cases rather than silently selecting a different sheet. Quasibound states use a different decaying-at-infinity boundary problem. No massive QNM roots are searched here.

## 6. Derived continued fraction for the GR solver

Use x=1−2/r and factor both endpoint behaviors exactly:

\[
\psi=e^{i\omega r}r^{2i\omega}x^{-2i\omega}A(x),\qquad
A(x)=\sum_{n=0}^{\infty}a_nx^n.
\tag{W1.16}
\]

At the horizon the nonanalytic factor is x^(−2i omega); at infinity the outgoing factor is exp(i omega r) r^(2i omega). Substitution in (W1.10) gives P A''+B A'+C A=0 with

\[
P=x(1-x)^2,
\]
\[
B=1-4i\omega+(-4+16i\omega)x+(3-8i\omega)x^2,
\]
\[
C=32\omega^2+8i\omega+3-L+(-16\omega^2-8i\omega-3)x.
\tag{W1.17}
\]

Matching powers of x produces alpha_n a_(n+1)+beta_n a_n+gamma_n a_(n−1)=0, with a_(−1)=0 and

\[
\alpha_n=(n+1)(n+1-4i\omega),
\]
\[
\beta_n=-2n^2+(-2+16i\omega)n+32\omega^2+8i\omega+3-L,
\quad\gamma_n=n^2-8i\omega n-16\omega^2-4.
\tag{W1.18}
\]

A minimal recurrence solution enforces the outgoing condition at x=1. With R_n=a_(n+1)/a_n, iterate R_n=−gamma_(n+1)/(beta_(n+1)+alpha_(n+1)R_(n+1)) downward from depth N and solve beta_0+alpha_0 R_0=0. The reference frequency is used only for comparison, never in the residual or coefficients. Two rough complex guesses seed a secant search. This is a targeted fundamental-mode solver, not a general overtone catalog.

There is no finite-radius boundary cutoff to vary: both endpoint phases are analytic factors. The numerical controls are recurrence depth and tail choice. Finite-radius reconstruction at r=2.1,3,6,10,20 checks the local equation independently of the recurrence residual. The saved code also checks the angular Einstein constraint using centered finite differences of the reconstructed metric, rather than canceling it algebraically by construction.

## 7. Reproduction commands and results

From the repository root, with Python 3 (standard library only for the numerical solver):

```sh
python3 scripts/weyl_ringdown/w1_gr_benchmark.py --output notes/validation/w1_gr_benchmark.json
python3 scripts/weyl_ringdown/check_w1_series.py > notes/validation/w1_boundary_series.json
```

Optional independent symbolic verification uses a temporary, pinned SymPy environment:

```sh
uv run --no-project --with sympy==1.14.0 python scripts/weyl_ringdown/verify_w1_algebra.py > notes/validation/w1_algebra.json
```

The GR result in units r_g=1 is

\[
\omega=0.3736716844180418-0.0889623156889357i.
\tag{W1.19}
\]

Comparison data: Emanuele Berti's `s2l2.dat`, n=0, gives 2M Re(omega)=0.7473433688360838 and a positive damping magnitude 2M |Im(omega)|=0.1779246313778714. Divide by two and use negative imaginary part for our time convention. Its quoted error is not a rigorous bound on our solver.

| Continued-fraction depth | Absolute difference from reference | Largest normalized RW residual over sampled radii |
|---:|---:|---:|
| 12 | 1.02e−5 | 3.35e−4 |
| 25 | 6.45e−8 | 9.54e−6 |
| 50 | 8.56e−11 | 3.58e−8 |
| 100 | 8.65e−15 | 2.71e−12 |
| 200 | 5.72e−17 | 2.60e−16 |
| 400 | 5.72e−17 | 2.60e−16 |

Differences at the last two depths reach floating-point roundoff; they are not a claim of certified sixteen-digit accuracy. Changing the terminal ratio among 0, 0.5 and 1 at depth 200, and varying the initial root guesses, preserves agreement below 1e−11. Centered-difference metric-constraint residuals decrease from 5.68e−5 to 8.87e−7 as the step decreases from 0.02 to 0.0025, consistent with second-order differencing. The local massive boundary series pass independent differential-operator residual checks at increasing expansion orders and changing evaluation radii; these checks use a generic frequency, not a massive eigenfrequency.

Machine-readable results are in `notes/validation/`. The symbolic checks confirm the trace, divergence, projected massive equations and remaining time-angular equation for the tested harmonics, plus the transformed GR ODE. This is stronger than comparing one frequency, but remains a working calculation requiring scientific review.

## 8. Scope and next gate

W1's acceptance criteria are met for the stated benchmark: covariant-to-axial reduction, exceptional multipoles, endpoint expansions, GR frequency convergence, and constraint checks. W2 still needs independent massive-mode root methods, named literature benchmark points, root tracking and the relevant stability analysis. The monopolar instability is not assessed by an axial calculation. No statement about rotating remnants, merger excitation, observed QQG corrections or a quantum ghost prescription follows from this result.

## References

1. Antoniou, Gualtieri and Pani, PRD **111**, 064059 (2025), [arXiv:2412.15037v2](https://arxiv.org/html/2412.15037v2), IV.1, appendix B and C: comparison of the axial conventions, constraints and spectra. Their physical-metric gauge differs from our homogeneous H choice.
2. Brito, Cardoso and Pani, PRD **88**, 023514 (2013), [arXiv:1304.6725](https://arxiv.org/abs/1304.6725): massive spin-two Schwarzschild perturbations and stability context.
3. Berti, [ringdown data and documentation](https://pages.jh.edu/eberti2/ringdown/), [s2l2.dat](https://pages.jh.edu/eberti2/ringdown/s2l2.dat), accessed 21 November 2026. The page links the associated references [arXiv:0905.2975](https://arxiv.org/abs/0905.2975) and [gr-qc/0512160](https://arxiv.org/abs/gr-qc/0512160). Only the n=0 row is used here.
