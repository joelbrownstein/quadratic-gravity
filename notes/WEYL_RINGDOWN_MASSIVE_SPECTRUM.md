# Gate W2 — selected massive Schwarzschild modes and stability checks

21 November 2026. Classical numerical reproduction following [W0](WEYL_RINGDOWN_CONVENTIONS.md) and [W1](WEYL_RINGDOWN_AXIAL_DERIVATION.md). Excluded from TeX. W2 completes a **selected fundamental-mode benchmark**, not an exhaustive spectrum, a stability theorem, a sourced waveform or a quantum coupling constraint.

## 1. Target, literature anchors and units

Retain the vacuum Einstein–Weyl theory and Schwarzschild branch. Set r_g=GM_BH=1, mu=r_g m2_pole, p=r_h m2_pole=2 mu and exp(−i omega t). Report omega in r_g units; multiply by two before comparing a paper's r_h omega axis.

The target is the l=2, n=0 tensor-led and vector-led branches at **p=0.8, 0.9, 1.0**. These are named points on the curves in Antoniou–Gualtieri–Pani (AGP), [arXiv:2412.15037v3](https://arxiv.org/abs/2412.15037v3), Figures 4 and 3 respectively, with the Appendix C recurrence as the coefficient-level literature benchmark. The points were fixed after exploratory continuation and before the production validation suite, not preregistered before all exploration.

The figures do not supply exact tabulated frequencies for these points. Accordingly, we do **not** invent digitized author values or treat visual agreement as a precision test. We separately implement the published recurrence to produce a reproducible literature-equation reference, compare it with our independently derived recurrence, and then compare both with direct integration of W1's differential equations. The independent numerical methods are the continued fraction and contour matching; the two recurrence implementations are an additional consistency check, not a third independent method.

For the separate radial test we import the polar-monopole master potential from Brito–Cardoso–Pani (BCP), [arXiv:1304.6725v3](https://arxiv.org/abs/1304.6725v3), Eq. (30). It is not derived from the axial system. The later QQG literature quotes the transition near p=0.876; compare also Held–Zhang, [arXiv:2209.01867](https://arxiv.org/abs/2209.01867). This is a classical perturbation benchmark, not a transfer of a graviton-mass observational bound from a different quantum theory.

## 2. Matrix continued fraction derived from W1

The two-field equation is U_** + [omega² I−F W] U=0, with F=1−2/r and W=mu² I+A/r²+B/r³. The matrices A,B are those independently checked in W1. Put

\[
x=1-2/r,\quad k^2=\omega^2-\mu^2,\quad
\chi=\frac{i(2\omega^2-\mu^2)}{k},\quad
U=e^{ikr}r^\chi x^{-2i\omega}\sum_{n=0}^{\infty}a_nx^n.
\tag{W2.1}
\]

The first factor enforces outgoing asymptotics; x^(−2i omega) enforces an ingoing horizon. Direct substitution yields

\[
x(1-x)^2 a''+(b_0+b_1x+b_2x^2)a'
+\left[(c_0+c_1x)I-A-\frac B2(1-x)\right]a=0,
\tag{W2.2}
\]

where

\[
\begin{split}
b_0&=1-4i\omega,\\
b_1&=6ik+8i\omega-4+2i\omega^2/k,\\
b_2&=-2ik-4i\omega+3-2i\omega^2/k,\\
c_0&=4k^2+12k\omega+3ik+12\omega^2+4i\omega
+4\omega^3/k+i\omega^2/k,\\
c_1&=-k^2-4k\omega-2ik-6\omega^2-4i\omega
-4\omega^3/k-2i\omega^2/k-\omega^4/k^2.
\end{split}
\tag{W2.3}
\]

Thus alpha_n a_(n+1)+beta_n a_n+gamma_n a_(n−1)=0, with

\[
\begin{split}
\alpha_n&=(n+1)(n+b_0)I,\\
\beta_n&=[-2n(n-1)+b_1n+c_0]I-A-B/2,\\
\gamma_n&=[(n-1)(n-2)+b_2(n-1)+c_1]I+B/2.
\end{split}
\tag{W2.4}
\]

No simultaneous diagonalization of A and B is assumed. Starting with R_N=t I at recurrence depth N, iterate

\[
R_n=-(\beta_{n+1}+\alpha_{n+1}R_{n+1})^{-1}\gamma_{n+1},
\qquad \det(\beta_0+\alpha_0R_0)=0.
\tag{W2.5}
\]

The default tail t=0 is varied to 0.5 and 1. The solver varies N from 100 to 1600. `verify_w2_recurrence.py` checks (W2.2) directly against the original differential equation as an exact symbolic identity for arbitrary scalar potential coefficients, then checks coefficient extraction. Since the potential enters linearly and the differential pieces multiply the identity, the same identity holds for the matrix coefficients. The published recurrence uses its own radial basis and q=−ik; we compare its roots, not assume that every matrix entry must equal (W2.4).

## 3. Independent direct integration and the resonance sheet

The second method solves the coupled differential equation in r with DOP853. Two ingoing basis columns start at r=2+epsilon using W1's horizon series. Two outgoing columns start at r_out=r_match+R exp(i theta) using W1's infinity series. They are integrated to a common real r_match. The determinant of the four matching columns (U,U') vanishes at a QNM. Column normalization reduces amplitude scaling effects; a singular-value check verifies that the smallest singular value is isolated from the next one.

The outgoing contour lies in a decay wedge Im[k exp(i theta)]>0. It avoids the singularities r=0,2; its positive real part fixes continuous logarithmic branches. Analytic continuation from this contour implements the outgoing resonance condition without trying to subtract exponentially disparate incoming and outgoing amplitudes at large real r. This is independent of the continued-fraction discretization, though both necessarily share the same physical differential equation and boundary conditions.

Production convergence varies all of the following:

| Setting | Coarse | Intermediate | Fine |
|---|---:|---:|---:|
| Matching radius | 5 | 5 | 6 |
| Contour length R | 35 | 55 | 70 |
| Contour angle theta (radians) | 1.40 | 1.40 | 1.45 |
| Horizon epsilon | 0.04 | 0.02 | 0.01 |
| Boundary-series order | 8 | 12 | 16 |
| Integration relative tolerance | 1e−8 | 2e−10 | 2e−12 |

These combined refinements give an empirical error check, not a rigorous independent error bound for each control. An additional perturbed starting frequency tests that the matching solver can recover the root rather than merely accept a continued-fraction seed. Roots are accepted only with small matching residuals, method agreement and the singular-value check.

**A necessary extension of W1:** its helper deliberately rejected Re(omega)<=mu. Some tensor-led resonances reach that region. In W2, k is continued from real omega>mu into the lower half-plane with Re(k)>0, Im(k)<0 along the tracked positive-frequency damped branch. A condition on Re(omega) alone does not turn a complex resonance into a quasibound state. Bound states require a different decaying-at-real-infinity sheet. We do not search that spectrum. The new `wave_number` argument lets W2 supply its chosen sheet to the existing boundary expansion, while retaining W1's default restriction. Neither k=0 nor crossings of the real-axis branch cut are included.

## 4. Mode identification, selected results and convergence

Forward continuation uses increments delta mu=0.025 from the limiting tensor-led GR frequency and vector-led electromagnetic frequency. The massless endpoint labels a branch mathematically; it is not interpreted as a finite-coupling QQG theory. Reversing the continuation gives maximum forward/reverse differences of approximately 1.3e−14 and 2.5e−13, respectively. This guards against accidental branch switching over the sampled interval; it is not a proof of globally unique mode labels.

The selected values below are rounded to nine decimals, in units **GM omega**:

| mu | p | Tensor-led fundamental | Vector-led fundamental | Radial background classification |
|---:|---:|---|---|---|
| 0.40 | 0.80 | 0.410399881 − 0.065864412 i | 0.540782087 − 0.070466380 i | Monopole unstable |
| 0.45 | 0.90 | 0.422007029 − 0.056373826 i | 0.561112681 − 0.065151916 i | Above the monopole threshold |
| 0.50 | 1.00 | 0.435090184 − 0.044425394 i | 0.583347597 − 0.059469216 i | Above the monopole threshold |

Every listed axial root is damped. Both independent methods agree comfortably within the roadmap's 1e−4 relative target; the production assertion is stricter, 1e−6. The largest fine-matching/continued-fraction relative discrepancy is about 1.55e−11, at the tensor-led mu=0.50 point. Its continued fraction converges more slowly: N=800 to 1600 changes the root by 5.74e−9 in absolute units, while the other selected points are much less sensitive. The note therefore does not infer accuracy from a small determinant alone.

A 50-decimal arithmetic calculation at fixed N=1600 agrees with the corresponding double-precision tensor-led mu=0.50 root to about 4.2e−17. This separates arithmetic roundoff from truncation error; it does **not** establish 40- or 50-digit accuracy of the physical eigenfrequency. Values agreeing at machine precision in other comparisons likewise do not constitute certified error bounds.

The published-coefficient recurrence reproduces the same branches to the validation tolerance. Exact data, settings, determinant residuals, singular values, tail tests, continuation records and version information are saved in `notes/validation/w2_massive_spectrum.json`. No precision comparison relies on reading coordinates from a plotted curve.

## 5. Separate polar-monopole stability check

For the imported BCP Eq. (30), use

\[
V_0=F\left[\frac2{r^3}+\mu^2+
\frac{24(1-r)\mu^2+6r^3(r-4)\mu^4}{(2+r^3\mu^2)^2}\right].
\tag{W2.6}
\]

A negative eigenvalue of −d²/dr_*²+V0 is omega²<0. It gives an instability omega=+i sqrt(−omega²): the time factor grows, while the eigenfunction decays towards both spatial ends. This is a different boundary-value problem from the outgoing damped axial QNMs.

We discretize this real self-adjoint operator with a centered tridiagonal difference matrix and Dirichlet boundaries, comparing tortoise-coordinate boxes from [−300,100] to [−700,200] and spacings from 0.1 to 0.025. The inverse tortoise map uses the real Lambert-W function; F is evaluated without subtractive cancellation near the horizon. The finest-grid growth rates are approximately

\[
GM\,\Omega(\mu=0.30)=0.0342797,\qquad
GM\,\Omega(\mu=0.40)=0.0112207.
\tag{W2.7}
\]

Here mu is already GM m2; the displayed GM on Omega restores the dimensionless frequency notation. Negative eigenvalues persist under domain and grid refinement. At mu=0.45 and 0.50 no negative eigenvalue is found in the tested boxes. Positive finite-box eigenvalues move toward the continuum edge as the box grows; they are **not** interpreted as a physical positive spectral gap or a proof of all-mode stability.

An independent zero-frequency shooting calculation imposes a regular constant solution at the horizon and a decaying solution at infinity, matching logarithmic derivatives. For omega=0 the large-r radial behavior is exp(−mu r) r^(−mu); this fixes the asymptotic initial derivative used by the solver. Varying the left/right cutoffs, matching location and integration tolerance gives

\[
\mu_{\rm crit}\simeq0.438080201,
\qquad p_{\rm crit}\simeq0.876160402.
\tag{W2.8}
\]

This reproduces the literature's approximately 0.876 horizon-scaled threshold. Across the tested zero-mode settings the p value changes by about 2.1e−9. That is empirical numerical convergence, not a rigorous uncertainty. The result demonstrates why the p=0.8 axial benchmark must be labeled a perturbation of a radially unstable background, even though its axial frequencies have negative imaginary parts.

## 6. Reproduction and acceptance

From the repository root:

```sh
uv run --no-project --with-requirements scripts/weyl_ringdown/requirements-w2.txt python scripts/weyl_ringdown/validate_w2.py --output notes/validation/w2_massive_spectrum.json
uv run --no-project --with sympy==1.14.0 python scripts/weyl_ringdown/verify_w2_recurrence.py > notes/validation/w2_recurrence_algebra.json
```

Pinned production dependencies: NumPy 2.4.6, SciPy 1.17.1, mpmath 1.3.0; symbolic verification uses SymPy 1.14.0 separately. The numerical suite writes its JSON after all acceptance assertions pass, and progress to stderr. Do not disable assertions with Python `-O`.

Acceptance passed: recurrence/operator identity; six selected roots with two independent numerical methods; depth, tail, endpoint and integration refinements; perturbed root seeds; forward/reverse branch tracking; arithmetic-precision check; negative radial eigenvalues below the threshold and a converged zero-mode transition. The reference frequencies are outputs of the literature recurrence, not hidden inputs to either solver's equations.

Limits remain explicit: only selected l=2 fundamentals are computed; no complete overtone survey, polar radiative spectrum, nonlinear stability test or quantum state prescription is supplied. The monopole potential is a cited input whose full polar reconstruction belongs to W3. No source amplitudes or detector signals follow from this spectrum alone. W3 is next: derive the polar system and source normalization before W4's sourced response.
