# W4 continuation — the sharp null-particle-passage front

21 November 2026. Working derivation and reproducible numerical benchmark, excluded from TeX. Read the [source derivation](WEYL_RINGDOWN_NULL_INFALL.md) and [threshold continuation](WEYL_RINGDOWN_THRESHOLD_RESPONSE.md) first. The result concerns the angular **auxiliary massive spin-two field**, not detector strain.

Subsequent work: the [time bridge](WEYL_RINGDOWN_TIME_BRIDGE.md) now validates the interval to the earlier late-time window, with an explicit front overlap and without extending this note's polynomial moment correction. The original local validation and saved data below are retained unchanged.

## 1. The question and conventions

A finite frequency band blurs a sharp particle-passage feature. Increasing that band alone is expensive, and a late-time inverse-frequency tail expansion is singular at passage. Here we derive the leading characteristic behavior, invert it analytically, and compute the remaining Schwarzschild response.

Retain the Einstein–Weyl benchmark, \(GM=1\), \(F=1-2/r\), \(x=r_*=r+2\log(r/2-1)\), \(\mu=GMm_2=0.5\), and \(r_{\rm obs}=200\). The particle has fixed null energy and trajectory \(t_p=-x\). Define \(v=t+x_{\rm obs}\). Passage occurs at \(v=0\), or \(t=-209.1902397002692\) at this observer. Every amplitude below is divided by \(\mathcal N=\kappa E Y_{00}^{*}\), with \(\kappa=8\pi G\). The inverse Fourier convention is \(e^{-i\omega t}\); no symmetric-transform factor is silently introduced.

The tensor amplitudes are \(q_{tt}=FaY_{00}\), \(q_{tr}=bY_{00}\), \(q_{rr}=cY_{00}/F\), and \(q_{AB}=r^2K\Omega_{AB}Y_{00}\). Null trace sourcing vanishes, so \(c=a-2K\). The source components divided by \(\mathcal N\) are \((1,F^{-1},F^{-2})\delta(v)/r^2\).

## 2. Deriving the characteristic behavior

Work at fixed exterior radius, with derivatives in \((v,r)\). Retarded support requires zero response at \(v<0\); there is no independently added incoming homogeneous field. Allow logarithmic characteristic data:

\[
a=b=c=\Theta(v)\,[A(r)\log v+B(r)]+\text{subleading terms},
\qquad
K=\Theta(v)v[C(r)\log v+D(r)]+o(v).
\tag{W4F.1}
\]

Equality of \(a,b,c\) here is a **leading-order statement**, not an ansatz for the full solution. Dimensionful logarithms mean their dimensionless \(GM\)-scaled arguments throughout.

The coefficient of the characteristic derivative in the W3 \(tt\) wave equation is

\[
\mathcal T[g]=2F g'+\frac{2g}{r}.
\tag{W4F.2}
\]

This is checked directly from the saved coordinate-derived W3 operator, after replacing a fixed-\(t\) radial derivative by \(\partial_r+F^{-1}\partial_v\). Distributional derivatives of \(\Theta(v)\log v\) are understood as finite-part distributions with the same reference scale as \(\log v\). Their coefficients must vanish; the ordinary step term supplies the prescribed delta source. Thus

\[
\mathcal T[A]=0,\qquad \mathcal T[B]=\frac1{r^2},
\qquad
A=\frac{a_\infty}{r-2},\quad
B=\frac{\tfrac12\log r+c_\infty}{r-2}.
\tag{W4F.3}
\]

Local transport alone does not determine the constants. Match to the asymptotically flat **retarded null-infall solution**, without an extra free characteristic disturbance. In flat space, the radial Klein–Gordon retarded integral gives, for \(0<v<2r\),

\[
a_{\rm flat}=-\frac1{2r}\int_{r-v/2}^{\infty}
 \frac{dr'}{r'}J_0\!\left(\mu\sqrt{2v(r'-r+v/2)}\right)
=-\frac{K_0(\mu\rho)}r,\qquad \rho^2=v(2r-v).
\tag{W4F.4}
\]

The oscillatory integral is defined by the retarded limiting prescription. The associated tensor amplitudes satisfy \(b=-ta/r\), \(c=t^2a/r^2\), and \(K=\rho^2a/(2r^2)\). Away from the source their equations reduce to \(f''+f'/\rho-\mu^2f=0\), with \(f=K_0(\mu\rho)\). Expanding at the front selects

\[
a_\infty=\frac12,\qquad
c_\infty=\log(\mu/\sqrt2)+\gamma_E.
\tag{W4F.5}
\]

This matching is a boundary-condition input, not a consequence of a local differential equation alone.

Add the two divergence constraints before taking the leading front limit. They give

\[
K_v=F a_r+(F'+2F/r)a+\text{subleading terms}.
\tag{W4F.6}
\]

Substitution of (W4F.3)–(W4F.5) yields

\[
\boxed{\displaystyle
K(v,r)=\frac{\Theta(v)v}{2r^2}
\left[\log\!\left(\frac{\mu^2rv}{2}\right)+2\gamma_E\right]+o(v).}
\tag{W4F.7}
\]

In particular, \(K(0)=0\), but its one-sided derivative diverges logarithmically. This is not a finite jump or a delta spike in \(K\). Other auxiliary components have logarithmic singularities. No conclusion about finite curvature invariants or quantum singularity resolution follows.

This expansion takes \(v\to0^+\) at fixed \(\mu>0\) and fixed exterior \(r\). A massless limit or a horizon limit cannot be interchanged with it without a separate analysis.

The cancellation of explicit factors of \(F\) in these first two coefficients is special to this characteristic calculation. It does **not** make the full curved-space response equal to its flat-space counterpart.

## 3. An analytic reference with an exact transform

Choose a causal reference that has both terms in (W4F.7) and decays for positive \(v\):

\[
P(v)=-\frac{\Theta(v)v}{r^2}K_0\!\left(\mu\sqrt{2rv}\right).
\tag{W4F.8}
\]

It is an analytic subtraction function, not the exact Schwarzschild solution or the full flat-space solution at all times. Let \(s=-i\omega\), \({\rm Re}\,s>0\), \(a_*=\mu^2r/2\), \(z=a_*/s\). Its Laplace transform is

\[
\widehat P(\omega)=
-\frac{(1+z)e^zE_1(z)-1}{2r^2s^2}.
\tag{W4F.9}
\]

The integral representation of \(K_0(2\sqrt{a_*v})\) gives the transform \(e^zE_1(z)/(2s)\); differentiating with respect to \(s\) supplies the factor \(v\). The validation script independently integrates (W4F.8) and compares (W4F.9) with 70-digit arithmetic.

At large \(|z|\), direct evaluation overflows or loses precision through cancellation. The implementation evaluates the combined numerator using

\[
(1+z)e^zE_1(z)-1
\sim\sum_{n=2}^{21}(-1)^n(n-1)(n-1)!\,z^{-n}.
\tag{W4F.10}
\]

The numerical switch is \(|z|=60\). This expansion is used only on the causal contour in its tested sector.

## 4. Reconstructing the curved remainder

The saved W4 response contains the trajectory phase. Remove it before subtracting:

\[
B(\omega)=e^{i\omega x_{\rm obs}}\widehat K(\omega),\qquad
\Delta B=B-\widehat P.
\]

For \(\omega=u+i\sigma\), the full inverse is

\[
K(v)=P(v)+\frac1\pi{\rm Re}\int_0^\infty
\Delta B(u+i\sigma)e^{(\sigma-iu)v}\,du.
\tag{W4F.11}
\]

The reference treats arbitrarily high frequencies analytically. The full radial solver supplies the curved remainder. Its source-tail boundary condition remains essential: this calculation does not replace the extended infall source by a compact source.

Because (W4F.7) fixes both the \(v\log v\) and \(v\) terms, the remainder obeys \(R(0)=R'(0)=0\). For a finite integration cutoff, we impose those two known moments:

\[
R_\Omega(v)=\frac1\pi{\rm Re}\int_0^\Omega
\Delta B(u+i\sigma)
\left[e^{sv}-1-sv\right]\,du,\qquad s=\sigma-iu.
\tag{W4F.12}
\]

We also impose retarded zero support for \(v\le0\). These are analytic constraints, **not independent numerical evidence of causality**. The saved waveform includes the uncorrected inverse; the ledger records the removed zeroth and first moments and prefront leakage. Frequency-cutoff, grid and contour changes test the finite-cutoff approximation. At fixed nonzero \(v\), the missing tail remains an error; moment subtraction does not prove it vanishes.

This method is local to the tested front window. Extending the polynomial correction to arbitrarily late times would be inappropriate. Use the earlier validated post-passage inversion for \(0\le t\le900\); do not splice the two windows without an overlap/convergence study.

## 5. Reproduction and acceptance

From the repository root, using the pinned W4 requirements:

~~~sh
uv run --no-project --with-requirements scripts/weyl_ringdown/requirements-w4.txt python scripts/weyl_ringdown/compute_w4_front.py --sigma .004 --output notes/validation/w4_front_high_sigma004.csv
uv run --no-project --with-requirements scripts/weyl_ringdown/requirements-w4.txt python scripts/weyl_ringdown/compute_w4_front.py --sigma .003 --output notes/validation/w4_front_high_sigma003.csv
uv run --no-project --with-requirements scripts/weyl_ringdown/requirements-w4.txt python scripts/weyl_ringdown/validate_w4_front.py --radial-only
uv run --no-project --with-requirements scripts/weyl_ringdown/requirements-w4.txt python scripts/weyl_ringdown/validate_w4_front.py
~~~

These commands reuse the saved low-frequency spectra and their solver hashes. The added spectra sample \(2\le u\le12\) with spacing 0.05. They use \(r_{\rm outer}=\max(800,200u)\), relative ODE tolerance \(2\times10^{-11}\), and the existing scaled Green solver. The independent radial refinements double the outer boundary and tighten tolerance and expansion orders at \(u=2,6,12\).

The local waveform grid spans \(-0.1\le v\le1\), with logarithmically spaced positive times down to \(10^{-9}\). The acceptance target for cutoff 4 and 8 versus 12, doubled grid spacing, and the independent \(\sigma=0.003\) contour is maximum difference below \(10^{-3}\) of the front peak. This is an **absolute, peak-normalized** convergence target, not a pointwise relative-error guarantee near a zero, nor a rigorous bound on all uncomputed frequencies. The analytical limiting behavior supplies the arbitrarily close-to-front statement.

The dominant front and the much smaller curved remainder have different relative accuracies. Passing the peak-normalized test does not establish a comparably precise fractional measurement of the remainder or of the weak signal at the far end of this window.

See [the validation ledger](validation/w4_front_validation.json), [waveform samples](validation/w4_null_front.csv), and [radial refinements](validation/w4_front_radial_checks.json) for measured results. The older w4_front_unresolved.csv remains a historical raw truncated inversion and must not be mistaken for the new result.

The saved run passes the declared comparisons:

| Check | Measured discrepancy |
|---|---:|
| Transform versus 70-digit evaluation | \(1.49\times10^{-12}\) relative or less |
| Transform versus direct Laplace quadrature | \(2.59\times10^{-13}\) relative or less |
| Refined radial boundaries/orders/tolerance, \(u=2,6,12\) | \(4.37\times10^{-9}\) relative or less |
| Cutoff 4 versus 12 | \(4.12\times10^{-4}\) of front peak |
| Cutoff 8 versus 12 | \(1.09\times10^{-4}\) of front peak |
| Doubled frequency spacing | \(2.13\times10^{-6}\) of front peak |
| Contour height 0.003 versus 0.004 | \(1.08\times10^{-5}\) of front peak |

The sampled peak magnitude is \(1.20378\times10^{-7}\) per \(\mathcal N\), near \(v=0.0239883\). Before imposing the analytic constraints, the maximum prefront leakage is \(4.34\times10^{-12}\); the removed moments are \(m_0=-4.34\times10^{-12}\) and \(m_1=-7.54\times10^{-11}\). The largest corrected remainder on the positive window is \(7.48\times10^{-11}\). Keeping these numbers visible matters: the first-moment correction is material to the small remainder even though it is tiny relative to the dominant front.

The [front figure](validation/w4_null_front.png) is generated from the saved waveform by:

~~~sh
uv run --no-project --with numpy==2.4.6 --with matplotlib==3.10.8 python scripts/weyl_ringdown/plot_w4_front.py
~~~

## 6. What this closes and what remains

The front's leading nonanalytic structure is now derived, and an analytic-subtraction solver makes its local profile testable without integrating the entire dominant high-frequency response numerically. This advances the null-infall auxiliary benchmark only. Full W4 still requires a precisely identified author-spectrum/waveform comparison, metric reconstruction, signed flux normalization and a convergence-controlled suppression scan. The gap between the local front window and the earlier late-time window is not newly validated here.

The PS scalar toy theory cannot supply this spin-two characteristic response. Nor can this calculation alone produce a bound on \(f_2\), a binary-merger waveform, or a quantum physical-state prescription.

### Literature and provenance

The tensor and source conventions are inherited from the primary references documented in [W3](WEYL_RINGDOWN_POLAR_DERIVATION.md) and [initial W4](WEYL_RINGDOWN_NULL_INFALL.md), including Brito, Cardoso and Pani, [arXiv:1304.6725](https://arxiv.org/abs/1304.6725), and the AGP2026 benchmark identified there. Equations (W4F.1)–(W4F.12) are this workbook's derivation, not a claim that those papers supplied this front subtraction. Special-function definitions follow [NIST DLMF, modified Bessel functions](https://dlmf.nist.gov/10.25) and [exponential integrals](https://dlmf.nist.gov/6.2).
