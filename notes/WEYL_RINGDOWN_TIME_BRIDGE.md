# W4 continuation — connecting the front and late-time windows

21 November 2026. Classical Einstein–Weyl working calculation, excluded from TeX. This continues the [sharp-front derivation](WEYL_RINGDOWN_NULL_FRONT.md) and [threshold-inclusive late-time inversion](WEYL_RINGDOWN_THRESHOLD_RESPONSE.md). The observable remains the angular auxiliary field \(K/\mathcal N\), where \(\mathcal N=\kappa E Y_{00}^{*}\). It is not detector strain.

## 1. The interval to be closed

The front calculation covered \(-0.1\le v\le1\), with \(v=t+r_*(r_{\rm obs})\). The earlier late-time calculation covered \(0\le t\le900\). At \(r_{\rm obs}=200\), \(r_*=209.1902397002692\), leaving \(1<v<209.1902397002692\) unvalidated.

The new bridge is tested on \(0.5\le v\le510\), or approximately \(-208.69024\le t\le300.80976\). It overlaps both existing windows. Retain \(GM=1\), \(\mu=GMm_2=0.5\), the same null trajectory, raw Fourier normalization, retarded boundary conditions and extended source-tail treatment.

The earlier front solver imposed two known moments by subtracting \(m_0+v m_1\) from a finite-band remainder. That is useful near passage, but extending the polynomial correction into the gap would introduce an artificial secular term. We therefore evaluate the ordinary residual inverse and estimate its high-frequency tail. We do not blend the two old curves or interpolate across missing time samples.

## 2. Subtraction, quadrature and the residual tail

Use the analytically invertible front reference already derived:

\[
P(v)=-\frac{\Theta(v)v}{r^2}K_0(\mu\sqrt{2rv}),\qquad
B(\omega)=e^{i\omega r_*}\widehat K(\omega),\qquad
\Delta B=B-\widehat P.
\]

With \(\omega=u+i\sigma\), the same retarded waveform is

\[
K(v)=P(v)+\frac{e^{\sigma v}}{\pi}{\rm Re}
\int_0^\infty \Delta B(u+i\sigma)e^{-iuv}\,du.
\tag{W4B.1}
\]

Removing the front reference makes the residual much smaller. Nevertheless, its Fourier integral is oscillatory: the upper-frequency samples have spacing 0.05, and at \(v=510\) the phase changes by 25.5 radians across one such interval. Direct low-order quadrature of those samples is inadequate.

The implementation constructs a complex quintic spline of \(\Delta B\), not of the rapidly oscillating product. For each frequency cell, it integrates the spline times the **exact oscillatory exponential** using a 16-node Gauss–Legendre rule:

\[
\int_{u_j}^{u_{j+1}}\!p_j(u)e^{-iuv}\,du
\simeq h_j\sum_{a=1}^{16}w_a
p_j(\bar u_j+h_j\xi_a)e^{-iv(\bar u_j+h_j\xi_a)},
\quad h_j=(u_{j+1}-u_j)/2.
\tag{W4B.2}
\]

Spline sampling error and integration error are separate. We test the former by coarsening the frequency grid and the latter by changing the Gauss order. The coarser-grid check uses 32 nodes to avoid confusing the two effects. No extra radial solutions are invented by the spline.

The initial cubic spline failed the \(2\times10^{-12}\) grid-refinement target in the stronger late overlap. Raising the spline degree to five resolved that failure; a separate seventh-degree calculation tests the remaining interpolation error. The original failed cubic-grid comparison is retained in the ledger. The acceptance target was not relaxed.

For the unsampled residual tail, repeated integration by parts gives

\[
\int_\Omega^\infty\Delta B(u+i\sigma)e^{-iuv}\,du
=e^{-i\Omega v}\sum_{n=0}^{N-1}
\frac{\partial_u^n\Delta B(\Omega+i\sigma)}{(iv)^{n+1}}
+\frac1{(iv)^N}\int_\Omega^\infty
\partial_u^N\Delta B(u+i\sigma)e^{-iuv}\,du.
\tag{W4B.3}
\]

Boundary terms at infinity vanish with the retarded prescription. The default estimate keeps \(N=3\); endpoint derivatives come from a degree-four polynomial fitted over the last unit of frequency. This is an asymptotic tail estimate, not a rigorous bound on its omitted remainder. We vary cutoff, retained orders, fit width and polynomial degree. Unlike the earlier late-time solver's discrete trapezoidal-tail formula, (W4B.3) corrects a continuous quadrature integral. Applying the old discrete correction here would double-count a different quadrature convention.

There is **no polynomial moment correction and no enforced time-domain zero** in the bridge solver. It accepts positive advanced times only; the validated domain begins at \(v=0.5\), not at the singular-slope front.

## 3. New data and independent checks

The main contour reuses the saved \(\sigma=0.004\) spectrum with spacing 0.001 below \(u=2\), and the sharp-front extension with spacing 0.05 from 2 to 12. The second contour now adds **1,000 new radial solutions** at the missing \(\sigma=0.003\) midpoints below 2. These use outer radius 800 and relative ODE tolerance \(2\times10^{-11}\); the older interleaved samples retain their recorded settings. The validation records all input hashes and solver versions.

The following acceptance targets apply to sampled times throughout the new interval:

- Required bridge comparisons: maximum absolute difference below \(2\times10^{-12}\) per \(\mathcal N\).
- Weak-response diagnostic interval \(10\le v\le400\): required differences below \(2\times10^{-13}\). This is an empirical numerical-floor criterion, not proof of an exactly empty signal interval.
- Front overlap \(0.5\le v\le1\): use the previous front calculation's declared tolerance, \(10^{-3}\) of its peak, approximately \(1.204\times10^{-10}\).
- Late overlap \(0\le t\le300.80976\): absolute difference below \(2\times10^{-12}\).
- Manufactured inverse: absolute error below \(10^{-13}\).

The manufactured test supplies the full transform

\[
B_{\rm test}=\widehat P+\frac{10^{-10}}{(s+0.7)^4},
\qquad s=-i\omega,
\]

whose exact inverse is \(P(v)+10^{-10}v^3e^{-0.7v}/6\). It checks phase, normalization, quadrature and tail signs without using the physical radial solution as its own reference.

Required variations include cutoffs 4 and 8 versus 12; one and two tail terms versus three; endpoint width 0.5 versus 1; endpoint polynomial degree three versus four; Gauss orders 12 and 20 versus 16; spline degree seven versus five; a coarser spectral grid; and the refined second contour. An eight-node rule and the old coarse second contour are retained as diagnostics, even if they fail a target. They are not used to certify the result.

The late overlap uses the independently formulated previous inversion of the **full** spectrum with its discrete endpoint correction. The front overlap uses the previous moment-corrected reconstruction. Neither comparison merely reruns the bridge function with another name. Their common radial data mean they are independent inversion checks, not independent derivations of the underlying field equations.

## 4. Reproduction and saved results

From the repository root:

~~~sh
uv run --no-project --with-requirements scripts/weyl_ringdown/requirements-w4.txt python scripts/weyl_ringdown/compute_w4_contour.py --sigma .003 --step .002 --offset .001 --maximum 2 --outer 800 --rtol 2e-11 --output notes/validation/w4_bridge_sigma003_midpoints.csv
uv run --no-project --with-requirements scripts/weyl_ringdown/requirements-w4.txt python scripts/weyl_ringdown/validate_w4_bridge.py
uv run --no-project --with numpy==2.4.6 --with matplotlib==3.10.8 python scripts/weyl_ringdown/plot_w4_bridge.py
~~~

The solver is [w4_bridge.py](../scripts/weyl_ringdown/w4_bridge.py). Outputs:

- [Validation ledger](validation/w4_bridge_validation.json), including input/source hashes and diagnostic failures.
- [Bridge waveform and numerical alternatives](validation/w4_bridge_waveform.csv).
- [Front overlap](validation/w4_bridge_front_overlap.csv).
- [Late overlap](validation/w4_bridge_late_overlap.csv).
- [Waveform and numerical-floor figure](validation/w4_time_bridge.png).

The saved run passes every required comparison:

| Check | Maximum absolute difference, per \(\mathcal N\) |
|---|---:|
| Cutoff 4 versus 12 | \(9.44\times10^{-13}\) |
| Cutoff 8 versus 12 | \(3.55\times10^{-13}\) |
| Coarser spectral grid, quintic spline | \(4.58\times10^{-13}\) |
| Refined second contour | \(1.20\times10^{-14}\) |
| Spline degree seven versus five | \(4.93\times10^{-15}\) |
| Gauss order 12 versus 16 | \(8.55\times10^{-16}\) |
| Independent late-time overlap | \(4.96\times10^{-15}\) |
| Manufactured inverse | \(9.05\times10^{-16}\) |

The front overlap differs by \(7.95\times10^{-11}\), or \(6.61\times10^{-4}\) of the previous front peak, within that calculation's \(10^{-3}\) tolerance. The bridge's stricter absolute target must not be retrospectively attributed to the older front reconstruction.

In \(10\le v\le400\), the largest required comparison is \(3.14\times10^{-14}\), while the main computed magnitude stays below \(9.69\times10^{-17}\). These values do not resolve a nonzero physical waveform in that interval. Quote the former as an empirical comparison scale, not the latter as a significant measurement or a rigorous error bound.

The retained failures are informative: cubic-grid refinement differed by \(1.35\times10^{-11}\), and an eight-node quadrature rule differed by \(4.77\times10^{-12}\). The original coarse second contour also failed the weak-interval target, with a discrepancy of \(7.21\times10^{-13}\). Quintic interpolation, higher-order quadrature and the new contour midpoints resolve those respective deficiencies.

## 5. Interpretation and remaining W4 work

The calculation tests whether one retarded auxiliary response connects the front and late-time windows. It does not assume that a small numerical value is an exact physical zero. In the weak middle interval, quote the observed comparison floor alongside the waveform; do not extract a damping law, exact arrival time or new suppression factor from values beneath that floor.

The accepted overlaps establish consistency at the documented absolute accuracy. They do not make the old front method as accurate as the bridge in its smaller late tail. In particular, the old front moment correction is material to its small curved remainder. The new solver supplies an alternative there without changing the historical saved front data.

A connected interval is not completion of W4. Remaining work includes a precisely identified author-spectrum/waveform comparison, metric reconstruction, signed flux normalization and a convergence-controlled mass/suppression scan. No detector strain, coupling bound, merger prediction or quantum-state prescription follows from the auxiliary waveform alone. This is still a Weyl spin-two calculation that the PS scalar toy theory cannot replace.

Literature and source conventions remain those cited in the [original null-infall benchmark](WEYL_RINGDOWN_NULL_INFALL.md). The bridge quadrature and overlap checks are this workbook's numerical continuation, not a new claim attributed to those references.
