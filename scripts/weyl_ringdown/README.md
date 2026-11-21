# Weyl Schwarzschild benchmarks

See [the derivation](../../notes/WEYL_RINGDOWN_AXIAL_DERIVATION.md) and [W0 conventions](../../notes/WEYL_RINGDOWN_CONVENTIONS.md). All calculations use `GM=1`, horizon radius `2`, and `exp(-i omega t)`.

From the repository root:

```sh
python3 scripts/weyl_ringdown/w1_gr_benchmark.py --output notes/validation/w1_gr_benchmark.json
python3 scripts/weyl_ringdown/check_w1_series.py > notes/validation/w1_boundary_series.json
uv run --no-project --with sympy==1.14.0 python scripts/weyl_ringdown/verify_w1_algebra.py > notes/validation/w1_algebra.json
```

The first two commands use only Python's standard library. The optional third command uses SymPy 1.14.0 in an isolated environment; it does not change system Python dependencies. Commands exit unsuccessfully if their acceptance assertions fail. Do not run them with Python assertion disabling (`-O`).

- `w1_gr_benchmark.py`: targeted l=2 GR fundamental root, continued-fraction convergence, reconstruction and finite-difference metric-constraint checks. The comparison reference is not used to calculate the root.
- `axial_system.py`: coupled massive axial potential and local horizon/infinity series for l>=2. Its outgoing-wave helper supports only a specified above-threshold sheet; it rejects unsupported frequencies.
- `check_w1_series.py`: local-series residuals evaluated in the differential operator at a generic frequency, not an eigenfrequency.
- `verify_w1_algebra.py`: constructs the connection/curvature and verifies projected massive equations and constraints for l=1,2,3; checks the GR continued-fraction ODE transformation. `tensor_angular_equation: false` for l=1 means that harmonic does not exist, not that a test failed.

The saved JSON files are reproducible validation records. Floating-point agreement is not a rigorous error bound.

## W2 massive-mode benchmark

See [the W2 derivation and stability analysis](../../notes/WEYL_RINGDOWN_MASSIVE_SPECTRUM.md).

```sh
uv run --no-project --with-requirements scripts/weyl_ringdown/requirements-w2.txt python scripts/weyl_ringdown/validate_w2.py --output notes/validation/w2_massive_spectrum.json
uv run --no-project --with sympy==1.14.0 python scripts/weyl_ringdown/verify_w2_recurrence.py > notes/validation/w2_recurrence_algebra.json
```

`w2_massive.py` implements a matrix continued fraction and independent complex-contour boundary matching. `w2_stability.py` tests the separately cited polar-monopole potential using finite-box eigenvalues and a zero-mode matching calculation. `validate_w2.py` covers six fixed l=2 fundamental-mode points, literature-recurrence consistency, forward/reverse continuation, numerical refinements and a 50-digit arithmetic check. It prints progress to stderr and writes JSON only after passing its assertions.

W2 supplies an explicit outgoing resonance-sheet wave number to the boundary helper; W1's default above-threshold restriction remains unchanged. Neither script is a general branch-sheet or quasibound-state solver. No sourced response, rotation or observational inference is implemented. Damped selected axial modes do not establish complete stability.

## W3 polar equations and infall source conventions

See [the polar derivation](../../notes/WEYL_RINGDOWN_POLAR_DERIVATION.md). Run these commands in order, from the repository root:

```sh
uv run --no-project --with sympy==1.14.0 python scripts/weyl_ringdown/w3_polar_system.py --output notes/validation/w3_polar_operators.json
uv run --no-project --with sympy==1.14.0 --with-requirements scripts/weyl_ringdown/requirements-w2.txt python scripts/weyl_ringdown/validate_w3.py --output notes/validation/w3_validation.json
```

- `w3_polar_system.py` derives the constrained polar tensor equations in Schwarzschild coordinates for ell=0 through 4 and extracts the general ell≥2 operator. Exact trigonometric reduction removes angular dependence without evaluating at a special angle.
- `w3_polar.py` reduces the generated equations, reconstructs eliminated fields/derivatives, and checks inconsistent/singular reductions. `PolarSystem(ell).ode(r,omega,mu,rhs)` returns a local first-order matrix, source vector and diagnostics. `reconstruction` returns every field's first three radial jets. Its separate monopole matching helper solves growing bound states using a published equation checked against this coordinate-derived reduction.
- `w3_infall_source.py` supplies raw covariant Fourier amplitudes, trace derivatives, and projected forcing for null or finite-energy radial infall. Supply the trajectory time and normalized harmonic factor explicitly; defaults isolate the radial normalization. Units, transform convention and the fixed-energy null limit are documented in the derivation.
- `validate_w3.py` checks the compact equations, horizon series, exact conservation, GR intertwining, reconstructed-jet derivative compatibility, sampled low/high multipoles and sources, and two unstable polar-monopole frequencies against an independently formulated finite-difference eigenproblem. JSON is written only after all assertions pass; it includes the operator artifact hash.

This provides local sourced equations, not a retarded response. The source is a test particle on Schwarzschild in Einstein–Weyl theory; finite-energy infall in generic R²+C² also needs the scalaron. W3's spectral comparison is the unstable branch of AGP2026 Appendix B.2. The damped polar branches of Appendix B.1 and regular boundary bases needed by a chosen W4 waveform still require validation. Do not use the local least-squares reducer through a singular coordinate/parameter limit without a separately checked limiting formulation.

## W4 initial null-infall monopole response

See [the derivation, source audit and limitations](../../notes/WEYL_RINGDOWN_NULL_INFALL.md). The initial benchmark is validated; the full W4 gate remains open.

```sh
uv run --no-project --with-requirements scripts/weyl_ringdown/requirements-w4.txt python scripts/weyl_ringdown/verify_w4_algebra.py > notes/validation/w4_source_algebra.json
uv run --no-project --with-requirements scripts/weyl_ringdown/requirements-w4.txt python scripts/weyl_ringdown/validate_w4.py --output-directory notes/validation
uv run --no-project --with numpy==2.4.6 --with matplotlib==3.10.8 python scripts/weyl_ringdown/plot_w4.py
```

- `verify_w4_algebra.py`: exact null-source elimination from W3, remaining-wave-equation check, Abel identity and radial/tortoise source conversion. In particular, it retains the derivative of 2K' in the radial constraint.
- `w4_monopole.py`: retarded monopole Green integration and independent sourced shooting, horizon and outgoing asymptotics, and forced-source tails. The validated interface requires omega>mu and mu≥0.44. `response` returns K, its tortoise derivative, the asymptotic outgoing amplitude and a Wronskian diagnostic. Amplitudes are per kappa E Y00* in the raw Fourier convention; they are not detector strain. Production settings are `outer=400, order=8, tail_order=14, left=-35, rtol=2e-11`.
- `validate_w4.py`: five named frequency/mass points with boundary/order/precision refinements, two methods, direct sourced-ODE residual, local series residuals, and 241 spectral samples. It saves CSV data and a JSON ledger after the checks pass. The waveform is a Tukey-tapered band over 0.52≤omega≤1 at mu=0.5 and r=200, with a frequency-quadrature refinement check. The signed flux diagnostic imports and convention-converts AGP Eq.(58); its stress-energy derivation is not independently validated here.
- `plot_w4.py`: creates PNG/PDF figures from the saved CSV files. It does not rerun the solver.

The analytic source tail matters: a sharp outer cutoff can dominate a weak outgoing signal even when two solvers agree. Both source-tail order and boundary radius must be refined. The code's `tail_order=0` is a diagnostic for this error, not the production boundary condition. The saved band waveform excludes subthreshold and omitted continuum frequencies; it is not a complete published-waveform reproduction. Threshold treatment, full frequency coverage, absolute literature normalization, a physical metric observable and a convergence-controlled suppression scan remain W4 work. No TeX changes are made by these scripts.
