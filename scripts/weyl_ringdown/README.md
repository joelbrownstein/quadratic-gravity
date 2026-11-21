# W1 axial Schwarzschild benchmark

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
