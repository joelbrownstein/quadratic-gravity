# From quadratic gravity to testable black-hole signals

An introduction for graduate students in computer science and astrophysics. This overview accompanies the working notes and calculations; it is not included in the manuscript's TeX build. Status: 21 November 2026.

This workbook explores a practical question: **if gravity contains additional degrees of freedom beyond general relativity, how do we turn them into a prediction that an observation could actually test?**

The broader project connects quadratic gravity, quantum field theory, black holes, and cosmology. The calculations completed so far concentrate on a controlled example: the response of a Schwarzschild black hole to an infalling particle.

Einstein's theory starts with an action proportional to the Ricci scalar $R$. Quadratic gravity adds terms containing two powers of curvature:

$$
S=\int d^4x\,\sqrt{-g}
\left[
\frac{M_{\rm Pl}^{2}}{2}R
+\alpha R^{2}
+\beta C_{\mu\nu\rho\sigma}C^{\mu\nu\rho\sigma}
\right].
$$

Here $C_{\mu\nu\rho\sigma}$ is the Weyl tensor. These additions change both the short-distance behavior of the theory and its propagating fields. Alongside the ordinary massless graviton, the linearized theory generally contains a massive scalar—the **scalaron**—and a massive spin-two field.

The quantum motivation is improved ultraviolet behavior. The difficulty is that the extra spin-two sector has a problematic kinetic sign. The workbook therefore distinguishes classical calculations from unresolved questions about the quantum theory. It does not assume that calculating a classical waveform settles quantum consistency.

A central strategy is to work through **controlled sectors**. A scalar toy theory can illuminate higher-derivative dynamics and perturbation theory. But its relationship to the physical scalaron must be demonstrated, and it cannot automatically reproduce the tensor modes responsible for black-hole signals. Even the “scalar-like” polarization of a massive spin-two field is distinct from a genuine spin-zero scalaron.

The black-hole calculations make that distinction concrete. They currently retain the Einstein–Weyl sector, isolating the additional massive spin-two field. In the workbook's conventions, $\beta=-1/(2f_2^2)$, and its mass is related to the Weyl coupling by

$$
m_2^2=\frac{f_2^2M_{\rm Pl}^2}{2}.
$$

An especially interesting feature is that a black hole can retain the Schwarzschild background geometry while having a richer perturbative response. Its equilibrium appearance can resemble GR even though its available oscillations differ.

For a CS student, the computational problem resembles constructing a constrained dynamical system and its transfer function:

- **Derive the state space.** Identify independent fields and enforce the constraints rather than evolving redundant variables indiscriminately.
- **Find the homogeneous response.** Determine oscillations, damping, and possible instabilities.
- **Specify the input.** Derive a conserved source with consistent normalization.
- **Compute the driven response.** Determine which modes are actually excited and with what amplitudes.
- **Validate the output.** Compare independent algorithms and vary numerical boundaries, resolution, and approximations.

The distinction between an available mode and an excited mode is crucial. Finding an extra frequency does not establish that a telescope or interferometer could detect it. The source might couple weakly to it, and propagation might further suppress or reshape the signal.

Our first driven example is **null radial infall**, focusing on the monopole. GR has no freely propagating monopolar gravitational wave; the additional massive spin-two sector provides a clean place to investigate different behavior. This is a test-particle calculation, not a binary-merger simulation.

The numerical work has already exposed two consequential issues. First, cutting the particle's trajectory off at a finite outer radius can overwhelm the weak outgoing signal. We therefore derive and include the exterior source tail. Second, the usual outgoing expansion fails at the mass threshold, where

$$
k=\sqrt{\omega^2-m_2^2}\rightarrow0.
$$

Below that threshold, the homogeneous field decays at infinity but can still contribute at a finite observation radius. A waveform calculation must retain that response.

The current workbook includes threshold-aware boundaries, scaled Green-function integration, an independent boundary-value solver, and causal-contour frequency inversion. These produce a convergence-tested auxiliary-field waveform over a specified post-passage time window. A subsequent [sharp-front calculation](WEYL_RINGDOWN_NULL_FRONT.md) derives a continuous angular auxiliary field with a logarithmically divergent slope at particle passage. It treats the dominant front analytically and tests the smaller curved-space remainder numerically in a separate local time window. The gap between these windows remains to be validated, and the result is not yet detector strain.

The motivating destination is a defensible chain:

$$
\text{action and couplings}
\;\longrightarrow\;
\text{degrees of freedom}
\;\longrightarrow\;
\text{excitation and propagation}
\;\longrightarrow\;
\text{observable prediction}.
$$

Black holes provide the present laboratory. Cosmological tests and AI-assisted exploration are later destinations. The workbook's value is that it makes the intermediate reasoning executable and reviewable—so an intriguing theory can gradually become a falsifiable calculation.

## Where to go next

- [Research roadmap](WEYL_RINGDOWN_ROADMAP.md): the staged questions, acceptance criteria, and remaining work.
- [Conventions and retained sectors](WEYL_RINGDOWN_CONVENTIONS.md): the action, coupling definitions, and classical scope.
- [Scalaron–PS matching roadmap](SCALARON_PS_MATCHING_ROADMAP.md): what must be established before identifying a scalar toy model with a sector of the parent theory.
- [Initial null-infall calculation](WEYL_RINGDOWN_NULL_INFALL.md): source derivation, Green function, and the importance of the exterior source tail.
- [Threshold-inclusive response](WEYL_RINGDOWN_THRESHOLD_RESPONSE.md): the latest finite-window waveform, numerical checks, and limitations.
- [Code and reproduction commands](../scripts/weyl_ringdown/README.md): how to inspect and rerun the calculations.

For the broader manuscript and its primary literature, return to the [repository README](../README.md) and [literature map](../LITERATURE_MAP.md). The linked derivation notes contain the detailed references and validation records behind this introduction.
