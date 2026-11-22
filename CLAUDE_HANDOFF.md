#    Handoff: give the scalar calculation its first numerical evidence

Prepared 20 November 2026. Joel is conserving his Codex weekly allowance and will continue in Claude. The working manuscript and this handoff are saved in `/software/github/joelbrownstein/quadratic-gravity/main` (resolved `/private/software/github/joelbrownstein/quadratic-gravity/main`). Changes are local and uncommitted. Preserve the existing working tree. No Overleaf or GitHub sync has been performed by Codex.

## Current scientific heartbeat

**Question:** Does the specified four-dimensional perfect-square scalar theory show controlled continuum scaling and renormalized running compatible with its perturbative benchmarks?

The manuscript defines a candidate periodic regulator, flat measure with the constant mode removed, an elementary finite-lattice integrability bound, the exact action derivative, free-field covariance, and a proposed finite-volume coupling and step-scaling protocol. The deterministic stencil checks pass. This is a research design with reproducible algebraic checks, not a completed nonperturbative result. No HMC sampler, interacting ensembles, continuum extrapolation, cosmological likelihood fit, or singularity-resolution proof has been produced.

## Files to read first

1. `LITERATURE_MAP.md`, `BOND_SEARCH_LOG.md`, `HARTLE_COSMOLOGY_READING_LIST.md`: evidence levels, actual Bond queries and primary-text checks. Bond was used through the existing Claude chat; it supplied metadata/abstracts and citation/version discovery, not full texts.
2. `MANUSCRIPT_DESIGN.md`: narrative and independent completion criteria.
3. `sections/model.tex`, `sections/lattice_action.tex`, `sections/observables.tex`, `sections/numerics.tex`, `sections/predictions.tex`.
4. `scripts/validate_scalar_regulator.py`: independent standard-library stencil checks. Run it before changing implementation.

The manuscript now has 68 bibliography entries. The 21 November AI pass adds `AI_RESEARCH_MAP.md` and expands the final section into a calibrated geometry–growth frontier experiment, optional corrected learned PS sampling, and matched-cost comparisons. No sampler or cosmological inference was run. Hartle's supervision and support are in the author's personal perspective. Cosmology now includes no-boundary histories, arrows of time and the nonrotating Schwarzschild/Kantowski–Sachs correspondence. The Hawking-visit anecdote is deliberately omitted. Hartle's classicality condition is |∇I_R| << |∇S| for exp[(-I_R+iS)/hbar]; Claude's discovery synthesis reversed it in one proposed numerical question. Use the primary paper, not that synthesis. A positive Euclidean weight is not a physical-history probability prescription.

## Next bounded implementation task

Implement a transparent CPU reference HMC sampler for the **existing** scalar regulator. Work at small fixed N first; GPU optimization follows validation. Keep the action, measure and field normalization unchanged unless a documented scientific reason requires a new regulator. Store scripts, exact parameters, seeds, environment, source hashes and raw scalar observable time series.

- Enforce zero field mean and project momenta onto the same subspace. Validate the action gradient against the independent reference stencil and finite differences; check reversibility and energy-error scaling for the integrator before interpreting samples.
- At lambda=0 compare HMC with independent Gaussian Fourier draws. The nonzero-mode covariance is 1/(p-hat squared)^2 with the manuscript's unitary FFT normalization, and E[S]=(N^4-1)/2. Assess agreement using correlation-aware errors, not selected visually favorable modes.
- Validate the accept/reject step and repeat at smaller step size. Use independent starts/seeds; retain all pilot attempts and failures. Estimate autocorrelation before reporting standard errors or effective samples.
- For a genuinely interacting fixed-N check, the scalar action decomposes into homogeneous pieces S=S2+S3+S4, where S2=1/2 sum(Delta s)^2, S3=lambda sum(Delta s)q, S4=lambda^2/2 sum q^2. Integration by parts on the V-1 dimensional zero-mean subspace gives E[s dot grad S]=V-1, hence E[2S2+3S3+4S4]=V-1. The finite-N tail bound justifies the boundary term. Independently verify this derivation and use it as a sampling check, not as evidence of a continuum Ward identity. The field-sign map also relates lambda to minus lambda and supplies a consistency check for even observables.
- Freeze a small pilot protocol before sampling. Report outcomes at fixed bare coupling and fixed regulator. Do not infer a beta function from bare-coupling comparisons. The lattice coupling, scheme matching and relevant-operator tuning are still required before a continuum running claim.

A temporary scientific Python environment was prepared at `/private/tmp/qqg-pilot-env` with Python 3.11.16, NumPy 2.4.6, SciPy 1.17.1 and Matplotlib 3.11.2. This provisioning ran no physics simulations and changed no system-wide packages. Create a portable dependency record for retained calculations; do not depend on a temporary path in published instructions.

## Manuscript and publication discipline

Do not replace the abstract's research-design status with results claims until real validated measurements exist. Label pilot findings as pilot findings. The scalar model does not establish the full QQG quantum limit, a Lorentzian probability interpretation, a black-hole interior, or evolving dark energy. The broad applications remain separately scoped future research.

Keep Joel as the responsible listed author; AI assistance is disclosed, not listed as coauthorship. Do not send correspondence, publish on Valency Hub, mint a DOI, commit/push or sync Overleaf merely to complete this handoff. Joel has said he will handle the sync. Correspondence remains on hold at Joel’s request pending the Redtail proposal pass.

## Completion criterion for the next milestone

A reproducible reference sampler with free-field validation and at least one statistically assessed interacting fixed-N pilot, including failed checks and cost/autocorrelation measurements. Only then update the paper with a small numerical evidence section and decide whether a larger allocation is justified.
