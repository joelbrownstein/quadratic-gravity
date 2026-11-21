# AI crescendo: literature and experiment design

21 November 2026. Discovery used Valency Bond through the existing [Claude chat](https://claude.ai/chat/89ceba0b-59c7-4646-9918-f601533b2e8e), followed by direct primary-source checks. Bond is not connected natively to this Codex session. Its completed memo disclosed abstract/metadata access, no full texts. This targeted search is not exhaustive.

## Selected evidence and consequences

“Text” means relevant primary sections inspected, not every derivation reproduced. “Abstract” means the primary abstract and record checked. Eleven new references were added to the manuscript, alongside its existing CosmoPower-JAX and simulation-based-calibration anchors.

| Primary source | Evidence | Design consequence and limit |
|---|---|---|
| [Fang, Hu & Lewis 2008](https://arxiv.org/abs/0808.3125), phantom-divide crossing | Abstract | Specify the dark-energy perturbation prescription. PPF is a phenomenological control, not PS matching. |
| [Albergo, Kanwar & Shanahan 2019](https://arxiv.org/html/1904.12072v3), flow-based lattice MCMC | Text, sec. II.1 | Metropolis correction with known density and full support on the target space. Demonstration in 2D phi-four does not establish transfer to PS. |
| [Del Debbio, Marsh Rossney & Wilson 2021](https://arxiv.org/abs/2105.12481), scalability | Abstract | Training costs can rise steeply near the continuum; measure amortized cost against HMC. |
| [Jense et al. 2024](https://arxiv.org/html/2405.07903v1), CosmoPower framework | Text, parameter ranges and sec. 4.3 | Package solver settings and training domains; use native-solver fallback. No generic QQG surrogate follows. |
| [Sáez-Casares, Rasera & Li 2023](https://arxiv.org/html/2303.08899v2), e-MANTIS | Text, secs. 3.2.1–3.2.2 | Restricted Hu–Sawicki boost emulator requires a baseline spectrum with its error. Abstract-rounded domain bounds differ from the precise text; manuscript avoids unnecessary numerical ranges. |
| [Linhart, Gramfort & Rodrigues 2023](https://arxiv.org/html/2306.03580v2), L-C2ST | Text, examples and sec. 6 | Local diagnostics complement global calibration. Finite classifier accuracy and power limit guarantees. |
| [Piras et al. 2024](https://arxiv.org/html/2405.12965v2), accelerated inference/model comparison | Text, secs. 3.4 and 4 | Validate evidence separately. Their simulated-shear benchmark and hardware timing are not a Redtail forecast. |
| [Lee 2026](https://arxiv.org/html/2511.16703v2), dark-energy null tests | Text, sec. V.5 and conclusions | Distinguish posterior displacement from model exclusion. Single-probe recovery does not calibrate a joint likelihood; retain convergence-selection outcomes. |
| [Rogers et al. 2019](https://arxiv.org/html/1812.04631v2), Bayesian emulator optimization | Text, acquisition and initialization | Adaptive simulation design has a cosmological precedent. Single-model posterior refinement does not establish model-discrimination savings. |
| [Ho, Bird & Shelton](https://arxiv.org/abs/2105.01081), multi-fidelity matter-power emulation | Abstract, 2021 preprint | Combine resolutions where validated. Their relative comparison does not establish converged absolute accuracy. |
| [Ye et al. 2025](https://arxiv.org/html/2510.24591v2), ReplicationBench | Text, tables 2–3, secs. 5.2 and 6 | Require executable outputs and reruns. Bond's “under 20%” summary obscures weighting: table 2 has a 0.22 unweighted score and 0.19 difficulty-weighted score for one configuration. Manuscript uses qualitative failure evidence, not an ambiguous headline score. Replication is not open-ended discovery. |

## Proposed contribution

A geometry-only versus geometry-plus-growth recovery-power map for LambdaCDM, specified w0–wa dark energy, and a viable f(R) control. Freeze decisions and priors, calibrate the joint false-alternative rate on development nulls, and confirm on hidden nulls/injections. Compare adaptive and space-filling designs at matched calibration and frontier accuracy, counting training, acquisition, failures and solver checks. Report inconclusive cases and Monte Carlo errors. No numerical results exist yet. The PS sampler comparison is separately costed and does not predict w(z).

## Discovery queries displayed in Claude's completed memo

All searches used arXiv. The memo reports twelve searches/exports but displays seven grouped operations below; this is a faithful grouped record, not an independently verified twelve-call audit. Exact strings are retained. Other retrieved leads were not automatically promoted to manuscript evidence.

| Operation / filter | Exact query or IDs | Representative returned IDs |
|---|---|---|
| export_papers_csv | 2305.06347, 2405.07903, 2405.12965 | Three abstracts |
| semantic_search_papers, start 2023 | emulator for modified gravity f(R) or dynamical dark energy nonlinear matter power spectrum and growth, training range limitations, accuracy versus N-body or Boltzmann solver | 2303.08899, 2405.05840, 2310.05786, 2607.29683, 2604.19613, 2511.13826 |
| semantic_search_papers, start 2022 | cosmology model misspecification in simulation-based or emulator-based inference: calibration of Bayesian evidence and model comparison, coverage tests with injected mock data, false detection of dark energy evolution | 2511.16703, 2605.00980, 2606.16248, 2601.11150, 2410.10616, 2507.22990 |
| semantic_search_papers, start 2019 | active learning or multi-fidelity design of cosmological simulation training sets for emulators; acquisition function for model discrimination between competing cosmological models rather than single-model posterior | 2305.19267, 2307.01138, 2105.01081, 2507.00514, 1905.07410 |
| semantic_search_papers, sort citations | normalizing flows for lattice field theory sampling with exact Metropolis-Hastings or reweighting correction; scaling of effective sample size with volume and critical slowing down; mode collapse; comparison with HMC cost | 2105.12481, 2211.07541, 2107.00734, 2302.14082, 2412.00200 |
| semantic_search_papers, start 2024, sort citations | benchmark evaluating LLM agents reproducing scientific papers or performing autonomous research tasks; computational reproducibility; failure rates and limitations of AI scientist systems | 2409.11363, 2510.24591, 2605.16616, 2608.05179, 2608.26753, 2610.39744 |
| semantic_search_papers, sort citations | parameterized post-Friedmann treatment of dark energy perturbations allowing phantom divide crossing w=-1 in Boltzmann codes | 0808.3125, astro-ph/0410680, astro-ph/0609040, 2508.01378 |

Claude's suggested numerical thresholds, exact-solver percentages and factor-two saving were proposals, not literature results; they were not adopted as established standards. Additional 2026 discovery leads remain unverified and uncited in this pass.
