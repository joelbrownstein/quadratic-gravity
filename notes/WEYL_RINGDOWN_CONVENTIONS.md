# Gate W0 — Weyl-ringdown conventions and retained sectors

21 November 2026. Codex analytical working derivation for review. This completes W0 at the classical, linearized benchmark level and supplies the conventions for W1. It is excluded from TeX. No mode spectrum, sourced waveform, stability result, quantum prescription or observational bound is newly computed here.

## 1. Scope and geometric conventions

Follow [the roadmap](WEYL_RINGDOWN_ROADMAP.md), [scalaron conventions](SCALARON_CONVENTIONS.md) and [packet B](SCALARON_MODE_DECOMPOSITION.md). Work in four dimensions, with hbar=c=1, signature (−+++), vanishing cosmological constant, and

\[
R^\rho{}_{\sigma\mu\nu}=\partial_\mu\Gamma^\rho_{\nu\sigma}
-\partial_\nu\Gamma^\rho_{\mu\sigma}
+\Gamma^\rho_{\mu\lambda}\Gamma^\lambda_{\nu\sigma}
-\Gamma^\rho_{\nu\lambda}\Gamma^\lambda_{\mu\sigma},
\quad R_{\sigma\nu}=R^\rho{}_{\sigma\rho\nu},\quad
\Box=g^{\mu\nu}\nabla_\mu\nabla_\nu.
\tag{W0.1}
\]

Start from

\[
S=\int d^4x\sqrt{-g}\left[\frac{M_{\rm Pl}^2}{2}R
+\alpha_LR^2+\beta_LC^2\right]+S_m[g,\psi],
\quad M_{\rm Pl}^2=(8\pi G)^{-1}>0.
\tag{W0.2}
\]

The W0/W1 model explicitly chooses **alpha_L=0**, beta_L<0, and finite nonzero M_Pl and spin-two mass. Matter couples minimally to the original metric g; define T_mu_nu=−2/sqrt(−g) delta S_m/delta g^mu_nu. The Schwarzschild background is vacuum. A later perturbative source must obey background stress-tensor conservation.

This is a classical Einstein–Weyl benchmark, not a claim that alpha_L=0 is preserved by quantum renormalization. We retain the fourth-order theory's extra massive mode rather than order-reducing it away as an EFT correction. Euler and integrated surface terms are discarded only in this local, fixed-topology bulk derivation; a finite-boundary variational principle is not established here.

## 2. Notation map and dimensions

Factor the Einstein normalization to obtain

\[
S_{\rm EW}=\frac{M_{\rm Pl}^2}{2}\int\sqrt{-g}
\left[R-a_WC^2\right]+S_m,
\quad
a_W=-\frac{2\beta_L}{M_{\rm Pl}^2}
=\frac{1}{f_2^2M_{\rm Pl}^2},\qquad
\beta_L=-\frac{1}{2f_2^2}.
\tag{W0.3}
\]

Choose f2>0 as a parametrization of this branch. The auxiliary-field normalization used below gives

\[
a_W=\frac{1}{2m_2^2},\qquad
m_2^2=-\frac{M_{\rm Pl}^2}{4\beta_L}
=\frac{f_2^2M_{\rm Pl}^2}{2}.
\tag{W0.4}
\]

| Our quantity | Meaning / source notation | Mass dimension |
|---|---|---:|
| alpha_L | Coefficient of R² in the unfactored action | 0 |
| beta_L | Coefficient of C² in the unfactored action | 0 |
| a_W | Papers' positive dimensionful Weyl coefficient alpha | −2 |
| b_R=2 alpha_L/M_Pl² | Papers' dimensionful R² coefficient beta, if restored | −2 |
| m2_pole = m_2 | Papers' spin-two mass mu | 1 |
| q_mu_nu | Papers' auxiliary tensor f_mu_nu; not the coupling f2 | 0 |
| M_BH | Physical remnant mass | 1 |
| r_g=G M_BH | Geometric mass; papers' M when G=c=1 | −1 |
| mu_RG | Reserved for a renormalization scale, unused at tree level | 1 |

The scalaron relation, if restored, is m0²=M_Pl²/(12 alpha_L)=1/(6 b_R). Do not interchange a_W and alpha_L. The unreduced Planck mass is m_Pl=G^(−1/2)=sqrt(8 pi) M_Pl. Frequencies and masses in this note are source-frame quantities.

## 3. First mass check: fourth-order metric action

Around Minkowski, set g=eta+h with no Planck factor in h. For a transverse-traceless massive-sector test field, linear Ricci curvature is −Box h_mu_nu/2 and linear R vanishes. Using

\[
C^2=E_4+2R_{\mu\nu}R^{\mu\nu}-\frac23R^2,
\]

the TT quadratic action, after integrations by parts, is

\[
S^{(2)}_{TT}=\int\left[
\frac{M_{\rm Pl}^2}{8}h^{\mu\nu}\Box h_{\mu\nu}
+\frac{\beta_L}{2}(\Box h_{\mu\nu})^2\right].
\tag{W0.5}
\]

Its operator is

\[
\frac{M_{\rm Pl}^2}{4}\Box+\beta_L\Box^2
=\beta_L\Box(\Box-m_2^2).
\tag{W0.6}
\]

For exp(−i omega t+i k.x), Box has eigenvalue s=omega²−|k|². The roots are s=0 and s=m2², with the latter giving the ordinary massive dispersion relation. The inverse TT factor is

\[
\left(\frac{M_{\rm Pl}^2s}{4}+\beta_Ls^2\right)^{-1}
=\frac4{M_{\rm Pl}^2}\left(\frac1s-\frac1{s-m_2^2}\right).
\tag{W0.7}
\]

The opposite massive residue is retained, not removed by choosing positive m2². This algebra fixes no quantum pole contour. At the massless pole physical counting requires gauge constraints, not a singular off-shell TT projector; packet B explains that distinction.

## 4. Independent mass check: auxiliary tensor and diagonalization

Use the bulk action

\[
S_{\rm aux}=\frac{M_{\rm Pl}^2}{2}\int\sqrt{-g}
\left[R+2q^{\mu\nu}G_{\mu\nu}
+m_2^2(q_{\mu\nu}q^{\mu\nu}-q^2)\right]+S_m[g,\psi],
\quad q=g^{\mu\nu}q_{\mu\nu}.
\tag{W0.8}
\]

Vary q to obtain G_mu_nu+m2²(q_mu_nu−q g_mu_nu)=0. Its trace gives q=−R/(3m2²), and hence

\[
q_{\mu\nu}=-\frac1{m_2^2}\left(R_{\mu\nu}-\frac16Rg_{\mu\nu}\right).
\tag{W0.9}
\]

Writing A=R_mu_nu R^mu_nu and B=R², direct contraction gives

\[
q^{\mu\nu}G_{\mu\nu}=-\frac{A-B/3}{m_2^2},\qquad
q_{\mu\nu}q^{\mu\nu}-q^2=\frac{A-B/3}{m_2^4}.
\]

Elimination therefore adds −(A−B/3)/m2² = −C²/(2m2²)+E4/(2m2²), proving (W0.3–4) modulo the stated Euler term.

Let G^(1) denote the self-adjoint linear Einstein operator on a Ricci-flat background. The vacuum quadratic action has the form

\[
S^{(2)}=M_{\rm Pl}^2\int\sqrt{-\bar g}
\left[-\frac14h^{\mu\nu}G^{(1)}_{\mu\nu}(h)
+q^{\mu\nu}G^{(1)}_{\mu\nu}(h)
+\frac{m_2^2}{2}(q_{\mu\nu}q^{\mu\nu}-q^2)\right].
\tag{W0.10}
\]

Here q now denotes its first-order perturbation about zero. Define H=h−2q. Completing the square leaves

\[
S^{(2)}=M_{\rm Pl}^2\int\sqrt{-\bar g}
\left[-\frac14H^{\mu\nu}G^{(1)}_{\mu\nu}(H)
+q^{\mu\nu}G^{(1)}_{\mu\nu}(q)
+\frac{m_2^2}{2}(q_{\mu\nu}q^{\mu\nu}-q^2)\right].
\tag{W0.11}
\]

Thus G^(1)(H)=0 and G^(1)(q)+(m2²/2)(q−bar g q)=0 in vacuum. In flat TT gauge the second equation is (Box−m2²)q_mu_nu=0. Both derivations agree in sign and normalization. The kinetic sign in (W0.11) also agrees with the opposite residue in (W0.7).

The physical metric perturbation is **h=H+2q**. For a conserved external perturbative stress tensor its coupling is (1/2) integral sqrt(−bar g) h_mu_nu T^mu_nu. The massive field cannot simply be excluded from detector response because the free action can be diagonalized. This is a linear field redefinition, not nonlinear decoupling.

## 5. Constraints on Schwarzschild and retained modes

Use bar q=0 and the Ricci-flat Schwarzschild background. The auxiliary equation and the Bianchi identity imply the exact bulk constraint

\[
\nabla^\mu q_{\mu\nu}=\nabla_\nu q.
\tag{W0.12}
\]

In vacuum Einstein–Weyl theory the metric trace equation gives R=0; consequently q=0. At linear order this yields four divergence constraints plus one trace constraint:

\[
\bar\nabla^\mu q_{\mu\nu}=0,\qquad q=0.
\tag{W0.13}
\]

Because bar q=0, its perturbation is invariant under infinitesimal coordinate changes; it has no independent gauge symmetry to remove additional components. Its ten symmetric components minus five constraints leave five massive polarizations. H retains linear diffeomorphism invariance and the usual Einstein constraints, leaving two massless helicities. A possible computational gauge is bar nabla^mu(H_mu_nu−bar g_mu_nu H/2)=0 with residual freedom handled separately; W1 may instead use a harmonic-adapted gauge.

Define explicitly, to avoid alternative Lichnerowicz signs,

\[
(\mathcal D q)_{\mu\nu}=\bar\Box q_{\mu\nu}
+2\bar R_{\mu\rho\nu\sigma}q^{\rho\sigma}.
\]

For vacuum TT q on a Ricci-flat background, G^(1)(q)=−mathcal D q/2, so

\[
(\mathcal D-m_2^2)q_{\mu\nu}=0.
\tag{W0.14}
\]

This is a covariant input for W1, not the yet-to-be-derived radial master system. The curvature term survives on Schwarzschild despite R_mu_nu=0. At generic l>=2 the massive multipoles divide into two axial and three polar degrees of freedom; exceptional l=0,1 must be reduced separately. The helicity-zero component of this massive spin-two field is not the independent spin-zero scalaron. Vacuum TT constraints must not be imposed unaltered at a matter source.

## 6. Why the scalaron is omitted, and when that fails

Derive the parent trace directly from (W0.2). Weyl-squared has zero bulk metric trace in four dimensions. Variation of R² contributes 6 alpha_L Box R to the trace before multiplying the full metric equation by two. With our matter normalization,

\[
-M_{\rm Pl}^2R+12\alpha_L\Box R=T.
\tag{W0.15}
\]

Therefore for alpha_L>0 and a Ricci-flat vacuum background,

\[
(\bar\Box-m_0^2)\delta R=\frac{\delta T}{12\alpha_L},
\qquad m_0^2=\frac{M_{\rm Pl}^2}{12\alpha_L}.
\tag{W0.16}
\]

Background R=0 does not remove this homogeneous scalar wave or its trace source. There are three distinct justifications one might use:

1. **Our benchmark:** set alpha_L=0 as an action choice. Then R=−T/M_Pl² is a constraint, not a propagating scalaron equation. Vacuum gives R=0, but matter need not.
2. **Linear axial restriction on a spherical background:** delta R is parity-even, so it is absent in the axial sector. This does not eliminate the polar scalaron of the parent theory.
3. **Selected scalar-free initial data with finite alpha_L:** delta R=0 and its initial time derivative remain zero at linear order only with compatible boundary conditions and delta T=0. A finite-mass point particle has a nonzero trace and generally fails this condition.

In the alpha_L=0 sourced model, (W0.9) and (W0.15) give

\[
q=\frac{T}{3m_2^2M_{\rm Pl}^2},\qquad
\nabla^\mu q_{\mu\nu}=\nabla_\nu q.
\tag{W0.17}
\]

The trace is algebraically fixed; it is not a new scalar degree of freedom. Away from source support the vacuum constraints return. As a check, varying (W0.11) plus the source coupling gives G^(1)(H)=T/M_Pl² and G^(1)(q)+(m2²/2)(q−bar g q)=−T/(2M_Pl²); their trace and divergence reproduce (W0.17).

For a future ultrarelativistic radial-infall test, specify m_particle→0 and gamma→infinity with E=gamma m_particle fixed. This trace-free null-particle limit is not the same as sending gamma→infinity at fixed particle mass and unbounded perturbation energy. Finite-energy timelike infall in generic QQG requires explicit scalar sourcing. Neither a field redefinition to an Einstein frame nor a vanishing background scalar licenses ignoring transformed matter couplings.

## 7. Coordinate, frequency and boundary conventions for W1

\[
ds^2=-Fdt^2+F^{-1}dr^2+r^2d\Omega^2,\quad
F=1-\frac{2r_g}{r},\quad r_g=GM_{\rm BH},\quad r_h=2r_g,
\]
\[
r_*=r+2r_g\ln\left(\frac{r}{2r_g}-1\right),\quad
x_2=r_gm_2,\quad p=r_hm_2=2x_2,\quad\hat\omega=r_g\omega.
\tag{W0.18}
\]

Use angular harmonics normalized to unit integral of |Y_lm|². Define Fourier reconstruction with exp(−i omega t) and measure d omega/(2 pi). Damped modes have Im(omega)<0. A source-free QNM is ingoing at the future horizon and outgoing at infinity; the retarded forced problem has no incoming radiation and causal source data. Near the horizon the oscillatory ingoing factor is exp(−i omega r_*). The massless outgoing factor is exp(+i omega r_*).

For a massive channel k²=omega²−m2². Fix the retarded branch by analytic continuation from Im(omega)>0, choosing k~omega there at large |omega| and Im(k)>0. For real positive omega>m2 take positive k; between zero and threshold take k=+i sqrt(m2²−omega²). Continue to resonance sheets consistently when searching for QNMs. The leading outgoing factor exp(+i k r_*) is not a complete asymptotic series: gravitational long-range terms can add powers/logarithmic phases, which W1 must derive. Quasibound modes with decaying spatial boundary conditions are a distinct spectral problem. Do not replace k by its absolute value or choose a fresh square-root sign at each root.

Retarded continuation specifies a classical response, not a Feynman prescription for a quantum ghost. If a background has growing modes, the causal inversion contour must initially lie above their poles; an unrestricted stable Fourier response cannot be presumed. Stability ranges are deferred to W2, not assumed from the positive flat-space m2².

At cosmological redshift z, omega_obs=omega_src/(1+z) and M_BH,det=(1+z)M_BH,src, so r_g,det omega_obs=r_g,src omega_src. The local pole mass in x2 remains a source-frame physical parameter; detector-frame mass must not be inserted into x2 without the corresponding redshift bookkeeping.

## 8. Acceptance checks and remaining boundaries

| W0 check | Result |
|---|---|
| Original-action massive pole | (W0.5–7) gives −M_Pl²/(4 beta_L) |
| Auxiliary elimination | (W0.9) recovers the curvature action modulo E4 |
| Independent auxiliary pole | (W0.10–14) gives the same mass and opposite kinetic sign |
| Units and normalization | a_W has length², x2 is dimensionless, p=2x2; Planck conventions separated |
| Retained modes | Two massless plus five massive; scalaron excluded by alpha_L=0 |
| Matter consistency | Trace and divergence constraints differ from vacuum as in (W0.17) |
| Boundary conventions | Time sign, retarded branch and spectral problem specified; radial series deferred |

These are analytical consistency checks, not peer review or a numerical solution. Exact-rational arithmetic spot checks also passed for the pole factorization, partial fractions, and auxiliary trace/elimination on representative diagonal Ricci tensors; they supplement rather than prove the general identities derived above. A symbolic-software check was not available because SymPy is not installed. W0 is complete for the stated classical model. W1 may derive and solve the axial Schwarzschild system. Rotation, the polar source normalization, numerical boundary expansions, stability domains, quantum contours and loop matching remain later work. PS ensembles cannot supply the omitted Weyl dynamics or determine f2 from this note.

## 9. Literature audit and provenance

- [Antoniou, Gualtieri and Pani (2025), arXiv:2412.15037v2](https://arxiv.org/html/2412.15037v2), section II and II.1, equations (4)–(9): auxiliary normalization and constraints. Their dimensionful alpha maps to a_W and mu to m2_pole. The derivations above explicitly check elimination and kinetic signs in our conventions. Their axial scope does not justify deleting the scalaron from a general polar sourced problem.
- [Antoniou, Gualtieri and Pani (2026), arXiv:2604.01631v2](https://arxiv.org/html/2604.01631v2), sections II.1, II.2 and II.4: comparison of field content, p normalization and infall setup. We define an explicit Einstein–Weyl benchmark, rather than importing an unrestricted scalar-decoupling claim for minimally coupled Jordan-frame matter.
- Local packet B supplies the independent metric-Hessian comparison. The trace derivation, source constraints and dimensional checks here are our working analysis. No quantum equivalence or generic QQG truncation theorem is inferred from either paper.
