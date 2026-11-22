"""check_scalaron_benchmark.py -- symbolic checks for packet A (labels match SCALARON_BENCHMARK_DERIVATION_CLAUDE.md) (scalaron benchmark), repo conventions:
signature (-,+,+,+), R^r_{s m n} = d_m G^r_{n s} - d_n G^r_{m s} + G^r_{m l} G^l_{n s} - G^r_{n l} G^l_{m s},
R_{s n} = R^r_{s r n}.  Algebra only; no simulation."""
import sympy as sp

def ricci_scalar(g, X):
    n = len(X); gi = g.inv()
    G = [[[sp.simplify(sum(gi[r, l]*(sp.diff(g[l, m], X[nn]) + sp.diff(g[l, nn], X[m]) - sp.diff(g[m, nn], X[l]))
            for l in range(n))/2) for nn in range(n)] for m in range(n)] for r in range(n)]
    def Riem(r, s, m, nn):
        e = sp.diff(G[r][nn][s], X[m]) - sp.diff(G[r][m][s], X[nn])
        e += sum(G[r][m][l]*G[l][nn][s] - G[r][nn][l]*G[l][m][s] for l in range(n))
        return e
    Ric = sp.Matrix(n, n, lambda s, nn: sum(Riem(r, s, r, nn) for r in range(n)))
    return sp.simplify(sum(gi[s, nn]*Ric[s, nn] for s in range(n) for nn in range(n)))

def box(f, g, X):
    gi = g.inv(); sg = sp.sqrt(-g.det()) if g.det().could_extract_minus_sign() else sp.sqrt(g.det())
    return sp.simplify(sum(sp.diff(sg*gi[m, nn]*sp.diff(f, X[nn]), X[m]) for m in range(4) for nn in range(4))/sg)

def grad2(f, g, X):
    gi = g.inv(); return sp.simplify(sum(gi[m, nn]*sp.diff(f, X[m])*sp.diff(f, X[nn]) for m in range(4) for nn in range(4)))

t, x, y, z = X = sp.symbols('t x y z', real=True)
out = {}

# 0. Sign convention: unit 2-sphere x R^2 Euclidean-like check -> R = +2 (positive for spheres)
th, ph = sp.symbols('theta phi', positive=True)
gS = sp.diag(1, sp.sin(th)**2, 1, 1); R0 = ricci_scalar(gS, [th, ph, y, z])
out['0 sphere_R_is_+2'] = sp.simplify(R0 - 2) == 0

# 1. Conformal map (A5): g = Omega^-2 gE, R[g] = Omega^2 [R_E + 6 BoxE u - 6 (gradE u)^2], u = ln Omega
a = sp.Function('a')(t); u = sp.Function('u')(t, x)
gE = sp.diag(-1, a**2, a**2, a**2)
g = sp.exp(-2*u)*gE
lhs = ricci_scalar(g, X)
rhs = sp.exp(2*u)*(ricci_scalar(gE, X) + 6*box(u, gE, X) - 6*grad2(u, gE, X))
out['1 A5_conformal_map'] = sp.simplify(lhs - rhs) == 0

# 2. Auxiliary identity (A4) and F
M = sp.symbols('M', positive=True); al, chi, R, J = sp.symbols('alpha chi R J', real=True)
aux = M**2*R/2 + 2*al*chi*R - al*chi**2
out['2a A4_identity'] = sp.expand(aux - (M**2*R/2 + al*R**2 - al*(chi - R)**2)) == 0
out['2b chi_eq_R'] = sp.solve(sp.diff(aux, chi), chi) == [R]
out['2c F_form'] = sp.expand(aux - ((M**2 + 4*al*chi)*R/2 - al*chi**2)) == 0

# 3. Einstein frame: kinetic normalization, potential, vacuum, mass (A6-A9)
phi = sp.symbols('phi', real=True); b = sp.sqrt(sp.Rational(2, 3))/M
F = M**2*sp.exp(b*phi); chi_phi = (F - M**2)/(4*al)
out['3a phi_inverts_F'] = sp.simplify(sp.expand_log(sp.sqrt(sp.Rational(3, 2))*M*sp.log(F/M**2), force=True) - phi) == 0
# kinetic: -3 M^2 (grad u)^2 with u = (1/2) ln(F/M^2), in units of (grad phi)^2
dudphi = sp.diff(sp.log(F/M**2)/2, phi)
out['3b kinetic_is_-1/2'] = sp.simplify(-3*M**2*dudphi**2 + sp.Rational(1, 2)) == 0
V = sp.simplify(M**4*al*chi_phi**2/F**2)
out['3c V_form'] = sp.simplify(V - M**4/(16*al)*(1 - sp.exp(-b*phi))**2) == 0
out['3d V0'] = sp.simplify(V.subs(phi, 0)) == 0
out['3e Vp0'] = sp.simplify(sp.diff(V, phi).subs(phi, 0)) == 0
m2 = sp.simplify(sp.diff(V, phi, 2).subs(phi, 0))
out['3f m0^2'] = sp.simplify(m2 - M**2/(12*al)) == 0
out['3g m0^2_in_f0'] = sp.simplify(m2.subs(al, 1/(6*sp.Symbol('f0')**2)) - sp.Symbol('f0')**2*M**2/2) == 0
out['3h dphi_dchi_nonzero_on_F>0'] = sp.simplify(sp.diff(sp.sqrt(sp.Rational(3, 2))*M*sp.log((M**2 + 4*al*chi)/M**2), chi)
                                               - 4*al*sp.sqrt(sp.Rational(3, 2))*M/(M**2 + 4*al*chi)) == 0

# 4. Trace equation algebra (A11, A13): L_R R - 2L + 3 Box L_R, with Box as symbol 'B' acting linearly
B = sp.symbols('B')  # stands for Box acting on what follows (linear order)
L = M**2*R/2 + al*R**2 + J*R
LR = sp.diff(L, R)
trace_alg = sp.expand(LR*R - 2*L)
out['4a trace_algebraic_part'] = sp.simplify(trace_alg - (-M**2*R/2 - J*R)) == 0
# linear order: 6 alpha B R - M^2 R/2 + 3 B J = 0
Rsol = sp.solve(sp.Eq(6*al*B*R - M**2*R/2 + 3*B*J, 0), R)[0]
out['4b A14_response'] = sp.simplify(Rsol - (-(1/(2*al))*B/(B - M**2/(12*al))*J)) == 0
out['4c A14_split'] = sp.simplify(Rsol - (-J/(2*al) - (M**2/(12*al))/(2*al)/(B - M**2/(12*al))*J)) == 0

# 5. Independent Einstein-frame OFF-SHELL response, contact term included.
#    Linear order: phi eq (B - m^2) phi = -3 b B J ; Einstein-frame trace with J R_E source: R_E = 6 B J / M^2;
#    Jordan curvature R = R_E + 3 b B phi (from A5 with u = b phi/2).
m0sq = M**2/(12*al)
phi_J = -3*b*B*J/(B - m0sq)
R_E = 6*B*J/M**2
R_jordan = R_E + 3*b*B*phi_J
out['5 einstein_frame_full_response_matches_A14'] = sp.simplify(R_jordan - Rsol) == 0

# 6. On-shell coupling c_R: R = c_R phi with chi = R at linear order
cR = sp.simplify(1/sp.diff(sp.sqrt(sp.Rational(3, 2))*M*sp.log((M**2 + 4*al*chi)/M**2), chi).subs(chi, 0))
out['6a c_R'] = sp.simplify(cR - M/(2*sp.sqrt(6)*al)) == 0
out['6b c_R^2=m0^2/(2alpha)'] = sp.simplify(cR**2 - m0sq/(2*al)) == 0
out['6c 3b_m0^2=c_R'] = sp.simplify(3*b*m0sq - cR) == 0

# 7. Restricted conformal mode g = exp(2 psi) eta (Lorentzian): exact R and quadratic dispersion
psi = sp.Function('psi')(t, x, y, z)
eta = sp.diag(-1, 1, 1, 1)
Rc = ricci_scalar(sp.exp(2*psi)*eta, X)
boxeta = lambda f: -sp.diff(f, t, 2) + sum(sp.diff(f, v, 2) for v in (x, y, z))
d2eta = lambda f: -sp.diff(f, t)**2 + sum(sp.diff(f, v)**2 for v in (x, y, z))
out['7a conformal_R_lorentzian'] = sp.simplify(Rc - (-6*sp.exp(-2*psi)*(boxeta(psi) + d2eta(psi)))) == 0
# quadratic Lagrangian (after integration by parts) L2 = 3 M^2 (d psi)^2 + 36 alpha (Box psi)^2 - 6 J Box psi
# EOM in momentum space with Box -> B (d psi)^2 -> -psi Box psi
psi_s = sp.symbols('psi_s')
EOM = -6*M**2*B*psi_s + 72*al*B**2*psi_s - 6*B*J   # delta/delta psi
psi_sol = sp.solve(sp.Eq(EOM, 0), psi_s)[0]
R_restricted = -6*B*psi_sol
out['7b restricted_conformal_response_matches_A14'] = sp.simplify(R_restricted - Rsol) == 0
# check the claimed quadratic sqrt(-g) R term: (M^2/2) sqrt(-g) R = 3 M^2 e^{2psi} (d psi)^2 + total derivative
eps = sp.symbols('epsilon')
lag = (M**2/2)*sp.exp(4*psi)*Rc
td = -3*M**2*sum(s*sp.diff(sp.exp(2*psi)*sp.diff(psi, v), v) for s, v in zip((-1, 1, 1, 1), X))
out['7c EH_conformal_reduction'] = sp.simplify(lag - (3*M**2*sp.exp(2*psi)*d2eta(psi) + td)) == 0

# 8. Euclidean PS identity (model.tex eq:conformal-r-squared) and 72 lambda^2 normalization
lam = sp.symbols('lambda', positive=True); sig = sp.Function('sigma')(t, x, y, z)
RE = ricci_scalar(sp.exp(2*lam*sig)*sp.eye(4), X)
lapE = lambda f: sum(sp.diff(f, v, 2) for v in X)
d2E = lambda f: sum(sp.diff(f, v)**2 for v in X)
P = lapE(sig) + lam*d2E(sig)
out['8a PS_R_identity'] = sp.simplify(RE - (-6*lam*sp.exp(-2*lam*sig)*P)) == 0
out['8b PS_sqrtg_R2'] = sp.simplify(sp.exp(4*lam*sig)*RE**2 - 36*lam**2*P**2) == 0
aE = sp.symbols('alpha_E', positive=True)
out['8c alphaE_norm'] = sp.solve(sp.Eq(36*aE*lam**2, sp.Rational(1, 2)), aE) == [1/(72*lam**2)]

# 9. Gaussian PS curvature correlator: <R R> = 36 lam^2 p^4 * (1/p^4) = 1/(2 alpha_E)
p2 = sp.symbols('p2', positive=True)
corr = 36*lam**2*p2**2*(1/p2**2)
out['9 PS_RR=36lam^2=1/(2alphaE)'] = sp.simplify(corr - 36*lam**2) == 0 and sp.simplify(corr - 1/(2/(72*lam**2))) == 0

# 10. Euclidean continuation of the benchmark response (Box -> -p^2): R/J = -(1/(2 alpha_L)) p^2/(p^2+m0^2)
RJ_E = sp.simplify((Rsol/J).subs(B, -p2))
out['10a euclid_response'] = sp.simplify(RJ_E - (-(1/(2*al))*p2/(p2 + m0sq))) == 0
out['10b UV_limit_contact'] = sp.limit(RJ_E, p2, sp.oo) == -1/(2*al)
out['10c M->0_limit'] = sp.simplify(sp.limit(RJ_E, M, 0) + 1/(2*al)) == 0

# 11. Euclidean auxiliary Gaussian: -[2 a chi R - a chi^2] = a (chi-R)^2 - a R^2  (convergent in chi iff a>0)
out['11 euclid_aux_split'] = sp.expand(-(2*al*chi*R - al*chi**2) - (al*(chi - R)**2 - al*R**2)) == 0

for k, v in out.items():
    print(f"{'PASS' if v else 'FAIL'}  {k}")
print('ALL PASS' if all(out.values()) else 'SOME FAIL')
