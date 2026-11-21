"""Reduced polar ODE and reconstruction from independently generated equations.

Load w3_polar_operators.json, produced by w3_polar_system.py. Elimination uses
all equations with row scaling; residuals are checked, never silently discarded.
This is a local ODE/reconstruction interface, not a retarded waveform solver.
"""
import json
from pathlib import Path
import numpy as np
import sympy as s
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

DEFAULT_OPERATORS=Path(__file__).resolve().parents[2]/'notes/validation/w3_polar_operators.json'


class PolarSystem:
    def __init__(self,ell,path=DEFAULT_OPERATORS):
        data=json.loads(Path(path).read_text())
        if not isinstance(ell,int) or ell<0:raise ValueError('ell must be a nonnegative integer')
        self.record=next((b for b in data['systems'] if b['ell']==ell),data.get('general_l_ge_2'))
        if self.record is None:raise ValueError('Regenerate the operator artifact for this multipole')
        self.ell=ell;self.names=self.record['fields'];self.n=len(self.names)
        self.retained=['K'] if ell==0 else ['K','e'] if ell==1 else ['K','e','G']
        cols=self.record['columns'];self.d=len(self.retained)
        self.known=[cols.index(n+'_'+str(j)) for j in (0,1) for n in self.retained]
        self.unknown=[i for i in range(len(cols)) if i not in self.known]
        self.acceleration=[cols.index(n+'_2') for n in self.retained]
        matrix=s.Matrix(self.record['matrix']).subs(s.Symbol('L'),ell*(ell+1))
        self.evaluate=s.lambdify(s.symbols('r omega mu'),matrix,'numpy',cse=True)

    def reconstruction(self,r,omega,mu,rhs=None):
        """Return jets=B@y+b, where y=(retained fields, their r derivatives)."""
        A=np.asarray(self.evaluate(r,omega,mu),dtype=complex)
        U=A[:,self.unknown];V=A[:,self.known]
        source=np.zeros(len(A),complex) if rhs is None else np.array([rhs[n] for n in self.record['rows']],complex)
        scale=np.maximum(np.linalg.norm(U,axis=1),1e-30)
        target=np.column_stack((-V,source))/scale[:,None]
        X,_,rank,singular=np.linalg.lstsq(U/scale[:,None],target,rcond=None)
        if rank!=len(self.unknown):raise ArithmeticError('Singular reduction: change variables or parameters')
        residual=np.linalg.norm(U@X-np.column_stack((-V,source)))/max(1.,np.linalg.norm(V),np.linalg.norm(source))
        if residual>1e-8:raise ArithmeticError(f'Inconsistent constrained source: {residual}')
        B=np.zeros((3*self.n,2*self.d),complex);b=np.zeros(3*self.n,complex)
        B[self.known]=np.eye(2*self.d);B[self.unknown]=X[:,:-1];b[self.unknown]=X[:,-1]
        return B,b,{'equation_residual':float(residual),'condition_number':float(singular[0]/singular[-1]),
                   'eliminated_rank':int(rank),'retained_second_order_channels':self.d}

    def ode(self,r,omega,mu,rhs=None):
        B,b,info=self.reconstruction(r,omega,mu,rhs)
        top=np.column_stack((np.zeros((self.d,self.d)),np.eye(self.d)))
        return np.vstack((top,B[self.acceleration])),np.r_[np.zeros(self.d),b[self.acceleration]],info


def monopole_coefficients(r,omega,mu):
    """Independent published K equation: AGP2026 Eqs. (32)-(33), GM=1.

    This reference is compared with coordinate-derived elimination, not used
    to generate the operator. K''+P K'+V K=0.
    """
    f=1-2/r
    D=r*(2*mu**2*f+mu**4*r*r+4*omega**2)
    P=2*(4*mu**2*r*f*f+mu**4*r*r*(r-1)+4*omega**2*(2*r-3))/(r*f*D)
    V=(-f*(mu**6*r**6+6*mu**2*r**4*omega**2)
       +4*omega**2*(-12+6*r+r**4*omega**2)
       +mu**4*r*r*(r**4*omega**2-4*(8-6*r+r*r)))/(D*f*f*r**3)
    return P,V


def horizon_regular_coefficient(omega,mu):
    """K=(r-2)^rho [1+b1(r-2)+...], rho=-2i omega."""
    rho=-2j*omega
    p0=(2*mu**4-mu**2+6*omega**2)/(4*(mu**4+omega**2))
    vm1=-(4*mu**6-8*mu**4*omega**2-4*mu**4+8*mu**2*omega**2-8*omega**4-3*omega**2)/(2*(mu**4+omega**2))
    return -(p0*rho+vm1)/(2*rho+1)


def growth_residual(growth,mu,epsilon=1e-6,outer=80.,match=4.,rtol=2e-10):
    """Monopole unstable bound state, omega=i growth; NOT a damped QNM.

    Uses the independently checked closed K equation for fast integration.
    First-order horizon / leading infinity data are converged by moving both
    boundaries; no finite-radius boundary is called an exact outgoing one.
    """
    omega=1j*growth;kappa=np.sqrt(mu*mu+growth*growth)
    power=-(mu*mu+2*growth*growth)/kappa-1
    def fun(r,y):
        P,V=monopole_coefficients(r,omega,mu)
        return [y[1],(-P*y[1]-V*y[0]).real]
    b1=horizon_regular_coefficient(omega,mu).real
    h=solve_ivp(fun,(2+epsilon,match),[1.,2*growth/epsilon+b1/(1+b1*epsilon)],method='DOP853',rtol=rtol,atol=rtol*.01)
    o=solve_ivp(fun,(outer,match),[1.,-kappa+power/outer],method='DOP853',rtol=rtol,atol=rtol*.01)
    if not h.success or not o.success:raise RuntimeError('Monopole integration failed')
    return h.y[1,-1]/h.y[0,-1]-o.y[1,-1]/o.y[0,-1]


def unstable_mode(mu,bracket,**settings):
    root=brentq(lambda z:growth_residual(z,mu,**settings),*bracket,xtol=2e-12)
    return {'mu':mu,'omega_rg':[0.,float(root)],'settings':settings,
            'matching_residual':float(abs(growth_residual(root,mu,**settings)))}
