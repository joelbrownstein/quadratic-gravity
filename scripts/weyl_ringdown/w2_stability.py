"""Limited radial stability checks using Brito et al. (2013), Eq. (30).

The polar master potential is imported explicitly, not derived from the axial
system. A finite-box positive eigenvalue is not a proof of full stability.
"""
import numpy as np
from scipy.special import lambertw
from scipy.linalg import eigh_tridiagonal
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


def radius(tortoise):
    u=lambertw(np.exp(np.asarray(tortoise)/2-1)).real
    # Compute F without cancellation near r=2.
    return 2*(1+u),u/(1+u)


def potential(tortoise,mu):
    r,f=radius(tortoise)
    return f*(2/r**3+mu*mu+(24*(1-r)*mu*mu+6*r**3*(r-4)*mu**4)/(2+r**3*mu*mu)**2)


def lowest_eigenvalue(mu,left=-500.,right=150.,spacing=.05):
    intervals=round((right-left)/spacing)
    x=np.linspace(left,right,intervals+1)[1:-1]
    h=(right-left)/intervals
    diagonal=2/h**2+potential(x,mu)
    eigen=eigh_tridiagonal(diagonal,np.full(len(x)-1,-1/h**2),
                          select='i',select_range=(0,0),check_finite=False,eigvals_only=True)[0]
    return {'mu':mu,'left':left,'right':right,'spacing':h,'omega_squared':float(eigen),
            'growth_rate':float(np.sqrt(-eigen)) if eigen<0 else None,
            'interpretation':'unstable bound state' if eigen<0 else 'no negative eigenvalue in this finite-box discretization'}


def zero_residual(mu,left=-60.,right=100.,match=0.,rtol=1e-10):
    def rhs(x,y): return [y[1],potential(x,mu)*y[0]]
    a=solve_ivp(rhs,(left,match),[1.,0.],method='DOP853',rtol=rtol,atol=rtol*.01)
    r,f=radius(right)
    b=solve_ivp(rhs,(right,match),[1.,f*(-mu-mu/r)],method='DOP853',rtol=rtol,atol=rtol*.01)
    if not a.success or not b.success: raise RuntimeError('Zero-mode integration failed')
    y=a.y[:,-1];z=b.y[:,-1]
    return y[1]/y[0]-z[1]/z[0]


def critical_mass(**kwargs):
    mu=brentq(lambda m:zero_residual(m,**kwargs),.42,.46,xtol=2e-12)
    return {'mu_critical':mu,'p_critical':2*mu,'settings':kwargs,
            'matching_residual':abs(zero_residual(mu,**kwargs))}
