"""W2 massive axial Schwarzschild QNMs: matrix CF and complex-contour matching.

Units GM=1; exp(-i omega t); l=2 by default. No sourced response.
Requires numpy and scipy. See W2 note for continuation and validation scope.
"""
import cmath
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from axial_system import matrices, massive_W, horizon_series, infinity_series, outgoing_k


def resonance_k(w,mu):
    """Outgoing resonance sheet reached from omega>mu into Im(omega)<0.

    Below the real mass threshold this is not the bound-state sheet. The
    benchmark continuation stays at Re(omega)>0, Im(omega)<0 and away from k=0.
    """
    if w.real<=0: raise ValueError('Only positive-frequency branch is implemented')
    k=cmath.sqrt(w*w-mu*mu)
    return k if k.real>0 else -k


def recurrence(n,w,mu,ell=2):
    k=resonance_k(w,mu)
    I=np.eye(2,dtype=complex)
    A,B=map(lambda a:np.array(a,dtype=complex),matrices(ell))
    b0=1-4j*w
    b1=6j*k+8j*w-4+2j*w*w/k
    b2=-2j*k-4j*w+3-2j*w*w/k
    c0=4*k*k+12*k*w+3j*k+12*w*w+4j*w+4*w**3/k+1j*w*w/k
    c1=-k*k-4*k*w-2j*k-6*w*w-4j*w-4*w**3/k-2j*w*w/k-w**4/k**2
    return ((n+1)*(n+b0)*I,
            (-2*n*(n-1)+b1*n+c0)*I-A-B/2,
            ((n-1)*(n-2)+b2*(n-1)+c1)*I+B/2)


def cf_matrix(w,mu,depth=200,ell=2,tail=0.):
    R=tail*np.eye(2,dtype=complex)
    for n in range(depth-1,-1,-1):
        a,b,c=recurrence(n+1,w,mu,ell)
        R=-np.linalg.solve(b+a@R,c)
    a,b,_=recurrence(0,w,mu,ell)
    return b+a@R


def complex_root(function,guess,tolerance=1e-10):
    def fun(x):
        z=function(complex(*x))
        return [z.real,z.imag]
    result=root(fun,[guess.real,guess.imag],tol=tolerance)
    w=complex(*result.x)
    residual=abs(function(w))
    if residual>1e-7 or not np.isfinite(residual):
        raise RuntimeError(f'root failed: {result.message}; residual={residual}')
    return w,residual,result.nfev


def cf_root(mu,guess,depth=200,tail=0.):
    return complex_root(lambda w:np.linalg.det(cf_matrix(w,mu,depth,tail=tail)),guess)


def series_value(coeff,x,inverse=False):
    v=np.zeros(2,dtype=complex); d=v.copy()
    for n,c in enumerate(coeff):
        p=-n if inverse else n
        v+=np.array(c)*x**p
        if p: d+=p*np.array(c)*x**(p-1)
    return v,d


def boundary_basis(w,mu,r,order,ingoing,ell=2):
    columns=[]
    for leading in ((1,0),(0,1)):
        if ingoing:
            c=horizon_series(w,mu,ell,leading,order)
            v,d=series_value(c,r-2)
            d+=(-1j*w/(1-2/r))*v
        else:
            k,chi,c=infinity_series(w,mu,ell,leading,order,wave_number=resonance_k(w,mu))
            v,d=series_value(c,r,True)
            d+=(1j*k+chi/r)*v
        columns.append(np.concatenate((v,d)))
    return np.column_stack(columns)


def integrate_basis(w,mu,start,end,basis,rtol=2e-10,ell=2):
    delta=end-start
    def rhs(t,flat):
        r=start+t*delta; f=1-2/r; fp=2/r**2
        y=flat.reshape(4,2)
        ypp=-(fp/f)*y[2:]-(w*w*np.eye(2)-f*np.array(massive_W(r,ell,mu)))@y[:2]/f**2
        return delta*np.vstack((y[2:],ypp)).ravel()
    sol=solve_ivp(rhs,(0.,1.),basis.ravel(),method='DOP853',rtol=rtol,atol=rtol*1e-2)
    if not sol.success: raise RuntimeError(sol.message)
    return sol.y[:,-1].reshape(4,2),sol.nfev


def matching_matrix(w,mu,rm=5.,radius=40.,angle=.5,epsilon=.02,order=10,rtol=2e-10):
    """Ingoing real segment and outgoing upper-half-r contour meet at rm.

    Outgoing infinity is approached within a decay wedge; endpoints avoid
    r=0,2 and never cross the logarithmic branch cut.
    """
    rin=2+epsilon
    rout=rm+radius*cmath.exp(1j*angle)
    k=resonance_k(w,mu)
    if (k*cmath.exp(1j*angle)).imag<=0:
        raise ValueError('Outgoing contour must lie in a decay wedge')
    incoming,_=integrate_basis(w,mu,rin,rm,boundary_basis(w,mu,rin,order,True),rtol)
    outgoing,_=integrate_basis(w,mu,rout,rm,boundary_basis(w,mu,rout,order,False),rtol)
    matrix=np.column_stack((incoming,-outgoing))
    return matrix/np.linalg.norm(matrix,axis=0)


def matching_root(mu,guess,**kwargs):
    return complex_root(lambda w:np.linalg.det(matching_matrix(w,mu,**kwargs)),guess)
