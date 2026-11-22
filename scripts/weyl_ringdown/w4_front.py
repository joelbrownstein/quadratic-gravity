"""Causal front subtraction for the null-infall auxiliary monopole.

The analytic reference has the Schwarzschild front's leading v log v and v
terms. The curved remainder is computed from the full W4 frequency response.
Here v=t+r_* at the observer; all amplitudes are per kappa E Y00*.
"""
import numpy as np
from scipy.special import exp1,k0
from scipy.integrate import simpson
from w4_monopole import tortoise


def analytic_front(v,r=200.,mu=.5):
    v=np.asarray(v,dtype=float)
    out=np.zeros_like(v);positive=v>0
    out[positive]=-v[positive]*k0(mu*np.sqrt(2*r*v[positive]))/r**2
    return out


def leading_front(v,r=200.,mu=.5):
    v=np.asarray(v,dtype=float);out=np.zeros_like(v);positive=v>0
    out[positive]=v[positive]*(np.log(mu*mu*r*v[positive]/2)+2*np.euler_gamma)/(2*r*r)
    return out


def front_transform(w,r=200.,mu=.5):
    """Laplace transform after removing exp(-i omega r_*); Re s > 0."""
    w=np.asarray(w,dtype=complex);s=-1j*w;a=mu*mu*r/2;z=a/s
    H=np.empty_like(z);regular=abs(z)<60
    H[regular]=(1+z[regular])*np.exp(z[regular])*exp1(z[regular])-1
    zz=z[~regular];value=np.zeros_like(zz);factorial=1
    # Stable expansion of (1+z) exp(z) E1(z)-1, avoiding cancellation.
    for n in range(2,22):
        factorial*=n-1
        value+=(-1)**n*(n-1)*factorial/zz**n
    H[~regular]=value
    return -H/(2*r*r*s*s)


def phase_removed(data,sigma,r=200.):
    u=data[:,0]
    return (data[:,1]+1j*data[:,2])*np.exp(1j*(u+1j*sigma)*tortoise(r))


def front_waveform(u,B,v,sigma=.004,r=200.,mu=.5,cutoff=None,moment_correct=True):
    """Front-local inverse with exact R(0)=R'(0)=0 moment conditions.

    Returns both the corrected result and the unprojected finite-cutoff
    remainder, so the causality/moment truncation error remains inspectable.
    This is a local v-window method, not a replacement for the late waveform.
    """
    u=np.asarray(u);B=np.asarray(B);v=np.asarray(v)
    if cutoff is not None:
        keep=u<=cutoff+1e-10;u=u[keep];B=B[keep]
    s=sigma-1j*u
    residual=B-front_transform(u+1j*sigma,r,mu)
    kernel=np.exp(np.outer(v,s))
    raw=simpson(np.real(kernel*residual[None,:]),x=u,axis=1)/np.pi
    m0=float(simpson(residual.real,x=u)/np.pi)
    m1=float(simpson(np.real(s*residual),x=u)/np.pi)
    corrected=raw-m0-v*m1 if moment_correct else raw
    answer=analytic_front(v,r,mu)+corrected
    if moment_correct:answer=np.where(v>0,answer,0.)
    return {'K':answer,'raw_K':analytic_front(v,r,mu)+raw,
            'remainder':corrected,'zeroth_moment':m0,'first_moment':m1}
