"""Post-passage Bromwich inversion with an estimated high-frequency tail.

The endpoint correction is matched to trapezoidal quadrature. It is a
large-(t+r_*) expansion, not a reconstruction of the null-front transient.
"""
import numpy as np
from numpy.polynomial import Polynomial
from w4_monopole import tortoise
from w4_threshold import invert_contour


def postpassage_waveform(u,K,times,sigma,extraction=200.,cutoff=None,terms=3,fit_width=.04,degree=3):
    u=np.asarray(u);K=np.asarray(K);times=np.asarray(times)
    if cutoff is not None:
        keep=u<=cutoff+1e-12;u=u[keep];K=K[keep]
    h=u[1]-u[0]
    if not np.allclose(np.diff(u),h,atol=1e-12,rtol=1e-8):raise ValueError('Uniform frequency grid required')
    if abs(u[0])>1e-12:raise ValueError('Include the imaginary-axis endpoint')
    if terms not in (1,2,3):raise ValueError('Supported tail orders are 1, 2, 3')
    x=tortoise(extraction);q=times+x
    if np.any(q<50) or np.any(abs(1-np.exp(-1j*q*h))<1e-6):
        raise ValueError('Tail estimate excludes the null front and quadrature alias points')
    raw=invert_contour(u,K,times,sigma)
    # Remove the known particle-trajectory phase before estimating endpoint derivatives.
    B=K*np.exp(1j*(u+1j*sigma)*x)
    count=max(degree+2,int(round(fit_width/h))+1)
    poly=Polynomial.fit(u[-count:]-u[-1],B[-count:],degree).convert()
    z=np.exp(-1j*q*h)
    weights=[h*(.5+z/(1-z)),h*h*z/(1-z)**2,h**3*z*(1+z)/(2*(1-z)**3)]
    correction=sum(poly.deriv(n)(0)*weights[n] for n in range(terms))
    tail=np.exp(sigma*q)*np.real(np.exp(-1j*u[-1]*q)*correction)/np.pi
    return raw+tail,raw,tail
