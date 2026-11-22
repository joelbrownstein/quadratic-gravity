"""Analytic-reference subtraction with oscillatory quadrature for the W4 bridge.

Unlike the front-local inversion, this does not subtract polynomial moments.
The small residual tail is estimated by integration by parts at the endpoint.
"""
import numpy as np
from scipy.interpolate import make_interp_spline
from numpy.polynomial import Polynomial
from numpy.polynomial.legendre import leggauss
from w4_front import analytic_front, front_transform


def bridge_waveform(u,B,v,sigma=.004,cutoff=None,tail_terms=3,fit_width=1.,
                    degree=4,quadrature_order=16,extraction=200.,mu=.5,spline_order=5):
    u=np.asarray(u);B=np.asarray(B);v=np.atleast_1d(v)
    if cutoff is not None:
        keep=u<=cutoff+1e-10;u=u[keep];B=B[keep]
    if u[0]!=0 or np.any(np.diff(u)<=0) or np.any(v<=0):
        raise ValueError('Ordered frequency grid from zero and positive advanced times required')
    if tail_terms not in (0,1,2,3):
        raise ValueError('Tail order must be 0, 1, 2 or 3')
    residual=B-front_transform(u+1j*sigma,extraction,mu)
    spline=make_interp_spline(u,residual,k=spline_order)
    nodes,weights=leggauss(quadrature_order)
    half=np.diff(u)/2
    points=((u[:-1]+half)[:,None]+half[:,None]*nodes).ravel()
    amplitudes=(half[:,None]*weights*spline(points.reshape(-1,quadrature_order))).ravel()
    raw=np.empty(v.shape)
    for start in range(0,len(v),32):
        vv=v[start:start+32]
        raw[start:start+32]=np.real(np.exp(-1j*np.outer(vv,points))@amplitudes)*np.exp(sigma*vv)/np.pi
    fit=u>=u[-1]-fit_width-1e-10
    if np.count_nonzero(fit)<degree+2:raise ValueError('Insufficient endpoint fit samples')
    poly=Polynomial.fit(u[fit]-u[-1],residual[fit],degree).convert()
    tail=np.zeros_like(v)
    if tail_terms:
        tail=np.real(np.exp((sigma-1j*u[-1])*v)*
                     sum(poly.deriv(n)(0)/(1j*v)**(n+1) for n in range(tail_terms)))/np.pi
    return {'K':analytic_front(v,extraction,mu)+raw+tail,'raw_K':analytic_front(v,extraction,mu)+raw,
            'tail':tail,'residual_integral':raw}
