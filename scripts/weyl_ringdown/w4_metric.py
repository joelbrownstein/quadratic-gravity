"""Sourced monopole metric and relational curvature, per N= kappa E Y00*.

Field coefficients still multiply Y00. Coordinates are Schwarzschild t,r;
the Einstein sector is the ingoing Vaidya perturbation written in this chart.
Tidal perturbations compare events at fixed areal radius and an instantaneous
radial observer velocity beta relative to the exterior Kodama frame.
"""
from functools import lru_cache
import numpy as np
import sympy as s
from w3_polar import monopole_coefficients
from w4_monopole import tortoise


@lru_cache(None)
def expressions():
    r,w,mu,K,J,Z=s.symbols('r omega mu K J Z')
    F=1-2/r;P,V=monopole_coefficients(r,w,mu)
    D=r*(2*mu**2*F+mu**4*r**2+4*w*w)
    S=2*(mu*mu*r+2*s.I*w)/(r*r*F*D)*Z
    Kpp=-P*J-V*K+S
    deriv=lambda e:s.diff(e,r)+s.diff(e,K)*J+s.diff(e,J)*Kpp-s.I*w*Z*s.diff(e,Z)/F
    linear=lambda e:sum(s.factor(s.diff(e,z))*z for z in (K,J,Z))
    a=linear(-r*r/(2*F)*(F*Kpp+2*(r-1)*J/r**2+
                        (w*w/F-mu*mu+(20-6*r)/r**3)*K))
    c=a-2*K
    mq=linear((3-r)*K+(r-2)*c-r*r*F*J)
    # The time Misner-Sharp constraint gives b without differentiating noisy
    # radial data. Full tensor/constraint reconstruction is checked separately.
    b=linear((2*s.I*w*mq-Z)/(mu*mu*r*r*F))
    return (r,w,mu,K,J,Z),{'a':a,'b':b,'c':c,'K':K,'mass_q':mq},deriv


@lru_cache(None)
def evaluator():
    args,expr,_=expressions()
    return s.lambdify(args,list(expr.values()),'numpy',cse=True)


def reconstruct(r,w,mu,K,K_x,beta=0.):
    if r<=2 or abs(beta)>=1 or np.any(np.asarray(w)==0):
        raise ValueError('Exterior radius, subluminal radial observer and nonzero contour frequency required')
    w=np.asarray(w,dtype=complex);F=1-2/r
    Z=np.exp(-1j*w*tortoise(r))
    a,b,c,k,mq=evaluator()(r,w,mu,K,np.asarray(K_x)/F,Z)
    # Retarded transform of Theta(t+r_*); the contour lies above omega=0.
    step=Z/(-1j*w)
    Htt=step/r
    mt=step/2+mq
    Ricci_uu=-mu*mu*(a+2*beta*b+beta*beta*c)/(1-beta*beta)
    Ricci_ang=-mu*mu*k
    Er_q=2*Ricci_ang-2*mq/r**3
    Et_q=Ricci_uu/2-Ricci_ang+mq/r**3
    return {'a':a,'b':b,'c':c,'K':k,
            'h_tt':Htt+2*F*a,'h_tr':Htt/F+2*b,'h_rr':Htt/F**2+2*c/F,
            'h_angular_fraction':2*k,'mass_q':mq,'mass_total':mt,
            'Ricci_uu':Ricci_uu,'Ricci_angular':Ricci_ang,
            'Ricci_null_out':-mu*mu*(a+2*b+c),
            'Ricci_null_in':-mu*mu*(a-2*b+c),
            'tidal_radial_q':Er_q,'tidal_transverse_q':Et_q,
            'tidal_radial_total':Er_q-step/r**3,
            'tidal_transverse_total':Et_q+step/(2*r**3)}
