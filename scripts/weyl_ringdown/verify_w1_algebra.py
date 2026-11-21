#!/usr/bin/env python3
"""Independent coordinate checks of W1. Requires sympy==1.14.0.

uv run --no-project --with sympy==1.14.0 python scripts/weyl_ringdown/verify_w1_algebra.py
No QNM package or literature coefficient table is used by this verifier.
"""
import json
from functools import lru_cache
import sympy as s

r,theta,w,mu=s.symbols('r theta omega mu', nonzero=True)
coords=s.symbols('t'),r,theta,s.symbols('phi')
f=1-2/r
g=s.diag(-f,1/f,r*r,r*r*s.sin(theta)**2)
gi=g.inv()


def clean(x):
    return s.factor(s.trigsimp(s.expand_trig(x)))


@lru_cache(None)
def gamma(a,b,c):
    return clean(sum(gi[a,d]*(s.diff(g[d,c],coords[b])+s.diff(g[d,b],coords[c])
                              -s.diff(g[b,c],coords[d]))/2 for d in range(4)))


@lru_cache(None)
def riemann(a,b,c,d):
    return clean(s.diff(gamma(a,d,b),coords[c])-s.diff(gamma(a,c,b),coords[d])
                 +sum(gamma(a,c,e)*gamma(e,d,b)-gamma(a,d,e)*gamma(e,c,b) for e in range(4)))


def harmonic_diff(expr,index):
    return -s.I*w*expr if index==0 else s.diff(expr,coords[index])


def projection(ell):
    # m=0 is sufficient for a coordinate check at each l by spherical symmetry.
    # We test l=1,2,3, including the exceptional tensor-free dipole.
    Y=s.legendre(ell,s.cos(theta))
    Sv=-s.sin(theta)*s.diff(Y,theta)
    Tv=s.sin(theta)*(s.diff(Y,theta,2)-s.cot(theta)*s.diff(Y,theta))
    Q,Z=s.Function('Q')(r),s.Function('Z')(r)
    L=ell*(ell+1)
    if ell==1:
        Z=s.S.Zero
    F1=Q/f
    F2=r*Z
    F0=s.I*f/w*(s.diff(Q,r)+2*Q/r+(L-2)*Z/r)
    q=s.zeros(4)
    q[0,3]=q[3,0]=F0*Sv
    q[1,3]=q[3,1]=F1*Sv
    q[2,3]=q[3,2]=F2*Tv

    @lru_cache(None)
    def nabla(a,b,c):
        return clean(harmonic_diff(q[b,c],a)-sum(gamma(e,a,b)*q[e,c]+gamma(e,a,c)*q[b,e] for e in range(4)))

    trace=clean(sum(gi[a,a]*q[a,a] for a in range(4)))
    divergence=clean(sum(gi[a,a]*nabla(a,a,3) for a in range(4)))
    assert trace==0 and divergence==0

    def operator(b,c):
        box=0
        for a in range(4):
            value=harmonic_diff(nabla(a,b,c),a)
            value-=sum(gamma(e,a,a)*nabla(e,b,c)+gamma(e,a,b)*nabla(a,e,c)
                       +gamma(e,a,c)*nabla(a,b,e) for e in range(4))
            box+=gi[a,a]*value
        curvature=sum(2*g[b,b]*riemann(b,a,c,d)*gi[a,a]*gi[d,d]*q[a,d]
                      for a in range(4) for d in range(4) if q[a,d]!=0)
        return clean(box+curvature-mu*mu*q[b,c])

    W11=mu*mu+(L+4)/r**2-16/r**3
    W12=2*(L-2)*(r-3)/r**3
    W21=2/r**2
    W22=mu*mu+(L-2)/r**2+2/r**3
    EQ=f*f*s.diff(Q,r,2)+f*s.diff(f,r)*s.diff(Q,r)+(w*w-f*W11)*Q-f*W12*Z
    EZ=f*f*s.diff(Z,r,2)+f*s.diff(f,r)*s.diff(Z,r)+(w*w-f*W22)*Z-f*W21*Q
    assert clean(operator(1,3)/Sv-EQ/f**2)==0
    if ell>=2:
        assert clean(operator(2,3)/Tv-r*EZ/f)==0
    # Check the remaining time-angular equation on the reduced ODE system.
    qpp=s.solve(EQ,s.diff(Q,r,2))[0]
    replacements={s.diff(Q,r,2):qpp}
    if ell>=2:
        replacements[s.diff(Z,r,2)]=s.solve(EZ,s.diff(Z,r,2))[0]
    qppp=s.diff(qpp,r).subs(replacements,simultaneous=True)
    remaining=(operator(0,3)/Sv).subs(s.diff(Q,r,3),qppp).subs(replacements,simultaneous=True)
    assert s.simplify(clean(remaining))==0
    return {'ell':ell,'trace':True,'divergence':True,'radial_angular_equation':True,
            'tensor_angular_equation':ell>=2,'time_angular_equation_on_shell':True}


def gr_recurrence():
    x,L=s.symbols('x L')
    rr=2/(1-x)
    dxdr=(1-x)**2/2
    A=s.Function('A')(x)
    lp=s.I*w+2*s.I*w/rr-2*s.I*w*dxdr/x
    yr=dxdr*s.diff(A,x)+lp*A
    yrr=dxdr*s.diff(yr,x)+lp*yr
    expr=x*x*yrr+x*(2/rr**2)*yr+(w*w-x*(L/rr**2-6/rr**3))*A
    P=x*(1-x)**2
    B=1-4*s.I*w+(-4+16*s.I*w)*x+(3-8*s.I*w)*x*x
    C=32*w*w+8*s.I*w+3-L+(-16*w*w-8*s.I*w-3)*x
    assert s.simplify(expr/(x*(1-x)**2/4)-P*s.diff(A,x,2)-B*s.diff(A,x)-C*A)==0
    return True


if __name__=='__main__':
    print(json.dumps({'sympy':s.__version__,'checks':'exact symbolic coordinate identities',
                      'massive_axial_projection':[projection(l) for l in (1,2,3)],
                      'GR_transformed_equation':gr_recurrence()},indent=2))
