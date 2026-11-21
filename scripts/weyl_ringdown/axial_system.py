"""W1 massive axial operator and local series; no massive QNM search.

Units GM=1. Fields U=(Q,Z), q_rA=Q/F S_A, q_AB=-2 r Z S_AB.
These helpers describe l>=2; the dipole is a separate scalar equation.
"""
import cmath
import math


def add(a,b):
    return tuple(x+y for x,y in zip(a,b))


def scale(c,a):
    return tuple(c*x for x in a)


def mv(a,v):
    return tuple(sum(x*y for x,y in zip(row,v)) for row in a)


def matrices(ell):
    if ell < 2:
        raise ValueError('Two-field system requires ell>=2; dipole has no Z field')
    L=ell*(ell+1)
    return ((L+4,2*(L-2)),(2,L-2)), ((-16,-6*(L-2)),(0,2))


def massive_W(r,ell,mu):
    A,B=matrices(ell)
    return tuple(tuple(mu*mu*(i==j)+A[i][j]/r**2+B[i][j]/r**3
                       for j in range(2)) for i in range(2))


def horizon_series(w,mu,ell,leading=(1+0j,0j),order=6):
    """U=e^(-i w r*) sum a_n (r-2)^n; arbitrary ingoing leading vector."""
    A,B=matrices(ell)
    def wp(j):
        c2=(-1)**j*(j+1)/2**(j+2)
        c3=(-1)**j*math.comb(j+2,2)/2**(j+3)
        return tuple(tuple(mu*mu*(i==k)*(j==0)+c2*A[i][k]+c3*B[i][k]
                           for k in range(2)) for i in range(2))
    coeff=[tuple(leading)]
    for n in range(1,order+1):
        rhs=(0j,0j)
        for j in range(n):
            rhs=add(rhs,mv(wp(j),coeff[n-1-j]))
        for m in range(1,n):
            j=n+1-m
            rhs=add(rhs,scale(-n*m*(-1)**(j-1)/2**j,coeff[m]))
        coeff.append(scale(1/(n*(n/2-2j*w)),rhs))
    return coeff


def outgoing_k(w,mu):
    """Continuation from positive real omega>mu to Re(w)>mu, Im(w)<=0.

    This deliberately rejects threshold/other sheets; it is not a global
    retarded-branch implementation or a quasibound-state prescription.
    """
    if mu < 0 or w.real <= mu:
        raise ValueError('This helper only implements the positive-frequency above-threshold sheet')
    k=cmath.sqrt(w*w-mu*mu)
    return k if k.real>0 else -k


def infinity_series(w,mu,ell,leading=(1+0j,0j),order=6,wave_number=None):
    """U=e^(i k r) r^chi sum b_n/r^n; local asymptotic expansion."""
    A,B=matrices(ell)
    k=outgoing_k(w,mu) if wave_number is None else wave_number
    chi=1j*(2*w*w-mu*mu)/k
    def mat(c,a=0,b=0):
        return tuple(tuple(c*(i==j)+a*A[i][j]+b*B[i][j] for j in range(2)) for i in range(2))
    C={2:mat(chi*chi-chi-8j*k*chi-4*k*k+2j*k,-1),
       3:mat(8j*k*chi-4*chi*chi+6*chi-4j*k,2,-1),
       4:mat(4*chi*chi-8*chi,0,2)}
    T=(1,-4,4)
    D={1:2*chi-8j*k,2:8j*k-8*chi+2,3:8*chi-4}
    coeff=[tuple(leading)]
    for n in range(1,order+1):
        rhs=(0j,0j)
        for j,t in enumerate(T):
            m=n-1-j
            if m>=0: rhs=add(rhs,scale(m*(m+1)*t,coeff[m]))
        for j,d in D.items():
            m=n-j
            if m>=0: rhs=add(rhs,scale(-m*d,coeff[m]))
        for j,c in C.items():
            m=n+1-j
            if m>=0: rhs=add(rhs,mv(c,coeff[m]))
        coeff.append(scale(1/(2j*k*n),rhs))
    return k,chi,coeff
