"""Finite-radius W4 response across the massive threshold.

Uses scaled Green integration and an independent collocation BVP. Complex
frequency is in the closed upper half-plane; no real-axis outgoing amplitude
or flux is assigned to subthreshold/threshold fields.
"""
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp,solve_bvp
from w3_polar import horizon_regular_coefficient,monopole_coefficients
from w4_monopole import (geometry,tortoise,coefficients,radial_series,
                         outgoing_series,source_tail)


def threshold_slope(r,mu,order=14):
    """Outgoing Riccati expansion for u=sqrt(r^5 F/D) K at omega=mu."""
    count=order+3
    P,V,*_=radial_series(mu,mu,count)
    Q=V-.25*np.convolve(P,P)[:count]
    for j in range(1,count):Q[j]+=.5*(j-1)*P[j-1]
    a=np.zeros(order+1,complex);a[1]=1j*np.sqrt(Q[1])
    for n in range(2,order+1):
        q=Q[(n+1)//2] if (n+1)%2==0 else 0.
        cross=sum(a[j]*a[n+1-j] for j in range(2,n))
        a[n]=(.5*(n-1)*a[n-1]-cross-q)/(2*a[1])
    p,_=monopole_coefficients(r,mu,mu)
    return (1-2/r)*(sum(a[n]*r**(-n/2) for n in range(1,order+1))-p/2)


def coulomb_slope(r,w,mu,dps=35):
    """Whittaker boundary retaining k^2+a/r+b/r^2 in the normal form."""
    k=np.sqrt(complex(w*w-mu*mu))
    a=2*(2*w*w-mu*mu);b=12*w*w-4*mu*mu-2-4*w*w/(mu*mu)
    with mp.workdps(dps):
        kk=mp.mpc(complex(k));kap=1j*mp.mpc(complex(a))/(2*kk)
        nu=mp.sqrt(mp.mpc(complex(.25-b)));z=-2j*kk*r
        fun=lambda zz:mp.whitw(kap,nu,zz)
        logu=complex(-2j*kk*mp.diff(fun,z)/fun(z))
    P,_=monopole_coefficients(r,w,mu)
    return (1-2/r)*(logu-P/2)


def boundary_slope(r,w,mu,order=10,threshold_order=14):
    k=np.sqrt(complex(w*w-mu*mu))
    if abs(k)<1e-10:return threshold_slope(r,mu,threshold_order),'threshold Riccati'
    if abs(k)<.18:return coulomb_slope(r,w,mu),'Whittaker Coulomb'
    k,p,a=outgoing_series(w,mu,order)
    n=np.arange(len(a));terms=a*r**(-n.astype(float));S=sum(terms)
    return (1-2/r)*(1j*k+p/r-sum(n*terms)/(r*S)),'inverse-radius outgoing/decaying'


def boundary_data(w,mu,left,outer,order=10,tail_order=16,threshold_order=14):
    r,F=geometry(left);eps=2*F/(1-F)
    b1=horizon_regular_coefficient(w,mu)
    left_slope=F*(-2j*w/eps+b1/(1+b1*eps))
    part=np.exp(-1j*w*left)*eps/(16*(mu*mu-1j*w)*(.25-1j*w))
    left_drive=(.5-1j*w-left_slope)*part
    right_slope,label=boundary_slope(outer,w,mu,order,threshold_order)
    part,px=source_tail(outer,w,mu,tail_order)
    return left_slope,left_drive,right_slope,px-right_slope*part,label


def scaled_propagate(fun,start,end,y,rtol,chunk):
    grid=np.linspace(start,end,max(1,int(np.ceil(abs(end-start)/chunk)))+1)
    y=np.array(y,complex);steps=0
    for a,b in zip(grid[:-1],grid[1:]):
        sol=solve_ivp(fun,(a,b),y,method='DOP853',rtol=rtol,atol=rtol*.01)
        if not sol.success:raise ArithmeticError(sol.message)
        y=sol.y[:,-1];scale=max(abs(y[0]),abs(y[1]))
        if not np.isfinite(scale) or scale==0:raise ArithmeticError('Degenerate scaled basis')
        y/=scale;steps+=sol.nfev
    return y,steps


def response(w,mu=.5,extraction=200.,left=-35.,outer=800.,rtol=2e-10,
             order=10,tail_order=16,threshold_order=14,chunk=40.,method='green'):
    w=complex(w)
    if w.real<0 or w.imag<0 or mu<.44:raise ValueError('Require Re omega, Im omega >=0 and mu>=0.44')
    if not 2<extraction<outer:raise ValueError('Require 2 < extraction < outer')
    xo=tortoise(extraction);xr=tortoise(outer)
    ls,ld,rs,rd,label=boundary_data(w,mu,left,outer,order,tail_order,threshold_order)
    if method=='green':
        def rhs(x,y):
            A,B,T,W=coefficients(x,w,mu)
            return [y[1],-A*y[1]-B*y[0],y[0]*T/W]
        h,nh=scaled_propagate(rhs,left,xo,[1,ls,ld/coefficients(left,w,mu)[3]],rtol,chunk)
        o,no=scaled_propagate(rhs,xr,xo,[1,rs,rd/coefficients(xr,w,mu)[3]],rtol,chunk)
        C=(h[0]*o[1]-h[1]*o[0])/coefficients(xo,w,mu)[3]
        K=(o[0]*h[2]-h[0]*o[2])/C
        Kx=(o[1]*h[2]-h[1]*o[2])/C
        diagnostic={'function_evaluations':nh+no,'scaled_matching_determinant':float(abs(C))}
    elif method=='bvp':
        mesh=np.linspace(left,xr,max(200,int((xr-left)/.5)))
        def real_matrix(J):
            return np.concatenate((np.concatenate((J.real,-J.imag),axis=1),
                                   np.concatenate((J.imag,J.real),axis=1)),axis=0)
        def fun(x,y):
            A,B,T,_=coefficients(x,w,mu)
            z=y[:2]+1j*y[2:]
            value=np.array([z[1],-A*z[1]-B*z[0]+T])
            return np.concatenate((value.real,value.imag))
        def jac(x,y):
            A,B,T,_=coefficients(x,w,mu)
            J=np.zeros((2,2,len(x)),complex);J[0,1]=1;J[1,0]=-B;J[1,1]=-A
            return real_matrix(J)
        def bc(a,b):
            a=a[:2]+1j*a[2:];b=b[:2]+1j*b[2:]
            z=np.array([a[1]-ls*a[0]-ld,b[1]-rs*b[0]-rd])
            return np.r_[z.real,z.imag]
        bcjac=lambda a,b:(real_matrix(np.array([[-ls,1],[0,0]],complex)),real_matrix(np.array([[0,0],[-rs,1]],complex)))
        sol=solve_bvp(fun,bc,mesh,np.zeros((4,len(mesh))),fun_jac=jac,bc_jac=bcjac,
                      tol=rtol,bc_tol=rtol*.01,max_nodes=60000)
        if not sol.success:raise ArithmeticError(sol.message)
        z=sol.sol(xo);K,Kx=z[:2]+1j*z[2:]
        diagnostic={'mesh_nodes':len(sol.x),'max_collocation_residual':float(max(sol.rms_residuals))}
    else:raise ValueError('Use green or bvp')
    return {'K':complex(K),'K_x':complex(Kx),'boundary':label,**diagnostic}


def invert_contour(frequencies,amplitudes,times,sigma):
    """Bromwich integral and conjugate half, for a real time-domain field."""
    return np.exp(sigma*np.asarray(times))*np.trapezoid(
        np.real(np.asarray(amplitudes)[None,:]*np.exp(-1j*np.outer(times,frequencies))),frequencies,axis=1)/np.pi
