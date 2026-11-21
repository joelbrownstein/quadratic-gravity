"""Null radial-infall monopole: finite-domain retarded Green-function solver.

GM=1, exp(-i omega t), t_particle=-r_*. Responses are per kappa E Y_00*.
Finite integration boundaries include asymptotic forced-tail data; moving
both boundaries is mandatory. It does not identify a band-limited waveform with the
complete retarded signal, nor treat unstable backgrounds as stationary ones.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import wrightomega
from w3_polar import horizon_regular_coefficient


def tortoise(r):
    return r+2*np.log(r/2-1)


def geometry(x):
    u=wrightomega(x/2-1).real
    return 2*(1+u),u/(1+u)


def coefficients(x,w,mu):
    r,F=geometry(x)
    D=r*(2*mu**2*F+mu**4*r*r+4*w*w)
    A=2*(4*mu**2*r*F*F+mu**4*r*r*(r-1)+4*w*w*(2*r-3))/(r*D)-2/r**2
    B=(-F*(mu**6*r**6+6*mu**2*r**4*w*w)+4*w*w*(-12+6*r+r**4*w*w)
       +mu**4*r*r*(r**4*w*w-4*(8-6*r+r*r)))/(D*r**3)
    T=2*F*(mu**2*r+2j*w)*np.exp(-1j*w*x)/(r*r*D)
    Wbase=D/r**5
    return A,B,T,Wbase


def divide_series(num,den,n):
    num=np.asarray(num,complex);den=np.asarray(den,complex);out=np.zeros(n,complex)
    for j in range(n):
        out[j]=((num[j] if j<len(num) else 0)-sum(den[k]*out[j-k] for k in range(1,min(j+1,len(den)))))/den[0]
    return out


def radial_series(w,mu,order):
    from numpy.polynomial import Polynomial as Poly
    z=Poly([0,1]);F=1-2*z
    D=Poly([mu**4,0,2*mu*mu+4*w*w,-4*mu*mu])
    Np=2*Poly([mu**4,-mu**4,4*mu*mu+8*w*w,-16*mu*mu-12*w*w,16*mu*mu])
    Nv=-F*(mu**6+6*mu*mu*w*w*z**2)+4*w*w*(w*w*z**2+6*z**5-12*z**6)+mu**4*(w*w-4*z**2+24*z**3-32*z**4)
    P=divide_series((z*Np).coef,(F*D).coef,order)
    V=divide_series(Nv.coef,(F*F*D).coef,order)
    return P,V,F,D,z


def outgoing_series(w,mu,order=6):
    """K=e^(ikr) r^p sum a_n/r^n; retarded positive-frequency branch."""
    k=np.sqrt(complex(w*w-mu*mu))
    if abs(k)<1e-10:raise ValueError('Threshold needs a separate uniform expansion')
    lam=1j*k;p=1j*(2*w*w-mu*mu)/k-1
    P,V,_,_,_=radial_series(w,mu,order+3)
    def contribution(n,j):
        out=(lam*P[j-n]+V[j-n]) if j>=n else 0j
        if j==n:out+=lam*lam
        if j>=n+1:out+=(p-n)*P[j-n-1]
        if j==n+1:out+=2*lam*(p-n)
        if j==n+2:out+=(p-n)*(p-n-1)
        return out
    a=[1.+0j]
    for n in range(1,order+1):
        a.append(-sum(a[j]*contribution(j,n+1) for j in range(n))/contribution(n,n+1))
    return k,p,np.array(a)


def outer_data(r,w,mu,order=6):
    k,p,a=outgoing_series(w,mu,order)
    powers=r**(-np.arange(len(a),dtype=float))
    series=np.dot(a,powers)
    derivative=np.dot(-np.arange(len(a))*a,powers)/r
    value=np.exp(1j*k*r+p*np.log(r))*series
    log_slope=(1-2/r)*(1j*k+p/r+derivative/series)
    return value,log_slope


def source_tail_series(w,mu,order=10):
    P,V,F,D,z=radial_series(w,mu,order+1)
    drift=P-2j*w*divide_series([1],F.coef,order+1)
    potential=(V+divide_series((-w*w+2j*w*z*z).coef,(F*F).coef,order+1)
               -1j*w*np.convolve(P,divide_series([1],F.coef,order+1))[:order+1])
    S=divide_series((2*z**4*(mu*mu+2j*w*z)).coef,(F*D).coef,order+1)
    b=np.zeros(order+1,complex)
    for n in range(order+1):
        known=sum(b[j]*(potential[n-j]-(j*drift[n-j-1] if n>=j+1 else 0)
                         +(j*(j+1) if n==j+2 else 0)) for j in range(n))
        b[n]=(S[n]-known)/potential[0]
    return b


def source_tail(r,w,mu,order=10):
    """Forced far-zone solution e^(-i omega r_*) sum b_n/r^n.

    It starts at n=4 and carries no freely specified incoming massive wave.
    Subtract it before extracting the homogeneous outgoing amplitude.
    """
    b=source_tail_series(w,mu,order)
    powers=r**(-np.arange(order+1,dtype=float));B=np.dot(b,powers)
    Br=np.dot(-np.arange(order+1)*b,powers)/r
    phase=np.exp(-1j*w*tortoise(r))
    return phase*B,phase*((1-2/r)*Br-1j*w*B)


def response(w,mu,extraction=25.,left=-35.,outer=400.,order=6,rtol=2e-10,method='green',tail_order=10):
    if w<=mu or mu<.44:raise ValueError('This validated stationary benchmark requires omega>mu and mu>=0.44')
    if not 2<extraction<outer:raise ValueError('Require 2 < extraction < outer')
    xl=left;xo=tortoise(extraction);xr=tortoise(outer)
    rl,fl=geometry(xl);eps=2*fl/(1-fl);rho=-2j*w
    b1=horizon_regular_coefficient(w,mu)
    slope=fl*(rho/eps+b1/(1+b1*eps))
    oval,oslope=outer_data(outer,w,mu,order)
    tail,tail_x=source_tail(outer,w,mu,tail_order) if tail_order else (0j,0j)
    hpart=np.exp(-1j*w*xl)*eps/(16*(mu*mu-1j*w)*(.25-1j*w))
    hpart_x=(.5-1j*w)*hpart
    options=dict(method='DOP853',rtol=rtol,atol=rtol*.01,dense_output=True)
    if method=='green':
        def rhs(x,y):
            A,B,T,W=coefficients(x,w,mu)
            return [y[1],-A*y[1]-B*y[0],y[0]*T/W]
        initial_integral=(hpart_x-slope*hpart)/coefficients(xl,w,mu)[3]
        h=solve_ivp(rhs,(xl,xr),[1.+0j,slope,initial_integral],**options)
        o=solve_ivp(rhs,(xr,xo),[1.+0j,oslope,0j],**options)
        if not h.success or not o.success:raise ArithmeticError('Green integration failed')
        hv=h.sol(xo);ov=o.y[:,-1];base=coefficients(xo,w,mu)[3]
        C=(hv[0]*ov[1]-hv[1]*ov[0])/base
        value=(ov[0]*hv[2]-hv[0]*ov[2])/C
        derivative=(ov[1]*hv[2]-hv[1]*ov[2])/C
        outgoing=h.y[2,-1]/C/oval if w>mu else None
        correction=(tail_x-oslope*tail)/(h.y[1,-1]-oslope*h.y[0,-1])
        value+=correction*hv[0];derivative+=correction*hv[1]
        if outgoing is not None:outgoing+=(correction*h.y[0,-1]-tail)/oval
        # Abel-Wronskian check at a second interior radius, independently
        # continuing the outer homogeneous basis to that point.
        xm=(xo+xr)/2
        hm=h.sol(xm);om=o.sol(xm)
        C2=(hm[0]*om[1]-hm[1]*om[0])/coefficients(xm,w,mu)[3]
        wr_error=abs(C2/C-1)
    elif method=='shoot':
        def rhs(x,y):
            A,B,T,_=coefficients(x,w,mu)
            return [y[1],-A*y[1]-B*y[0],y[3],-A*y[3]-B*y[2]+T]
        sol=solve_ivp(rhs,(xl,xr),[1.+0j,slope,hpart,hpart_x],**options)
        if not sol.success:raise ArithmeticError('Sourced shooting integration failed')
        end=sol.y[:,-1];c=(tail_x-oslope*tail-end[3]+oslope*end[2])/(end[1]-oslope*end[0])
        y=sol.sol(xo);value=y[2]+c*y[0];derivative=y[3]+c*y[1]
        outgoing=(end[2]+c*end[0]-tail)/oval if w>mu else None
        wr_error=None
    else:raise ValueError('Unknown method')
    return {'K':value,'K_x':derivative,'outgoing':outgoing,'wronskian_relative_error':wr_error}


def band_waveform(frequencies,amplitudes,times):
    """Real band-limited K, inverse d omega/(2pi); includes conjugate band."""
    return np.trapezoid(np.real(np.asarray(amplitudes)[None,:]*np.exp(-1j*np.outer(times,frequencies))),frequencies,axis=1)/np.pi


def signed_flux(w,mu,outgoing,energy=1.):
    """AGP Eq.(58), converted to our raw Fourier convention and G=1.

    outgoing is per kappa E Y00*. The multiplication restores kappa=8pi and
    Y00=1/sqrt(4pi). Negative sign is retained, not relabeled positive energy.
    """
    if w<=mu:return 0.
    physical=4*np.sqrt(np.pi)*energy*outgoing/np.sqrt(2*np.pi)
    return -8*w*(w*w-mu*mu)**1.5/mu**2*abs(physical)**2
