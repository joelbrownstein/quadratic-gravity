"""Reproducible W3 algebra, source, reconstruction and spectral checks."""
import argparse
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
import scipy
import sympy as s
from scipy.integrate import quad
from w3_polar import PolarSystem,monopole_coefficients,unstable_mode,DEFAULT_OPERATORS,horizon_regular_coefficient
from w3_infall_source import source_data,forcing,trajectory_derivative
from w2_stability import lowest_eigenvalue


def source_checks():
    r,w,E,g=s.symbols('r omega E gamma',positive=True)
    f=1-2/r;v=s.sqrt(g*g-f);tp=-g/(f*v)
    # Every amplitude multiplies exp(i omega t_p). Differentiate that phase.
    D=lambda a:s.diff(a,r)+s.I*w*tp*a
    tt=E*g/(r*r*v);tr=E/(r*r*f);rr=E*v/(g*r*r*f*f)
    tau=-E/(g*r*r*v)
    a,b,c=tt/f,tr,f*rr
    divt=s.I*w*a+f*D(b)+(s.diff(f,r)+2*f/r)*b
    divr=s.I*w*b/f+D(c)+s.diff(f,r)*(a+c)/(2*f)+2*c/r
    assert s.simplify(divt)==0 and s.simplify(divr)==0
    assert s.simplify(-tt/f+f*rr-tau)==0
    assert s.simplify(s.limit(tt,g,s.oo)-E/r**2)==0
    assert s.simplify(s.limit(rr,g,s.oo)-E/(r*r*f*f))==0
    assert s.limit(tau,g,s.oo)==0
    reference=s.lambdify((r,w,E,g),[tt,tr,rr,tau,D(tau),D(D(tau)),tp],'numpy')
    errors=[]
    for rv in (2.05,3.,8.,30.):
        for gv in (1.,2.,10.):
            actual=source_data(rv,.3-.04j,energy=.001,gamma=gv)
            aa=np.array([actual[k] for k in ('tt','tr','rr','trace','trace_prime','trace_second','time_prime')])
            bb=np.array(reference(rv,.3-.04j,.001,gv))
            errors.append(float(np.max(abs(aa-bb)/np.maximum(abs(bb),1e-14))))
    assert max(errors)<1e-10
    # GR polar/axial intertwining, separate from the massive sector.
    lam=s.symbols('lambda',positive=True)
    C=lam*(lam+1)/3
    W=C+3*f/(r*(lam*r+3))
    rw=f*((2*lam+2)/r**2-6/r**3)
    z=2*f*(lam**2*(lam+1)*r**3+3*lam**2*r*r+9*lam*r+9)/(r**3*(lam*r+3)**2)
    assert s.factor(W*W-f*s.diff(W,r)-C*C-rw)==0
    assert s.factor(W*W+f*s.diff(W,r)-C*C-z)==0
    return {'stress_conservation':'exact symbolic zero (time, radial); angular components identically zero',
            'trace_and_fixed_energy_null_limit':'exact symbolic identities',
            'source_implementation_max_relative_error':max(errors),
            'GR_Zerilli_Regge_Wheeler_intertwining':'exact symbolic identity; applies to massless sector only'}


def compact_equation_checks():
    """Check every compact equation printed in the derivation against coordinates."""
    data=json.loads(DEFAULT_OPERATORS.read_text())['general_l_ge_2']
    r,w,mu,L=s.symbols('r omega mu L');F=1-2/r
    a,b,c,d,e,K,G=s.symbols('a_0 b_0 c_0 d_0 e_0 K_0 G_0')
    ap,bp,cp,dp,ep,Kp,Gp=s.symbols('a_1 b_1 c_1 d_1 e_1 K_1 G_1')
    app,bpp,cpp,dpp,epp,Kpp,Gpp=s.symbols('a_2 b_2 c_2 d_2 e_2 K_2 G_2')
    B=w*w/F-mu*mu-L/r**2
    equations=[
        F*F*app+2*F*(r-1)/r**2*ap+(F*B-2/r**4)*a-4*s.I*w*b/r**2-2*(2*r-3)*c/r**4+4*F*K/r**3,
        F*bpp+2*(r-1)*bp/r**2+(B-2/r**2-4/(r**4*F))*b-2*s.I*w*(a+c)/(r*r*F)+2*L*d/r**3,
        cpp+2*(r-1)*cp/(r*r*F)+(B/F-4/r**2-2/(r**4*F*F))*c-2*(2*r-3)*a/(r**4*F*F)-4*s.I*w*b/(r*r*F*F)+4*(r-3)*K/(r**3*F)+4*L*e/r**3,
        F*dpp+(B+4/r**3)*d+2*F*b/r-2*s.I*w*e/r**2,
        F*epp+4*ep/r**2+(B-4*F/r**2)*e+2*(c-K)/r+2*(L-2)*G/r-2*s.I*w*d/(r*r*F*F),
        F*Kpp+2*(r-1)*Kp/r**2+(B-2*(r-4)/r**3)*K+2*a/r**3+2*(r-3)*c/r**3-2*L*F*e/r**3,
        F*Gpp+2*(r-1)*Gp/r**2+(B+2/r**2)*G+2*F*e/r**3]
    actual=s.Matrix(data['matrix'])[:7,:]*s.Matrix(s.symbols(' '.join(data['columns'])))
    assert all(s.cancel(a-b)==0 for a,b in zip(actual,equations))
    P,V=monopole_coefficients(r,w,mu)
    pscaled=s.cancel((r-2)*P);vscaled=s.cancel((r-2)**2*V)
    assert s.cancel(pscaled.subs(r,2)-1)==0
    assert s.cancel(vscaled.subs(r,2)-4*w*w)==0
    rho=-2*s.I*w
    b1=-(rho*s.diff(pscaled,r).subs(r,2)+s.diff(vscaled,r).subs(r,2))/(2*rho+1)
    assert s.cancel(b1-horizon_regular_coefficient(w,mu))==0
    return 'all seven compact equations agree exactly with the coordinate-derived operator'


def local_checks():
    records=[]
    for ell in (0,1,2,3,4,7):
        system=PolarSystem(ell)
        for rv in (2.3,4.,11.):
            omega=.31-.06j;mu=.4
            for source in ('vacuum','null','massive'):
                def mapping(x):
                    rhs=None
                    if source!='vacuum':
                        time=quad(lambda z:-1/(1-2/z) if source=='null' else trajectory_derivative(z,2.),4.,x,epsabs=1e-12)[0]
                        rhs=forcing(ell,x,omega,mu,source_data(x,omega,energy=.001,gamma=2.,time=time,null=source=='null'))
                    B,b,info=system.reconstruction(x,omega,mu,rhs)
                    return B,b,info
                B,b,info=mapping(rv);n=system.n;d=system.d
                A=np.vstack((np.column_stack((np.zeros((d,d)),np.eye(d))),B[system.acceleration]))
                z=np.r_[np.zeros(d),b[system.acceleration]]
                h=2e-4
                bm2,zm2,_=mapping(rv-2*h);bm1,zm1,_=mapping(rv-h)
                bp1,zp1,_=mapping(rv+h);bp2,zp2,_=mapping(rv+2*h)
                Bp=(bm2-8*bm1+8*bp1-bp2)/(12*h)
                zp=(zm2-8*zm1+8*zp1-zp2)/(12*h)
                # Derivatives of reconstructed fields AND first derivatives
                # must agree with the next jet, not merely with algebraic rows.
                err=np.linalg.norm(Bp[:2*n]+B[:2*n]@A-B[n:])/max(1.,np.linalg.norm(B[n:]))
                errb=np.linalg.norm(zp[:2*n]+B[:2*n]@z-b[n:])/max(1.,np.linalg.norm(b[n:]))
                assert err<2e-6 and errb<2e-6,(ell,rv,source,err,errb)
                item={'ell':ell,'r':rv,'source':source,**info,
                      'jet_derivative_relative_residual':float(err),'source_jet_derivative_scaled_residual':float(errb)}
                if ell==0:
                    P,V=monopole_coefficients(rv,omega,mu)
                    error=np.linalg.norm(A[-1]-[-V,-P])/max(1.,abs(P),abs(V))
                    assert error<1e-10
                    item['AGP_monopole_coefficient_relative_error']=float(error)
                if ell>=2:
                    f=1-2/rv;L=ell*(ell+1)
                    expected=np.array([0.,-2/rv**3,-(omega**2-f*(mu**2+(L-2)/rv**2))/f**2,0.,0.,-2*(rv-1)/(rv**2*f)])
                    error=np.linalg.norm(A[-1]-expected)/max(1.,np.linalg.norm(expected))
                    assert error<1e-10
                    item['AGP_G_row_relative_error']=float(error)
                records.append(item)
    return records


def spectrum_checks():
    records=[]
    for mu,bracket in ((.3,(.025,.045)),(.4,(.007,.017))):
        modes=[unstable_mode(mu,bracket,epsilon=eps,outer=outer,match=match,rtol=tol)
               for eps,outer,match,tol in ((1e-5,60.,4.,1e-9),(1e-6,80.,4.,2e-10),(1e-7,100.,5.,5e-11))]
        fd=[lowest_eigenvalue(mu,left=-900.,right=150.,spacing=h) for h in (.1,.05,.025)]
        boundary_check=lowest_eigenvalue(mu,left=-700.,right=150.,spacing=.025)
        reference=(4*fd[-1]['growth_rate']-fd[-2]['growth_rate'])/3
        growth=modes[-1]['omega_rg'][1]
        difference=abs(growth-reference)
        assert abs(growth-modes[-2]['omega_rg'][1])<2e-8
        assert difference<2e-8
        records.append({'mu':mu,'family':'polar monopole unstable bound state',
                        'literature_anchor':'AGP2026 Appendix B.2 Fig.7; BCP2013 Eq.(30)',
                        'reference_type':'independent published-potential computation, not digitized author frequencies',
                        'K_equation_matching':modes,'independent_master_potential_finite_difference':fd,
                        'finite_difference_inner_boundary_check':boundary_check,
                        'Richardson_growth_reference':reference,'absolute_frequency_difference':difference})
    return records


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result={'gate':'W3','units':'GM=1, exp(-i omega t)',
            'operator_sha256':hashlib.sha256(DEFAULT_OPERATORS.read_bytes()).hexdigest(),
            'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'sympy':s.__version__}}
    print('Checking exact source and GR identities',flush=True);result['source']=source_checks()
    result['compact_equations']=compact_equation_checks()
    print('Checking constrained reduction and sourced jet compatibility',flush=True);result['local_checks']=local_checks()
    print('Checking independent polar-monopole frequencies',flush=True);result['spectral_checks']=spectrum_checks()
    result['passed']=True
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print('All W3 checks passed',flush=True)
