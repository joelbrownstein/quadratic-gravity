"""Reproduce the characteristic-front algebra and spectral convergence ledger."""
import hashlib
import json
from pathlib import Path
import platform
import sys
import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad
import sympy as sp
from w4_front import analytic_front, leading_front, front_transform, phase_removed, front_waveform
from w4_monopole import tortoise

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'notes/validation'


def algebra():
    # Recover characteristic transport from the independently generated W3 tt row.
    record=json.loads((DATA/'w3_polar_operators.json').read_text())['systems'][0]
    r,w,mu=sp.symbols('r omega mu');F=1-2/r
    A,Ap,App=sp.symbols('A Ap App')
    eq=(sp.Matrix(record['matrix'])*sp.Matrix(sp.symbols(' '.join(record['columns']))))[0]
    subs={}
    for name in ('a','b','c'):
        subs[sp.Symbol(name+'_0')]=A
        subs[sp.Symbol(name+'_1')]=Ap-sp.I*w*A/F
        subs[sp.Symbol(name+'_2')]=App-2*sp.I*w*Ap/F+sp.I*w*sp.diff(F,r)*A/F**2-w*w*A/F**2
    subs.update({sp.Symbol('K_'+str(j)):0 for j in range(3)})
    transport=sp.expand(eq.subs(subs)).coeff(w)/(-sp.I)
    assert sp.simplify(transport-(2*F*Ap+2*A/r))==0
    C=sp.symbols('C')
    A0=1/(2*(r-2));B0=(sp.log(r)/2+C)/(r-2)
    assert sp.simplify(2*F*sp.diff(A0,r)+2*A0/r)==0
    assert sp.simplify(2*F*sp.diff(B0,r)+2*B0/r-1/r**2)==0
    H=sp.diff(F,r)+2*F/r
    c1=sp.simplify(F*sp.diff(A0,r)+H*A0)
    d1=sp.simplify(F*sp.diff(B0,r)+H*B0-c1)
    assert sp.simplify(c1-1/(2*r*r))==0
    assert sp.simplify(d1-B0*(r-2)/r**2)==0
    # Flat-space radial Klein-Gordon equation away from the null source.
    t,rho=sp.symbols('t rho',positive=True);fun=sp.Function('f')
    z=sp.sqrt(r*r-t*t)
    a=-fun(z)/r
    flat=sp.diff(a,r,2)+2*sp.diff(a,r)/r-sp.diff(a,t,2)-mu**2*a
    expected=-(sp.Subs(sp.diff(fun(rho),rho,2),rho,z)+
               sp.Subs(sp.diff(fun(rho),rho),rho,z)/z-mu**2*fun(z))/r
    assert sp.simplify(flat-expected)==0
    b=-t*a/r;c=t*t*a/r**2;K=(a-c)/2
    wave=lambda y:sp.diff(y,r,2)+2*sp.diff(y,r)/r-sp.diff(y,t,2)-mu**2*y
    ode=-r*expected
    for eq,field in ((wave(b)-2*b/r**2,b),
                     (wave(c)-4*(c-K)/r**2,c),
                     (wave(K)+2*(c-K)/r**2,K)):
        assert sp.simplify(eq-field/fun(z)*ode)==0
    assert sp.simplify(-sp.diff(a,t)+sp.diff(b,r)+2*b/r)==0
    assert sp.simplify(-sp.diff(b,t)+sp.diff(c,r)+2*(c-K)/r)==0
    return {'tt_characteristic_transport':'2 F A_prime + 2 A/r',
            'log_coefficient_A':'1/[2(r-2)]',
            'K_vlogv_coefficient':'1/(2r^2)',
            'K_v_coefficient':'[log(mu/sqrt(2))+EulerGamma+log(r)/2]/r^2',
            'flat_all_tensor_wave_equations_and_divergences':True,'passed':True}


def transforms():
    mp.mp.dps=70
    errors=[]
    for u in (0,.1,.4,.42,1,12,100,10000):
        w=u+.004j;s=-1j*mp.mpc(w);z=25/s
        ref=-((1+z)*mp.exp(z)*mp.e1(z)-1)/(80000*s*s)
        got=complex(front_transform(w))
        errors.append(float(abs(mp.mpc(got)-ref)/abs(ref)))
    quadrature=[]
    for s in (.1,1,10):
        val,err=quad(lambda v:float(analytic_front(v))*np.exp(-s*v),0,np.inf,
                     epsabs=1e-19,epsrel=1e-10)
        quadrature.append(float(abs(val-front_transform(1j*s))/abs(val)))
    v=np.geomspace(1e-12,1e-8,20)
    asymptotic=float(np.max(abs((analytic_front(v)-leading_front(v))/analytic_front(v))))
    assert max(errors)<1e-9
    assert max(quadrature)<1e-8
    assert asymptotic<1e-6
    return {'max_transform_relative_error_70_digit_reference':max(errors),
            'max_transform_relative_error_direct_Laplace_quadrature':max(quadrature),
            'leading_front_max_relative_error_v_1e_12_to_1e_8':asymptotic}


def load_low(name):
    meta=json.loads((DATA/(name+'.json')).read_text())
    digest=hashlib.sha256((ROOT/'scripts/weyl_ringdown/w4_threshold.py').read_bytes()).hexdigest()
    assert meta['solver_sha256']==digest
    return np.loadtxt(DATA/(name+'.csv'),delimiter=',',skiprows=1)


def radial_checks():
    from w4_threshold import response
    rows=[]
    for u in (2.,6.,12.):
        w=u+.004j;radius=max(800.,200*u)
        base=response(w,outer=radius,rtol=2e-11)['K']
        refined=response(w,outer=2*radius,rtol=4e-12,order=12,tail_order=20)['K']
        phase=np.exp(1j*w*tortoise(200.))
        rows.append({'u':u,'outer':radius,'refined_outer':2*radius,
                     'relative_difference':float(abs(base-refined)/abs(refined)),
                     'absolute_phase_removed_difference':float(abs((base-refined)*phase)),
                     'B_real':float((refined*phase).real),'B_imag':float((refined*phase).imag),
                     'reference_residual':float(abs(refined*phase-front_transform(w)))})
        print(rows[-1],flush=True)
    (DATA/'w4_front_radial_checks.json').write_text(json.dumps(rows,indent=2)+'\n')


def waveform():
    coarse=load_low('w4_contour_sigma004')
    middle=load_low('w4_contour_sigma004_midpoints')
    other=load_low('w4_contour_sigma003')
    fine=np.vstack((coarse,middle));fine=fine[np.argsort(fine[:,0])]
    high={}
    for name in ('004','003'):
        meta=json.loads((DATA/('w4_front_high_sigma'+name+'.json')).read_text())
        assert meta['threshold_solver_sha256']==hashlib.sha256((ROOT/'scripts/weyl_ringdown/w4_threshold.py').read_bytes()).hexdigest()
        high[name]=np.loadtxt(DATA/('w4_front_high_sigma'+name+'.csv'),delimiter=',',skiprows=1)
    def combined(low,hi,sigma):
        hi=hi[hi[:,0]>2+1e-10]
        return np.r_[low[:,0],hi[:,0]],np.r_[phase_removed(low,sigma),hi[:,3]+1j*hi[:,4]]
    u,B=combined(fine,high['004'],.004)
    v=np.unique(np.r_[np.linspace(-.1,0,101),np.geomspace(1e-9,1,901)])
    best=front_waveform(u,B,v)
    alternatives={}
    for cut in (2,4,8):
        alternatives['cutoff_'+str(cut)]=front_waveform(u,B,v,cutoff=cut)
    uc,Bc=combined(coarse,high['004'][::2],.004)
    alternatives['coarser_grid']=front_waveform(uc,Bc,v)
    uo,Bo=combined(other,high['003'],.003)
    alternatives['contour_sigma003']=front_waveform(uo,Bo,v,sigma=.003)
    peak=float(np.max(abs(best['K'])))
    comparisons={}
    for name,result in alternatives.items():
        comparisons[name]={'max_absolute':float(max(abs(result['K']-best['K']))),
                           'max_relative_to_peak':float(max(abs(result['K']-best['K']))/peak)}
    # Fixed target for this local window: 0.1% of the resolved front peak.
    # Cutoff 2 is a diagnostic, not an acceptance comparison.
    for name in ('cutoff_4','cutoff_8','coarser_grid','contour_sigma003'):
        assert comparisons[name]['max_relative_to_peak']<1e-3,(name,comparisons[name])
    radial=json.loads((DATA/'w4_front_radial_checks.json').read_text())
    assert max(row['relative_difference'] for row in radial)<1e-5
    np.savetxt(DATA/'w4_null_front.csv',
               np.column_stack((v,v-tortoise(200.),best['K'],analytic_front(v),
                                best['remainder'],best['raw_K'],
                                *(o['K'] for o in alternatives.values()))),
               delimiter=',',comments='',
               header='v,t,K_front,analytic_reference,corrected_remainder,unprojected_K,'+','.join(alternatives))
    return {'mu':.5,'extraction':200.,'v_window':[-.1,1.],
            'frequency_maximum':float(max(u)),'sigma':[.004,.003],
            'peak_abs_K_per_kappa_E_Y00':peak,'v_at_peak':float(v[np.argmax(abs(best['K']))]),
            'comparison_tolerance_relative_to_peak':1e-3,
            'comparisons':comparisons,
            'unprojected_max_prefront_leakage':float(max(abs(best['raw_K'][v<=0]))),
            'zeroth_residual_moment':best['zeroth_moment'],
            'first_residual_moment':best['first_moment'],
            'max_corrected_remainder_postfront':float(max(abs(best['remainder'][v>0]))),
            'causality_imposed_from_characteristics_not_an_independent_test':True,
            'radial_checks':radial,
            'scope':'analytic leading front plus convergence-tested moment-corrected local auxiliary waveform; no metric strain or flux',
            'full_W4_complete':False}


if __name__=='__main__':
    if '--radial-only' in sys.argv:
        radial_checks()
        sys.exit(0)
    result={'stage':'W4 null-particle-passage front','algebra':algebra(),'transforms':transforms()}
    if '--algebra-only' not in sys.argv:
        result['waveform']=waveform()
        result['versions']={'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__}
        result['source_sha256']={name:hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()
                                 for name in ('w4_front.py','compute_w4_front.py','validate_w4_front.py','w4_threshold.py')}
        result['checks_passed']=True
        (DATA/'w4_front_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
