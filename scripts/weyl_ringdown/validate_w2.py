#!/usr/bin/env python3
"""Reproduce selected l=2, n=0 branches; see W2 note for scope/limitations.

Run in the pinned environment in requirements-w2.txt. Writes JSON only when
--output is given. Progress goes to stderr. No plots are digitized or invented.
"""
import argparse
import json
import platform
import sys
from pathlib import Path
import numpy as np
import scipy
from w2_massive import cf_root,cf_matrix,matching_root,matching_matrix,resonance_k,complex_root
from w2_stability import lowest_eigenvalue,critical_mass


def pair(w): return [float(w.real),float(w.imag)]


def literature_matrix(w,mu,depth=800):
    """AGP 2025 Appendix C Eqs. (88)-(91), q=-ik, GM=1, l=2.

    This is a published-coefficient consistency check, NOT the independent
    differential-equation integration method.
    """
    q=-1j*resonance_k(w,mu);I=np.eye(2,dtype=complex);L=6
    def abc(n):
        a=(n+1)*(n+1-4j*w)
        b=2-L-2*(n*n+n-1)+w*w*(2*n-4j*w+1)/q-3*q*(2*n-4j*w+1)+4j*(2*n+1)*w-4*q*q+12*w*w
        c=(q*q*(n*n-4j*n*w-6*w*w-9)+q**4+w**4+2*q**3*(n-2j*w)-2*q*w*w*(n-2j*w))/(q*q)
        return a*I,np.array([[b,L-2],[-2,b-3]]),np.array([[c,6-3*L],[0,c+9]])
    R=np.zeros((2,2),complex)
    for n in range(depth-1,-1,-1):
        a,b,c=abc(n+1);R=-np.linalg.solve(b+a@R,c)
    a,b,_=abc(0)
    return b+a@R


def high_precision_cf(mu,guess,depth=1600):
    import mpmath as mp
    mp.mp.dps=50
    mu=mp.mpf(str(mu));I=mp.eye(2)
    A=mp.matrix([[10,8],[2,4]]);B=mp.matrix([[-16,-24],[0,2]])
    def determinant(w):
        k=mp.sqrt(w*w-mu*mu)
        if mp.re(k)<0:k=-k
        b0=1-4j*w;b1=6j*k+8j*w-4+2j*w*w/k;b2=-2j*k-4j*w+3-2j*w*w/k
        c0=4*k*k+12*k*w+3j*k+12*w*w+4j*w+4*w**3/k+1j*w*w/k
        c1=-k*k-4*k*w-2j*k-6*w*w-4j*w-4*w**3/k-2j*w*w/k-w**4/k**2
        R=mp.zeros(2)
        for n in range(depth,0,-1):
            a=(n+1)*(n+b0)*I
            b=(-2*n*(n-1)+b1*n+c0)*I-A-B/2
            c=((n-1)*(n-2)+b2*(n-1)+c1)*I+B/2
            R=-mp.inverse(b+a*R)*c
        return mp.det(c0*I-A-B/2+b0*R)
    seed=mp.mpc(str(guess.real),str(guess.imag))
    w=mp.findroot(determinant,(seed,seed+mp.mpf('1e-7')),tol=mp.mpf('1e-38'),maxsteps=20)
    return {'precision_digits':50,'depth':depth,'real':mp.nstr(w.real,40),'imag':mp.nstr(w.imag,40),
            'residual':mp.nstr(abs(determinant(w)),6),'difference_from_double':abs(complex(w)-guess)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    # Fixed before production: points on AGP Figs. 3 and 4, p=.8,.9,1.0.
    # There are no tabulated exact values supplied by those figures.
    points=(.4,.45,.5)
    configs=[dict(rm=5.,radius=35.,angle=1.4,epsilon=.04,order=8,rtol=1e-8),
             dict(rm=5.,radius=55.,angle=1.4,epsilon=.02,order=12,rtol=2e-10),
             dict(rm=6.,radius=70.,angle=1.45,epsilon=.01,order=16,rtol=2e-12)]
    rows=[];tracks=[]
    for label,seed,figure in [('tensor',.37367-.08896j,4),('vector',.4576-.095j,3)]:
        w=seed;track=[]
        for step in range(21):
            mu=step*.025
            w,res,_=cf_root(mu,w,depth=800)
            track.append({'mu':mu,'p':2*mu,'omega_rg':pair(w),'residual':res})
        tracks.append({'family':label,'method':'stepwise continuation from the mu=0 limiting spectrum',
                       'massless_endpoint_is_not_finite_coupling_QQG':True,'points':track})
        for mu in points:
            seed=complex(*track[round(mu/.025)]['omega_rg'])
            conv=[]
            for depth in (100,200,400,800,1600):
                w,res,_=cf_root(mu,seed,depth=depth)
                conv.append({'depth':depth,'omega_rg':pair(w),'determinant_residual':res})
            reference=complex(*conv[-1]['omega_rg'])
            di=[]
            for config in configs:
                wd,res,nfev=matching_root(mu,reference,**config)
                singular=np.linalg.svd(matching_matrix(wd,mu,**config),compute_uv=False)
                di.append({'settings':config,'omega_rg':pair(wd),'difference_from_CF':abs(wd-reference),
                           'relative_difference_from_CF':abs(wd-reference)/abs(reference),
                           'determinant_residual':res,'singular_values':singular.tolist(),'root_evaluations':nfev})
            # Independent root solve with perturbed starting frequency, same ODE method.
            repeat,_,_=matching_root(mu,reference+(.0005+.0003j),**configs[1])
            lit,lr,_=complex_root(lambda z:np.linalg.det(literature_matrix(z,mu,1600)),reference)
            tails=[]
            for tail in (.5,1.):
                wt,_,_=cf_root(mu,reference,depth=1600,tail=tail)
                tails.append({'tail':tail,'difference':abs(wt-reference)})
            change=abs(complex(*conv[-2]['omega_rg'])-reference)
            assert change/abs(reference)<1e-6
            assert max(t['difference'] for t in tails)/abs(reference)<1e-6
            assert di[-1]['relative_difference_from_CF']<1e-6
            assert abs(lit-reference)/abs(reference)<1e-7
            assert abs(repeat-complex(*di[1]['omega_rg']))/abs(reference)<1e-6
            assert reference.imag<0 and all(d['omega_rg'][1]<0 for d in di)
            assert di[-1]['singular_values'][-1] < 1e-7*di[-1]['singular_values'][-2]
            row={'family':label,'ell':2,'overtone':0,'mu':mu,'p':2*mu,
                 'literature_anchor':f'AGP 2025 arXiv:2412.15037v3 Fig. {figure}, p={2*mu:g}; Appendix C recurrence',
                 'reference_type':'published recurrence reproduced; not a tabulated or digitized author data point',
                 'CF_convergence':conv,'matching_convergence':di,'tail_checks':tails,
                 'published_recurrence_omega_rg':pair(lit),'published_recurrence_difference':abs(lit-reference),
                 'perturbed_matching_seed_result':pair(repeat),'CF_last_depth_change':change,
                 'radial_classification':'below radial stability threshold' if 2*mu<.87616 else 'above radial stability threshold'}
            rows.append(row)
            print(label,mu,reference,'DI relative error',di[-1]['relative_difference_from_CF'],file=sys.stderr,flush=True)
    # Reverse continuation checks that labels have not switched branches.
    reverse=[]
    for track in tracks:
        w=complex(*track['points'][-1]['omega_rg']);errors=[]
        for point in reversed(track['points'][:-1]):
            w,_,_=cf_root(point['mu'],w,depth=800)
            errors.append(abs(w-complex(*point['omega_rg'])))
        assert max(errors)<1e-8
        reverse.append({'family':track['family'],'max_forward_reverse_difference':max(errors)})
    print('Checking radial spectrum and threshold',file=sys.stderr,flush=True)
    radial=[]
    for mu in (.3,.4,.45,.5):
        for left,right,h in [(-300.,100.,.1),(-500.,150.,.1),(-500.,150.,.05),(-700.,200.,.025)]:
            radial.append(lowest_eigenvalue(mu,left,right,h))
    thresholds=[critical_mass(left=l,right=r,match=m,rtol=t) for l,r,m,t in
                [(-40.,60.,0.,1e-9),(-60.,100.,0.,1e-10),(-70.,150.,2.,1e-12)]]
    assert max(t['p_critical'] for t in thresholds)-min(t['p_critical'] for t in thresholds)<1e-6
    assert abs(thresholds[-1]['p_critical']-.876)<.001
    for mu in (.3,.4):
        assert all(r['omega_squared']<0 for r in radial if r['mu']==mu)
    for mu in (.45,.5):
        assert all(r['omega_squared']>0 for r in radial if r['mu']==mu)
    print('Checking 50-digit arithmetic at the hardest tensor point',file=sys.stderr,flush=True)
    hardest=next(r for r in rows if r['family']=='tensor' and r['mu']==.5)
    precision=high_precision_cf(.5,complex(*hardest['CF_convergence'][-1]['omega_rg']))
    assert precision['difference_from_double']<1e-10
    result={'gate':'W2','units':'GM=1; p=2 GM m2; exp(-i omega t)',
            'scope':'selected fundamental axial branches plus separate imported polar-monopole stability equation',
            'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
            'benchmarks':rows,'continuation':tracks,'reverse_continuation':reverse,
            'radial_finite_boxes':radial,'radial_zero_mode_threshold':thresholds,
            'high_precision_check':precision,'checks_passed':True,
            'not_claimed':['complete spectrum','all-mode or nonlinear stability','sourced waveform','quantum coupling bound'],
            'sources':['https://arxiv.org/abs/2412.15037v3','https://arxiv.org/abs/1304.6725v3']}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
