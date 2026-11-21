"""Validate W4 threshold boundaries and the extended post-passage waveform."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import numpy as np
import scipy
import sympy as s
from w3_polar import monopole_coefficients
from w4_monopole import response as original_response
from w4_threshold import response,threshold_slope
from w4_inversion import postpassage_waveform


def pair(z):return [float(z.real),float(z.imag)]


def algebra():
    z,w,mu=s.symbols('z omega mu')
    P,V=monopole_coefficients(1/z,w,mu)
    Q=V+z*z*s.diff(P,z)/2-P*P/4
    series=s.series(Q,z,0,3).removeO().expand()
    expected=[w*w-mu*mu,2*(2*w*w-mu*mu),12*w*w-4*mu*mu-2-4*w*w/(mu*mu)]
    assert all(s.simplify(series.coeff(z,j)-v)==0 for j,v in enumerate(expected))
    residuals=[]
    for r in (200.,400.,800.):
        def logu(x):
            p,_=monopole_coefficients(x,.5,.5)
            return threshold_slope(x,.5)/(1-2/x)+p/2
        h=.02
        Sp=(logu(r-2*h)-8*logu(r-h)+8*logu(r+h)-logu(r+2*h))/(12*h)
        p,v=monopole_coefficients(r,.5,.5)
        pp=(monopole_coefficients(r-2*h,.5,.5)[0]-8*monopole_coefficients(r-h,.5,.5)[0]
            +8*monopole_coefficients(r+h,.5,.5)[0]-monopole_coefficients(r+2*h,.5,.5)[0])/(12*h)
        q=v-pp/2-p*p/4
        err=float(abs(Sp+logu(r)**2+q)/abs(q));assert err<2e-8
        residuals.append({'r':r,'Riccati_relative_residual':err})
    return {'normal_form_coefficients':'exact symbolic identities through 1/r^2',
            'threshold_boundary_residuals':residuals}


def points():
    cases=[(.0,[400.,800.]),(.2,[400.,800.]),(.49,[400.,800.]),
           (.499,[800.,1600.]),(.4999,[3200.,6400.,12800.]),
           (.5,[300.,400.,800.]),(.501,[1600.,3200.,6400.]),
           (.5+.004j,[400.,800.]),(2+.004j,[400.,800.])]
    rows=[]
    for w,radii in cases:
        values=[response(w,outer=R) for R in radii]
        change=abs(values[-1]['K']-values[-2]['K'])/abs(values[-1]['K'])
        assert change<2e-5,(w,change)
        rows.append({'omega':pair(complex(w)),'radii':radii,'K':[pair(v['K']) for v in values],
                     'boundary':values[-1]['boundary'],'last_relative_change':float(change)})
        print('Boundary check',w,flush=True)
    comparisons=[]
    for w in (.2,.49,.5,.5+.004j,.7+.004j):
        a=response(w,outer=400.);b=response(w,outer=400.,method='bvp',rtol=2e-12)
        err=abs(a['K']-b['K'])/abs(a['K']);assert err<2e-5,(w,err)
        comparisons.append({'omega':pair(complex(w)),'relative_difference':float(err),
                             'green_K':pair(a['K']),'collocation_K':pair(b['K']),
                             'collocation_residual':b['max_collocation_residual'],'collocation_tolerance':2e-12,'mesh_nodes':b['mesh_nodes']})
    overlaps=[]
    for w in (.55,.7,1.):
        a=response(w,outer=400.)['K']
        b=original_response(w,.5,extraction=200.,outer=400.,order=10,tail_order=16,rtol=2e-10)['K']
        err=abs(a-b)/abs(a);assert err<2e-6
        overlaps.append({'omega':w,'relative_difference':float(err)})
    a=response(.5+.004j,chunk=20.)['K'];b=response(.5+.004j,chunk=80.)['K']
    scalecheck=float(abs(a-b)/abs(a));assert scalecheck<2e-7
    return {'boundary_refinements':rows,'independent_collocation':comparisons,
            'legacy_above_threshold_overlap':overlaps,'segment_rescaling_relative_difference':scalecheck}


def waveform(directory):
    def load(name):
        info=json.loads((directory/(name+'.json')).read_text())
        assert info['solver_sha256']==hashlib.sha256(Path(__file__).with_name('w4_threshold.py').read_bytes()).hexdigest()
        data=np.loadtxt(directory/(name+'.csv'),delimiter=',',skiprows=1)
        return data,info
    coarse,ci=load('w4_contour_sigma004');mid,mi=load('w4_contour_sigma004_midpoints')
    other,oi=load('w4_contour_sigma003')
    fine=np.vstack((coarse,mid));fine=fine[np.argsort(fine[:,0])]
    t=np.arange(0.,901.)
    def inv(data,sigma,**kwargs):
        return postpassage_waveform(data[:,0],data[:,1]+1j*data[:,2],t,sigma,**kwargs)
    best,raw,tail=inv(fine,.004)
    reference_norm=np.linalg.norm(best);peak=max(abs(best))
    alternatives={'frequency_step_002':inv(coarse,.004)[0],
                  'contour_sigma003':inv(other,.003)[0],
                  'cutoff_15':inv(fine,.004,cutoff=1.5)[0],
                  'cutoff_12':inv(fine,.004,cutoff=1.2)[0],
                  'two_tail_terms':inv(fine,.004,terms=2)[0],
                  'one_tail_term':inv(fine,.004,terms=1)[0],
                  'endpoint_fit_width008':inv(fine,.004,fit_width=.08,degree=4)[0]}
    errors={name:{'relative_L2':float(np.linalg.norm(y-best)/reference_norm),
                  'max_absolute':float(max(abs(y-best))),'max_relative_to_peak':float(max(abs(y-best))/peak)}
            for name,y in alternatives.items()}
    for name in ('frequency_step_002','contour_sigma003','cutoff_15','two_tail_terms','endpoint_fit_width008'):
        assert errors[name]['relative_L2']<1e-3,(name,errors[name])
        assert errors[name]['max_relative_to_peak']<2e-3,(name,errors[name])
    table=np.column_stack((t,best,raw,tail,*alternatives.values()))
    np.savetxt(directory/'w4_extended_waveform.csv',table,delimiter=',',comments='',
               header='t,K_postpassage,K_raw_cutoff2,tail_estimate,'+','.join(alternatives))
    # No extrapolation to the null front; retain the raw integral separately.
    from w4_threshold import invert_contour
    front_times=np.arange(-260.,1.)
    front=invert_contour(fine[:,0],fine[:,1]+1j*fine[:,2],front_times,.004)
    np.savetxt(directory/'w4_front_unresolved.csv',np.c_[front_times,front],delimiter=',',comments='',
               header='t,raw_finite_cutoff_K_NOT_converged_front')
    return {'mu':.5,'extraction':200.,'validated_time_window':[0.,900.],
            'frequency_interval':[0.,2.],'fine_step':.001,'coarse_step':.002,
            'contour_heights':[.004,.003],'high_frequency_tail_terms':3,
            'tail_model':'post-passage, phase-aware endpoint expansion matched to trapezoidal quadrature',
            'peak_abs_K_per_kappa_E_Y00':float(peak),'peak_time':float(t[np.argmax(abs(best))]),
            'comparisons':errors,'maximum_tail_correction':float(max(abs(tail))),
            'null_front_resolved':False,'full_W4_complete':False,
            'scope':'finite-radius auxiliary monopole; all frequency sectors included by causal contour, with an estimated high-frequency tail on the stated time window'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--directory',type=Path,required=True)
    parser.add_argument('--points-only',action='store_true');args=parser.parse_args()
    result={'stage':'W4 threshold and extended inversion','versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
            'threshold_solver_sha256':hashlib.sha256(Path(__file__).with_name('w4_threshold.py').read_bytes()).hexdigest(),
            'inversion_sha256':hashlib.sha256(Path(__file__).with_name('w4_inversion.py').read_bytes()).hexdigest()}
    result['algebra']=algebra();result['points']=points()
    if not args.points_only:result['waveform']=waveform(args.directory)
    result['checks_passed']=True
    name='w4_threshold_points.json' if args.points_only else 'w4_threshold_validation.json'
    (args.directory/name).write_text(json.dumps(result,indent=2)+'\n')
    print('Saved',name,flush=True)
