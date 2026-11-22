"""Independent metric geometry, sourced reconstruction and waveform checks."""
import json
import hashlib
from pathlib import Path
import platform
import numpy as np
import scipy
import sympy as s
from w4_metric import expressions,reconstruct
from w3_polar import PolarSystem
from w3_infall_source import forcing,source_data
from w4_monopole import tortoise
from w4_inversion import postpassage_waveform

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'notes/validation'


def algebra_and_geometry():
    args,e,d=expressions();r,w,mu,K,J,Z=args;F=1-2/r
    mass=e['mass_q']+Z/(-2*s.I*w)
    constraints={'div_time':s.I*w*e['a']+F*d(e['b'])+(s.diff(F,r)+2*F/r)*e['b'],
                 'div_radial':s.I*w*e['b']/F+d(e['c'])+s.diff(F,r)*(e['a']+e['c'])/(2*F)+2*(e['c']-K)/r,
                 'mass_radial':d(mass)+mu*mu*r*r*e['a']/2}
    for label,eq in constraints.items():
        assert all(s.simplify(s.factor(s.diff(eq,z)))==0 for z in (K,J,Z)),label
    fields=[e[n] for n in ('a','b','c','K')]
    qjets=s.lambdify(args,fields+[d(x) for x in fields]+[d(d(x)) for x in fields],'numpy',cse=True)
    step=Z/(-s.I*w)
    h=[step/r+2*F*e['a'],step/(r*F)+2*e['b'],
       step/(r*F**2)+2*e['c']/F,2*K]
    hjets=s.lambdify(args,h+[d(x) for x in h]+[d(d(x)) for x in h],'numpy',cse=True)
    record=json.loads((DATA/'w4_metric_geometry.json').read_text())
    geom=s.lambdify((r,w),s.Matrix(record['matrix']),'numpy',cse=True)
    polar=PolarSystem(0)
    operator_errors=[];geometry_errors=[];independent_errors=[]
    for rv in (2.7,8.,30.,200.):
        for omega in (.3+.004j,.5+.004j,.8+.003j,1.5+.004j):
            f=1-2/rv;phase=np.exp(-1j*omega*tortoise(rv))
            for k,kp in ((.1+.2j,-.03+.04j),(0j,0j)):
                params=(rv,omega,.5,k,kp,phase)
                qj=np.asarray(qjets(*params),complex)
                src=forcing(0,rv,omega,.5,source_data(rv,omega,energy=1.,time=-tortoise(rv),null=True))
                rhs=np.array([src[n] for n in polar.record['rows']])
                A=np.asarray(polar.evaluate(rv,omega,.5),complex)
                error=A@qj-rhs
                scales=np.sum(abs(A*qj[None,:]),axis=1)+abs(rhs)+1e-20
                operator_errors.append(float(max(abs(error)/scales)))
                hh=np.asarray(hjets(*params),complex)
                actual=geom(rv,omega)@hh
                out=reconstruct(rv,omega,.5,k,f*kp)
                expected=np.array([-.25*f*out['a'],-.25*out['b'],-.25*out['c']/f,
                                   -.25*k,0,out['mass_total'],out['tidal_radial_total'],out['tidal_transverse_total']])
                geometry_errors.append(float(max(abs(actual-expected)/np.maximum(1.,abs(expected)))))
                if rv<=30:
                    B,b,info=polar.reconstruction(rv,omega,.5,src)
                    independent=B@np.array([k,kp])+b
                    independent_errors.append(float(max(abs(independent-qj))/max(1.,max(abs(qj)))))
    assert max(operator_errors)<1e-9,max(operator_errors)
    assert max(geometry_errors)<1e-8,max(geometry_errors)
    assert max(independent_errors)<1e-8,max(independent_errors)
    # General spherical Riemann decomposition, including radial boosts.
    a,b,c,k,m,beta=s.symbols('a b c k m beta')
    Ru=-mu**2*(a+2*beta*b+beta**2*c)/(1-beta**2)
    Ra=-mu**2*k
    Er=2*Ra-2*m/r**3;Et=Ru/2-Ra+m/r**3
    assert s.simplify(Er+2*Et-Ru)==0
    assert s.simplify((Er-Et)- (2*Ra-2*m/r**3-(Ru/2-Ra+m/r**3)))==0
    return {'exact_divergence_and_mass_constraints':True,
            'coordinate_geometry_checks':record['checks'],
            'full_sourced_W3_rowwise_relative_residual':max(operator_errors),
            'independent_full_metric_curvature_max_scaled_error':max(geometry_errors),
            'independent_overdetermined_reconstruction_max_scaled_error':max(independent_errors),
            'radial_boost_tidal_trace_identity':True,
            'sample_radii':[2.7,8.,30.,200.],'sample_complex_frequencies':[[.3,.004],[.5,.004],[.8,.003],[1.5,.004]]}


def waveforms():
    inputs=[]
    digest=hashlib.sha256((Path(__file__).parent/'w4_threshold.py').read_bytes()).hexdigest()
    def load(name):
        meta=json.loads((DATA/(name+'.json')).read_text());assert meta['solver_sha256']==digest
        inputs.extend((name+'.csv',name+'.json'))
        return np.loadtxt(DATA/(name+'.csv'),delimiter=',',skiprows=1)
    def merge(a,b):
        z=np.vstack((a,b));return z[np.argsort(z[:,0])]
    coarse=load('w4_contour_sigma004')
    fine=merge(coarse,load('w4_contour_sigma004_midpoints'))
    other=merge(load('w4_contour_sigma003'),load('w4_bridge_sigma003_midpoints'))
    r=200.;F=1-2/r;times=np.arange(901.)
    names=('a','b','c','K','mass_q','Ricci_uu','Ricci_angular',
           'Ricci_null_out','Ricci_null_in','tidal_radial_q','tidal_transverse_q')
    def invert(data,sigma,**options):
        u=data[:,0];R=reconstruct(r,u+1j*sigma,.5,data[:,1]+1j*data[:,2],data[:,3]+1j*data[:,4])
        return {key:postpassage_waveform(u,R[key],times,sigma,**options)[0] for key in names}
    best=invert(fine,.004)
    alt={'grid_002':invert(coarse,.004),'contour_003':invert(other,.003),
         'cutoff_15':invert(fine,.004,cutoff=1.5),
         'two_tail_terms':invert(fine,.004,terms=2),
         'endpoint_width_08':invert(fine,.004,fit_width=.08,degree=4)}
    checks={}
    for key in names:
        peak=float(max(abs(best[key])))
        checks[key]={'peak_abs':peak,'comparisons':{}}
        for name,result in alt.items():
            err=float(max(abs(result[key]-best[key])))
            checks[key]['comparisons'][name]={'max_absolute':err,'relative_to_peak':err/peak}
            assert err<1e-3*peak,(key,name,err,peak)
    # Add the Einstein-sector step analytically, not as a truncated pole integral.
    columns={'t':times,'v':times+tortoise(r),**best,
             'h_tt_q':2*F*best['a'],'h_tr_q':2*best['b'],'h_rr_q':2*best['c']/F,
             'h_tt_total':1/r+2*F*best['a'],
             'h_tr_total':1/(r*F)+2*best['b'],
             'h_rr_total':1/(r*F**2)+2*best['c']/F,
             'h_angular_fraction_total':2*best['K'],'mass_total':.5+best['mass_q'],
             'tidal_radial_total':-1/r**3+best['tidal_radial_q'],
             'tidal_transverse_total':1/(2*r**3)+best['tidal_transverse_q']}
    trace_error=float(max(abs(columns['tidal_radial_total']+2*columns['tidal_transverse_total']-best['Ricci_uu'])))
    assert trace_error<1e-17
    np.savetxt(DATA/'w4_metric_waveform.csv',np.column_stack(list(columns.values())),
               delimiter=',',comments='',header=','.join(columns))
    # Carry the angular Ricci observable through the already validated bridge.
    bridge=np.genfromtxt(DATA/'w4_bridge_waveform.csv',delimiter=',',names=True)
    np.savetxt(DATA/'w4_metric_bridge_angular_curvature.csv',
               np.c_[bridge['v'],bridge['t'],-.25*bridge['K_bridge']],
               delimiter=',',comments='',header='v,t,Ricci_angular')
    inputs.extend(('w4_bridge_waveform.csv','w4_bridge_validation.json'))
    return {'t_window':[0.,900.],'mu':.5,'extraction':r,
            'comparison_target_relative_to_each_massive_quantity_peak':1e-3,
            'massless_sector_added_analytically':True,'quantities':checks,
            'tidal_trace_max_absolute_error':trace_error,
            'angular_Ricci_bridge':'-mu^2 K; inherited bridge errors scale by mu^2',
            'limits':'No full-component metric front or complete bridge inversion; no detector transfer function or flux',
            'input_sha256':{n:hashlib.sha256((DATA/n).read_bytes()).hexdigest() for n in inputs}}


if __name__=='__main__':
    result={'stage':'W4 monopole metric and observer curvature',
            'algebra_and_geometry':algebra_and_geometry(),'waveforms':waveforms(),
            'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'sympy':s.__version__},
            'checks_passed':True,'full_W4_complete':False}
    result['source_sha256']={n:hashlib.sha256((Path(__file__).parent/n).read_bytes()).hexdigest()
                            for n in ('w4_metric.py','derive_w4_metric_geometry.py','validate_w4_metric.py','w4_inversion.py')}
    result['geometry_sha256']=hashlib.sha256((DATA/'w4_metric_geometry.json').read_bytes()).hexdigest()
    (DATA/'w4_metric_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
