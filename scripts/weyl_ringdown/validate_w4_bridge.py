"""Validate the interval connecting the front and late-time W4 inversions."""
import hashlib
import json
from pathlib import Path
import platform
import numpy as np
import scipy
from w4_bridge import bridge_waveform
from w4_front import analytic_front,front_transform,phase_removed,front_waveform
from w4_inversion import postpassage_waveform
from w4_monopole import tortoise

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'notes/validation'
INPUTS=[]


def load(name,high=False):
    meta=json.loads((DATA/(name+'.json')).read_text())
    digest=hashlib.sha256((Path(__file__).parent/'w4_threshold.py').read_bytes()).hexdigest()
    assert meta['threshold_solver_sha256' if high else 'solver_sha256']==digest
    INPUTS.extend((name+'.csv',name+'.json'))
    return np.loadtxt(DATA/(name+'.csv'),delimiter=',',skiprows=1)


def merged(a,b):
    out=np.vstack((a,b))
    return out[np.argsort(out[:,0])]


def combined(low,high,sigma):
    high=high[high[:,0]>2+1e-10]
    return np.r_[low[:,0],high[:,0]],np.r_[phase_removed(low,sigma),high[:,3]+1j*high[:,4]]


def run():
    coarse=load('w4_contour_sigma004')
    fine=merged(coarse,load('w4_contour_sigma004_midpoints'))
    other_coarse=load('w4_contour_sigma003')
    other=merged(other_coarse,load('w4_bridge_sigma003_midpoints'))
    high=load('w4_front_high_sigma004',True);other_high=load('w4_front_high_sigma003',True)
    u,B=combined(fine,high,.004)
    v=np.unique(np.r_[np.linspace(.5,1,101),np.arange(1,510.01,.25),
                       tortoise(200.)+np.arange(301.)])
    best=bridge_waveform(u,B,v)
    quiet=(v>=10)&(v<=400)
    comparisons={};alternatives={}
    def compare(name,result,required=True):
        delta=result['K']-best['K']
        comparisons[name]={'max_absolute':float(max(abs(delta))),
                           'quiet_interval_max_absolute':float(max(abs(delta[quiet]))),
                           'required':required}
        alternatives[name]=result['K']
        if required:
            assert max(abs(delta))<2e-12,(name,comparisons[name])
            assert max(abs(delta[quiet]))<2e-13,(name,comparisons[name])
    for cut in (4,8):
        compare('cutoff_'+str(cut),bridge_waveform(u,B,v,cutoff=cut))
    for terms in (1,2):
        compare('tail_terms_'+str(terms),bridge_waveform(u,B,v,tail_terms=terms))
    compare('endpoint_width_05',bridge_waveform(u,B,v,fit_width=.5))
    compare('endpoint_degree_3',bridge_waveform(u,B,v,degree=3))
    for order in (12,20):
        compare('gauss_order_'+str(order),bridge_waveform(u,B,v,quadrature_order=order))
    compare('spline_order_7',bridge_waveform(u,B,v,spline_order=7))
    # Document why a small Gauss rule is insufficient on the wide high-u cells.
    compare('gauss_order_8_diagnostic',bridge_waveform(u,B,v,quadrature_order=8),False)
    uc,Bc=combined(coarse,high[::2],.004)
    compare('coarser_spectral_grid',bridge_waveform(uc,Bc,v,quadrature_order=32))
    cubic_fine=bridge_waveform(u,B,v,spline_order=3)['K']
    cubic_coarse=bridge_waveform(uc,Bc,v,spline_order=3,quadrature_order=32)['K']
    cubic_error=float(max(abs(cubic_fine-cubic_coarse)))
    uo,Bo=combined(other,other_high,.003)
    compare('refined_contour_sigma003',bridge_waveform(uo,Bo,v,sigma=.003))
    uoc,Boc=combined(other_coarse,other_high,.003)
    compare('old_coarse_sigma003_diagnostic',bridge_waveform(uoc,Boc,v,sigma=.003),False)
    # Independent previous inversion uses the full K spectrum and a discrete
    # endpoint-tail formula, rather than spline integration of the residual.
    late_mask=v>=tortoise(200.)
    times=v[late_mask]-tortoise(200.)
    late=postpassage_waveform(fine[:,0],fine[:,1]+1j*fine[:,2],times,.004)[0]
    late_error=float(max(abs(late-best['K'][late_mask])))
    assert late_error<2e-12,late_error
    front_mask=v<=1
    front=front_waveform(u,B,v[front_mask])['K']
    front_error=float(max(abs(front-best['K'][front_mask])))
    front_peak=json.loads((DATA/'w4_front_validation.json').read_text())['waveform']['peak_abs_K_per_kappa_E_Y00']
    INPUTS.append('w4_front_validation.json')
    assert front_error<1e-3*front_peak
    # Independent manufactured transform: L[v^3 exp(-a v)/6]=(s+a)^-4.
    s=.004-1j*u
    test_B=front_transform(u+.004j)+1e-10/(s+.7)**4
    manufactured=bridge_waveform(u,test_B,v)['K']
    exact=analytic_front(v)+1e-10*v**3*np.exp(-.7*v)/6
    manufactured_error=float(max(abs(manufactured-exact)))
    assert manufactured_error<1e-13
    np.savetxt(DATA/'w4_bridge_waveform.csv',
               np.column_stack((v,v-tortoise(200.),best['K'],best['raw_K'],best['tail'],
                                *(alternatives[k] for k in alternatives))),
               delimiter=',',comments='',
               header='v,t,K_bridge,K_without_residual_tail,residual_tail,'+','.join(alternatives))
    np.savetxt(DATA/'w4_bridge_front_overlap.csv',
               np.c_[v[front_mask],best['K'][front_mask],front,best['K'][front_mask]-front],
               delimiter=',',comments='',header='v,K_bridge,K_previous_front,difference')
    np.savetxt(DATA/'w4_bridge_late_overlap.csv',
               np.c_[times,best['K'][late_mask],late,best['K'][late_mask]-late],
               delimiter=',',comments='',header='t,K_bridge,K_previous_late,difference')
    result={'stage':'W4 connecting time interval','mu':.5,'extraction':200.,
            'amplitude_units':'per kappa E Y00*, auxiliary angular K',
            'validated_v_window':[.5,510.],'validated_t_window':[.5-tortoise(200.),510-tortoise(200.)],
            'quiet_diagnostic_v_window':[10.,400.],
            'absolute_comparison_target':2e-12,'quiet_absolute_comparison_target':2e-13,
            'comparisons':comparisons,
            'front_overlap':{'v_window':[.5,1.],'max_absolute':front_error,
                             'relative_to_previous_front_peak':front_error/front_peak,
                             'acceptance':'previous front tolerance, 1e-3 of front peak'},
            'late_overlap':{'t_window':[float(min(times)),float(max(times))],
                            'max_absolute':late_error,'absolute_target':2e-12},
            'manufactured_inverse_max_absolute_error':manufactured_error,
            'rejected_initial_cubic_grid_check':{
                'max_absolute':cubic_error,'target':2e-12,
                'passed':cubic_error<2e-12,
                'resolution':'quintic interpolation; independent seventh-order and spectral-grid checks'},
            'quiet_max_abs_K':float(max(abs(best['K'][quiet]))),
            'quiet_max_required_comparison':max(x['quiet_interval_max_absolute'] for x in comparisons.values() if x['required']),
            'quiet_interpretation':'values below the empirical comparison floor are not evidence of exact vanishing',
            'max_abs_residual_tail':float(max(abs(best['tail']))),
            'full_W4_complete':False,'checks_passed':True,
            'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
            'source_sha256':{name:hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()
                             for name in ('w4_bridge.py','validate_w4_bridge.py','w4_threshold.py','w4_front.py','w4_inversion.py')},
            'input_sha256':{name:hashlib.sha256((DATA/name).read_bytes()).hexdigest() for name in sorted(set(INPUTS))}}
    (DATA/'w4_bridge_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    run()
