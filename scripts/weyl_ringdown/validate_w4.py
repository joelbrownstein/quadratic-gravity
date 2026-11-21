"""W4 initial benchmark; explicitly not certification of the entire gate.

Writes an above-threshold spectrum and band-limited waveform only after
method, boundary, source-tail, and Fourier-quadrature checks pass.
"""
import argparse
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
import scipy
from scipy.signal.windows import tukey
from w3_polar import monopole_coefficients
from w4_monopole import (response,coefficients,geometry,tortoise,outgoing_series,
                         source_tail_series,outer_data,band_waveform,signed_flux)
from verify_w4_algebra import verify


def pair(z):return [float(np.real(z)),float(np.imag(z))]


def series_checks():
    rows=[]
    for w in (.55,.7,1.):
        mu=.5
        for r in (200.,400.,800.):
            P,V=monopole_coefficients(r,w,mu);F=1-2/r
            k,p,a=outgoing_series(w,mu,8);n=np.arange(len(a));powers=r**(-n.astype(float))
            B=sum(a*powers);Bp=sum(-n*a*powers)/r;Bpp=sum(n*(n+1)*a*powers)/r**2
            log=1j*k+p/r
            out=Bpp+2*log*Bp+(log*log-p/r**2)*B+P*(Bp+log*B)+V*B
            b=source_tail_series(w,mu,14);n=np.arange(len(b));powers=r**(-n.astype(float))
            C=sum(b*powers);Cp=sum(-n*b*powers)/r;Cpp=sum(n*(n+1)*b*powers)/r**2
            log=-1j*w/F
            forced=Cpp+2*log*Cp+(log*log+2j*w/(r*r*F*F))*C+P*(Cp+log*C)+V*C
            D=r*(2*mu*mu*F+mu**4*r*r+4*w*w)
            S=2*(mu*mu*r+2j*w)/(r*r*F*D)
            row={'omega':w,'r':r,'homogeneous_scaled_residual':float(abs(out)/max(1.,abs(V*B))),
                 'forced_tail_relative_residual':float(abs(forced/S-1))}
            if r>=400:assert row['homogeneous_scaled_residual']<1e-9 and row['forced_tail_relative_residual']<1e-7,row
            rows.append(row)
    return rows


def point_checks():
    rows=[]
    for mu,w in ((.44,.484),(.5,.55),(.5,.7),(.5,1.),(.6,.84)):
        settings=[dict(outer=200.,order=6,tail_order=10,left=-25.,rtol=2e-10),
                  dict(outer=400.,order=8,tail_order=14,left=-35.,rtol=2e-11),
                  dict(outer=800.,order=10,tail_order=16,left=-40.,rtol=5e-12)]
        computed=[response(w,mu,**setting) for setting in settings]
        shoot=response(w,mu,method='shoot',**settings[1])
        best=computed[-1];mid=computed[-2]
        refinement=abs(mid['outgoing']-best['outgoing'])/abs(best['outgoing'])
        agreement=abs(mid['outgoing']-shoot['outgoing'])/abs(mid['outgoing'])
        assert refinement<2e-5 and agreement<2e-5,(mu,w,refinement,agreement)
        assert max(c['wronskian_relative_error'] for c in computed)<2e-6
        rows.append({'mu':mu,'omega':w,'runs':[{'settings':st,'outgoing':pair(c['outgoing']),
                      'K_at_r25':pair(c['K']),'Wronskian_relative_error':c['wronskian_relative_error']}
                      for st,c in zip(settings,computed)],'shooting_outgoing':pair(shoot['outgoing']),
                     'final_refinement_relative_change':float(refinement),'independent_method_relative_difference':float(agreement)})
    # Check the sourced ODE using derivatives of independently evaluated K_x.
    w=.7;mu=.5;x=tortoise(25.);h=.002
    opts=dict(outer=400.,order=8,tail_order=14,rtol=2e-11)
    mid=response(w,mu,**opts)
    minus=response(w,mu,extraction=float(geometry(x-h)[0]),**opts)
    plus=response(w,mu,extraction=float(geometry(x+h)[0]),**opts)
    A,B,T,_=coefficients(x,w,mu)
    residual=((plus['K_x']-minus['K_x'])/(2*h)+A*mid['K_x']+B*mid['K']-T)
    scaled=float(abs(residual)/abs(T));assert scaled<2e-5,scaled
    # An explicit counterexample to trusting a sharp source cutoff.
    cut=response(1.,.5,outer=200.,tail_order=0,rtol=2e-11)
    corrected=response(1.,.5,outer=400.,tail_order=14,rtol=2e-11)
    return {'points':rows,'sourced_ODE_relative_residual':scaled,
            'sharp_cutoff_r200_relative_amplitude_error_at_omega1':float(abs(cut['outgoing']-corrected['outgoing'])/abs(corrected['outgoing']))}


def make_band(directory):
    mu=.5;extraction=200.;frequencies=np.linspace(.52,1.,241)
    values=[];fields=[];wr=[]
    for j,w in enumerate(frequencies):
        out=response(float(w),mu,extraction=extraction,outer=400.,order=8,tail_order=14,rtol=2e-11)
        values.append(out['outgoing']);fields.append(out['K']);wr.append(out['wronskian_relative_error'])
        if j%40==0:print(f'Spectrum {j+1}/{len(frequencies)}',flush=True)
    values=np.array(values);fields=np.array(fields)
    radiative=np.array([a*outer_data(extraction,float(w),mu,8)[0] for w,a in zip(frequencies,values)])
    flux=np.array([signed_flux(w,mu,a) for w,a in zip(frequencies,values)])
    assert np.all(flux<0)
    times=np.linspace(-200,1200,1401);window=tukey(len(frequencies),.25)
    wave=band_waveform(frequencies,radiative*window,times)
    coarse=band_waveform(frequencies[::2],radiative[::2]*window[::2],times)
    full=band_waveform(frequencies,fields*window,times)
    error=float(np.linalg.norm(wave-coarse)/np.linalg.norm(wave))
    assert error<.02,error
    # Save unwindowed frequency amplitudes as well as the chosen tapered band.
    table=np.column_stack((frequencies,values.real,values.imag,fields.real,fields.imag,flux,window))
    np.savetxt(directory/'w4_monopole_spectrum.csv',table,delimiter=',',
               header='omega,outgoing_real,outgoing_imag,K_r200_real,K_r200_imag,signed_flux_dE_domega_over_E2,window',comments='')
    np.savetxt(directory/'w4_monopole_waveform.csv',np.column_stack((times,wave,full,coarse)),delimiter=',',
               header='t,radiative_K_band,full_K_band,coarse_radiative_K_band',comments='')
    return {'mu':mu,'extraction_radius':extraction,'frequency_band':[.52,1.],
            'frequency_step':.002,'coarse_frequency_step':.004,'tukey_alpha':.25,
            'time_range':[-200.,1200.],'waveform_quadrature_relative_L2_difference':error,
            'signed_energy_in_sampled_band_over_E2':float(np.trapezoid(flux,frequencies)),
            'max_Wronskian_relative_error':float(max(wr)),
            'outgoing_amplitude_units':'per kappa E Y00*; raw Fourier transform',
            'waveform_scope':'above-threshold tapered band only; excludes subthreshold response and omitted continuum bands',
            'flux_scope':'AGP Eq.(58) convention conversion; not an independent stress-energy derivation'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output-directory',type=Path,required=True)
    args=parser.parse_args();args.output_directory.mkdir(parents=True,exist_ok=True)
    result={'gate':'W4 initial null-infall monopole benchmark','full_W4_complete':False,
            'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
            'solver_sha256':hashlib.sha256(Path(__file__).with_name('w4_monopole.py').read_bytes()).hexdigest()}
    print('Checking exact source, boundary series, and independent methods',flush=True)
    result['algebra']=verify();result['series']=series_checks();result['numerical']=point_checks()
    result['band_benchmark']=make_band(args.output_directory)
    result['initial_benchmark_passed']=True
    (args.output_directory/'w4_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Initial W4 benchmark passed; full gate remains open',flush=True)
