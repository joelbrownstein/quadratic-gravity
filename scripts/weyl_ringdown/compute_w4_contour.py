"""Compute a reproducible finite-radius Bromwich-contour spectrum."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from w4_threshold import response


def sample(u,sigma,outer,extraction,rtol):
    result=response(complex(u,sigma),outer=outer,extraction=extraction,rtol=rtol)
    return [u,result['K'].real,result['K'].imag,result['K_x'].real,result['K_x'].imag]


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--sigma',type=float,required=True);p.add_argument('--step',type=float,default=.002)
    p.add_argument('--offset',type=float,default=0.)
    p.add_argument('--maximum',type=float,default=2.);p.add_argument('--outer',type=float,default=400.)
    p.add_argument('--extraction',type=float,default=200.);p.add_argument('--rtol',type=float,default=2e-10)
    p.add_argument('--workers',type=int,default=2);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.sigma<=0:raise ValueError('Bromwich contour must be above the real axis')
    count=int(np.floor((args.maximum-args.offset)/args.step+1e-9))+1
    u=args.offset+np.arange(count)*args.step
    fun=partial(sample,sigma=args.sigma,outer=args.outer,extraction=args.extraction,rtol=args.rtol)
    rows=[];start=time.monotonic()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for j,row in enumerate(pool.map(fun,u,chunksize=8)):
            rows.append(row)
            if j%100==0:print(f'sigma={args.sigma}: {j+1}/{len(u)}',flush=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    np.savetxt(args.output,np.array(rows),delimiter=',',header='Re_omega,K_real,K_imag,K_x_real,K_x_imag',comments='')
    info={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()}
    info.update({'mu':.5,'left':-35.,'homogeneous_order':10,'source_tail_order':16,
                 'threshold_order':14,'elapsed_seconds':time.monotonic()-start,
                 'solver_sha256':hashlib.sha256(Path(__file__).with_name('w4_threshold.py').read_bytes()).hexdigest(),
                 'amplitude_units':'per kappa E Y00*, raw Fourier transform'})
    args.output.with_suffix('.json').write_text(json.dumps(info,indent=2)+'\n')
    print('Saved',args.output,flush=True)
