"""Supplement the saved low-frequency spectra for front reconstruction."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import hashlib
import json
from pathlib import Path
import numpy as np
from w4_threshold import response
from w4_monopole import tortoise


def sample(u,sigma):
    radius=max(800.,200.*u)
    value=response(complex(u,sigma),outer=radius,rtol=2e-11)['K']
    B=value*np.exp(1j*(u+1j*sigma)*tortoise(200.))
    return [u,value.real,value.imag,B.real,B.imag,radius]


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--sigma',type=float,default=.004)
    p.add_argument('--maximum',type=float,default=12.)
    p.add_argument('--step',type=float,default=.05)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--workers',type=int,default=2)
    args=p.parse_args();u=2+np.arange(round((args.maximum-2)/args.step)+1)*args.step
    rows=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for j,row in enumerate(pool.map(partial(sample,sigma=args.sigma),u,chunksize=2)):
            rows.append(row)
            if j%20==0:print(f'Front spectrum sigma={args.sigma}: {j+1}/{len(u)}',flush=True)
    np.savetxt(args.output,rows,delimiter=',',comments='',
               header='Re_omega,K_real,K_imag,B_real,B_imag,outer_radius')
    metadata={'mu':.5,'extraction':200.,'sigma':args.sigma,'frequency_step':args.step,
              'maximum':args.maximum,'outer_rule':'max(800,200*Re_omega)','rtol':2e-11,
              'threshold_solver_sha256':hashlib.sha256(Path(__file__).with_name('w4_threshold.py').read_bytes()).hexdigest()}
    args.output.with_suffix('.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('Saved',args.output,flush=True)
