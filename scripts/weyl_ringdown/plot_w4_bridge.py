"""Plot the saved connecting interval and empirical numerical floor."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=Path(__file__).resolve().parents[2]/'notes/validation'
d=np.genfromtxt(p/'w4_bridge_waveform.csv',delimiter=',',names=True)
j=json.loads((p/'w4_bridge_validation.json').read_text())
v=d['v'];K=d['K_bridge']
fig,axes=plt.subplots(2,2,figsize=(11,7.5),constrained_layout=True)
ax=axes[0,0];early=v<=10
ax.plot(v[early],K[early]);ax.set_yscale('symlog',linthresh=2e-13)
ax.set(xlabel=r'$v=t+r_*$',ylabel=r'$K/\mathcal{N}$',title='Decay away from the sharp front')
ax=axes[0,1];ax.plot(v,1e9*K)
ax.axvline(209.1902397002692,color='0.5',ls='--',label='Previous late window begins')
ax.set(xlabel=r'$v$',ylabel=r'$10^9 K/\mathcal{N}$',title='Bridge through the onset of the late waveform')
ax.legend(fontsize=8)
ax=axes[1,0];ax.semilogy(v,np.maximum(abs(K),1e-22),label='Computed magnitude')
floor=j['quiet_max_required_comparison']
ax.hlines(floor,10,400,color='tab:orange',label='Largest required comparison in weak interval')
ax.axvspan(10,400,color='0.5',alpha=.1)
ax.set(xlabel=r'$v$',ylabel=r'$|K|/\mathcal{N}$',title='Small values must be compared with the numerical floor')
ax.legend(fontsize=8)
ax=axes[1,1]
for name in ('cutoff_4','cutoff_8','coarser_spectral_grid','refined_contour_sigma003','gauss_order_12'):
    ax.semilogy(v,np.maximum(abs(d[name]-K),1e-22),label=name.replace('_',' '))
ax.axhline(2e-12,color='k',ls=':',lw=1,label='Bridge acceptance target')
ax.set(xlabel=r'$v$',ylabel=r'$|\Delta K|/\mathcal{N}$',title='Convergence checks across the connecting interval')
ax.legend(fontsize=7)
fig.suptitle(r'W4 time bridge: $\mu=0.5$, $r=200$ — auxiliary field, not detector strain')
for ext in ('png','pdf'):fig.savefig(p/('w4_time_bridge.'+ext),dpi=180)
