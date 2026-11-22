"""Plot the saved front and its convergence diagnostics; no solver rerun."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

data=Path(__file__).resolve().parents[2]/'notes/validation'
d=np.genfromtxt(data/'w4_null_front.csv',delimiter=',',names=True)
v=d['v'];positive=v>0
fig,axes=plt.subplots(2,2,figsize=(11,7.5),constrained_layout=True)
ax=axes[0,0]
ax.plot(v,d['K_front']*1e7,label='Schwarzschild + analytic front')
ax.plot(v,d['analytic_reference']*1e7,'--',label='Analytic reference')
ax.axvline(0,color='0.6',lw=.7);ax.set_xlim(-.05,.3)
ax.set(xlabel=r'$v=t+r_*$',ylabel=r'$10^7 K/\mathcal{N}$',title='Continuous passage, sharp initial slope')
ax.legend(fontsize=8)
ax=axes[0,1]
ax.semilogx(v[positive],(2*200**2*d['K_front'][positive]/v[positive]))
ax.semilogx(v[positive],np.log(.5**2*200*v[positive]/2)+2*np.euler_gamma,'--')
ax.set(xlabel=r'$v$',ylabel=r'$2r^2K/(\mathcal{N}v)$',title=r'Approach to $\log(\mu^2rv/2)+2\gamma_E$')
ax=axes[1,0]
for name in ('cutoff_4','cutoff_8','coarser_grid','contour_sigma003'):
    ax.plot(v[positive],(d[name]-d['K_front'])[positive]*1e12,label=name.replace('_',' '))
ax.set(xlabel=r'$v$',ylabel=r'$10^{12}\Delta K/\mathcal{N}$',title='Numerical convergence differences')
ax.legend(fontsize=8)
ax=axes[1,1]
ax.plot(v,(d['unprojected_K']-d['analytic_reference'])*1e12,label='Raw finite-band remainder')
ax.plot(v[positive],d['corrected_remainder'][positive]*1e12,label='Moment-corrected remainder')
ax.axvline(0,color='0.6',lw=.7)
ax.set(xlabel=r'$v$',ylabel=r'$10^{12}R/\mathcal{N}$',title='Analytic constraints and truncation error')
ax.legend(fontsize=8)
fig.suptitle(r'W4 null-infall auxiliary monopole: $\mu=0.5$, $r=200$ (not detector strain)')
for suffix in ('png','pdf'):
    fig.savefig(data/('w4_null_front.'+suffix),dpi=180)
