"""Render saved threshold-inclusive W4 waveform and convergence evidence."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


if __name__=='__main__':
    directory=Path(__file__).resolve().parents[2]/'notes/validation'
    wave=np.genfromtxt(directory/'w4_extended_waveform.csv',delimiter=',',names=True)
    spectrum=np.genfromtxt(directory/'w4_contour_sigma004.csv',delimiter=',',names=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(11,7),layout='constrained');grid=fig.add_gridspec(2,2)
    top=fig.add_subplot(grid[0,:]);left=fig.add_subplot(grid[1,0]);right=fig.add_subplot(grid[1,1])
    top.plot(wave['t'],wave['K_postpassage'],lw=1.1,color='#235d89')
    top.set(xlabel=r'$t/(GM)$',ylabel=r'$K/(\kappa E Y_{00}^{*})$',title='Retarded finite-radius monopole, post-passage window')
    left.semilogy(spectrum['Re_omega'],np.hypot(spectrum['K_real'],spectrum['K_imag']),color='#235d89')
    left.axvline(.5,color='#b15c28',ls='--',label=r'$\mathrm{Re}\,\omega=m_2$')
    left.set(xlabel=r'$GM\,\mathrm{Re}\,\omega$',ylabel=r'$|\widehat K/(\kappa E Y_{00}^{*})|$',
             title=r'Causal contour: $GM\,\mathrm{Im}\,\omega=0.004$')
    left.legend(frameon=False)
    peak=max(abs(wave['K_postpassage']))
    for name,label,color in [('frequency_step_002','Frequency spacing','#235d89'),
                              ('contour_sigma003','Contour height','#b15c28'),
                              ('cutoff_15','Upper frequency limit','#48845b')]:
        error=abs(wave[name]-wave['K_postpassage'])/peak
        right.semilogy(wave['t'],np.maximum(error,1e-12),label=label,color=color,lw=.9)
    right.set(xlabel=r'$t/(GM)$',ylabel='Difference / peak waveform',title='Independent numerical refinements',ylim=(1e-10,3e-3))
    right.legend(frameon=False)
    for ax in (top,left,right):ax.grid(alpha=.18)
    fig.suptitle(r'W4 extension: $GMm_2=0.5$, $r/(GM)=200$'+'\nThreshold and subthreshold response included; estimated high-frequency tail',fontsize=12)
    fig.savefig(directory/'w4_threshold_waveform.png',dpi=170)
    fig.savefig(directory/'w4_threshold_waveform.pdf')
    print('Saved threshold-inclusive waveform figures')
