"""Plot saved W4 benchmark data; no solver is run by this script."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


if __name__=='__main__':
    directory=Path(__file__).resolve().parents[2]/'notes/validation'
    spectrum=np.genfromtxt(directory/'w4_monopole_spectrum.csv',delimiter=',',names=True)
    wave=np.genfromtxt(directory/'w4_monopole_waveform.csv',delimiter=',',names=True)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
    axes[0].semilogy(spectrum['omega'],abs(spectrum['signed_flux_dE_domega_over_E2']),color='#235d89',lw=2)
    axes[0].set(xlabel=r'$GM\omega$',ylabel=r'$|dE_{\rm ghost}/d\omega|/E_{\rm particle}^{2}$',
                title='Sampled outgoing spectrum')
    axes[0].text(.97,.97,'Signed flux is negative',transform=axes[0].transAxes,ha='right',va='top',fontsize=10)
    norm=max(abs(wave['radiative_K_band']))
    axes[1].plot(wave['t'],wave['radiative_K_band']/norm,color='#235d89',lw=1.2,label=r'$\Delta\omega=0.002$')
    axes[1].plot(wave['t'][::8],wave['coarse_radiative_K_band'][::8]/norm,'.',ms=2.2,color='#ce7434',label=r'$\Delta\omega=0.004$')
    axes[1].set(xlabel=r'$t/(GM)$',ylabel='Normalized radiative auxiliary K',
                title=r'Tapered band waveform, $r/(GM)=200$')
    axes[1].legend(frameon=False,loc='upper right')
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle(r'Null-infall monopole: $GMm_2=0.5$, frequency band $0.52\leq GM\omega\leq1.00$'+'\nInitial W4 benchmark — not a complete waveform or detector strain',fontsize=12)
    fig.savefig(directory/'w4_monopole_benchmark.png',dpi=170)
    fig.savefig(directory/'w4_monopole_benchmark.pdf')
    print('Saved benchmark PNG and PDF; matplotlib',matplotlib.__version__)
