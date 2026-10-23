#!/usr/bin/env python3
"""Static scientific summary of the declared screens, not exclusion contours."""
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]


def main():
    with (ROOT/'data/scan.csv').open() as f:
        rows=list(csv.DictReader(f))
    rows=[r for r in rows if int(r['Q'])==39944 and float(r['lambda_phi'])==.01]
    mass=np.array([float(r['m_s_GeV']) for r in rows])
    prompt=np.array([float(r['kappa_for_lifetime_1s_central']) for r in rows])
    band_lo=np.array([float(r['kappa_for_lifetime_1s_band_fast']) for r in rows])
    band_hi=np.array([float(r['kappa_for_lifetime_1s_band_slow']) for r in rows])
    lz=np.array([float(r['kappa_LZ_same_mass_count_ESTIMATE']) for r in rows])
    fixed=json.loads((ROOT/'FIXED_YIELD_RESULTS.json').read_text())['rows']
    plt.rcParams.update({'font.size':11,'axes.titlesize':13,'axes.spines.top':False,
                         'axes.spines.right':False,'figure.facecolor':'white'})
    fig,axs=plt.subplots(1,2,figsize=(12.3,5.1),gridspec_kw={'width_ratios':[1.2,1]})
    ax=axs[0]
    ax.fill_between(mass,band_lo,band_hi,color='#d46a42',alpha=.18,label='Stated width-model envelope')
    ax.loglog(mass,prompt,color='#bf522e',lw=2,label=r'Nominal $\tau_s=1$ s')
    ax.loglog(mass,lz,color='#225f94',lw=2,label='LZ same-mass count estimate')
    ax.set(xlabel=r'$m_s$ [GeV], varying $v$ at fixed $\lambda_\phi=.01$',
           ylabel=r'Portal coupling $\kappa$',title='Original branch: lifetime and recoil tension')
    ax.grid(alpha=.16,which='both');ax.legend(fontsize=9,loc='upper right')
    ax.text(.04,.06,'One-second cut is not a BBN fit.\nLZ mapping is approximate.',
            transform=ax.transAxes,fontsize=9,color='#444444')
    ax=axs[1]
    for i,r in enumerate(fixed):
        lo=r['overlap_lower_kappa']/1e-9;hi=r['overlap_upper_kappa']/1e-9
        ax.plot([lo,hi],[i,i],color='#2a806b',lw=9,solid_capstyle='butt')
        ax.plot(r['joint_screen_kappa']/1e-9,i,'o',color='#122f29',ms=7)
        ax.annotate(f"{100*r['object_halo_fraction']:.2f}% object matter\n"
                    +rf"$\tau_s={r['joint_screen_population']['scalar_lifetime_seconds']:.2g}$ s",
                    (.02,i+.13),xycoords=('axes fraction','data'),fontsize=10)
    ax.set(yticks=[0,1],yticklabels=[r'$Q=39\,944,\ \lambda_\phi=.021$',r'$Q=200\,000,\ \lambda_\phi=.03$'],
           xlabel=r'$\kappa/10^{-9}$',title='Two fixed-charge-yield examples',ylim=(-.45,1.55),xlim=(.8,1.10))
    ax.grid(axis='x',alpha=.2)
    ax.text(.02,.02,'Bar: historical median reach to lower LZ estimate.\n'
            'Supplied Q; calibrated v; nonreacting daughters.',transform=ax.transAxes,fontsize=9,color='#444444')
    fig.suptitle('A conditional decay/detection overlap remains; formation and relic evolution are unresolved',fontsize=14,y=.98)
    fig.tight_layout(rect=[0,0,.99,.92],w_pad=2.3)
    fig.savefig(ROOT/'overview.png',dpi=160,bbox_inches='tight')


if __name__=='__main__':main()
