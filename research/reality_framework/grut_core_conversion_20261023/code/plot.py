#!/usr/bin/env python3
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scattering import ROOT


def main():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(1,3,figsize=(15,4.7),layout='constrained')
    c=json.loads((ROOT/'CONVERGENCE_RESULTS.json').read_text());by={r['name']:r for r in c['runs']}
    for names,label,col in [(['H025','H0125','H00625','H003125'],'2nd order space / 4th time','#c56c3b'),
                            (['N025','N0125','N00625'],'4th order space / 4th time','#246b8e')]:
        axs[0].plot([by[n]['parameters']['h'] for n in names],[100*by[n]['core60_late_mean_fraction'] for n in names],
                    'o-',color=col,label=label)
    axs[0].set_xscale('log',base=2);axs[0].set_xlabel('Radial spacing h (1/mX)')
    axs[0].set_ylabel('Late charge inside radius 60 (%)');axs[0].set_title('Prepared spherical patch')
    axs[0].legend(fontsize=8);axs[0].grid(alpha=.18)
    rate=json.loads((ROOT/'RATE_RESULTS.json').read_text())['rows'];u=[r['p_over_m'] for r in rate]
    axs[1].loglog(u,[r['sigma_v_conversion_cm3_s'] for r in rate],color='#246b8e',label='Hidden conversion')
    axs[1].loglog(u,[r['sigma_v_visible_photon10_cm3_s'] for r in rate],color='#ae5167',label='Leading visible (+ photon sensitivity)')
    axs[1].set_xlabel('Incoming anti-X momentum / mX');axs[1].set_ylabel('sigma v (cm³/s)')
    axs[1].set_title('Fixed Q=200000 ground object');axs[1].legend(fontsize=8);axs[1].grid(alpha=.18)
    trans=json.loads((ROOT/'TRANSPORT_RESULTS.json').read_text())['rows']
    for z,col in [(1e6,'#c56c3b'),(1e8,'#246b8e'),(1.2e11,'#75538e')]:
        row=next(r for r in trans if r['formation_z']==z and r['nodes']==1601)
        hist=np.array(row['history']);elapsed=hist[:,0]/(1/(1+z))
        axs[2].loglog(elapsed,np.maximum(hist[:,2],1e-24),color=col,label=f'Supplied z={z:.2g}')
    axs[2].set_ylim(1e-24,2);axs[2].set_xlim(1,1e5)
    axs[2].set_xlabel('Scale factor / supplied source scale factor')
    axs[2].set_ylabel('Anti-X remaining / initial anti-X')
    axs[2].set_title('Conditional source; inverse conversion included')
    axs[2].legend(fontsize=8);axs[2].grid(alpha=.18)
    fig.suptitle('C1-HP-03/v1: numerical progress and a hidden energy outlet',fontsize=14)
    fig.savefig(ROOT/'construction_results.png',dpi=180)


if __name__=='__main__':main()
