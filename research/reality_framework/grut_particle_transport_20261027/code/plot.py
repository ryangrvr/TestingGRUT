"""Standalone scientific figures; all particle-like quantities retain action labels."""
from pt_common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    result=read('RESULTS.json');fig,ax=plt.subplots(2,2,figsize=(12,8.2),constrained_layout=True)
    windows=[read(f'WINDOW_base48_{a}_{b}.json') for a,b in [(0,6),(6,12),(14,24),(24,34)]]
    x=np.array([np.mean(w['window_s']) for w in windows]);coh=[w['coherent_negative_frequency_action'] for w in windows]
    wave=[w['wave_negative_frequency_action'] for w in windows];total=[w['total_negative_frequency_action'] for w in windows]
    a=ax[0,0];a.plot(x,total,'o-',color='black',label='Total');a.plot(x,coh,'o-',label='Coherent mean');a.plot(x,wave,'o-',label='Connected waves')
    a.set(xlabel='Conformal window midpoint s',ylabel='Classical negative-frequency action / net current',title='A stable total; a changing split',ylim=(0,1.92));a.legend(frameon=False)
    a=ax[0,1];labels=[]
    for j,r in enumerate(result['rows']):
        labels.append(r['case']['tag'].replace('base48','base').replace('IR12_','L12 ').replace('IR24_','L24 ').replace('_',' '))
        a.plot([j-.12,j+.12],[w['total_negative_action'] for w in r['windows']],color='.7',lw=.8)
        for k,c in enumerate(['#126c9a','#bd6126']):a.scatter(j+(-.12 if k==0 else .12),r['windows'][k]['total_negative_action'],color=c,s=25,label=['s=14–24','s=24–34'][k] if j==0 else None)
    a.set(xticks=np.arange(len(labels)),xticklabels=labels,ylabel='Classical total action / net current',title='IR, width, volume and numerical controls',ylim=(1.737,1.756));a.tick_params(axis='x',rotation=45);a.legend(frameon=False)
    old=np.load(PARENT/'FIELD_CURVES.npz');t=old['N48_tau'];y=old['N48_energy_outside_spatial_mean_fraction'];sel=(t>=700)&(t<=800+1e-9)
    mean=next(r for r in read('LEGACY_ORBIT_AVERAGES.json')['rows'] if r['tag']=='N48')['statistics']['energy_outside_spatial_mean_fraction']['cosmic_time_mean']
    a=ax[1,0];a.plot(t[sel],100*y[sel],color='#334d66',lw=1.3);a.axhline(100*mean,color='#bd6126',ls='--',label=f'Time mean {100*mean:.2f}%');a.scatter([t[sel][-1]],[100*y[sel][-1]],color='black',s=25,label='Old final-phase sample')
    a.set(xlabel='Cosmic time label tau',ylabel='Energy outside spatial mean (%)',title='The orbit envelope matters',ylim=(40,80));a.legend(frameon=False)
    a=ax[1,1];w=read('WINDOW_base48_14_24.json');m=[0,7,12,20,26];v=[w['equal_time'][f'ref_{n}']['hann2_mean'] for n in m]
    a.plot(m,v,'o-',color='#126c9a',label='Time-averaged fixed-mass references');a.axhline(w['wave_negative_frequency_action'],color='black',ls='--',label='Measured wave frequency action')
    a.set(xlabel='Reference mass squared (H_i units)',ylabel='Nonzero-mode reference / net current',title='Massive reference insensitivity after averaging');a.legend(frameon=False,fontsize=8)
    fig.suptitle('Measured classical field budget — reduced X/phi model',fontsize=15)
    fig.savefig(ROOT/'BUDGET.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(12,4.6),constrained_layout=True);co=read('COHERENT_MATCHING.json');a=ax[0]
    for r in co['rows']:a.scatter(max(r['duration'],.04),r['free_nminus_per_charge'],s=18,color='#126c9a',alpha=.65)
    a.axhline(co['radial_action_per_charge'],color='black',ls='--',label=f"Radial action {co['radial_action_per_charge']:.4f}")
    a.axhline(co['lab_negative_frequency_action_per_charge'],color='#bd6126',ls=':',label=f"Initial lab-frequency action {co['lab_negative_frequency_action_per_charge']:.4f}")
    a.set_xscale('log');a.set(xticks=[.04,.2,2,20,80],xticklabels=['instant','0.2','2','20','80'],xlabel='Supplied matching duration (conformal units)',ylabel='Free out antiparticles / coherent charge',title='A coherent field has a dynamical out-state',ylim=(-.2,19));a.legend(frameon=False,fontsize=8)
    a=ax[1]
    for f in sorted(ROOT.glob('RESPONSE_s*.json')):
        r=read(f.name);a.plot([v['k_center'] for v in r['rows']],[v['anomalous_over_normal_RMS'] for v in r['rows']],'o-',ms=3,label=f"s0={r['s0']:g}, probe {r['probe_seed']}")
    a.axvspan(0,1.6,color='#126c9a',alpha=.08,label='Count-dominant infrared');a.axhline(0,color='.6',lw=.8)
    a.set(xlabel='Continuum shell momentum k/H_i',ylabel='Conjugate-channel / normal response RMS',title='The medium response has substantial mixing',ylim=(0,1.9),xlim=(0,8.2));a.legend(frameon=False,fontsize=8)
    fig.suptitle('Particle matching controls — no cosmological transition inferred',fontsize=14)
    fig.savefig(ROOT/'MATCHING.png',dpi=180);plt.close(fig)
    print('Rendered BUDGET.png and MATCHING.png')

if __name__=='__main__':main()
