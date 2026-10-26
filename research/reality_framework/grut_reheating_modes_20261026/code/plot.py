import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mode_common import ROOT

def main():
    coarse=np.load(ROOT/'MODE_PROFILES.npz');fine=np.load(ROOT/'REFINED_PROFILES.npz');phase=np.load(ROOT/'PHASE_PROFILES.npz')
    results=json.loads((ROOT/'MODE_RESULTS.json').read_text());rows={r['tag']:r for r in results['rows']};t=coarse['times']
    orbit=json.loads((ROOT/'QUARTIC_ORBIT.json').read_text())
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    colors=['#146B8A','#C76A28'];tags=['b3_g0.5','b5_g1']
    for tag,col,lab in zip(tags,colors,[r'$(m_I,\Gamma_I)/H_i=(3,0.5)$',r'$(5,1)$']):
        a=ax[0,0];a.plot(t,coarse[tag][:,0,5]/rows[tag]['late_NX'],color=col,label=lab)
        a=ax[0,1];ks=phase[tag+'_ks'];rr=phase[tag][:,ks>=.01]
        a.plot(t,3*np.max(abs(rr[:,:,0]),axis=1),color=col,label=lab+' amplitude')
        a.plot(t,3*np.max(abs(rr[:,:,1]),axis=1),color=col,ls='--',label=lab+' local charge')
        a=ax[1,0];d=fine[tag]
        a.plot(fine['ks'],d[:,:,0].max(axis=0),color=col,label=lab+r' through $H_it=80$')
        a.plot(fine['ks'],d[t<=3,:,0].max(axis=0),color=col,ls=':',label=lab+' through 3')
    ax[0,0].axvspan(0,3,color='.8',alpha=.25)
    ax[0,0].set(xlim=(0,10),xlabel=r'$H_i t$',ylabel=r'$N_X(t)/\overline{N_X}$',title='A  Charge is established early')
    ax[0,0].legend(fontsize=9)
    ax[0,1].axvline(3,color='.5',ls=':')
    ax[0,1].set(xlabel=r'$H_i t$',ylabel='linear modulation amplitude (%)',title='B  Supplied 0.03-radian phase ripple')
    ax[0,1].legend(fontsize=8)
    ax[1,0].set(xlabel=r'comoving $k/H_i$',ylabel='relative-amplitude response coefficient',yscale='log',title='C  Later response is wavelength dependent')
    ax[1,0].legend(fontsize=8)
    ax[1,1].plot([r['k'] for r in orbit['rows']],[r['mu_conformal'] for r in orbit['rows']],color='#773B83',lw=2)
    ax[1,1].set(xlabel=r'comoving $k/H_i$',ylabel=r'$\mu$ per conformal time $H_i\eta$',title='D  Isolated quartic orbit has a growing band')
    ax[1,1].text(.98,.92,'Auxiliary matched orbit; not the full cosmology\nCircular control: no resolved growth',ha='right',va='top',transform=ax[1,1].transAxes,fontsize=8.5)
    fig.suptitle('Coupled linear response: early charge creation and later self-resonance',fontsize=14)
    fig.savefig(ROOT/'mode_response.png',dpi=160)

if __name__=='__main__':main()
