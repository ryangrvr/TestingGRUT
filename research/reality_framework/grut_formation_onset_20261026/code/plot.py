"""Plot limited onset calculations with their different assumptions visible."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from common import *

def main():
    c=read('CIRCULAR_RESULTS.json');e=read('ELLIPTIC_RESULTS.json');r=read('ROLL_RESULTS.json')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(1,3,figsize=(15.5,4.6),layout='constrained')
    rows=[x for x in c['rows'] if x['a']<=1]
    ax[0].plot([x['a'] for x in rows],[x['adiabatic_growth'] for x in rows],label='Adiabatic elimination',color='#b75018',ls='--')
    ax[0].plot([x['a'] for x in rows],[x['full_growth'] for x in rows],label='Binding field retained',color='#2463a8')
    ax[0].scatter([1],[0],color='#2463a8',zorder=5)
    ax[0].set(xlabel='Circular amplitude |X| / x_c',ylabel='Maximum growth exponent / m_X',title='Circular vacuum diagnostic')
    ax[0].legend(fontsize=8,loc='upper left')
    ax[0].text(.45,.04,'No exponential growth at the exact\nendpoint in the full circular problem',transform=ax[0].transAxes,fontsize=8)
    rows=e['rows']
    ax[1].plot([x['k'] for x in rows],[x['mu_phi'] for x in rows],label='Neutral binding field',color='#b75018')
    ax[1].plot([x['k'] for x in rows],[x['mu_X'] for x in rows],label='Charged field',color='#2463a8')
    ax[1].set(xlabel='Momentum k / m_X',ylabel='Floquet exponent / m_X',title='Matched cold elliptic orbit',xlim=(0,2.5))
    ax[1].legend(fontsize=8,loc='upper right')
    ax[1].text(.48,.55,'μ_phi ≈ 0.859 at k=0\nμ_X ≈ 0.103 at k≈0.341\nNeutral growth has no linear X charge',transform=ax[1].transAxes,fontsize=8)
    for i,row in enumerate(r['rows']):
        ls='-' if i==0 else '--';tag='10⁻⁶' if i==0 else '10⁻¹²'
        ax[2].plot(row['times'],np.abs(row['phi_over_v']),ls=ls,color='#b75018',label='|phi|/v, seed '+tag)
        ax[2].plot(row['times'],row['X_amplitude_for_supplied_seed'],ls=ls,color='#2463a8',label='X response, phi seed '+tag)
    ax[2].axhline(.1,color='.4',ls=':',lw=1)
    ax[2].set(yscale='log',ylim=(1e-12,.3),xlabel='Local time m_X t',ylabel='Relative field amplitude',title='Supplied seeds; stop at |phi|=0.1v')
    ax[2].legend(fontsize=7,loc='lower right')
    fig.suptitle('GRUT formation onset — vacuum controls, no plasma or final object population',fontsize=13)
    fig.savefig(ROOT/'ONSET.png',dpi=170);plt.close(fig)

if __name__=='__main__':main()
