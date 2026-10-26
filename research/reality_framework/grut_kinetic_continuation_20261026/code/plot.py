"""Scientific figure of conditional controls; all stopping rules explicit."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bootstrap import ROOT
import ks_common as K

def main():
    cases=json.loads((ROOT/'MOMENT_CONTROLS.json').read_text())['cases']
    follow=json.loads((ROOT/'FOLLOWUP_CONTROLS.json').read_text())['cases']
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(1,3,figsize=(15,4.7),layout='constrained')
    colors=['#2463a8','#b75018','#33865f','#8a4fa0']
    for (key,label),color in zip([('e1e7_band_gas','Gas, e=10⁻⁷'),('e1e7_band_instant','Instant objects, e=10⁻⁷'),
            ('e1e6_band_gas','Gas, e=10⁻⁶'),('e1e6_band_instant','Instant objects, e=10⁻⁶')],colors):
        r=cases[key];tt=r['trajectories'];x=[t['x'] for t in tt];y=[max(t['f_at_stop'],1e-8) for t in tt]
        ax[0].plot(x,y,label=label,color=color)
        marker='x' if r['stop'].startswith('kinetic') else 'o'
        ax[0].scatter(x[-1],y[-1],marker=marker,color=color,zorder=5)
    ax[0].axhline(.01,color='.4',ls='--',lw=1)
    ax[0].set(yscale='log',ylim=(1e-8,30),xlabel='ln(a / a_ref)',ylabel='Symmetric mass / fixed DM budget',title='Forced-equilibrium moment recipes')
    ax[0].legend(fontsize=8,loc='center right',bbox_to_anchor=(.99,.48))
    ax[0].text(.02,.05,'×  relaxation diagnostic reaches H\n○  supplied stop; values below 10⁻⁸ clipped',transform=ax[0].transAxes,fontsize=8)
    scan=sorted((r['R'],r['final']['f_at_stop']) for k,r in follow.items() if k.startswith('gas_scan'))
    ax[1].plot(*np.array(scan).T,'o-',color=colors[0],label='Gas-only diagnostic')
    for k in ['e1e7_band_gas','e1e7_broad_gas','e1e6_band_gas']:
        r=cases[k];ax[1].scatter(r['R'],r['final']['f_at_stop'],marker='s',s=40,color=colors[1],zorder=4)
    ax[1].axhline(.01,color='.4',ls='--',lw=1,label='1% allowance')
    ax[1].set(xscale='log',yscale='log',xlabel='Total extra energy ratio R',ylabel='Conditional ratio at T=0.002 m_X',title='A recipe scan, not a kinetic threshold')
    ax[1].legend(fontsize=8)
    obj=K.Objects();ax[2].plot(obj.Q,obj.EQ,color=colors[2],lw=2)
    ax[2].axhline(1,color='.4',ls='--',lw=1)
    for q,c in [(25000.,colors[1]),(K.Q_MIN,'.45'),(K.THM['objects']['E_over_Q_equals_m_at']['Q'],colors[0])]:
        ax[2].axvline(q,color=c,ls=':',lw=1.4)
    ax[2].text(27000,.88,'Fastest fragments\nQ ≈ 25,000',color=colors[1],fontsize=8)
    ax[2].text(80000,1.012,'Binding crossing\nQ ≈ 65,620',color=colors[0],fontsize=8)
    ax[2].text(53000,.845,'Local branch floor ≈ 51,964',color='.35',fontsize=8)
    ax[2].set(xscale='log',xlim=(20000,1e6),ylim=(.83,1.04),xlabel='Object charge Q',ylabel='E(Q) / (Q m_X)',title='Merger outcome remains unprovided')
    fig.suptitle('GRUT kinetic continuation — conditional calculations, no final cosmological abundance',fontsize=13)
    fig.savefig(ROOT/'CONTROLS.png',dpi=170);plt.close(fig)

if __name__=='__main__':main()
