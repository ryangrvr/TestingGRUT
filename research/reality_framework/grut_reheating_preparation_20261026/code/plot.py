"""Standalone scientific figure from the saved numerical results."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import ROOT


def main():
    re=json.loads((ROOT/'REHEATING_RESULTS.json').read_text())
    rows={r['tag']:r for r in re['rows']};cols=re['columns']
    boundary=json.loads((ROOT/'BOUNDARY_RESULTS.json').read_text())
    portals=json.loads((ROOT/'PORTAL_RESULTS.json').read_text())
    profiles=np.load(ROOT/'REHEATING_PROFILES.npz')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(2,2,figsize=(12,8.3),layout='constrained')
    a=axs[0,0]
    for tag,label,color in [('b2_g0.1',r'$m_I/H_i=2$','#146B8A'),('b10_g0.1',r'$m_I/H_i=10$','#C76A28')]:
        data=profiles[tag];mask=data[:,0]<=10
        a.plot(data[mask,0],data[mask,cols.index('R')]/data[mask,cols.index('H')]**2,label=label,color=color,lw=1.4)
    a.axhline(3,color='.55',ls=':',label='matter average')
    a.set(xlabel=r'$H_i t$',ylabel=r'$R/H^2$',title='A  Curvature evolves through release')
    a.legend(fontsize=9,loc='lower right');a.set_ylim(-7,13)
    a=axs[0,1]
    for beta,color in [(2,'#146B8A'),(10,'#C76A28')]:
        dec=np.array([.03,.1,.3])
        yy=[abs(rows[f'b{beta}_g{d:g}']['yield_over_YQ']) for d in dec]
        a.plot(dec,yy,'o-',label=fr'$m_I/H_i={beta}$',color=color)
    a.axhline(1,color='.3',ls='--',label='required yield')
    a.set(xscale='log',xlabel=r'$\Gamma_I/H_i$',ylabel=r'$|Y_X|/Y_Q$',ylim=(0,1.09),
          title='B  All six unchanged-scale yields are low')
    a.set_xticks([.03,.1,.3],['0.03','0.1','0.3']);a.legend(fontsize=9)
    a=axs[1,0]
    fine=boundary['models'][1]['rows'];rr=[fine[2],fine[0],fine[1],fine[3]]
    a.bar(np.arange(4),[2*r['light_charge_torque'] for r in rr],color=['#146B8A','.65','#C76A28','.65'])
    a.axhline(0,color='.3',lw=.8)
    a.set_xticks(range(4),[r'$-\pi/2$','0',r'$+\pi/2$',r'$\pi$'])
    a.set(xlabel='supplied relative phase',ylabel='signed torque, toy units',ylim=(-3.5,4.2),
          title='C  Equal listed initial properties leave a sign choice')
    a.text(.04,.94,'Opposite phases: equal energy, occupations,\nentropies and initial hidden/SM independence',
           transform=a.transAxes,va='top',fontsize=8.5)
    a=axs[1,1];sr=portals['source_rows'][1:]
    ep=[r['epsilon_BX'] for r in sr]
    a.loglog(ep,[r['R_total'] for r in sr],color='#773B83',label='total raw source',lw=2)
    a.loglog(ep,[r['R_primordial'] for r in sr],color='#146B8A',label='primordial hidden branch',ls='--')
    a.loglog(ep,[r['R_thermal'] for r in sr],color='#C76A28',label='thermal inverse-decay source',ls=':')
    a.axhline(portals['analytic_bounds']['minimum_lower_bound_R'],color='.4',ls='-.',label='no-mixing primordial band')
    a.set(xlabel=r'$\epsilon_{BX}$',ylabel='conformal energy ratio R',ylim=(1e-3,2e6),
          title='D  A visible decay portal also produces vectors')
    a.legend(fontsize=8,loc='lower right')
    fig.suptitle('Conditional GRUT calculations: reheating, preparation and portal sources',fontsize=14)
    fig.savefig(ROOT/'reheating_preparation.png',dpi=160)


if __name__=='__main__':main()
