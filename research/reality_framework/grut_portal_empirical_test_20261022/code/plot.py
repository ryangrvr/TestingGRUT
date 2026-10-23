"""Standalone scientific figure, assembled from the archived calculation outputs."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from model import ROOT


def run():
 det=json.loads((ROOT/'DETECTOR_RESULTS.json').read_text());shape=json.loads((ROOT/'SHAPE_RESULTS.json').read_text())
 photo=json.loads((ROOT/'PHOTODISINTEGRATION_RESULTS.json').read_text())
 rows=[r for r in det['cases'] if r['Q']==39944]
 v=np.array([r['v_phi_GeV'] for r in rows]);p=np.array([r['xenonnt_projected_kappa90_percentiles'] for r in rows])*1e9
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':170})
 fig,ax=plt.subplots(2,2,figsize=(12,8.5),layout='constrained')
 a=ax[0,0]
 a.fill_between(v,p[:,0],p[:,2],color='#3586bd',alpha=.18,label='Historical projection: 68% toy spread')
 a.plot(v,p[:,1],color='#16689c',marker='o',label='Historical XENONnT 20 ty median')
 a.plot(v,[r['xenon1t_kappa90_recommended']*1e9 for r in rows],color='#323b45',marker='s',label='XENON1T approximate observed limit')
 a.plot(v,[r['lz2024_count_reweighted_kappa90_ESTIMATE']*1e9 for r in rows],color='#b35b22',ls='--',label='LZ same-mass count estimate')
 a.scatter([1],[shape['LZ_alternative_shape_matched_ESTIMATE']['kappa90_ESTIMATE']*1e9],color='#b35b22',marker='x',s=70,label='LZ shape-matched estimate (reference)')
 a.set(xlabel=r'$v_\phi$ [GeV]',ylabel=r'Coupling $\kappa_{90}$ [$10^{-9}$]',title='Detector reach: Q = 39,944')
 a.legend(fontsize=7,loc='upper left');a.set_ylim(0,4.2)
 a=ax[0,1];data=np.loadtxt(ROOT/'data/reference_shape_comparison.csv',delimiter=',',skiprows=1)
 a.plot(data[:,0],data[:,1],color='#16689c',lw=2,label='3.64 TeV extended object')
 a.plot(data[:,0],data[:,2],color='#b35b22',ls='--',label='Best point WIMP: 25.5 GeV')
 a.set(xlim=(0,35),xlabel='Reconstructed recoil energy [keV]',ylabel='Signal probability per bin',title='Xenon shapes are nearly degenerate')
 a.text(.97,.78,'7.06 expected signals at reference coupling\nIdeal shape-information ratio: 0.00084',ha='right',va='top',transform=a.transAxes,fontsize=8)
 a.legend(fontsize=8)
 a=ax[1,0];hot=np.loadtxt(ROOT/'data/reference_hot_history.csv',delimiter=',',skiprows=1)
 a.plot(hot[:,1],hot[:,5],color='#b35b22',label='Assumed rapid chemical equilibrium')
 a.plot(hot[:,1],hot[:,6],color='#16689c',label='Collisionless self-curvature diagnostic')
 a.axhline(1,color='#333333',ls=':',label='Positive thermal curvature exceeds vacuum scale')
 a.set(xscale='log',xlim=(80,.5),xlabel='SM bath temperature [GeV]',ylabel=r'Induced curvature / $\lambda_\phi v_\phi^2$',title='Hot history: chemical evolution is unresolved')
 a.legend(fontsize=8)
 a=ax[1,1];names=['Nominal','Enlarged source','Long-lived photon\nenvelope','Finer-grid\nenvelope']
 d=np.array([r['relative_changes']['D_over_H_relative'] for r in photo['cases']])*100
 he=np.array([r['relative_changes']['He3_over_D_relative'] for r in photo['cases']])*100
 x=np.arange(len(names));a.bar(x-.18,d,.36,color='#16689c',label='D/H change')
 a.bar(x+.18,he,.36,color='#b35b22',label='He-3 / D change')
 a.set_xticks(x,names,fontsize=8);a.set(ylabel='Increment relative to supplied BBN [%]',title='Cold 20 MeV history: ACROPOLIS reference case')
 a.legend(fontsize=8)
 fig.suptitle('GRUT proposed Higgs portal: conditional empirical test',fontsize=16)
 fig.savefig(ROOT/'data/empirical_summary.png');fig.savefig(ROOT/'data/empirical_summary.pdf')
 plt.close(fig)


if __name__=='__main__':run()
