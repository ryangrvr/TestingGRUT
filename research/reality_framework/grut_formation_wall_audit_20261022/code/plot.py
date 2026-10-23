"""Scientific summary of conditional onset and energy accounting."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from calculate import ROOT


def run():
 result=json.loads((ROOT/'RESULTS.json').read_text())
 ref=next(r for r in result['cases'] if r['Q']==39944 and r['v_phi_GeV']==1.)
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':170})
 fig,ax=plt.subplots(2,2,figsize=(11.7,8.1),layout='constrained')
 d=np.loadtxt(ROOT/'data/reference_dispersion.csv',delimiter=',',skiprows=1)
 a=ax[0,0];g=ref['linear_condensate']['maximum_growth_GeV']
 a.plot(np.sqrt(d[:,0]),d[:,1]/g,color='#176c96',lw=2,label='Full coupled dispersion')
 a.plot(np.sqrt(d[:,0]),d[:,2]/g,color='#b15b28',ls='--',label='Dilute nonrelativistic limit')
 a.set(xlabel=r'$k / \sqrt{a f_0^2}$',ylabel=r'Growth / maximum growth',title='A supplied condensate has a growing band')
 a.text(.34,.12,r'$a=g^4/\lambda_\phi-\lambda_X=0.009$'+'\n'+r'$\Gamma_{\max}/H=4.60\times10^8$',transform=a.transAxes)
 a.legend(fontsize=8)
 a=ax[0,1];rows=[r for r in result['cases'] if r['Q']==39944]
 v=np.array([r['v_phi_GeV'] for r in rows]);charges=[r['linear_condensate']['charge_in_wavelength_cube_DIAGNOSTIC'] for r in rows]
 a.plot(v,charges,color='#176c96',marker='o',label='Charge in preferred wavelength cube')
 a.axhline(39944,color='#b15b28',ls='--',label='Target object charge')
 a.set(yscale='log',xlabel=r'$v_\phi$ [GeV]',ylabel='Charge',title='Linear length does not select the target objects')
 a.legend(fontsize=8);a.text(.04,.32,'The curve is a cell-size diagnostic,\nnot a predicted final charge distribution.',transform=a.transAxes,fontsize=9)
 a=ax[1,0];budgets=[ref['empty_cold_portal_max_DM_fraction_nominal'],ref['empty_cold_portal_max_DM_fraction_enlarged_response'],ref['binding_release_over_object_DM_energy']]
 a.bar(['Cold portal\nnominal','Cold portal\nenlarged source','Free carriers to\nreference objects'],budgets,color=['#83adc4','#176c96','#b15b28'])
 a.axhline(1,color='#444444',ls=':',label='Supplied final DM energy')
 a.set(yscale='log',ylim=(1e-10,3),ylabel='Energy per entropy / final DM input',title='Formation has a much larger release to track')
 a.legend(fontsize=8)
 a=ax[1,1];v=np.linspace(.5,1.25,100)
 tension=(2*np.sqrt(.02)*v**3/3)**(1/3)*1000
 a.plot(v,tension,color='#b15b28',label='Neutral tree-level wall')
 a.axhline(.81,color='#444444',ls='--',label='2026 stable-network bound')
 a.set(yscale='log',xlabel=r'$v_\phi$ [GeV]',ylabel=r'$\sigma^{1/3}$ [MeV]',title='A surviving standard wall network fails')
 a.legend(fontsize=8);a.text(.03,.39,'Large tension does not prove\nrestoration or network survival.',transform=a.transAxes,fontsize=9)
 fig.suptitle('GRUT formation prerequisites: supplied charge, unresolved final state',fontsize=15)
 fig.savefig(ROOT/'data/formation_summary.png');plt.close(fig)


if __name__=='__main__':run()
