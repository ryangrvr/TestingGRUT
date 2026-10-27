"""Scientific plots of actual local energy and motion diagnostics."""
from cw_common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    final=read('FINAL_AUDIT.json')['cases'];history=read('HISTORY_AUDIT.json')['cases']
    fig,axes=plt.subplots(2,2,figsize=(11,7.6),layout='constrained')
    main=[c for c in final['mixed_N96']['candidates'] if 'spheres' in c]
    for c in main:
        rr=c['spheres'];r=[x['radius'] for x in rr]
        label='Candidate '+str(c['candidate']['id'])
        axes[0,0].plot(r,[x['E_over_Q'] for x in rr],'o-',label=label)
        axes[0,1].plot(r,[x['E_X_over_Q'] for x in rr],'o-',label=label)
    for ax,title in zip(axes[0],['Total energy in fixed spheres','Charged energy; neutral sector booked separately']):
        ax.axhline(1,color='k',ls='--',lw=.8);ax.set(xlabel=r'Radius $m_X R$',ylabel=r'$E/(m_X Q)$',title=title);ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8)
    for tag,rr in history.items():
        t=[];kin=[];scalar=[]
        for row in rr:
            c=row['candidates'];K=sum(x['mask']['kinetic_X'] for x in c)
            excess=sum(x['mask']['kinetic_radial']+x['mask']['kinetic_angular_dispersion']+x['mask']['kinetic_at_exact_nodes'] for x in c)
            t.append(row['t']);kin.append(excess/K if K else np.nan)
            scalar.append(row['whole_box']['E_phi']/row['whole_box']['E_total'])
        name='Mixed seed' if tag=='mixed_N96' else 'Mostly one-sign seed'
        axes[1,0].plot(t,kin,'o-',label=name)
        axes[1,1].plot(t,scalar,'o-',label=name)
    axes[1,0].set(xlabel=r'Time $m_X t$',ylabel='Excess / charged kinetic energy',title='Motion beyond common rotation')
    axes[1,1].set(xlabel=r'Time $m_X t$',ylabel=r'$E_\phi/E_{\rm total}$',title='Scalar Hamiltonian reservoir; not a particle count')
    for ax in axes[1]:ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('Existing cold 3D fields: local energy and common-rotation audit',fontsize=13)
    fig.savefig(ROOT/'CORE_LEDGER.png',dpi=155);plt.close(fig)

if __name__=='__main__':main()
