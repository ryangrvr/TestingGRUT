"""Scientific figures and small, explicitly downsampled visualization arrays."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from analyze import *

def blockmean(a,n=48):
    N=a.shape[-1];assert N%n==0;f=N//n
    return a.reshape(n,f,n,f,n,f).mean(axis=(1,3,5))

def main():
    runs=json.loads((ROOT/'RUNS.json').read_text())['cases'];analysis=json.loads((ROOT/'ANALYSIS.json').read_text())['cases'];profiles={}
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(11,9),layout='constrained')
    for col,tag in enumerate(['mixed_N96','onesign_N96']):
        r,q,p,j,labs=analyze_snapshot(runs[tag]['snapshots'][-1],details=True);N=r['N'];L=r['L']
        candidates=r['catalogs'][1]['candidates'];z=int(round(candidates[0]['centre'][2]/L*N))%N if candidates else N//2
        charge=j[:,:,z]/j.mean();ph=q[2,:,:,z]/C.NU
        im=axes[0,col].imshow(charge.T,origin='lower',extent=(0,L,0,L),cmap='RdBu_r',norm=SymLogNorm(linthresh=1,vmin=-20,vmax=100),interpolation='nearest')
        axes[0,col].contour((labs[:,:,z]>0).T,levels=[.5],colors='black',linewidths=.45,extent=(0,L,0,L))
        axes[0,col].set(title=('Zero-mean phi seed' if col==0 else 'Mostly one-sign phi seed')+' — charge current',xlabel='m_X x',ylabel='m_X y')
        fig.colorbar(im,ax=axes[0,col],shrink=.8,label='j0 / mean j0')
        im=axes[1,col].imshow(ph.T,origin='lower',extent=(0,L,0,L),cmap='coolwarm',vmin=-1.5,vmax=1.5,interpolation='nearest')
        axes[1,col].set(title=f'Binding field; slice m_X z={z*L/N:.1f}',xlabel='m_X x',ylabel='m_X y')
        fig.colorbar(im,ax=axes[1,col],shrink=.8,label='phi / v')
        for name,a in [('charge',j),('phi',q[2]),('rho',local_energy(q,p,L))]:profiles[tag+'_'+name]=blockmean(a).astype(np.float32)
        profiles[tag+'_L']=np.array(L);profiles[tag+'_t']=np.array(r['t'])
    fig.suptitle('Cold 3D pilot at m_X t=240 — slices through largest candidate; masks are not bound-object certificates',fontsize=12)
    fig.savefig(ROOT/'STRUCTURES.png',dpi=165);plt.close(fig)
    np.savez_compressed(ROOT/'VISUALIZATION_FIELDS.npz',**profiles)

    fig,ax=plt.subplots(1,3,figsize=(15,4.7),layout='constrained')
    styles={'mixed_N64':('#8a4fa0','--'),'mixed_N96':('#2463a8','-'),'mixed_N128':('#33865f','-'),
        'mixed_dt_half':('#2463a8',':'),'mixed_seed2':('#777777','--'),'mixed_box64':('#ad7d12','-'),'onesign_N96':('#b75018','-')}
    for tag,case in analysis.items():
        color,ls=styles[tag];rr=case['rows'];t=[r['t'] for r in rr]
        ax[0].plot(t,[r['charge_contrast_RMS'] for r in rr],label=tag,color=color,ls=ls)
        ax[1].plot(t,[r['catalogs'][1]['Q_fraction_in_masks'] for r in rr],color=color,ls=ls)
        if tag in ['mixed_N64','mixed_N96','mixed_N128']:
            ax[2].plot(t,[r['spectra']['boundary_shell_energy_fractions']['0.75'] for r in rr],color=color,ls=ls,label=tag)
    ax[0].set(xlabel='m_X t',ylabel='RMS current contrast / mean current',title='Nonlinear charge modulation')
    ax[0].legend(fontsize=7)
    ax[1].set(xlabel='m_X t',ylabel='Q in threshold-4 masks / net Q',title='Candidate-mask statistic; no relic fraction')
    ax[2].set(xlabel='m_X t',ylabel='Fraction of kinetic + gradient energy',title='Grid-boundary shell: max |k_i| > 0.75 Nyquist')
    ax[2].legend(fontsize=8)
    fig.suptitle('Resolution, time-step, seed and box controls — vacuum classical fields, no injected vector heat',fontsize=12)
    fig.savefig(ROOT/'CONTROLS.png',dpi=165);plt.close(fig)

if __name__=='__main__':main()
