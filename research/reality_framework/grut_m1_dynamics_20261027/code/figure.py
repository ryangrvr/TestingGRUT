"""Standalone scientific plots of measured finite-time M1 controls."""
from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analysis_io import read_analysis

def main():
    data=read_analysis()['cases']
    fig,axs=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    colors=['#155b89','#d45b2b','#34805b']
    for tag,title,color in zip(('object64','object64_perturbed','object96_perturbed'),('Unperturbed 64³','5% disturbance 64³','5% disturbance 96³'),colors):
        rows=data[tag]['rows'];t=[r['t'] for r in rows]
        axs[0,0].plot(t,[r['central_spheres'][1]['Q']/200000 for r in rows],label=title,color=color)
        axs[0,1].plot(t,[r['central_spheres'][1]['E_over_Q'] for r in rows],label=title,color=color)
    axs[0,0].set(title='Supplied object retains its charge',ylabel=r'$Q(r<12/m_X)/Q_0$')
    axs[0,1].axhline(1,color='black',ls=':',label='Free-particle rest energy')
    axs[0,1].set(title='Core energy stays below free-particle threshold',ylabel=r'$E(r<12/m_X)/(m_X Q)$')
    axs[0,0].legend(fontsize=9);axs[0,1].legend(fontsize=8)
    for tag,title,color in zip(('cold64','dilute64','ellipse64'),('Circular F=1','Circular F=0.5','Elliptic F=1'),colors):
        rows=data[tag]['rows'];t=[r['t'] for r in rows]
        axs[1,0].plot(t,[r['charge_contrast'] for r in rows],color=color,label=title)
        axs[1,1].plot(t,[r['catalogs'][1]['mask_charge_fraction'] for r in rows],color=color,label=title)
    axs[1,0].set(title='Nonlinear density structure develops',ylabel='Charge density standard deviation / mean')
    mask_max=max(r['catalogs'][1]['mask_charge_fraction'] for tag in ('cold64','dilute64','ellipse64') for r in data[tag]['rows'])
    axs[1,1].set(title='Dense-mask charge is not final object abundance',ylabel='Signed net-charge fraction in catalog masks',ylim=(0,max(1,1.05*mask_max)))
    axs[1,0].legend(fontsize=9);axs[1,1].legend(fontsize=9)
    for ax in axs.flat:ax.set_xlabel(r'$m_X t$');ax.grid(alpha=.18)
    fig.suptitle('GRUT-M1: conservative local scalar stress tests',fontsize=15)
    fig.savefig(ROOT/'DYNAMICS.png',dpi=160);plt.close(fig)
    v=np.load(ROOT/'FINAL_PROJECTIONS.npz');fig,axs=plt.subplots(3,2,figsize=(11,12),layout='constrained')
    for row,(tag,title) in enumerate([('cold64','Circular F=1'),('dilute64','Circular F=0.5'),('ellipse64','Elliptic F=1')]):
        _,L,t=v[tag+'_geometry'];mean=data[tag]['rows'][-1]['Q']/L**3
        a=v[tag+'_charge_max_projection'];b=v[tag+'_phi_min_projection_over_v']
        im=axs[row,0].imshow(np.log10(1+np.maximum(a,0)/mean).T,origin='lower',extent=(0,L,0,L),cmap='magma',vmin=0,vmax=1.7)
        fig.colorbar(im,ax=axs[row,0],label=r'$\log_{10}(1+\max_z j^0/\langle j^0\rangle)$')
        im=axs[row,1].imshow(b.T,origin='lower',extent=(0,L,0,L),cmap='viridis',vmin=-1 if tag=='ellipse64' else 0,vmax=1)
        fig.colorbar(im,ax=axs[row,1],label=r'$\min_z\phi/v$')
        axs[row,0].set_title(f'{title}; charge, $m_X t={t:g}$');axs[row,1].set_title('Binding-field depletion')
        if tag=='ellipse64':axs[row,1].set_title('N64 hot stress control; spatially unresolved')
        for ax in axs[row]:ax.set(xlabel=r'$m_X x$',ylabel=r'$m_X y$')
    fig.suptitle('Clumps and waves depend on the supplied preparation',fontsize=15)
    fig.savefig(ROOT/'STRUCTURE.png',dpi=150);plt.close(fig)
    print('Rendered DYNAMICS.png and STRUCTURE.png')

if __name__=='__main__':main()
