"""Measured nonlinear fields, control comparisons and basis dependence."""
from en_common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    curves=np.load(ROOT/'FIELD_CURVES.npz');result=read('RESULTS.json');control=read('LINEAR_PHI_CONTROL.json')
    fig,ax=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for b,color in [(3,'#245c9c'),(5,'#bc651b')]:
        d=np.load(ROOT/f'dense_b{b}.npz');mask=d['times']<=160
        ax[0,0].semilogy(d['times'][mask],.03*d['maximum_amplitude_per_radian'][mask],color=color,label=f'({b},'+('.5)' if b==3 else '1)'))
    ax[0,0].axhline(.1,color='grey',ls=':');ax[0,0].axhline(1,color='k',ls='--')
    ax[0,0].set(xlabel=r'$\tau=H_i t$',ylabel='Linear radial modulation / current radius',title='Dense preparation-phase control (.03 rad)',ylim=(1e-4,3));ax[0,0].legend(fontsize=8)
    for tag,style in [('N32',':'),('N48','-'),('N64','--')]:
        ax[0,1].plot(curves[tag+'_tau'],curves[tag+'_energy_outside_spatial_mean_fraction'],style,label=tag)
    ax[0,1].set(xlabel=r'$\tau$',ylabel='Fraction of active field energy',title='Energy outside the coherent spatial mean',ylim=(0,1));ax[0,1].legend(fontsize=8)
    for tag in ['N48','N64']:
        ax[1,0].semilogy(curves[tag+'_tau'],curves[tag+'_phi_plus_interaction_energy_fraction'],label=tag+' nonlinear')
    ax[1,0].semilogy([r['tau'] for r in control['rows']],[r['phi_plus_interaction_raw_fraction'] for r in control['rows']],'k--',label='Fixed X, linear phi extrapolation')
    ax[1,0].axhline(1,color='grey',ls=':');ax[1,0].set(xlabel=r'$\tau$',ylabel=r'$(E_\phi+E_{\rm mix})/E_{\rm initial}$',title='Same phi seed; different dynamical continuation',ylim=(1e-4,3));ax[1,0].legend(fontsize=7)
    d=np.load(ROOT/'VIEW_FIELDS.npz');field=d['charge'];z=int(np.unravel_index(np.argmax(abs(field)),field.shape)[2]);limit=np.quantile(abs(field[:,:,z]),.97)
    im=ax[1,1].imshow(field[:,:,z]/field.mean(),origin='lower',cmap='RdBu_r',vmin=-limit/field.mean(),vmax=limit/field.mean(),extent=(0,12,0,12))
    ax[1,1].set(xlabel=r'$H_i x$',ylabel=r'$H_i y$',title='Block-averaged final classical current / mean');fig.colorbar(im,ax=ax[1,1],shrink=.75)
    for a in ax.flat:
        if a is not ax[1,1]:a.grid(alpha=.2)
    fig.suptitle('Matched, regulated X/phi scalar-limit pilot — not a full quantum source',fontsize=13)
    fig.savefig(ROOT/'NONLINEAR.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(11,4.1),layout='constrained');names=['N32','N48','N64','seed2_48','cut10_48','box16_64'];xx=np.arange(len(names));w=.36
    a=[result['final'][tag]['raw_nonzero_reference_nminus_over_net'][0] for tag in names];b=[result['final'][tag]['raw_nonzero_reference_nminus_over_net'][1] for tag in names]
    ax[0].bar(xx-w/2,a,w,label='Massless reference');ax[0].bar(xx+w/2,b,w,label='Isotropic Hartree reference')
    ax[0].set(xticks=xx,xticklabels=names,ylabel='Raw reference count / net current',title='Reference projections; not physical antiparticle counts',ylim=(0,3));ax[0].tick_params(axis='x',rotation=30);ax[0].legend(fontsize=8)
    f=read('FIELD_RESULTS.json')['cases']
    for tag in ['N32','N48','N64']:
        sp=f[tag]['spectra'];edges=np.array(sp['kinetic_gradient_k_edges']);cent=(edges[:-1]+edges[1:])/2;vals=np.array(sp['kinetic_gradient_spectrum'])
        ax[1].semilogy(cent[1:],vals[1:]/np.diff(edges)[1:],label=tag)
    ax[1].set(xlabel=r'$k/H_i$',ylabel='Kinetic + gradient energy per k',title='Final spectrum; zero-mode bin omitted');ax[1].legend(fontsize=8);ax[1].grid(alpha=.2)
    fig.savefig(ROOT/'REFERENCES.png',dpi=150);plt.close(fig)

if __name__=='__main__':main()
