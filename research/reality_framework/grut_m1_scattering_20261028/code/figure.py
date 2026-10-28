"""Scientific plots of cold scattering and explicitly conditional rates."""
from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    d=np.load(ROOT/'PARTIAL_WAVES.npz');k=d['k'];threshold=np.sqrt((np.sqrt(6)-OMEGA)**2-1)
    fig,ax=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    ax[0,0].plot(k,d['P_antiparticle'][:,0],label='Reflected antiparticle',color='#155b89')
    ax[0,0].plot(k,d['P_X'][:,0],label='Fast X',color='#d45b2b')
    ax[0,0].plot(k,d['P_phi'][:,0],label='Neutral scalar',color='#34805b')
    ax[0,0].axvline(threshold,color='black',ls=':',label='Neutral channel opens')
    ax[0,0].set(xlim=(.1,3.3),ylim=(0,1),title='Monochromatic monopole channels',xlabel='Incoming momentum / mX',ylabel='Outgoing number-flux probability')
    ax[0,0].legend(fontsize=8)
    resolved=np.load(ROOT/'RESONANCE_DATA.npz');kr=resolved['k'];outside=(k<kr[0])|(k>kr[-1])
    kk=np.r_[k[outside],kr];ss=np.r_[d['sigma_X_over_mXminus2'][outside],resolved['h005_sigma']];order=np.argsort(kk)
    ax[0,1].semilogy(kk[order],ss[order],label='Forward antiparticle to X conversion',color='#d45b2b')
    positive=d['sigma_phi_over_mXminus2']>0
    ax[0,1].semilogy(k[positive],d['sigma_phi_over_mXminus2'][positive],label='Neutral scalar: pair returned after decay',color='#34805b')
    ax[0,1].set(title='Angular sum through l = 32',xlabel='Incoming momentum / mX',ylabel='Cross section × mX²')
    ax[0,1].legend(fontsize=8)
    src=np.load(ROOT/'SOURCE_SPECTRA.npz')
    for tag,label,color in [('radial_fine','5% pulse','#155b89'),('pulse10','10% pulse','#d45b2b'),('pulse20','20% pulse','#34805b')]:
        q=src[tag+'_k'];w=np.maximum(src[tag+'_negative_number_weight'],0)
        keep=q<2;density=w/np.maximum(np.gradient(q),1e-30)/w.sum()
        ax[1,0].plot(q[keep],density[keep],label=label,color=color)
    ax[1,0].axvline(threshold,color='black',ls=':');ax[1,0].set(title='Measured outgoing antiparticle spectra',xlabel='Momentum / mX',ylabel='Normalized number density in momentum',xlim=(.2,1.6));ax[1,0].legend(fontsize=8)
    p=json.loads((ROOT/'PACKET_RESULTS.json').read_text())['rows']
    chosen=[next(r for r in p if r['tag']==tag) for tag in ('packet04','packet10_fine','packet15')]
    x=np.arange(3)
    ax[1,1].bar(x-.15,[r['measured_P_X_anti_phi'][0] for r in chosen],width=.3,label='Nonlinear packet',color='#d45b2b')
    ax[1,1].bar(x+.15,[r['stationary_packet_weighted_P_X_anti_phi'][0] for r in chosen],width=.3,label='Weighted linear prediction',color='#155b89')
    ax[1,1].set_xticks(x,['k=0.4','k=1.0 (finer)','k=1.5']);ax[1,1].set(title='Separate nonlinear evolution checks conversion',ylabel='Fast-X probability');ax[1,1].legend(fontsize=8)
    for a in ax.flat:a.grid(alpha=.15)
    fig.suptitle('Fixed GRUT-M1: forward channels; inverse reactions also required',fontsize=15)
    fig.savefig(ROOT/'SCATTERING.png',dpi=160);plt.close(fig)
    print('Rendered SCATTERING.png')

if __name__=='__main__':main()
