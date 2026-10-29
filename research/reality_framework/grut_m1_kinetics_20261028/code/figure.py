"""Scientific figure of the measured bounded transport result and its heat."""
from setup import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'svg.hashsalt':'grut-m1-kinetics-20261028'})
    fig,ax=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    colors=['#64748b','#b85c00','#155f9a']
    cases=[('translated_core_64','Core reactions only'),('translated_resonant_64','+ resonant redistribution'),
           ('translated_full_96','+ ordinary elastic scattering')]
    for (tag,label),color in zip(cases,colors):
        d=np.load(ROOT/(tag+'.npz'));n=len(d['E']);t=d['tau'];y=d['Y'];ok=t>0
        nm=y[n:2*n].sum(axis=0);ax[0,0].loglog(t[ok],nm[ok],label=label,color=color,lw=2)
    ax[0,0].set(xlabel=r'Optical time $\tau=n_{\rm core}t/m_X^2$',ylabel='Antiparticles per reference core',title='A. Antiparticle reduction is conditional')
    ax[0,0].legend(fontsize=8,loc='lower left');ax[0,0].grid(alpha=.15)
    d=np.load(ROOT/'translated_full_96.npz');n=len(d['E']);t=d['tau'];y=d['Y'];ok=t>0
    np_=y[:n].sum(axis=0);nm=y[n:2*n].sum(axis=0);q=QREF-np_+nm
    ax[0,1].semilogx(t[ok],q[ok]/QREF*100,color='#155f9a',lw=2,label='Core charge / original charge')
    ax[0,1].set(xlabel=r'Optical time $\tau$',ylabel='Percent',title='B. Additional core charge becomes free X')
    ax[0,1].grid(alpha=.15);ax[0,1].legend(fontsize=8)
    e=d['E'];gas=e@(y[:n]+y[n:2*n]);kin=gas-np_-nm;cost=(1-OMEGA)*np_+(1+OMEGA)*nm
    ax[1,0].semilogx(t[ok],kin[ok],label='Gas kinetic energy',color='#b85c00',lw=2)
    ax[1,0].semilogx(t[ok],cost[ok],label='Gas rest energy minus core energy released',color='#155f9a',lw=2)
    ax[1,0].axhline((kin+cost)[0],color='#111827',ls='--',label='Conserved sum',lw=1.4)
    ax[1,0].set(xlabel=r'Optical time $\tau$',ylabel=r'Energy per reference core / $m_X$',ylim=(0,4250),title='C. Redistribution does not remove heat')
    ax[1,0].legend(fontsize=8);ax[1,0].grid(alpha=.15)
    for target,label,color in [(0,'Initial','#64748b'),(1,'After one optical time','#b85c00'),(1000,'At the stopping time','#155f9a')]:
        j=int(np.argmin(abs(t-target)));w=np.empty(n);w[1:-1]=(e[2:]-e[:-2])/2;w[0]=(e[1]-e[0])/2;w[-1]=(e[-1]-e[-2])/2
        ax[1,1].semilogy(e,y[:n,j]/w,label=label,color=color,lw=1.6)
    ax[1,1].axvline(1+2*OMEGA,color='#9f1239',ls='--',lw=1.2,label='Reverse threshold')
    ax[1,1].set(xlim=(1,4.2),ylim=(1e-3,3e5),xlabel=r'Positive-particle energy $E/m_X$',ylabel=r'$dN_X/d(E/m_X)$',title='D. The positive-particle spectrum remains warm')
    ax[1,1].legend(fontsize=8);ax[1,1].grid(alpha=.15)
    fig.suptitle('GRUT-M1: supplied 10% pulse, fixed-core dilute transport',fontsize=15,fontweight='bold')
    fig.savefig(ROOT/'KINETICS.png',dpi=180);fig.savefig(ROOT/'KINETICS.svg',metadata={'Date':None})
    svg=ROOT/'KINETICS.svg';svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')

if __name__=='__main__':main()
