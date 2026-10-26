"""Scientific figure from recorded spectra; no new physics computation."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
from background import ROOT


def main():
    v=json.loads((ROOT/'VECTOR_RESULTS.json').read_text())
    c=json.loads((ROOT/'COUPLED_BAND.json').read_text())
    b=json.loads((ROOT/'BUDGET_RESULTS.json').read_text())
    unit=b['rows'][0]['deposit_R_per_vector']
    fine=v['modes'][:97]
    xx=np.array([r['k_over_kstar'] for r in fine]);yy=np.array([r['dN_dlnk'] for r in fine])*unit
    cx=np.array([r['k_over_kstar'] for r in c['points']]);cy=np.array([r['dN_dlnk'] for r in c['points']])*unit
    integral=cumulative_trapezoid(yy,np.log(xx),initial=0)
    fig,axs=plt.subplots(2,1,figsize=(8,6.7),sharex=True)
    axs[0].plot(xx,yy,color='#12616c',lw=2.2,label='Actual-norm Proca, 97 momenta')
    axs[0].plot(cx,cy,'o',ms=3.3,mfc='none',color='#b35634',label='Coupled common mode, 49 momenta')
    axs[0].set_ylabel(r'$dR_A/d\ln k$')
    axs[0].set_title('A primordial-vector band exceeds the corrected budget',loc='left',fontsize=12)
    axs[0].legend(frameon=False,fontsize=9)
    axs[1].plot(xx,integral,color='#12616c',lw=2.2,label='Cumulative vector band')
    axs[1].axhline(8.2,color='#9c4945',ls='--',lw=1.5,label='ke threshold, about 8.2')
    axs[1].axhline(20,color='#8e793e',ls=':',lw=1.7,label='ke + scalar route, upper bracket 20')
    axs[1].set_ylabel('Cumulative injection R')
    axs[1].set_xlabel(r'Momentum $k/k_*$')
    axs[1].set_ylim(0,45)
    axs[1].legend(frameon=False,fontsize=9,loc='upper left')
    for ax in axs:
        ax.set_xscale('log');ax.grid(True,alpha=.2);ax.set_xlim(.2,5)
    fig.text(.12,.017,r'$e=10^{-7}$; standard inflationary common-mode vacuum. Thresholds are equilibrium proxies.',fontsize=9)
    fig.tight_layout(rect=[0,.045,1,1])
    fig.savefig(ROOT/'vector_budget.png',dpi=180)
    plt.close(fig)


if __name__=='__main__':main()
