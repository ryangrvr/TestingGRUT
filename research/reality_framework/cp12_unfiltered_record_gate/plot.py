"""Scientific plots from saved controlled intervals and exact capacity bounds."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=Path(__file__).resolve().parent


def errorbar(ax,x,y,interval,color,label=None,marker='o',scale=1):
    lo,hi=interval
    ax.errorbar(x*scale,y,xerr=np.array([[max(0,x-lo)],[max(0,hi-x)]])*scale,
                fmt=marker,color=color,capsize=4,ms=6,label=label,lw=1.5)


def run():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'svg.fonttype':'none','savefig.facecolor':'white'})
    fig,axes=plt.subplots(2,2,figsize=(12.4,8.7),gridspec_kw={'hspace':.50,'wspace':.40})
    blue='#145c92';orange='#bd5d19';gray='#6b7280';green='#18775e'
    main4=json.loads((P/'UNFILTERED_N4_Q4_B256_W10_L0.5.json').read_text())
    main5=json.loads((P/'UNFILTERED_N5_Q5_B128_W10_L0.5.json').read_text())
    refs=json.loads((P/'ENERGY_REFERENCE_RESPONSE.json').read_text())
    a=axes[0,0]
    cases=[main4['no_initial_filter'],main5['no_initial_filter']]+[r['response'] for r in refs['rows']]
    labels=['N4 neutral evolution','N5 neutral evolution','N4 spectral reference','N4 max-entropy reference']
    for j,r in enumerate(cases):
        d=r['susceptibilities']['density'];c=r['density_response_inside_final_circulation']
        errorbar(a,d['value'],3-j+.12,d['interval'],blue,'Plain density' if j==0 else None,scale=1e3)
        errorbar(a,c['value'],3-j-.12,c['interval'],orange,'Within final circulation' if j==0 else None,marker='s',scale=1e3)
    a.set_yticks(range(4),labels[::-1]);a.set_xlim(0,10.6);a.set_ylim(-.5,3.7)
    a.set_xlabel('Susceptibility × 10³');a.set_title('A   Positive response; no initial loop filter',loc='left',fontweight='bold',pad=13)
    a.axvline(0,color=gray,lw=.8);a.legend(frameon=False,loc='upper left',bbox_to_anchor=(.02,1.02),fontsize=9)
    a.grid(axis='x',alpha=.15)
    a.text(.02,-.30,'Wait 10, lag 0.5 for main states. Final conditioning is a readout.',transform=a.transAxes,fontsize=8.7,color=gray)
    b=axes[0,1];co=json.loads((P/'SHORT_LAG_CONTROLS.json').read_text())['main_state_coefficients']
    for j,r in enumerate(co):
        value=r['projected_circulation_weighted_density_T2_coefficient']
        errorbar(b,value,2-j,r['full_law_circulation_weighted_density_T2_coefficient_interval'],green,scale=1e4)
    b.plot([0,0],[0,-1],'x',color=gray,ms=7)
    b.set_yticks([-1,0,1,2],['N4 max-entropy reference','N4 spectral reference','N5 neutral evolution','N4 neutral evolution'])
    b.set_xlim(-.12,2.75);b.set_ylim(-1.5,2.6);b.axvline(0,color=gray,lw=.8);b.grid(axis='x',alpha=.15)
    b.set_xlabel('Unforced phase–circulation coefficient α × 10⁴')
    b.set_title('B   Unforced phase–circulation correlation',loc='left',fontweight='bold',pad=13)
    b.text(.02,-.30,'Real references: α = 0. No dimension is inferred.',transform=b.transAxes,fontsize=8.7,color=gray)
    c=axes[1,0]
    for j,(m,col,off) in enumerate([(main4,blue,.11),(main5,orange,-.11)]):
        for k,r in enumerate(m['record_comparisons']):
            z=r['response_interference']['density']
            errorbar(c,z['value'],1-k+off,z['interval'],col,'N4' if j==0 and k==0 else 'N5' if j==1 and k==0 else None,scale=1e5)
    c.set_yticks([0,1],['Circulation record','Occupied-flux record']);c.set_ylim(-.45,1.55)
    c.axvline(0,color=gray,lw=.8);c.set_xlim(-12,12);c.grid(axis='x',alpha=.15)
    c.set_xlabel('Coherent − recorded density response × 10⁵')
    c.set_title('C   Recording alters the N4 response',loc='left',fontweight='bold',pad=13);c.legend(frameon=False,loc='upper right')
    c.text(.02,-.30,'N5 intervals include zero. Recording is a separate intervention.',transform=c.transAxes,fontsize=8.7,color=gray)
    d=axes[1,1];n=np.arange(4,17)
    d.plot(n,2*n,color=blue,lw=2,label='Spanning circulation: B ≥ 2N')
    d.plot(n,4*n,color=green,lw=2,label='Full degree-four circulation: B ≥ 4N')
    d.scatter([6],[10],marker='x',s=65,color=orange,zorder=5,label='Existing dense N6 budget 10')
    d.set_xlabel('Vertices N');d.set_ylabel('Necessary electric-square budget B')
    d.set_ylim(0,70);d.set_xlim(3.6,16.4);d.set_xticks([4,6,8,12,16]);d.grid(alpha=.15)
    d.set_title('D   Electric capacity for spanning loops',loc='left',fontweight='bold',pad=13)
    d.legend(frameon=False,fontsize=8.6,loc='upper left')
    d.text(.02,-.30,'Neutral resource bounds only; matter raises the requirements.',transform=d.transAxes,fontsize=8.7,color=gray)
    fig.suptitle('CP-12: unfiltered relational response and the next formation gate',x=.045,ha='left',y=.99,fontsize=17,fontweight='bold')
    fig.text(.045,.943,'Unchanged Hamiltonian and supplied neutral input • analytic electric-error bounds evaluated in floating point',fontsize=10.5,color=gray)
    fig.text(.045,.014,'Finite microscopic results. Autonomous extended geometry, dimensional growth and a causal cone remain unestablished.',fontsize=10,color=gray)
    fig.subplots_adjust(top=.87,bottom=.13,left=.185,right=.98)
    for ext in ['png','svg']:fig.savefig(P/f'unfiltered_record_gate.{ext}',dpi=180)
    plt.close(fig)


if __name__=='__main__':run()
