"""Static scientific figure: rigorous sector bounds and supplied-host controls."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bounds import vacuum_bound,low_matter_bound

P=Path(__file__).resolve().parent


def run():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'svg.fonttype':'none','savefig.facecolor':'white'})
    fig,(a,b)=plt.subplots(1,2,figsize=(12.6,5.5));blue='#145c92';orange='#b75b18';green='#18775e'
    n=np.unique(np.geomspace(32,1048576,150).astype(int))
    a.loglog(n,[vacuum_bound(v,2*v,2*v) for v in n],color=blue,lw=2,label='Full circulation, vertex vacuum')
    a.loglog(n,[low_matter_bound(v,2*v,2*v,1) for v in n],color=orange,lw=2,label='Full circulation, nᵥ ≤ 1')
    a.loglog(n,[low_matter_bound(v,2*v,2*v,int(np.sqrt(v))) for v in n],color=green,lw=2,label='Full circulation, nᵥ ≤ √N')
    a.set_ylim(4e-5,1.35);a.set_xlim(32,1048576);a.grid(alpha=.16,which='both')
    a.set_xlabel('Vertices N, with Q = 2N');a.set_ylabel('Uniform-in-time probability upper bound')
    a.set_title('A   Extensive circulation needs extensive matter',loc='left',fontweight='bold',fontsize=11)
    a.legend(frameon=False,fontsize=8.8,loc='lower left')
    a.text(0,-.24,'Analytic upper bounds, not dynamical probabilities.\nCharged sectors with nᵥ of order N remain allowed.',transform=a.transAxes,fontsize=9,color='#666666')
    data=json.loads((P/'MATCHED_HOSTS.json').read_text())['cases']
    labels={'one_dimensional_ring_power':'Supplied 1D ring power','two_dimensional_square':'Supplied 2D square torus','three_dimensional_diamond':'Supplied 3D diamond host'}
    colors={'one_dimensional_ring_power':blue,'two_dimensional_square':orange,'three_dimensional_diamond':green}
    for r in data:
        if r['N']!=11664:continue
        y=r['ball_volume_from_representative_root'];x=np.arange(len(y))
        b.loglog(x[1:],y[1:],marker='o',ms=3,lw=1.8,label=labels[r['supplied_host']],color=colors[r['supplied_host']])
    b.set_xlabel('Intrinsic graph radius r');b.set_ylabel('Vertices in the ball V(r)')
    b.grid(alpha=.16,which='both');b.set_xticks([1,2,4,8,12],['1','2','4','8','12']);b.legend(frameon=False,fontsize=8.8,loc='upper left')
    b.set_title('B   Matched energy resources; different ball growth',loc='left',fontweight='bold',fontsize=11)
    b.text(0,-.24,'All three admit pure states with ⟨H⟩ = N and Var(H) = 4N.\nGraphs and phase programs are supplied, not dynamically selected.',transform=b.transAxes,fontsize=9,color='#666666')
    fig.suptitle('CP-12: the extensive-circulation probability and energy gate',x=.06,ha='left',y=.98,fontsize=16,fontweight='bold')
    fig.text(.06,.905,'Unchanged frozen law • exact analytic theorems with independent finite controls',fontsize=10.5,color='#666666')
    fig.subplots_adjust(left=.075,right=.98,top=.82,bottom=.26,wspace=.30)
    fig.text(.06,.025,'No large-H evolution, autonomous geometry formation or causal-cone claim is made.',fontsize=10,color='#666666')
    for ext in ['png','svg']:fig.savefig(P/f'extensive_circulation_bounds.{ext}',dpi=180)
    plt.close(fig)


if __name__=='__main__':run()
