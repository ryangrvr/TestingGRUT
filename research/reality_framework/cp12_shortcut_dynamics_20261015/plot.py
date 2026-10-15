"""Static scientific figure from the saved full-law error-controlled results."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,FixedFormatter,NullFormatter

P=Path(__file__).resolve().parent


def main():
    d=json.loads((P/'TRANSPORT_N4_Q5_B256.json').read_text())
    yes=[r for r in d['results'] if r['shortcut_present_at_record']]
    no=[r for r in d['results'] if not r['shortcut_present_at_record']]
    ts=np.array([r['lag'] for r in yes])
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','svg.hashsalt':'GRUT_CP12_shortcut_20261015'})
    fig,ax=plt.subplots(1,3,figsize=(12.8,4.4),layout='constrained')
    for rows,model,label,color,marker in [
        (yes,'full_frozen','Full CP-12, direct diagonal','#b9403d','o'),
        (yes,'recorded_ring_hopping_only','Same state, off-ring hopping removed','#245b91','s'),
        (no,'full_frozen','Full CP-12, other diagonal','#757575','^')]:
        rs=[r['models'][model]['receiver_response'] for r in rows]
        vals=np.array([r['conditional_value'] for r in rs]);I=np.array([r['conditional_interval'] for r in rs])
        ax[0].errorbar(ts,vals,yerr=np.maximum(0.,np.array([vals-I[:,0],I[:,1]-vals])),
            marker=marker,color=color,label=label,lw=1.5,capsize=2,ms=4)
    ax[0].set_yscale('symlog',linthresh=1e-7);ax[0].axhline(0,color='#b8bec5',lw=.6)
    ax[0].set_ylabel('Opposite-vertex susceptibility');ax[0].set_title('Off-ring hopping changes response',fontsize=10)
    ax[0].legend(fontsize=7.5,loc='upper left',bbox_to_anchor=(0,-.22),frameon=False)
    for model,label,color in [('full_frozen','Full CP-12','#b9403d'),
        ('recorded_ring_hopping_only','Off-ring hopping removed','#245b91')]:
        vals=np.array([r['models'][model]['initial_ring_survival_conditional'] for r in yes])*100
        I=np.array([r['models'][model]['initial_ring_survival_conditional_interval'] for r in yes])*100
        ax[1].plot(ts,vals,'o-',color=color,label=label,ms=4)
        ax[1].fill_between(ts,I[:,0],I[:,1],color=color,alpha=.18)
    ax[1].set_ylim(55,101);ax[1].set_ylabel('Original-ring survival (%)')
    ax[1].set_title('The ring persists during the effect',fontsize=10)
    ax[1].legend(fontsize=8,loc='upper left',bbox_to_anchor=(0,-.22),frameon=False)
    vals=np.array([r['models']['full_frozen']['new_chord_and_ring_conditional_probability'] for r in yes])
    I=np.array([r['models']['full_frozen']['new_chord_and_ring_conditional_interval'] for r in yes])
    ax[2].errorbar(ts,vals,yerr=np.maximum(0.,np.array([vals-I[:,0],I[:,1]-vals])),
        color='#5e368b',marker='o',label='Full CP-12',ms=4,capsize=2)
    short=np.linspace(.125,.32,60);coef=yes[0]['new_chord_and_ring_probability_T4_coefficient_conditional']
    ax[2].plot(short,coef*short**4,'--',color='#333333',label='Exact local onset: 0.131409 τ⁴')
    ax[2].set_xscale('log');ax[2].set_yscale('log');ax[2].set_ylabel('New diagonal flow AND original ring')
    ax[2].set_title('Direct charge transfer is allowed',fontsize=10)
    ax[2].xaxis.set_major_locator(FixedLocator(ts));ax[2].xaxis.set_major_formatter(FixedFormatter(['0.125','0.25','0.5','1']))
    ax[2].xaxis.set_minor_formatter(NullFormatter())
    ax[2].legend(fontsize=8,loc='upper left',bbox_to_anchor=(0,-.22),frameon=False)
    for a in ax:
        a.set_xlabel('Probe duration τ');a.grid(alpha=.15);a.spines[['top','right']].set_visible(False)
    fig.suptitle('Frozen CP-12, N=4, Q=5 • neutral wait 5 • explicit initial ring record',fontsize=12)
    target=P/'shortcut_transport.svg'
    fig.savefig(target,metadata={'Date':None})
    target.write_text('\n'.join(line.rstrip() for line in target.read_text().splitlines())+'\n')
    plt.close(fig)


if __name__=='__main__':main()
