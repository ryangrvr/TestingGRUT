"""Actual computed response, auxiliary band and source-clock chronology."""
from timing_common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    result=read('RESULTS.json');clocks=read('SOURCE_CLOCKS.json');scan=read('QUARTIC_SCAN.json')
    fig=plt.figure(figsize=(11,8),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[1,1.15])
    ax=fig.add_subplot(gs[0,0]);bx=fig.add_subplot(gs[0,1]);cx=fig.add_subplot(gs[1,:])
    ax.plot([r['k'] for r in scan['rows']],[r['mu_conformal'] for r in scan['rows']],color='k',lw=1.6)
    ax.plot([0,.6],[0,.6*.4398],'--',lw=1,label=r'IR slope $\mu/k\simeq0.440$')
    ax.set(xlabel=r'$k/H_i$',ylabel=r'$\mu$ per $H_i\eta$',title='Auxiliary isolated quartic orbit',xlim=(0,6.3),ylim=(-.01,.68));ax.legend(fontsize=8);ax.grid(alpha=.2)
    colors=['#2462a7','#b65c18']
    for i,row in enumerate(result['rows']):
        z=np.load(ROOT/(row['tag']+'.npz'));time=z['times'];amp=.03*np.max(abs(z['amplitude_per_radian']),axis=1)
        stop=int(np.flatnonzero(amp>=1)[0]);label=f"({row['beta']:g},{row['Gamma_over_H']:g})"
        bx.semilogy(time[1:stop+1],amp[1:stop+1],color=colors[i],label=label)
    bx.axhline(.1,color='grey',ls=':',lw=1);bx.axhline(1,color='k',ls='--',lw=1)
    bx.set(xlabel=r'$\tau=H_i t$',ylabel='Maximum amplitude modulation',title='Supplied .03-rad phase; stop at unity',ylim=(1e-4,2));bx.grid(alpha=.2);bx.legend(title=r'$(m_I,\Gamma_I)/H_i$',fontsize=8,title_fontsize=8)
    labels=['X response reaches1 (.03 rad)','Y_L decay','A decay, e=1e-6','A decay, e=1e-7','Reference charge density']
    for i,row in enumerate(clocks['histories']):
        shift=(-.1 if i==0 else .1);color=colors[i]
        source=next(r for r in result['rows'] if r['tag']==row['tag']);event=source['supplied_phase_crossings'][0]['amplitude_unity']
        cx.plot(event['a_bracket'][1],shift,'o',color=color)
        for j,key in enumerate(['YL','A_e1e-6','A_e1e-7'],1):
            quant=row['sources'][key]['number_quantiles'];lo,mid,hi=quant[1]['a'],quant[2]['a'],quant[3]['a']
            cx.plot([lo,hi],[j+shift,j+shift],color=color,lw=3);cx.plot(mid,j+shift,'o',color=color)
        cx.plot(row['a_at_old_reference_charge_density'],4+shift,'D',color=color,fillstyle='none')
    cx.set(xscale='log',xlim=(10,2e13),yticks=range(5),yticklabels=labels,xlabel='Scale factor a, with a=1 at the supplied reheating start',title='Time ordering on each unperturbed quiet history')
    cx.invert_yaxis();cx.grid(axis='x',alpha=.25)
    cx.text(.01,-.25,'Decay bars: 10–90% by number; dots: medians. Reference density is not a formation prediction.',transform=cx.transAxes,fontsize=9)
    fig.suptitle('Conditional nonlinear-response warning precedes prescribed decay heating',fontsize=13)
    fig.savefig(ROOT/'TIMING.png',dpi=150);plt.close(fig)

if __name__=='__main__':main()
