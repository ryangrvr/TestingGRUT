#!/usr/bin/env python3
"""Scientific figures and concise metrics from preserved execution records."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent


def main():
    raw=json.loads((HERE/'results_initial.json').read_text())
    challenge=json.loads((HERE/'challenge_results.json').read_text())
    rows=[]
    for name,x in raw['runs'].items():
        ts=np.array(x['retention']['times']);ps=np.array(x['retention']['success'])
        rows.append(dict(run=name,L=x['L'],minimum_memory_t10=float(ps[ts<=10].min()),
              memory_t40=float(ps[-1]),best_readout=x['readout']['best_success'],
              readout_time=x['readout']['best_registered_time'],memory_at_readout=x['readout']['memory_success_at_best']))
    summary=dict(core=rows,registered_core_conditions=raw['checks'],
          adversarial_conditions=challenge['checks'],
          source_of_core_times='dimensionless Hamiltonian time, hbar/J=1',
          lyapunov_exponent=None,gravitational_prediction=None)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(2,2,figsize=(12,8.5),constrained_layout=True)
    for name,label,color in [('candidate_L4','M1, L=4','#156a8a'),('candidate_L5','M1, L=5','#54a3b0'),('unfacilitated_L4','Saturation removed','#b44f49')]:
        x=raw['runs'][name]['retention'];axs[0,0].plot(x['times'],x['success'],label=label,color=color)
    axs[0,0].axhline(.9,color='#777',ls=':',lw=1)
    axs[0,0].set(title='Prepared repetition bit: local read succeeds',xlabel='Time (hbar/J)',ylabel='Equal-prior guessing probability',ylim=(.48,1.02))
    axs[0,0].legend(fontsize=9)
    for name,label,color in [('candidate_L4','Feedback present','#156a8a'),('no_feedback_L4','Feedback removed','#777')]:
        x=raw['runs'][name]['readout'];axs[0,1].plot(x['times'],x['success'],label=label,color=color)
    axs[0,1].set(title='A matter-spin measurement reads the memory',xlabel='Time (hbar/J)',ylabel='Guessing probability',ylim=(.48,.84))
    axs[0,1].legend(fontsize=9)
    for L,color in [('4','#156a8a'),('5','#54a3b0')]:
        x=challenge['local_bit'][L];axs[1,0].plot(x['times'],x['success'],label='One-cell bit, L='+L,color=color)
    x=challenge['autonomous_write']['4'];axs[1,0].plot(x['times'],x['success'],label='Bare matter-to-memory write, L=4',color='#b44f49')
    axs[1,0].axhline(.9,color='#777',ls=':',lw=1)
    axs[1,0].set(title='Stronger claims fail in the core model',xlabel='Time (hbar/J)',ylabel='Guessing probability',ylim=(.48,1.02))
    axs[1,0].legend(fontsize=8)
    eq=raw['runs']['candidate_L4']['equilibrium']
    axs[1,1].plot(eq['memory']['times'],eq['memory']['autocorrelation'],label='Thermal memory (energy removed)',color='#156a8a')
    o=eq['otoc'][0];axs[1,1].plot(o['times'],o['squared_commutator'],label='Matter commutator, adjacent sites',color='#b57728')
    axs[1,1].set(title='Memory and spreading are separate observables',xlabel='Time (hbar/J)',ylabel='Normalized correlation / commutator',xlim=(0,4))
    axs[1,1].legend(fontsize=8)
    fig.suptitle('GRUT-M1: a constructed closed memory model — finite-system results',fontsize=15)
    fig.savefig(HERE/'core_summary.png',dpi=160)
    plt.close(fig)
    writer_path=HERE/'writer_results_refined.json'
    if writer_path.exists():
        writer=json.loads(writer_path.read_text())
        if 'checks' in writer:
            x=writer['runs'][-1]
            summary['writer_checks']=writer['checks']
            summary['writer_refinement']=writer['refinement']
            summary['writer_minimum_memory_t10']=x['minimum_all_cell_success_through_t10']
            summary['writer_best_readout']=x['readout']['best_success']
            fig,axs=plt.subplots(1,2,figsize=(11,4.3),constrained_layout=True)
            for i,line in enumerate(x['memory_success']):
                axs[0].plot(x['minus']['storage']['times'],line,label='Memory cell '+str(i))
            axs[0].axhline(.9,color='#777',ls=':',lw=1)
            axs[0].set(title='After the internal writer: unconditional record',xlabel='Storage time after scattering',ylabel='Guessing probability',ylim=(.85,1.005))
            axs[0].legend(fontsize=8)
            axs[1].plot(x['readout']['times'],x['readout']['success'],color='#156a8a')
            axs[1].set(title='Matter response reads the newly written bit',xlabel='Storage time after scattering',ylabel='Guessing probability',ylim=(.48,.84))
            fig.suptitle('GRUT-M1-W: scattering, followed by storage with residual coupling neglected',fontsize=12)
            fig.savefig(HERE/'writer_summary.png',dpi=160)
            plt.close(fig)
    (HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
