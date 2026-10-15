"""Standalone scientific figure from published numerical interval records."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=Path(__file__).resolve().parent
BLUE='#225ea8';ORANGE='#d95f0e';GREEN='#238b45';GRAY='#555555'


def load(name):return json.loads((P/name).read_text())


def interval_errorbar(ax,x,values,intervals,**kwargs):
    values=np.asarray(values);intervals=np.asarray(intervals,float)
    # A central estimate can be outside an interval only for an unresolved
    # denominator. Do not invent an errorbar for that case.
    assert np.all(intervals[:,0]<=values) and np.all(values<=intervals[:,1])
    return ax.errorbar(x,values,yerr=np.array([values-intervals[:,0],intervals[:,1]-values]),
        capsize=4,elinewidth=1.5,**kwargs)


def main():
    n5=load('N5_B50.json');n6=load('N6_SHELL2.json');rare=load('N5_RARE_PAULI.json')
    old=load('../cp12_shortcut_dynamics_20261015/SPLIT_N6_B12.json')
    adapt=load('N6_ADAPT225K_PILOT.json')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
        'axes.spines.right':False,'svg.hashsalt':'grut-dynamical-suppression-20261015'})
    fig,axs=plt.subplots(2,2,figsize=(12.4,9.0))
    fig.suptitle('CP-12: full-law dynamical suppression tests',fontsize=17,fontweight='bold',y=.98)
    fig.text(.5,.94,'Original neutral preparation · no hopping deletion · intervals include omitted-state error',
        ha='center',fontsize=10,color=GRAY)
    ax=axs[0,0]
    for data,offset,color in [(n5,-.025,BLUE),(n6,.025,ORANGE)]:
        waits=[];central=[];ints=[]
        for row in data['rows']:
            if row['wait'] not in [.5,1.]:continue
            probe=next(z for z in row['probes'] if z['lag']==.25)
            if probe['eta_interval'][1] is None:continue
            waits.append(row['wait']+offset);central.append(probe['eta']);ints.append(probe['eta_interval'])
        interval_errorbar(ax,waits,central,ints,fmt='o',color=color,label=f"N={data['N']}, Q={data['Q']}")
    ax.axhline(1,color=GRAY,lw=1,ls='--');ax.set_xticks([.5,1.]);ax.set_yscale('log')
    ax.set(xlabel='Neutral waiting time',ylabel=r'Action ratio $\eta$ at delay 0.25',
        title='A. Matched early-time finite-size test')
    ax.legend(frameon=False);ax.set_xlim(.3,1.2)
    ax.text(.02,.02,'Log scale · later N6 comparison has a wide error interval',
        transform=ax.transAxes,fontsize=8,color=GRAY)
    ax=axs[0,1]
    row=next(r for r in n5['rows'] if r['wait']==2.)
    probe=next(r for r in row['probes'] if r['lag']==.25)
    for offset,kind,color in [(-.17,'coherent',BLUE),(.17,'diagonal reference',ORANGE)]:
        names=['outside','inside'] if kind=='coherent' else ['dephased_outside_action','dephased_inside_action']
        central=[probe['joint_integrals'][k]['value']/row['formation_probability'] for k in names]
        ints=[probe['conditional_integrals'][k] for k in names]
        ax.bar(np.arange(2)+offset,central,width=.3,color=color,label=kind,alpha=.9)
        interval_errorbar(ax,np.arange(2)+offset,central,ints,fmt='none',ecolor='black')
    ax.set_xticks([0,1],['Outside','Inside']);ax.set(ylabel='Conditional action integral',
        title='B. Constructive interference at N=5')
    ax.legend(frameon=False,fontsize=9);ax.text(.02,.98,'Wait 2 · delay 0.25 · formation ≈5.963%',
        transform=ax.transAxes,va='top',fontsize=9,color=GRAY)
    ax.set_ylim(0,5.7)
    ax=axs[1,0]
    rs=[r for r in rare['rows'] if r['lag']<=.5]
    ts=[r['lag'] for r in rs]
    interval_errorbar(ax,ts,[r['eta'] for r in rs],[r['eta_interval'] for r in rs],
        fmt='o-',color=GREEN,label='Rare filled-ring branch')
    small=np.linspace(0,.25,100)
    ax.plot(small,2*small**2/3,ls='--',color=GRAY,label=r'Leading small-time term $2T^2/3$')
    ax.set(xlabel='Delay after formation',ylabel=r'Action ratio $\eta$',
        title='C. Real Pauli suppression, rare formation')
    ax.legend(frameon=False,fontsize=8)
    ax.text(.03,.97,'Formation ≈7×10⁻⁸; immediate density phase response =0',
        transform=ax.transAxes,va='top',fontsize=8,color=GRAY)
    ax.set_ylim(0,.14)
    ax=axs[1,1]
    for wait,color in [(.5,BLUE),(1.,ORANGE)]:
        before=next(r['full_state_error'] for r in old['trajectory'] if r['time']==wait)
        middle=next(r['waiting_state_error'] for r in adapt['rows'] if r['wait']==wait)
        after=next(r['waiting_state_error'] for r in n6['rows'] if r['wait']==wait)
        ax.plot(np.array([old['dimension'],adapt['dimension'],n6['dimension']])/1e6,
            [before,middle,after],'o-',color=color,label=f'Wait {wait:g}')
    ax.set_yscale('log');ax.set(xlabel='Retained six-site orbits (millions)',
        ylabel='Full-state error upper bound',title='D. Six-site outgoing-boundary improvement')
    ax.legend(frameon=False);ax.grid(axis='y',alpha=.2)
    fig.subplots_adjust(left=.075,right=.98,bottom=.09,top=.865,hspace=.42,wspace=.28)
    fig.text(.5,.025,'A–B: primary quarter-threshold event. C: exploratory postselection. D: same original neutral ray and full H.',
        ha='center',fontsize=9,color=GRAY)
    target=P/'dynamical_suppression.svg'
    fig.savefig(target,metadata={'Date':None,'Creator':'GRUT CP-12 numerical research'})
    target.write_text('\n'.join(line.rstrip() for line in target.read_text().splitlines())+'\n')
    print(str(target))


if __name__=='__main__':main()
