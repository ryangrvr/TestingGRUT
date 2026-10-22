#!/usr/bin/env python3
"""Render the scientific figure for the explicit metastable version."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from metastable import run,potential,ROOT,sc,cfg,m,vev,delta

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                     'axes.spines.right':False,'svg.fonttype':'none','svg.hashsalt':'grut-q-stage2-20261021'})
row,sol=run(cfg['initial_s_values'][1],True)
t=sol.t; N,x,v,h,s,w=sol.y
rho_s=.5*w*w+potential(s)[0]
rd=sc['rho_d_initial']*np.exp(-3*N); rr=sc['rho_r_initial']*np.exp(-4*N)
z=np.linspace(-1.45,1.45,800); curve=potential(z*vev)[0]
fig,axes=plt.subplots(1,2,figsize=(11.4,4.3),gridspec_kw={'width_ratios':[1,1.2]})
ax=axes[0]
ax.plot(z,curve,color='#285c8c',lw=2)
ax.scatter([-1,1],[delta,0],s=50,c=['#bc7035','#285c8c'],zorder=3)
ax.annotate('Occupied false phase\nV = supplied gap',(-1,delta),xytext=(-1.3,.09),
            arrowprops={'arrowstyle':'-','color':'#666'},fontsize=9)
ax.annotate('Eventual true phase\nV = 0 by matching',(1,0),xytext=(.33,.095),
            arrowprops={'arrowstyle':'-','color':'#666'},fontsize=9)
ax.set(xlabel='Scalar field s / v',ylabel='Potential energy',ylim=(-.006,.145),
       title='QV1-M: one potential, two phases')
ax.grid(axis='y',alpha=.15)
ax=axes[1]
for values,label,color in [(rr,'Radiation','#407da8'),(rd,'Dust','#738c5a'),
                            (rho_s,'Scalar energy','#bc7035'),(.75*m*m*x*x,'q potential piece','#8b76a6')]:
    ax.semilogy(t,values,label=label,lw=1.8,color=color)
ax.axhline(delta,color='#bc7035',ls='--',lw=1,label='Supplied gap')
ax.set(xlabel='Jordan proper time',ylabel='Energy density',title='Conserved evolution near the false phase',ylim=(1e-9,1))
ax.grid(axis='y',alpha=.15); ax.legend(frameon=True,facecolor='white',edgecolor='none',framealpha=.95,fontsize=8,loc='upper right')
fig.suptitle('GRUT vacuum-sector construction: explicit inputs and common gravitational response',fontsize=12,y=.99)
fig.text(.02,.012,'Illustrative compressed units; c_q = 0. The gap and preparation are supplied. No tunneling lifetime or observational fit.',fontsize=8,color='#555')
fig.tight_layout(rect=(0,.055,1,.94))
fig.savefig(ROOT/'FIGURE.svg',bbox_inches='tight',metadata={'Date':None})
svg_path=ROOT/'FIGURE.svg'
svg_path.write_text('\n'.join(line.rstrip() for line in svg_path.read_text().splitlines())+'\n')
fig.savefig('/tmp/grut_stage2_figure.png',dpi=160,bbox_inches='tight')
print('Rendered FIGURE.svg and a scratch PNG for visual inspection.')
