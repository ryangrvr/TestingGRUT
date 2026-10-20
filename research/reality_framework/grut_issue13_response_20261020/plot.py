#!/usr/bin/env python3
"""Scientific figure from the six full-grid per-step diagnostic histories."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

P=Path(__file__).resolve().parent
plt.rcParams.update({'font.size':10,'svg.fonttype':'none','svg.hashsalt':'grut-wall-diagnostics-v1'})
fig,axes=plt.subplots(1,2,figsize=(10,3.9),layout='constrained')
for prefix,ax,hs in [('pair_generic',axes[0],[.75,.5,.375]),('production',axes[1],[1,2/3,.5])]:
    for tag,h in zip(['coarse','medium','fine'],hs):
        a=np.genfromtxt(P/'runs'/f'{prefix}_{tag}.csv.walls.csv',delimiter=',',names=True)
        ax.plot(a['t'],a['minS'],label=f'h = {h:.3g}',lw=1)
    ax.set_xlabel(r'Time $\tau=m_X t$')
    ax.set_ylabel(r'Full-grid minimum of $\phi/v$')
    ax.grid(alpha=.2)
    ax.legend(frameon=False,loc='upper right')
axes[0].set_title('Same collision, three spatial grids')
axes[0].set_ylim(0,.061)
axes[0].axhline(0,color='black',lw=.8)
axes[1].set_title('Same nonthermal production seed')
axes[1].set_ylim(-1.08,1.02)
for value,label in [(0,'sign change'),(-.5,'loose opposite-branch threshold'),(-1,'negative vacuum')]:
    axes[1].axhline(value,color='black',lw=.8,ls='--',alpha=.65)
    axes[1].text(4,value+(.025 if value==-1 else -.025),label,fontsize=8,va='bottom' if value==-1 else 'top')
fig.savefig(P/'WALL_DIAGNOSTICS.svg',metadata={'Date':None})
svg=P/'WALL_DIAGNOSTICS.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
fig.savefig('/tmp/grut_issue13_wall_diagnostics.png',dpi=170)
