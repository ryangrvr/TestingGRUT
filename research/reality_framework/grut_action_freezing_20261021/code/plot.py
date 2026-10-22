"""Optional figure from archived, executed trajectories; matplotlib 3.10.8."""
from pathlib import Path
import io,tarfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
p=Path(__file__).resolve().parents[1]
with tarfile.open(p/'RAW_DATA.tar.gz','r:gz') as archive:
 data={name:np.genfromtxt(io.BytesIO(archive.extractfile(name).read()),delimiter=',',names=True) for name in archive.getnames()}
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(1,2,figsize=(11.5,4.4),layout='constrained')
for name,label,col in [('original','KP0: Higgs-controlled','#176a94'),('completed','KP1: product derivative','#be5033')]:
 d=data[name+'.csv'];ax[0].semilogy(d['N'],d['h_over_v'],label=label,color=col,lw=2)
ax[0].set(xlabel='Expansion since release, N',ylabel='Higgs field h/v',ylim=(5e-7,2),title='Full action, compressed zero-temperature example');ax[0].legend(loc='center right',frameon=False)
for b,c,col,style in [(3,1,'#176a94','--'),(1,2.5,'#be5033','-'),(3,2.5,'#7856a4','-'),(10,2.5,'#27875b','-')]:
 d=data['late_b%g_c%g.csv'%(b,c)];ax[1].semilogx(d['a'],d['w_phi'],label='c=%g, b=%g'%(c,b),color=col,ls=style,lw=2)
ax[1].set(xlabel='Scale factor a',ylabel='Scalar equation of state w',xlim=(1e-3,1),ylim=(-1.025,-.52),title='Late solutions assume the remnant was prepared');ax[1].axhline(-1,color='#777777',lw=.7);ax[1].legend(loc='upper left',frameon=False)
for a in ax:a.grid(alpha=.15)
fig.savefig(p/'overview.png',dpi=160);plt.close(fig)
