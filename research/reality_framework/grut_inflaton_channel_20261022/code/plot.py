from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parents[1]
d=json.loads((P/'RESULTS.json').read_text())
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(1,3,figsize=(12.6,3.7),layout='constrained')
r=np.loadtxt(P/d['pulse'][0]['file'],delimiter=',',skiprows=1)
ax[0].plot(r[:,0],r[:,1],label='Mode evolution',color='#006d77')
ax[0].plot(r[:,0],r[:,2],'--',label='First Born',color='#d98424')
ax[0].set(xlabel=r'$k/m_I$',ylabel='Produced occupation',title='Derivative wave production')
ax[0].set_xlim(.43,.57);ax[0].legend(frameon=False)
r=np.loadtxt(P/d['receiver_refined']['file'],delimiter=',',skiprows=1)
total=r[:,3]+r[:,4]+r[:,5]
ax[1].plot(r[:,1],r[:,3]/total,label='Inflaton',color='#725193')
ax[1].plot(r[:,1],r[:,4]/total,label='Remnant waves',color='#006d77')
ax[1].set(xlabel='Expansion since declared start (e-folds)',ylabel='Energy fraction',title='Conserved receiver, no visible bath')
ax[1].set_ylim(0,1.04);ax[1].legend(frameon=False)
x=np.geomspace(.5,100,300);a=d['dark_radiation']['coefficient']
ax[2].loglog(x,a/x,color='#006d77')
ax[2].axhline(.2,color='#ad5b35',ls='--',label='Illustrative tolerance 0.2')
ax[2].set(xlabel=r'Supplied $\Gamma_{vis}/\Gamma_\chi$',ylabel=r'Conditional $\Delta N_{eff}$',title='Conditional radiation relation')
ax[2].legend(frameon=False,fontsize=8)
fig.suptitle('One action: a decay rate, an energy receiver, and a conditional radiation relation',fontsize=13)
fig.savefig(P/'overview.png',dpi=160)
