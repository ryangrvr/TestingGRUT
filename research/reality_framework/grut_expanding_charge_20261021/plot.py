from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
RAW=HERE/'raw'

def main():
    r=json.loads((HERE/'RESULTS.json').read_text())['runs']
    fig=plt.figure(figsize=(13,8.3),layout='constrained')
    grid=fig.add_gridspec(2,6,height_ratios=[1,1.1])
    ax=fig.add_subplot(grid[0,:3]); ar=fig.add_subplot(grid[0,3:])
    names=['main_s101','main_s202','main_s303','cold_s101','gaussian_s101','gaussian_s202','gaussian_large_s101']
    colors=['#a11b31','#c96a20','#2f7d76','#707780','#2749a6','#548fc9','#482b7a']
    for name,color in zip(names,colors):
        if name not in r:continue
        d=np.atleast_1d(np.genfromtxt(RAW/(name+'.csv'),delimiter=',',names=True))
        den=d['fplus_05']+d['fminus_05'];fraction=np.divide(d['fplus_05'],den,out=np.full_like(den,np.nan),where=den>.1)
        label=name.replace('main_','64 modes, ').replace('cold_s101','cold control').replace('gaussian_large_','Gaussian, larger box ').replace('gaussian_','Gaussian, ')
        style='--' if 'gaussian' in name else '-'
        ax.plot(d['a'],fraction,label=label,color=color,ls=style,lw=1.8)
        top=[json.loads(x) for x in (RAW/(name+'.topology.jsonl')).read_text().splitlines()]
        ta=np.array([x['a'] for x in top]);area=np.array([x['basin_link_area_comoving_over_V'] for x in top])
        ar.plot(ta,np.where(area>0,area,np.nan),color=color,ls=style,lw=1.8)
    ax.axhline(.5,color='#b0b0b0',lw=.8);ax.set(xlabel='Scale factor a',ylabel='Positive fraction among |S| > 0.5 cells',ylim=(-.03,1.03),title='The supplied spatial preparation matters')
    ax.legend(fontsize=8,loc='center left',frameon=True,facecolor='white',edgecolor='none',framealpha=.97)
    ar.set(xlabel='Scale factor a',ylabel='Basin interface area / comoving volume',yscale='log',title='Interface estimate from genuine-vacuum seeds')
    ar.grid(alpha=.15)
    panels=[('main_s101','64 modes, seed 101'),('main_s202','64 modes, seed 202'),('gaussian_s101','Gaussian spectrum, seed 101')]
    for i,(name,title) in enumerate(panels):
        bx=fig.add_subplot(grid[1,2*i:2*i+2])
        if name not in r:bx.set_axis_off();continue
        N=r[name]['N'];L=r[name]['L'];p=RAW/(name+'.slice.f32')
        if not p.exists():bx.text(.5,.5,'Rerun solver to regenerate slice',ha='center');continue
        data=np.fromfile(p,dtype='<f4').reshape(N,N,3)
        im=bx.imshow(data[:,:,0].T,origin='lower',extent=[0,L,0,L],cmap='RdBu_r',vmin=-1,vmax=1,interpolation='nearest')
        bx.set(title=title,xlabel='Comoving x / ell',ylabel='Comoving y / ell')
        if i==2:fig.colorbar(im,ax=bx,label='S = phi/v',fraction=.045,pad=.03)
    fig.suptitle('Candidate 1: independent finite-time 3D preparations\nD=0, prescribed radiation expansion; these runs do not establish cosmic persistence',fontsize=13)
    fig.savefig(HERE/'FINITE_3D.png',dpi=150)
    plt.close(fig)

if __name__=='__main__':main()
