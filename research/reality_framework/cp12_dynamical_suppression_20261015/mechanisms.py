"""Bounded, gauge-invariant diagnostics for Pauli availability and detuning.

The dephased references below are diagonal parts of L^dagger L in the full
labelled physical basis, not in the permutation-orbit basis. They diagnose
interference in the actual coherent state without applying a dephasing channel.
"""
from itertools import combinations
from pathlib import Path
import hashlib,json,sys,time
import numpy as np
from scipy.sparse import coo_matrix
from boundary import cache_root,RF
sys.path.insert(0,str(RF/'cp12_shortcut_dynamics_20261015'))
from operators import flow


def source_statistics(n,q,s):
    es=list(combinations(range(n),2));M=len(es);g,bits,E=s
    cm,reach=flow(n,s);nv=(bits&((1<<n)-1)).bit_count()
    din=dout=large4=large8=cap2=neutral=available=0.
    correction=.1*(M-q)/M
    for e,(i,j) in enumerate(es):
        if not(g>>e&1):continue
        x,y=E[2*e:2*e+2];nb=(bits>>(n+e))&1
        if x==y==0:
            assert nb==1
            neutral+=1;available+=2-((bits>>i)&1)-((bits>>j)&1)
        for side,v in [(0,i),(1,j)]:
            if nb==((bits>>v)&1):continue
            shift=1 if nb else -1
            xx=x+(shift if side==0 else 0);yy=y+(shift if side==1 else 0)
            target=False
            if xx*yy<0:
                a,b=(i,j) if xx>0 else (j,i);target=bool(reach[b]>>a&1)
            wi=.5*(bool(cm>>e&1)+target);wo=1-wi
            din+=wi*wi;dout+=wo*wo
            removed_before=nb and x==y==0
            removed_after=(not nb) and xx==yy==0
            leg=x if side==0 else y
            delta=2*shift+(2*leg*shift+1)/8+correction*(removed_after-removed_before)
            d=abs(delta);weight=wo*wo
            large4+=weight*(d>=4);large8+=weight*(d>=8);cap2+=weight*min(d,8.)**2
    return np.array([dout,din,neutral,available,large4,large8,cap2,
        cm.bit_count()/q,nv/n,float(cm.bit_count()>=int(np.ceil(q/4)) and nv>=int(np.ceil(n/4)))])


NAMES=['dephased_outside_action','dephased_inside_action','neutral_edges','neutral_edge_vacancies',
    'outside_detuning_ge4','outside_detuning_ge8','outside_detuning_square_capped8',
    'circulation_fraction','matter_fraction','charged_event']


def statistics(m,n,q,name):
    digest=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:12]
    path=cache_root()/f'statistics_{name}_{digest}.npz'
    bounds=dict(zip(NAMES,[2*q,2*q,q,2*q,2*q,2*q,128*q,1,1,1]))
    if path.exists():
        data=np.load(path);return {name:data[name] for name in NAMES},bounds
    start=time.monotonic();last=start;values=np.empty((len(NAMES),len(m['states'])))
    for i,s in enumerate(m['states']):
        values[:,i]=source_statistics(n,q,s)
        if time.monotonic()-last>15:
            print(json.dumps(dict(stage='mechanism_statistics',processed=i,dimension=len(m['states']),
                seconds=time.monotonic()-start)),flush=True);last=time.monotonic()
    result=dict(zip(NAMES,values));path.parent.mkdir(parents=True,exist_ok=True);np.savez(path,**result)
    return result,bounds


def action_operators(m,n):
    """V and K_inside, including every outgoing action, excluding rewiring.

    In an adaptive basis, rewiring can leave the retained set at unchanged
    electric energy. Thus B is not generally the outgoing hopping matrix.
    """
    nv=np.array([(s[1]&((1<<n)-1)).bit_count() for s in m['states']])
    a=m['H'].tocoo();keep=nv[a.row]!=nv[a.col]
    K=coo_matrix((a.data[keep],(a.row[keep],a.col[keep])),shape=a.shape).tocsr()
    if 'outer_states' in m:
        no=np.array([(s[1]&((1<<n)-1)).bit_count() for s in m['outer_states']])
        b=m['B'].tocoo();use=no[b.row]!=nv[b.col]
        BK=coo_matrix((b.data[use],(b.row[use],b.col[use])),shape=b.shape).tocsr()
    else:
        # Complete electric-square-budget closure contains all rewiring targets.
        BK=m['B']
    return {'outside':(m['V'],m['BV']),'inside':(K-m['V'],BK-m['BV'])}
