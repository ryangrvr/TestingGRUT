"""State-selected Galerkin extension, never a physical formation filter.

Every selected state is a complete signed permutation orbit. All other
outgoing states remain in the residual matrix, without an electric cutoff.
The frozen Hamiltonian and all its escape diagonals are unchanged.
"""
from pathlib import Path
from math import lcm,sqrt
import json,os,pickle,sys,tempfile,time
import numpy as np
from scipy.sparse import bmat,coo_matrix,csr_matrix,vstack,hstack

P=Path(__file__).resolve().parent
RF=P.parent
sys.path.insert(0,str(RF/'cp12_shortcut_dynamics_20261015'))
from operators import quotient_model,Triplets,flow,outside_weight,hopping_edge
sys.path.insert(0,str(RF/'cp12_ground_connectivity'))
from exact_destination import transitions,gauss


def cache_root():
    return Path(os.environ.get('GRUT_SUPPRESSION_CACHE',str(Path(tempfile.gettempdir())/'grut_suppression_cache')))


def load_model(n,q,budget,name=None):
    if name:
        with (cache_root()/f'{name}.pkl').open('rb') as f:return pickle.load(f)
    return quotient_model(n,q,budget)


def recover_outer(m,n,q,budget):
    """Recover row identities in the legacy B construction's exact order."""
    cache=cache_root();name=cache/f'outer_n{n}_q{q}_b{budget}.pkl'
    if name.exists():
        with name.open('rb') as f:return pickle.load(f)
    started=time.monotonic();sub=m['sub'];outer={};vv=Triplets();sc=lcm(8,10*sub.M)
    last=started;cols=np.unique(m['B'].indices)
    for count,col in enumerate(cols):
        s=m['states'][col]
        for target,c in transitions(n,s,sc):
            if sum(x*x for x in target[2])<=budget:continue
            key,sign,size=sub.canonical(target)
            if not sign:continue
            row=outer.setdefault(key,len(outer))
            vv.add(row,int(col),c/sc*sign*sqrt(m['sizes'][col]/size))
        if time.monotonic()-last>15:
            print(json.dumps(dict(stage='recover_boundary',processed=count,columns=len(cols),
                outside=len(outer),seconds=time.monotonic()-started)),flush=True);last=time.monotonic()
    assert len(outer)==m['B'].shape[0]
    new=vv.csr(m['B'].shape);err=float(np.linalg.norm((new-m['B']).data));assert err<1e-10
    result=dict(states=list(outer),row_identity_matrix_error=err,seconds=time.monotonic()-started)
    sub.cache.clear();cache.mkdir(parents=True,exist_ok=True)
    with name.open('wb') as f:pickle.dump(result,f)
    return result


def extend(m,n,q,selected,name):
    """Extend by chosen outgoing row IDs; retain every new outgoing row."""
    cache=cache_root();target=cache/f'{name}.pkl'
    if target.exists():
        with target.open('rb') as f:return pickle.load(f)
    start=time.monotonic();sub=m['sub'];old=m['states'];outer=m['outer_states']
    selected=np.sort(np.asarray(selected,dtype=np.int64));take=set(int(x) for x in selected)
    extra=[outer[i] for i in selected];states=old+extra;nb=len(old);ne=len(extra)
    index={s:i for i,s in enumerate(states)}
    remaining=[i for i in range(len(outer)) if i not in take]
    outside=[outer[i] for i in remaining];outside_index={s:i for i,s in enumerate(outside)}
    sizes=np.r_[m['sizes'],np.asarray([sub.canonical(s)[2] for s in extra])]
    hh=Triplets();vv=Triplets();bb=Triplets();bv=Triplets();sc=lcm(8,10*sub.M)
    cross=Triplets();vcross=Triplets();last=start
    for col,s in enumerate(extra):
        gauss(n,s);info=flow(n,s)
        for t,c in transitions(n,s,sc):
            e=hopping_edge(n,s,t);weight=outside_weight(n,s,t,e,info) if e is not None else 0.
            key,sign,size=sub.canonical(t)
            if not sign:continue
            value=c/sc*sign*sqrt(sizes[nb+col]/size)
            row=index.get(key)
            if row is not None:
                if row<nb:cross.add(row,col,value);vcross.add(row,col,weight*value)
                else:hh.add(row-nb,col,value);vv.add(row-nb,col,weight*value)
            else:
                row=outside_index.get(key)
                if row is None:row=len(outside);outside_index[key]=row;outside.append(key)
                bb.add(row,col,value);bv.add(row,col,weight*value)
        if time.monotonic()-last>15:
            print(json.dumps(dict(stage='extend_boundary',name=name,processed=col,total=ne,
                outside=len(outside),seconds=time.monotonic()-start)),flush=True);last=time.monotonic()
    C=m['B'][selected,:].T.tocsr();CV=m['BV'][selected,:].T.tocsr()
    ce=float(np.linalg.norm((cross.csr((nb,ne))-C).data))
    ve=float(np.linalg.norm((vcross.csr((nb,ne))-CV).data))
    assert max(ce,ve)<1e-9
    H=bmat([[m['H'],C],[C.T,hh.csr((ne,ne))]],format='csr')
    V=bmat([[m['V'],CV],[CV.T,vv.csr((ne,ne))]],format='csr')
    B0=vstack([m['B'][remaining,:],csr_matrix((len(outside)-len(remaining),nb))],format='csr')
    BV0=vstack([m['BV'][remaining,:],csr_matrix((len(outside)-len(remaining),nb))],format='csr')
    B=hstack([B0,bb.csr((len(outside),ne))],format='csr')
    BV=hstack([BV0,bv.csr((len(outside),ne))],format='csr')
    herror=float(np.linalg.norm((H-H.T).data));verror=float(np.linalg.norm((V-V.T).data))
    assert max(herror,verror)<1e-8
    psi=np.r_[m['psi'],np.zeros(ne,dtype=complex)]
    sub.cache.clear()
    info=dict(base_dimension=nb,added_dimension=ne,dimension=H.shape[0],outside=B.shape[0],
        H_nnz=H.nnz,base_cross_H_error=ce,base_cross_V_error=ve,H_hermiticity=herror,
        V_hermiticity=verror,seconds=time.monotonic()-start,
        minimum_added_electric_square=min(sum(x*x for x in s[2]) for s in extra),
        maximum_added_electric_square=max(sum(x*x for x in s[2]) for s in extra))
    result=dict(sub=sub,H=H,V=V,B=B,BV=BV,states=states,outer_states=outside,sizes=sizes,
        psi=psi,extension=info,model_checks=m.get('checks',m.get('model_checks')))
    cache.mkdir(parents=True,exist_ok=True)
    with target.open('wb') as f:pickle.dump(result,f)
    print(json.dumps(dict(stage='extended_matrix',name=name,**info)),flush=True)
    return result


def main():
    import argparse
    a=argparse.ArgumentParser();a.add_argument('--n',type=int,default=6);a.add_argument('--q',type=int,default=12)
    a.add_argument('--budget',type=int,default=12);a.add_argument('--count',type=int,default=225000)
    a.add_argument('--priority');a.add_argument('--name',required=True)
    a.add_argument('--parent-model');a.add_argument('--first',type=int);a.add_argument('--all',action='store_true')
    args=a.parse_args()
    if args.parent_model:
        with (cache_root()/f'{args.parent_model}.pkl').open('rb') as f:m=pickle.load(f)
        outer=dict(row_identity_matrix_error=0.)
    else:
        m=quotient_model(args.n,args.q,args.budget)
        outer=recover_outer(m,args.n,args.q,args.budget);m['outer_states']=outer['states']
    if args.all:selected=np.arange(len(m['outer_states']));coverage=1.
    elif args.first is not None:selected=np.arange(args.first);coverage=None
    else:
        if not args.priority:a.error('Supply --priority, --first or --all')
        priority=np.load(args.priority)['priority'];assert len(priority)==len(m['outer_states'])
        selected=np.argsort(-priority,kind='stable')[:args.count]
        coverage=float(priority[selected].sum()/priority.sum())
    out=extend(m,args.n,args.q,selected,args.name)
    selection=('All outgoing rows of the parent basis.' if args.all else
        'First outgoing rows in the deterministic stored order.' if args.first is not None else
        'Largest outgoing fluence proxy from the base neutral trajectory and normalized charged trajectory.')
    info=dict(name=args.name,parent_model=args.parent_model,selected_count=len(selected),priority_coverage=coverage,
        recovered_boundary=outer['row_identity_matrix_error'],extension=out['extension'],
        selection=selection+' Numerical basis selection only; original neutral amplitudes and physical H unchanged.')
    (P/f'{args.name}_BUILD.json').write_text(json.dumps(info,indent=2)+'\n')


if __name__=='__main__':main()
