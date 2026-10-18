"""Frozen B2 experiment; exact invariant sectors, no removed returns.

Shared ancestry: scalar catalog/diagonal from read-only B1. Colored hopping
and conversion are rebuilt directly below, not assumed to tensorize.
"""
import hashlib, importlib.util, json, math, os, pathlib, resource, time
from collections import deque
from fractions import Fraction as F
import numpy as np
from scipy.linalg import eigh, norm
from scipy.sparse import coo_matrix, diags, eye, kron

HERE=pathlib.Path(__file__).resolve().parent
ROOT=HERE.parents[2]
B1=HERE.parent/'grut_b1_neutral_escape_20261018'
spec=importlib.util.spec_from_file_location('pinned_b1_model',B1/'model.py')
b1=importlib.util.module_from_spec(spec);spec.loader.exec_module(b1)


def gamma(n):
    u=2.**-53
    return n*u/(1-n*u)


class Propagator:
    def __init__(self,h,source):
        self.h=h.tocsr();self.n=h.shape[0];self.source=source
        dense=h.toarray();self.e,self.v=eigh(dense)
        self.coeff=self.v[source,:].conj()
        vf=float(norm(self.v,'fro'));hf=float(norm(dense,'fro'))
        self.orth_calc=float(norm(self.v.T@self.v-np.eye(self.n),'fro'))
        self.res_calc=float(norm(h@self.v-self.v*self.e[None,:],'fro'))
        nnz=int(np.max(np.diff(h.indptr)))
        # Positive conservative error factors under the declared binary64
        # contract; no assumption on the eigensolver's accuracy is needed.
        self.orth=self.orth_calc+2*gamma(4*self.n+16)*(vf*vf+math.sqrt(self.n))
        self.rho=self.res_calc+2*gamma(4*nnz+16)*(hf+max(abs(self.e)))*vf
        rowsum=float(np.max(np.array(abs(h).sum(axis=1)).ravel()))
        self.repr=128*2.**-53*(1+rowsum)
        self.vf=vf

    def error(self,t):
        phase=32*2.**-53*(1+max(abs(self.e))*abs(t))
        matvec=2*gamma(4*self.n+16)*self.vf*math.sqrt(1+self.orth)
        return (self.orth+abs(t)*(self.rho*(1+self.orth)+self.repr)
                +phase*(1+self.orth)+matvec)

    def state(self,t):
        return self.v@(np.exp(-1j*self.e*t)*self.coeff)

    def account(self):
        return {'residual_frobenius_computed':self.res_calc,
                'residual_with_arithmetic_account':self.rho,
                'orthogonality_computed':self.orth_calc,
                'orthogonality_with_arithmetic_account':self.orth,
                'intended_rational_H_representation_bound':self.repr,
                'end_time_state_error_bound':self.error(6/b1.JPAIR),
                'scope':'Evaluated enclosure under ARITHMETIC.md contract, not outward-rounded certification.'}


def interval(prob,error):
    s=math.sqrt(max(0,prob))
    return [max(0,s-error)**2,min(1,s+error)**2]


def color_coordinates(model):
    coords=[]
    for x in model.sites:
        coords.append(-1 if x==0 else model.length if x==model.length-1 else x+.25)
    return coords


class Colored:
    def __init__(self,model):
        assert model.resource==2
        self.model=model;self.states=[];self.scalar=[]
        fcoords=color_coordinates(model)
        # Catalog is scalar-major then input order 00,01,10,11.
        for si,(p,m,f) in enumerate(model.states):
            carriers=[(x,'p',a) for a,x in enumerate(p)]
            carriers += [(fcoords[a],'f',a) for a in range(len(model.sites)) if f>>a&1]
            carriers.sort();assert len(carriers)==2
            for code in range(4):
                colors=(code//2,code%2);pc=[-1]*len(p);fc=[-1]*len(model.sites)
                for color,(_,kind,a) in zip(colors,carriers):
                    (pc if kind=='p' else fc)[a]=color
                self.states.append((p,m,f,tuple(pc),tuple(fc)));self.scalar.append(si)
        self.index={s:i for i,s in enumerate(self.states)}
        entries=[]
        for col,(p,m,f,pc,fc) in enumerate(self.states):
            entries.append((col,col,float(model.D[self.scalar[col],self.scalar[col]])))
            for species,slots,coupling in [(0,p,b1.JP),(1,m,b1.JM)]:
                for a,old in enumerate(slots):
                    for new in (old-1,old+1):
                        if new<0 or new>=model.length or new in slots:continue
                        if species==0:
                            pairs=list(zip(p,pc));pairs[a]=(new,pc[a]);pairs.sort()
                            pp=tuple(x for x,c in pairs);pc2=tuple(c for x,c in pairs);mm=m
                        else:
                            pp=p;pc2=pc;mm=tuple(sorted(new if x==old else x for x in m))
                        if max(map(abs,b1.field(pp,mm,model.length)),default=0)>1:continue
                        target=(pp,mm,f,pc2,fc)
                        entries.append((self.index[target],col,-coupling))
            for a,x in enumerate(model.sites):
                if f>>a&1:
                    if x in p or x in m:continue
                    pairs=sorted(list(zip(p,pc))+[(x,fc[a])])
                    pp=tuple(y for y,c in pairs);pc2=tuple(c for y,c in pairs)
                    mm=tuple(sorted(m+(x,)));fc2=list(fc);fc2[a]=-1
                    target=(pp,mm,f^(1<<a),pc2,tuple(fc2))
                else:
                    if x not in p or x not in m:continue
                    pos=p.index(x);fc2=list(fc);fc2[a]=pc[pos]
                    pp=tuple(y for y in p if y!=x);mm=tuple(y for y in m if y!=x)
                    pc2=tuple(c for aa,c in enumerate(pc) if aa!=pos)
                    target=(pp,mm,f^(1<<a),pc2,tuple(fc2))
                entries.append((self.index[target],col,-b1.JPAIR))
        rows,cols,vals=zip(*entries)
        self.h=coo_matrix((vals,(rows,cols)),shape=(len(self.states),len(self.states))).tocsr()
        assert np.max(abs((self.h-self.h.T).data),initial=0)<1e-15


def spectral_one():
    model=b1.PathModel(16);src=int(np.argmax(abs(model.initial())))
    prop=Propagator(model.H,src);weights=abs(prop.coeff)**2
    high_edge=2+5/4-5/96
    high=prop.e>=high_edge-1e-12
    gap=F(2)+F(5,4)-F(5,96)-(F(2)-F(13,23040))
    bound=(F(1,3840)/gap)**2
    assert bound==F(36,27613**2)
    observed=float(np.sum(weights[high]))
    assert observed<float(bound)
    recon=[]
    for theta in [2,4,6]:
        t=theta/b1.JPAIR
        amp=np.sum(weights*np.exp(-1j*prop.e*t))
        direct=prop.state(t)[src]
        assert abs(amp-direct)<1e-13
        recon.append({'theta':theta,'fuel_return_probability':float(abs(amp)**2)})
    return {'dimension':model.dim,'low_cluster_dimension':int(np.sum(~high)),
            'high_cluster_fuel_weight':observed,'high_weight_upper':float(bound),
            'high_weight_exact_fraction':str(bound),'return_reconstruction':recon,
            'all_eigenvalues':prop.e.tolist(),'all_initial_fuel_residues':weights.tolist(),
            'scope':'Complete finite spectrum. Not infinite-volume trapped-state classification.'}


def branch_word():
    model=b1.PathModel(3,sites=(1,),resource=2);colored=Colored(model)
    state=((0,2),(0,2),0,(0,1),(-1,));initial=state
    steps=[('p',0,1),('m',0,1),('absorb',1,1),('p',2,1),('m',2,1),
           ('p',1,0),('m',1,0),('emit',1,1),('p',1,2),('m',1,2)]
    record=[];product=F(1)
    for kind,old,new in steps:
        p,m,f,pc,fc=state;prev=colored.index[state]
        if kind=='p':
            pairs=[(new if x==old else x,c) for x,c in zip(p,pc)];pairs.sort()
            p=tuple(x for x,c in pairs);pc=tuple(c for x,c in pairs);v=-F(1,64)
        elif kind=='m':m=tuple(sorted(new if x==old else x for x in m));v=-F(1,96)
        elif kind=='absorb':
            a=p.index(1);fc=(pc[a],);pc=tuple(c for i,c in enumerate(pc) if i!=a)
            p=tuple(x for x in p if x!=1);m=tuple(x for x in m if x!=1);f=1;v=-F(1,3840)
        else:
            pairs=sorted(list(zip(p,pc))+[(1,fc[0])]);p=tuple(x for x,c in pairs)
            pc=tuple(c for x,c in pairs);m=tuple(sorted(m+(1,)));fc=(-1,);f=0;v=-F(1,3840)
        state=(p,m,f,pc,fc);now=colored.index[state]
        assert abs(colored.h[now,prev]-float(v))<1e-18
        fields=b1.field(p,m,3);assert max(map(abs,fields),default=0)<=1
        product*=v
        record.append({'operation':[kind,old,new],'state':state,'fields':fields,'matrix_element':str(v)})
    target=((0,2),(0,2),0,(1,0),(-1,));assert state==target
    start=colored.index[initial];end=colored.index[target]
    graph=colored.h.copy();graph.setdiag(0);graph.eliminate_zeros()
    dist={start:0};queue=deque([start])
    while queue:
        x=queue.popleft()
        if x==end:break
        for y in graph.indices[graph.indptr[x]:graph.indptr[x+1]]:
            if y not in dist:dist[y]=dist[x]+1;queue.append(int(y))
    # All primitive off-diagonal signs are negative: shortest words cannot
    # cancel at their leading time-series order.
    assert np.all(graph.data<0)
    diff=colored.h-kron(model.H,eye(4),format='csr')
    assert np.max(abs(diff.data),initial=0)>0
    return {'scalar_dimension':model.dim,'colored_dimension':len(colored.states),
            'ten_step_word':record,'word_product':str(product),'shortest_primitive_word':dist[end],
            'ordered_color_factorization_entry_defect':float(np.max(abs(diff.data),initial=0)),
            'scope':'Local allowed exchange. No quantitative growing-geometry arrival probability.'}


def main():
    begin=time.perf_counter()
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())['sources']
    unchanged={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==sha for p,sha in pins.items()}
    assert all(unchanged.values())
    model=b1.PathModel(8,sites=(0,7),resource=2)
    colored=Colored(model)
    defect=colored.h-kron(model.H,eye(4),format='csr')
    entry_defect=float(np.max(abs(defect.data),initial=0));assert entry_defect<1e-15
    src=int(np.argmax(abs(model.initial())));physical=Propagator(model.H,src)
    contact=np.array([.5*sum(x in set(p)&set(m) and x+1 in set(p)&set(m) for x in range(7))
                      for p,m,f in model.states])
    reference=Propagator(model.H-diags(contact),src)
    outgoing=np.array([f==0 and p==m and len(p)==2 and p[0] in (1,2) and p[1] in (5,6)
                       for p,m,f in model.states])
    encounter=np.array([f==0 and p==m and p in ((2,4),(3,5)) for p,m,f in model.states])
    count=model.pairs;fuel_left=np.array([bool(f&1) for p,m,f in model.states])
    fuel_right=np.array([bool(f&2) for p,m,f in model.states])
    rows=[];pR=None;pRef=None;enc=None
    e_init=2*b1.OMEGA
    for theta in range(7):
        t=theta/b1.JPAIR;psi=physical.state(t);prob=abs(psi)**2
        e=physical.error(t);assert e<1e-4
        drift=abs(float(np.vdot(psi,model.H@psi).real)-e_init)
        assert drift<1e-9 and abs(np.sum(prob)-1)<1e-10
        pr=float(np.sum(prob[outgoing]));pc=float(np.sum(prob[encounter]))
        rows.append({'theta':theta,'time':t,'two_pairs_probability':float(np.sum(prob[count==2])),
                     'at_least_one_pair_probability':float(np.sum(prob[count>=1])),
                     'fuel_left_expectation':float(prob@fuel_left),
                     'fuel_right_expectation':float(prob@fuel_right),
                     'mean_contact_energy':float(prob@contact),'encounter_probability':pc,
                     'outgoing_probability':pr,'evaluated_state_error_upper':e,'energy_drift':drift})
        if theta==3:enc=interval(pc,e)
        if theta==6:pR=interval(pr,e)
    t=6/b1.JPAIR;refpsi=reference.state(t)
    pref=float(np.sum(abs(refpsi[outgoing])**2));pRef=interval(pref,reference.error(t))
    difference_lower=max(0,pR[0]-pRef[1],pRef[0]-pR[1])
    difference_upper=max(abs(pR[0]-pRef[1]),abs(pR[1]-pRef[0]))
    def verdict(lo,hi,threshold):
        return 'PASS' if lo>=threshold else 'FAIL' if hi<threshold else 'UNRESOLVED'
    gates={'outgoing':verdict(*pR,.1),'encounter':verdict(*enc,.25),
           'contact_sensitivity':verdict(difference_lower,difference_upper,.01),
           'conditional_information':'PASS by exact ordered-carrier factorization (error=0)'}
    result={'baseline_pins_unchanged':unchanged,'one_composite_spectrum':spectral_one(),
            'two_composites':{'scalar_dimension':model.dim,'colored_dimension':len(colored.states),
                'colored_matrix_entry_defect':entry_defect,'time_points':rows,
                'outgoing_probability_enclosure':pR,'encounter_probability_enclosure':enc,
                'contact_reference_outgoing_probability':pref,
                'contact_reference_outgoing_enclosure':pRef,
                'contact_sensitivity_enclosure':[difference_lower,difference_upper],
                'physical_arithmetic_account':physical.account(),
                'reference_arithmetic_account':reference.account(),'gates':gates},
            'interior_fuel_exchange':branch_word()}
    result['elapsed_seconds']=time.perf_counter()-begin
    result['peak_rss_MiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
    assert result['elapsed_seconds']<60 and result['peak_rss_MiB']<1024
    result['thread_settings']={k:os.environ.get(k) for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']}
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
