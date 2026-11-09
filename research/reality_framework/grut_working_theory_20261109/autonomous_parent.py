"""W1 closed parent with an explicit parity-dephasing ancilla.

Finite theoretical consistency check; no device or physical measurements.
Independent bit-index assembly, distinct from the four-qubit density solver.
"""
import argparse,json,math,time
from pathlib import Path
import numpy as np
from scipy.sparse import bmat,csr_matrix,eye
from scipy.sparse.linalg import expm_multiply

HERE=Path(__file__).resolve().parent
REGISTERS=('S','M1','M2','A','E','D1','D2');D=2**len(REGISTERS);L=9
BITS=np.array([[(i>>(6-j))&1 for j in range(7)] for i in range(D)])
S,M1,M2,A,E,D1,D2=BITS.T;PARITY=M1^M2;READOUT=4*S+2*PARITY+A

def diagonal(values):return csr_matrix(np.diag(values))
def flip(target,mask):
    row=np.arange(D)^((np.asarray(mask,dtype=int))<<(6-target))
    return csr_matrix((np.ones(D),(row,np.arange(D))),shape=(D,D),dtype=complex)
def pulse(theta,selector):
    I=eye(D,dtype=complex,format='csr');X=flip(2,np.ones(D,dtype=int));P=diagonal(selector)
    return I+P@((math.cos(theta/2)-1)*I-1j*math.sin(theta/2)*X)
def links(theta1,theta2,beta):
    return [flip(2,S),diagonal(np.exp(-1j*np.pi*(1-2*M1)/4)),pulse(theta1,D1),
      diagonal((-1.)**(E*PARITY)),diagonal(np.exp(-.5j*beta*(1-2*PARITY))),
      pulse(theta2,D2),flip(3,PARITY),eye(D,dtype=complex,format='csr')]
def preparation(eta,program):
    p=np.zeros(D,complex)
    for s in (0,1):
        for memory in (0,1):
            for environment in (0,1):
                bits=(s,memory,memory,0,environment,*program);idx=sum(b<<(6-j) for j,b in enumerate(bits))
                p[idx]=.5*math.sqrt((1+eta)/2 if environment==0 else (1-eta)/2)
    return p
def readout(p):return np.bincount(READOUT,weights=abs(p)**2,minlength=8)
def formula(theta1,theta2,beta,eta,program):
    a=theta1*program[0];b=theta2*program[1]
    return .5*(1+math.cos(a)*math.cos(b)-eta*math.sin(a)*math.sin(b)*math.cos(beta))

def run():
    cfg=json.loads((HERE/'parent_protocol.json').read_text());start=time.monotonic();records=[]
    if (cfg['dimension'],cfg['clock_levels'],cfg['omega'])!=(D*L,L,1):raise ValueError('Unregistered parent geometry')
    if cfg['programs']!=[[0,0],[1,0],[0,1],[1,1]]:raise ValueError('Registered selector coverage changed')
    if not cfg['cases'] or 1 not in cfg['times_over_pi'] or not 0<cfg['clock_probability_floor']<1:raise ValueError('Incomplete registered coverage')
    if any(not 0<=case['eta']<=1 for case in cfg['cases']):raise ValueError('Invalid coherence preparation')
    maximum={'hermiticity':0.,'link_unitarity':0.,'total_norm':0.,'clock_distribution':0.,'conditional_tv':0.,'formula_error':0.,'final_clock_transfer':0.,'energy_mean_drift':0.}
    for case in cfg['cases']:
        theta1,theta2,beta,eta=[case[k]*np.pi if k.endswith('_over_pi') else case[k] for k in ('theta1_over_pi','theta2_over_pi','beta_over_pi','eta')]
        U=links(theta1,theta2,beta);blocks=[[None]*L for _ in range(L)]
        for k,gate in enumerate(U):
            weight=math.sqrt((k+1)*(8-k))/2;blocks[k+1][k]=weight*gate;blocks[k][k+1]=weight*gate.conj().T
            error=gate.conj().T@gate-eye(D,format='csr');maximum['link_unitarity']=max(maximum['link_unitarity'],float(np.max(abs(error.data))) if error.nnz else 0.)
        H=bmat(blocks,format='csr');difference=H-H.conj().T
        maximum['hermiticity']=max(maximum['hermiticity'],float(np.max(abs(difference.data))) if difference.nnz else 0.)
        for program in ((0,0),(1,0),(0,1),(1,1)):
            state=preparation(eta,program);ordered=[state]
            for gate in U:ordered.append(gate@ordered[-1])
            expected=[readout(p) for p in ordered]
            target=formula(theta1,theta2,beta,eta,program)
            selected=expected[-1];agreement=float(sum(selected[i] for i in range(8) if ((i>>2)&1)==(i&1)))
            maximum['formula_error']=max(maximum['formula_error'],abs(agreement-target))
            initial=np.zeros(D*L,complex);initial[:D]=state;energy=float(np.vdot(initial,H@initial).real)
            for scale in cfg['times_over_pi']:
                t=scale*np.pi;psi=expm_multiply(-1j*t*H,initial,traceA=0.);by_clock=psi.reshape(L,D);weights=np.sum(abs(by_clock)**2,axis=1)
                b=math.sin(t/2)**2;theory=np.array([math.comb(8,k)*b**k*(1-b)**(8-k) for k in range(L)])
                maximum['total_norm']=max(maximum['total_norm'],abs(float(weights.sum())-1))
                maximum['clock_distribution']=max(maximum['clock_distribution'],float(np.max(abs(weights-theory))))
                maximum['energy_mean_drift']=max(maximum['energy_mean_drift'],abs(float(np.vdot(psi,H@psi).real)-energy))
                if scale==1:maximum['final_clock_transfer']=max(maximum['final_clock_transfer'],abs(float(weights[-1])-1))
                tv={}
                for k,weight in enumerate(weights):
                    if weight>=cfg['clock_probability_floor']:
                        observed=readout(by_clock[k])/weight;tv[str(k)]=float(np.sum(abs(observed-expected[k]))/2)
                        maximum['conditional_tv']=max(maximum['conditional_tv'],tv[str(k)])
                records.append({'case':case['name'],'program':list(program),'t_over_pi':scale,'clock_probabilities':weights.tolist(),'scored_tick_tv':tv,'unscored_ticks':[k for k,w in enumerate(weights) if w<cfg['clock_probability_floor']],
                  'ordered_final_agreement':agreement,'analytic_final_agreement':target})
    return {'scope':'Theoretical autonomous parent with nine clock levels and seven qubits. No hardware or empirical data; does not reclassify the historical512-dimensional parent.',
      'dimension':D*L,'registers':['C']+list(REGISTERS),'maxima':maximum,'tolerance':cfg['absolute_tolerance'],
      'passed':all(v<=cfg['absolute_tolerance'] for v in maximum.values()),'records':records,'seconds':time.monotonic()-start}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    result=run()
    with args.output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({'passed':result['passed'],'maxima':result['maxima'],'seconds':result['seconds']}))
    raise SystemExit(0 if result['passed'] else 1)
