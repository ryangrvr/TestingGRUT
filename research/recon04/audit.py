"""Standard mixed-state locality and one central-spin interaction control.

No GRUT candidate, fit, selector or optimization occurs in this script.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import numpy as np
import scipy
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parent
SEED, TOL = 20261010, 1e-10
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1, -1]).astype(complex)
PLUS = np.ones((2, 2), complex)/2
OMEGA = 0.3


def enc(a):
    a = np.asarray(a)
    return {"real": a.real.tolist(), "imag": a.imag.tolist()}


def kron_all(arrays):
    out = np.array([[1]], complex)
    for a in arrays:
        out = np.kron(out, a)
    return out


def embed(a, i, n):
    return np.kron(np.eye(2**i), np.kron(a, np.eye(2**(n-i-1))))


def pair_embed(a, i, n):
    return np.kron(np.eye(2**i), np.kron(a, np.eye(2**(n-i-2))))


def partial(rho, sites, n):
    sites = list(sites)
    rest = [i for i in range(n) if i not in sites]
    perm = sites + rest + [i+n for i in sites] + [i+n for i in rest]
    a = rho.reshape([2]*(2*n)).transpose(perm)
    a = a.reshape(2**len(sites), 2**len(rest), 2**len(sites), 2**len(rest))
    return np.einsum('abcb->ac', a)


def distance(a, b):
    return float(np.abs(np.linalg.eigvalsh(a-b)).sum()/2)


def helstrom_observable(a, b):
    eig, vec = np.linalg.eigh(a-b)
    sign = np.where(abs(eig)<1e-12, 0, np.sign(eig))
    return (vec*sign)@vec.conj().T


def packing(intervals, n, threshold):
    dp = [0]*(n+1)
    for end in range(n):
        dp[end+1] = dp[end]
        for start in range(end+1):
            if intervals[start,end] >= threshold-TOL:
                dp[end+1] = max(dp[end+1], dp[start]+1)
    return dp[n]


def initial_pair(kind, n, rng):
    dim = 2**n
    if kind == 'product_0p8':
        return [kron_all([(I+s*0.8*Z)/2]*n) for s in (1,-1)]
    if kind == 'shared_entangled_noise':
        a,b = np.zeros(dim), np.zeros(dim)
        a[0],b[-1] = 1,1
        ghz = (a+b)/np.sqrt(2)
        return [0.7*np.outer(v,v)+0.3*np.outer(ghz,ghz) for v in (a,b)]
    if kind == 'identical_maximally_mixed':
        return [np.eye(dim)/dim, np.eye(dim)/dim]
    if kind == 'wishart':
        out=[]
        for _ in range(2):
            a=rng.normal(size=(dim,dim))+1j*rng.normal(size=(dim,dim))
            rho=a@a.conj().T
            out.append(rho/np.trace(rho))
        return out
    raise ValueError(kind)


def locality_controls(rng):
    fixtures=[]
    counts={'sites':0,'intervals':0,'transported_measurements':0,'fixed_blocks':0}
    maxerr=0.0
    for n in (4,6):
        for d in (0,1,2,3):
            for kind in ('product_0p8','shared_entangled_noise',
                         'identical_maximally_mixed','wishart'):
                states=initial_pair(kind,n,rng)
                layers=[]
                u=np.eye(2**n,dtype=complex)
                for depth in range(d):
                    layer=[]
                    for i in range(depth%2,n-1,2):
                        a=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
                        h=(a+a.conj().T)/2
                        h*=np.pi*rng.uniform(0.15,1)/np.linalg.norm(h,2)
                        g=expm(-1j*h)
                        u=pair_embed(g,i,n)@u
                        layer.append({'edge':i,'u':g,'h':h})
                    layers.append(layer)
                final=[u@rho@u.conj().T for rho in states]
                sites=[]
                for i in range(n):
                    initial=[partial(rho,[i],n) for rho in states]
                    di=distance(*initial)
                    a0=helstrom_observable(*initial)
                    support={i}
                    cone=np.eye(2**n,dtype=complex)
                    for layer in layers:
                        old=support.copy()
                        for g in layer:
                            edge={g['edge'],g['edge']+1}
                            if old&edge:
                                support|=edge
                                cone=pair_embed(g['u'],g['edge'],n)@cone
                    block=range(max(0,i-d),min(n,i+d+1))
                    assert set(support)<=set(block)
                    ai=embed(a0,i,n)
                    transported=u@ai@u.conj().T
                    cone_error=float(np.max(abs(transported-cone@ai@cone.conj().T)))
                    expected=[float(np.real(np.trace(rho@ai))) for rho in states]
                    actual=[float(np.real(np.trace(rho@transported))) for rho in final]
                    err=max(cone_error,*(abs(a-b) for a,b in zip(actual,expected)),
                            abs((actual[0]-actual[1])/2-di))
                    db=distance(*[partial(rho,block,n) for rho in final])
                    assert db+TOL>=di and err<TOL
                    maxerr=max(maxerr,err)
                    sites.append({'site':i,'initial_D':di,'block':list(block),
                                  'final_block_D':db,'transported_expectations':actual,
                                  'probability_or_cone_error':err})
                    counts['sites']+=1
                    counts['transported_measurements']+=1
                intervals={}
                for a in range(n):
                    for b in range(a,n):
                        intervals[a,b]=distance(*[partial(rho,range(a,b+1),n) for rho in final])
                        counts['intervals']+=1
                blockrows=[]
                length=2*d+1
                for a in range(0,n-length+1,length):
                    center=a+d
                    value=intervals[a,a+length-1]
                    assert value+TOL>=sites[center]['initial_D']
                    blockrows.append({'start':a,'end':a+length-1,'center_initial_D':sites[center]['initial_D'],
                                      'final_D':value})
                    counts['fixed_blocks']+=1
                deltas=[]
                for delta in (0.1,0.15):
                    threshold=1-2*delta
                    guaranteed=sum(row['center_initial_D']>=threshold-TOL for row in blockrows)
                    cap=packing(intervals,n,threshold)
                    assert cap>=guaranteed
                    deltas.append({'delta':delta,'certified_block_count':guaranteed,
                                   'measured_packing':cap})
                fixtures.append({'N':n,'depth':d,'initial_kind':kind,
                                 'initial_state_hashes':[hashlib.sha256(np.ascontiguousarray(r,dtype='<c16').tobytes()).hexdigest() for r in states],
                                 'initial_state_recipe':'initial_pair(kind,N,rng); fixed global seed and disclosed loop order',
                                 'layers':[[{'edge':g['edge'],'u':enc(g['u']),'h':enc(g['h'])}
                                            for g in layer] for layer in layers],
                                 'sites':sites,'fixed_blocks':blockrows,'packing':deltas})
    return {'fixture_count':len(fixtures),'counts':counts,'maximum_error':maxerr,'fixtures':fixtures}


def rho_memory(r):
    return (I-r*Z)/2


def local_formula(h,g,r,t):
    scale=np.hypot(h,g)
    q=1-2*g*g/(scale*scale)*np.sin(scale*t)**2
    return float(q),float(r*np.sqrt(max(0,1-q*q))),float(h*r*(1-q))


def local_unitaries(h,g,t):
    return [expm(-1j*t*(h*Z+s*g*X)) for s in (1,-1)]


def model_controls():
    times=(0,0.2,np.pi/(2*np.sqrt(2)),1.7,np.pi/np.sqrt(2))
    locals_=[]
    maxerr=0.0
    for h,g in ((1,1),(0.5,0.7),(1,0.3),(1,0)):
        for r in (0,0.4,0.8,1):
            for t in times:
                rho=rho_memory(r)
                up,um=local_unitaries(h,g,t)
                states=[v@rho@v.conj().T for v in (up,um)]
                q,d,energy=local_formula(h,g,r,t)
                actualq=np.trace(up@rho@um.conj().T)
                actuald=distance(*states)
                energies=[float(np.real(np.trace(h*Z@(a-rho)))) for a in states]
                restriction=actuald**2-2*r*energies[0]/h+(energies[0]/h)**2
                error=max(abs(actualq-q),abs(actuald-d),*(abs(e-energy) for e in energies),abs(restriction))
                assert error<TOL
                maxerr=max(maxerr,float(error))
                locals_.append({'h':h,'g':g,'r':r,'time':t,'q':q,'D':actuald,
                                'memory_energy_change':energies[0],
                                'energy_record_residual':restriction,'maximum_error':float(error)})

    profiles=[((1,1,1,1),(0.8,0.8)),((1,1,1,1),(0,0)),
              ((1,0.5,0.3,0.7),(0.4,1)),((0.5,1,0.7,1),(1,0.4)),
              ((1,1,0,0),(0.8,0.8))]
    globals_=[]
    general_phase=[]
    channels=0
    for (hf,hr,gf,gr),(rf,rr) in profiles:
        initial_e=np.kron(rho_memory(rf),rho_memory(rr))
        total_h=OMEGA/2*embed(Z,0,3)+hf*embed(Z,1,3)+hr*embed(Z,2,3)
        total_h+=gf*embed(Z,0,3)@embed(X,1,3)+gr*embed(Z,0,3)@embed(X,2,3)
        assert np.linalg.norm(total_h,2)<=8+TOL
        for t in times:
            u=expm(-1j*t*total_h)
            initial=np.kron(PLUS,initial_e)
            evolved=u@initial@u.conj().T
            conservation_error=float(abs(np.trace(total_h@(evolved-initial))))
            qf,df,ef=local_formula(hf,gf,rf,t)
            qr,dr,er=local_formula(hr,gr,rr,t)
            q=np.exp(-1j*OMEGA*t)*qf*qr
            errors=[conservation_error]
            for b in range(2):
                for c in range(2):
                    unit=np.zeros((2,2),complex);unit[b,c]=1
                    actual=partial(u@np.kron(unit,initial_e)@u.conj().T,[0],3)
                    expected=unit*(1 if b==c else (q if b==0 else q.conjugate()))
                    errors.append(float(np.max(abs(actual-expected))))
                    channels+=1
            cond=[]
            for b in range(2):
                unit=np.zeros((2,2));unit[b,b]=1
                cond.append(partial(u@np.kron(unit,initial_e)@u.conj().T,[1,2],3))
            actualdf=distance(*[partial(a,[0],2) for a in cond])
            actualdr=distance(*[partial(a,[1],2) for a in cond])
            errors.extend([abs(actualdf-df),abs(actualdr-dr)])
            assert max(errors)<TOL
            maxerr=max(maxerr,max(errors))
            globals_.append({'h_F':hf,'h_R':hr,'g_F':gf,'g_R':gr,'r_F':rf,'r_R':rr,
                             'time':t,'channel_coherence':enc(q),'D_F':actualdf,'D_R':actualdr,
                             'energy_F':ef,'energy_R':er,'total_energy_conservation_error':conservation_error,
                             'maximum_error':max(errors)})
            for phi in (0,np.pi/3,np.pi/2,np.pi):
                phase_u=expm(-1j*(hf*np.kron(Z,I)+hr*np.kron(I,Z)))@np.diag([1,1,1,np.exp(-1j*phi)])
                out=[phase_u@a@phase_u.conj().T for a in cond]
                finaldf=distance(*[partial(a,[0],2) for a in out])
                finaldr=distance(*[partial(a,[1],2) for a in out])
                kf=np.sqrt(max(0,1-(1-rr*rr*qr*qr)*np.sin(phi/2)**2))
                kr=np.sqrt(max(0,1-(1-rf*rf*qf*qf)*np.sin(phi/2)**2))
                energy=[float(np.real(np.trace(op@((out[0]+out[1])/2-initial_e))))
                        for op in (hf*np.kron(Z,I),hr*np.kron(I,Z))]
                pe=max(abs(finaldf-df*kf),abs(finaldr-dr*kr),abs(energy[0]-ef),abs(energy[1]-er))
                assert pe<TOL
                maxerr=max(maxerr,float(pe))
                general_phase.append({'h_F':hf,'h_R':hr,'g_F':gf,'g_R':gr,'r_F':rf,'r_R':rr,
                                      'time':t,'phi':phi,'D_F':finaldf,'D_R':finaldr,
                                      'energy_F':energy[0],'energy_R':energy[1],
                                      'maximum_error':float(pe)})

    phase=[]
    tstar=np.pi/(2*np.sqrt(2))
    p11=np.diag([0,0,0,1])
    bare_e=np.kron(Z,I)+np.kron(I,Z)
    drift=expm(-1j*bare_e)  # Same tau=1 in every comparison.
    for r in (0,0.4,0.8,1):
        rho=rho_memory(r)
        ups=local_unitaries(1,1,tstar)
        initial_e=np.kron(rho,rho)
        conditional=[np.kron(v,v)@initial_e@np.kron(v,v).conj().T for v in ups]
        quench_work=-sum((-1)**b*np.real(np.trace((np.kron(X,I)+np.kron(I,X))@conditional[b]))/2 for b in (0,1))
        assert abs(quench_work-2*r)<TOL
        qwrite=np.trace(np.kron(ups[0],ups[0])@initial_e@np.kron(ups[1],ups[1]).conj().T)
        for phi in (0,np.pi/3,np.pi/2,np.pi):
            control=phi*p11
            assert np.linalg.norm(control,2)<=np.pi+TOL
            assert np.max(abs(control@np.kron(Z,I)-np.kron(Z,I)@control))<TOL
            assert np.max(abs(control@np.kron(I,Z)-np.kron(I,Z)@control))<TOL
            assert np.linalg.norm(control+bare_e,2)+OMEGA/2<=8+TOL
            v=drift@expm(-1j*control)
            out=[v@a@v.conj().T for a in conditional]
            df=distance(*[partial(a,[0],2) for a in out])
            dr=distance(*[partial(a,[1],2) for a in out])
            predicted=r*abs(np.cos(phi/2))
            expected_e=r
            energies=[float(np.real(np.trace(embed(Z,i,2)@((out[0]+out[1])/2-initial_e)))) for i in (0,1)]
            errors=[abs(df-predicted),abs(dr-predicted),*(abs(e-expected_e) for e in energies)]
            # Cross-branch operator determines the whole channel, not just S=|+>.
            cross=np.kron(ups[0],ups[0])@initial_e@np.kron(ups[1],ups[1]).conj().T
            qfinal=np.trace(v@cross@v.conj().T)
            errors.append(abs(qfinal-qwrite))
            for b in range(2):
                for c in range(2):
                    cb=np.kron(ups[b],ups[b]); cc=np.kron(ups[c],ups[c])
                    actual=np.trace(v@cb@initial_e@cc.conj().T@v.conj().T)
                    target=1 if b==c else (qwrite if b==0 else qwrite.conjugate())
                    errors.append(abs(actual-target));channels+=1
            conservation=[float(abs(np.trace(embed(Z,i,2)@(out[b]-conditional[b]))))
                          for i in (0,1) for b in (0,1)]
            errors.extend(conservation)
            work_on=float(np.real(np.trace(control@((conditional[0]+conditional[1])/2))))
            work_off=-float(np.real(np.trace(control@((out[0]+out[1])/2))))
            assert abs(work_on+work_off)<TOL
            assert max(errors)<TOL
            maxerr=max(maxerr,float(max(errors)))
            residual=df**2-2*r*energies[0]+energies[0]**2
            expected_residual=-r*r*np.sin(phi/2)**2
            assert abs(residual-expected_residual)<TOL
            phase.append({'r':r,'phi':phi,'D_F':df,'D_R':dr,'D_FR':distance(*out),
                          'system_visibility':float(abs(qfinal)),
                          'energy_F':energies[0],'energy_R':energies[1],
                          'energy_record_residual':residual,
                          'unchanged_initial_resources':True,'common_quench_work':float(quench_work),
                          'ideal_phase_switch_work_on':work_on,'ideal_phase_switch_work_off':work_off,
                          'ideal_net_phase_switch_work':work_on+work_off,
                          'phase_control_norm':float(np.linalg.norm(control,2)),
                          'energy_conservation_errors':conservation,'maximum_error':float(max(errors))})

    rates=[]
    for alpha in (np.pi/4,3*np.pi/4):
        t=alpha/np.sqrt(2)
        q=np.cos(alpha)**4
        gamma=2*np.sqrt(2)*np.tan(alpha)
        choi=np.array([[1,0,0,q],[0,0,0,0],[0,0,0,0],[q,0,0,1]])/2
        assert min(np.linalg.eigvalsh(choi))>=-TOL
        rates.append({'time':t,'visibility':float(q),'gamma':float(gamma),
                      'normalized_choi_eigenvalues':np.linalg.eigvalsh(choi).tolist()})
    assert rates[0]['gamma']>0 and rates[1]['gamma']<0
    commutant=[]
    hint=np.kron(Z,np.kron(X,I)+np.kron(I,X))
    for name,a in (('I',I),('X',X),('Y',Y),('Z',Z)):
        a=embed(a,0,3)
        commutant.append({'Pauli':name,'commutator_norm':float(np.linalg.norm(a@hint-hint@a,2))})
    assert commutant[0]['commutator_norm']<TOL and commutant[3]['commutator_norm']<TOL
    assert commutant[1]['commutator_norm']>0 and commutant[2]['commutator_norm']>0
    bell=np.array([1,0,0,1],complex)/np.sqrt(2)
    correlated=[]
    ups=local_unitaries(1,1,tstar)
    ub=[np.kron(v,v) for v in ups]
    for name,rho in (('product_I_over_4',np.eye(4)/4),('correlated_Bell',np.outer(bell,bell.conj()))):
        assert max(float(np.max(abs(partial(rho,[i],2)-I/2))) for i in (0,1))<TOL
        q=np.trace(ub[0]@rho@ub[1].conj().T)
        expected=0 if name=='product_I_over_4' else 1
        assert abs(q-expected)<TOL
        correlated.append({'preparation':name,'q':enc(q),'local_r_F':0,'local_r_R':0,
                           'same_local_marginals':True,'inside_product_baseline':name=='product_I_over_4'})
    return {'single_memory_fixtures':locals_,'two_memory_fixtures':globals_,
            'phase_fixtures':phase,'general_phase_formula_fixtures':general_phase,
            'channel_matrix_units':channels,
            'rate_controls':rates,'correlated_preparation_control':correlated,
            'pointer_algebra_controls':commutant,'maximum_error':maxerr}


def main():
    locality=locality_controls(np.random.default_rng(SEED))
    model=model_controls()
    result={'status':'STANDARD_MIXED_RECORD_AND_INTERACTION_CONTROLS_PASS_NO_GRUT_LAW',
            'date':'2026-10-10','seed':SEED,'numpy':np.__version__,'scipy':scipy.__version__,
            'floating_tolerance':TOL,'interval_or_experimental_certificate':False,
            'budget':{'cards_consumed':1,'cards_remaining':2,'Card2':'UNOPENED'},
            'locality':locality,'model':model}
    (ROOT/'results.json').write_text(json.dumps(result,separators=(',',':'))+'\n')
    lines=['# Generated interaction controls','',
           'Synthetic diagnostics generated by audit.py from results.json.','',
           '| Same r=0.8 preparation and writing H | Phase phi/pi | Visibility | D_F | D_R | Energy F | Energy R | Energy-record residual |',
           '|---|---:|---:|---:|---:|---:|---:|---:|']
    for row in model['phase_fixtures']:
        if row['r']==0.8:
            lines.append('| Common one-slot ceiling | '+f"{row['phi']/np.pi:.8g} | {row['system_visibility']:.8g} | {row['D_F']:.8g} | {row['D_R']:.8g} | {row['energy_F']:.8g} | {row['energy_R']:.8g} | {row['energy_record_residual']:.8g} |")
    lines+=['','Energy units E0; hbar=1. Near-zero visibility/residual roundoff is not a physical detection threshold.','',
            f"- Mixed-record circuit fixtures: {locality['fixture_count']}"]
    lines.extend(f'- Mixed-record {k}: {v}' for k,v in locality['counts'].items())
    lines.extend([f"- Single-memory exact-formula controls: {len(model['single_memory_fixtures'])}",
                  f"- Two-memory complete-H controls: {len(model['two_memory_fixtures'])}",
                  f"- Matched phase controls: {len(model['phase_fixtures'])}",
                  f"- General phase-formula controls: {len(model['general_phase_formula_fixtures'])}",
                  f"- Interaction/channel matrix-unit checks: {model['channel_matrix_units']}",
                  f"- Maximum mixed-record discrepancy: {locality['maximum_error']:.16g}",
                  f"- Maximum interaction/control discrepancy: {model['maximum_error']:.16g}",
                  '', 'All declared fixtures pass. No interval certificate, experiment or independent review is claimed.',''])
    (ROOT/'TABLES.md').write_text('\n'.join(lines))
    print(result['status'])
    print('\n'.join(lines[-15:]))


if __name__=='__main__':
    main()
