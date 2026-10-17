"""Bounded exact controls. Ordinary double precision, no formation trajectory."""
import hashlib
import importlib.util
import json
import math
from fractions import Fraction
from itertools import product
from pathlib import Path
import resource
import time
import numpy as np
from scipy.linalg import expm
from scipy.sparse import diags
from model import (Model, JP, JM, KAPPA, chi, edges, evolve, fields,
                   local_state, product_prepared)

HERE = Path(__file__).resolve().parent


def maximum(a):
    return float(np.max(np.abs(a.data), initial=0))


def expect(v, h):
    return float(np.vdot(v, h@v).real)


def rabi(d, a, t=math.pi/2):
    z = math.sqrt(a*a+d*d/4)
    return a*a/(z*z)*math.sin(t*z)**2


def tail(k, x):
    # Positive series; sufficient for these small diagnostic arguments.
    term = x**k/math.factorial(k)
    total = term
    for n in range(k+1, max(k+500, 1000)):
        term *= x/n
        total += term
        if term < max(1e-300, total*1e-17):
            break
    return total


def operation_bound(a, length, t):
    c, v = 2/(1-math.exp(-1)), 8*(math.e-1)
    return min(2., tail(a+1, 1.1*t)
               +10*math.sqrt(c)/v*math.exp(-length/2)*math.expm1(v*t/2))


def tensor_control():
    # Independent ambient tensor: +0,-0,+1,-1,E,f. No solved Gauss fields.
    annih = np.array([[0, 1], [0, 0]], complex)
    ids = [np.eye(2)]*4+[np.eye(3), np.eye(2)]
    def embed(operators):
        out = np.array([[1]], complex)
        for k, identity in enumerate(ids):
            out = np.kron(out, operators.get(k, identity))
        return out
    c = [embed({i: annih}) for i in range(4)]
    nn = [v.conj().T@v for v in c]
    es = np.diag([-1, 0, 1]); us = np.array([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
    ee, u, f = embed({4: es}), embed({4: us}), embed({5: annih})
    nf = f.conj().T@f
    gg = [ee-nn[0]+nn[1], -ee-nn[2]+nn[3]]
    conv = f@u.conj().T@c[0].conj().T@c[3].conj().T
    conv += conv.conj().T
    hp = c[2].conj().T@u@c[0]; hm = c[3].conj().T@u.conj().T@c[1]
    hopping = -JP*(hp+hp.conj().T)-JM*(hm+hm.conj().T)
    density = KAPPA*(nn[0]@nn[3]+nn[1]@nn[2])
    base = sum(nn)+ee@ee+3.3*nf
    loadrot = base-3*nf+.97*(f+f.conj().T)
    full = [base, base+1.04*conv, base+hopping+density, loadrot]
    comm = max(float(np.max(abs(g@h-h@g))) for g in gg for h in full)
    assert comm < 1e-12
    assert np.array_equal(es@us-us@es, -us)
    physical = np.flatnonzero(np.diag(gg[0]@gg[0]+gg[1]@gg[1]) == 0)
    mod = Model(2, fixed_word=(0, 1), modes=1, fuel_delta=[.3], eta=[.97], lam=[1.04])
    ix = []
    for w, q, rec, p, m, fuel in mod.states:
        val = fields(2, w, p, m)[0]
        ix.append(np.ravel_multi_index((int(0 in p), int(0 in m), int(1 in p),
                                       int(1 in m), val+1, fuel), (2,2,2,2,3,2)))
    assert set(ix) == set(physical)
    err = max(float(np.max(abs(h[np.ix_(ix, ix)]-hh.toarray())))
              for h, hh in zip(full, [mod.H0, mod.HC, mod.Hop, mod.HLrot]))
    assert err < 1e-12
    assert abs(conv[0, 0]) == 0  # no diagonal production energy.
    return {'ambient_dimension':96, 'physical_dimension':len(ix),
            'constraint_commutator_max':comm, 'matrix_entry_max':err,
            'partial_shift_identity_exact':True}


def pulse_control():
    bound = Fraction(324,349)*Fraction(2280,2401)
    assert bound > Fraction(7,8)
    worst, error, energy_error = 1., 0., 0.
    for d, a, b in product([-.5, 0, .5], [.9, 1., 1.1], [.9, 1., 1.1]):
        v = local_state(d, a, b)
        pl, pc = rabi(d, a), rabi(d, b)
        p = abs(v[2])**2
        worst = min(worst, p)
        error = max(error, abs(p-pl*pc))
        diag = (3+d)*abs(v[1])**2+3*p
        hc = np.array([[0,0,0],[0,0,b],[0,b,0]], complex)
        coupling = np.vdot(v,hc@v).real
        energy_error = max(energy_error, abs(coupling-d*p),
                           abs(diag-((3+d)*pl-d*p)))
    assert worst > 49/64 and error < 1e-12 and energy_error < 1e-12
    v = local_state(0, 1, 1)
    hc = np.array([[0,0,0],[0,3,1],[0,1,3]], complex)
    returned = expm(-1j*(math.pi/2)*hc)@v
    assert abs(returned[2])**2 < 1e-25
    assert abs(returned[1])**2 > 1-1e-12
    return {'parameter_combinations':27,'pair_probability_corner_min':worst,
            'uniform_analytic_pair_lower':49/64,'rabi_max_error':error,
            'energy_ledger_max_error':energy_error,
            'single_stage_rational_lower':str(bound),
            'reverse_pulse_matter_probability':float(abs(returned[2])**2)}


def fixed_control():
    mod = Model(4, matching=((0,1),(2,3)), fixed_word=(0,1,2,3), modes=1,
                fuel_delta=[.3,-.4], eta=[.97,1.06], lam=[1.04,.92])
    psi = mod.initial()
    loaded = mod.load(psi, math.pi/2)
    final = evolve(mod.HC, loaded, math.pi/2)
    target = product_prepared(mod, psi)
    error = float(np.linalg.norm(final-target))
    assert error < 1e-10 and abs(np.linalg.norm(final)-1) < 1e-12
    probabilities = [rabi(d,a)*rabi(d,b) for d,a,b in zip(mod.fdelta,mod.eta,mod.lam)]
    energy_load = expect(loaded,mod.H0)
    energy_off = expect(final,mod.H0)
    prediction = sum((3+d)*rabi(d,a)-d*p for d,a,p in
                     zip(mod.fdelta,mod.eta,probabilities))
    assert abs(energy_off-prediction) < 1e-10
    assert abs(expect(final,mod.C)-sum(d*p for d,p in zip(mod.fdelta,probabilities))) < 1e-10
    n = mod.particles
    p_macro = float(np.vdot(final,(n>=2)*final).real)
    assert abs(p_macro-(1-np.prod(1-np.array(probabilities)))) < 1e-10
    rows = []
    for species, rate, receiver in [(0,JP,1),(1,JM,0)]:
        src = mod.source(species,0,1); obs = mod.number_at(species,receiver)
        variance = expect(final,src@src)-expect(final,src)**2
        assert abs(variance-probabilities[0]/4) < 1e-10
        t0 = 1e-5
        response = chi(mod.Hop,final,src,obs,t0)
        lower = probabilities[0]*rate*t0-64*t0*t0
        assert -response >= lower-1e-13 and lower > 0
        rows.append({'species': '+' if species==0 else '-', 'response':response,
                     'variance':variance,'taylor_lower_magnitude':lower})
    # A single-particle hopping expectation must vanish even though
    # coherent pair/vacuum and cross-pair sectors have NOT been removed.
    assert abs(expect(final,mod.T)) < 1e-12
    # At least one vacuum-pair off-diagonal is nonzero for these detunings.
    assert np.count_nonzero(abs(final)>1e-12) > 4
    work_density = expect(final,mod.V)
    measured_energy = expect(final,diags(mod.field))
    assert abs(measured_energy-sum(probabilities)) < 1e-10
    assert work_density < 1
    return {'dimension':mod.dim,'full_state_error':error,'macro_P_at_least_2':p_macro,
            'pair_probabilities':probabilities,'load_energy':energy_load,
            'post_conversion_energy':energy_off,'field_energy':measured_energy,
            'operating_switch_work':work_density,'near_responses':rows}


def geometry_control():
    n=4
    mod = Model(n, matching=((0,1),), modes=2,
                delta=.23+np.linspace(-.25,.25,n),fuel_delta=[.3],eta=[.97],lam=[1.04])
    psi = mod.initial()
    prefix=(0,1,2)
    # Short, nontrivial boundary control; not the primary preparation time.
    h=.02
    actual = evolve(mod.HC,mod.load(psi,h),h)
    ref = evolve(mod.reference(mod.HC,prefix),mod.load(psi,h,prefix),h)
    ep = float(np.linalg.norm(actual-ref))
    ep_upper = min(2.,tail(2,1.1*2*h))
    assert ep <= ep_upper+1e-10
    # Exact constraints, field cap and new empty-head requirement.
    constraint_max=0
    neutral_head=0
    for index,(w,q,r,p,m,f) in enumerate(mod.states):
        gauss=np.zeros(n,int)
        for (i,j),e in zip(edges(n),mod.ef[index]):
            gauss[i]+=e; gauss[j]-=e
            assert e==0 or (i,j) in {tuple(sorted(z)) for z in zip(w,w[1:])}
        for i in p: gauss[i]-=1
        for i in m: gauss[i]+=1
        constraint_max=max(constraint_max,int(np.max(abs(gauss))))
        assert q.bit_count()-1==sum(r) and max(map(abs,mod.ef[index]),default=0)<=1
        if len(w)>1 and w[-1] in p and w[-1] in m:
            # Without the empty-head gate a zero-flux pop could leave
            # this neutral pair off graph. All actual word-changing
            # entries out of this state must retain its occupied head.
            for row in mod.G[:,index].nonzero()[0]:
                ww=mod.states[row][0]
                assert ww==w or len(ww)>=len(w)
            neutral_head+=1
    assert constraint_max==0 and neutral_head>0

    # Independently paid local preparation, with actual coherent receiver,
    # then short operating comparison. This is not a primary N=4 success.
    initial = product_prepared(mod,psi)
    t=.02
    aa=evolve(mod.Hop,initial,t)
    rr=evolve(mod.reference(mod.Hop,prefix),initial,t)
    error=float(np.linalg.norm(aa-rr))
    upper=operation_bound(1,0,t)
    assert error<=upper+1e-10
    prefix_lost=float(np.vdot(aa,(~mod.retained(prefix))*aa).real)
    inconsistent=np.array([q!=sum(1<<x for x in w) for w,q,r,p,m,f in mod.states])
    storage_error=float(np.vdot(aa,inconsistent*aa).real)
    assert 0<storage_error<=upper**2+1e-10

    # Frozen D1 vacuum block, including all storage-inconsistent states.
    source=HERE.parent/'wm2_many_carrier_20261017'/'model.py'
    spec=importlib.util.spec_from_file_location('d1_frozen_model',source)
    old=importlib.util.module_from_spec(spec); spec.loader.exec_module(old)
    baseline=old.KernelModel(n,0,0)
    indices=[]
    for w,q,n0,a,b in baseline.states:
        rec=(n0,q.bit_count()-1-n0)
        indices.append(mod.index[(w,q,rec,(),(),0)])
    difference=maximum(mod.H0[indices,:][:,indices]-baseline.H)
    assert difference<1e-12
    return {'events':n,'receiver_modes':2,'physical_receiver_modes_would_be':64*n*n,
            'dimension':mod.dim,'operation_nonzeros':mod.Hop.nnz,
            'short_pulses_time_each':h,'preparation_state_error':ep,
            'preparation_bound':ep_upper,'operation_time':t,
            'operation_state_error':error,'operation_bound':upper,
            'prefix_loss_probability':prefix_lost,'storage_error_probability':storage_error,
            'constraint_max':constraint_max,'neutral_head_blocked_controls':neutral_head,
            'vacuum_matrix_entry_error':difference,
            'all_physical_number_sectors_retained':True}


def formulas():
    rows=[]
    for n in [256,512,1024]:
        k=n//32
        failed=math.exp(-289*n/65536)
        ep=math.exp(-n/4)
        d=operation_bound(n//4,3*n/16-1,n/65536)
        near=operation_bound(n//4,3*n/16-1,1e-5)
        coarse=math.exp(-n/5)+3*math.exp(-n/12)
        assert d<=coarse and 2*(ep+near)<2e-14
        assert (49/64)*JM*1e-5-64e-10-2*(ep+near)>7.3e-8
        rows.append({'N':n,'conditional_matter_lower':1-failed-2*ep,
                     'operating_error_at_window_end':d,'coarse_operating_error':coarse,
                     'near_response_correction':2*(ep+near)})
    cubic=[]
    for k in [1,2,4,8,16]:
        fs=list(range(1,k+1))+list(range(k-1,0,-1))
        cost=sum(x*x for x in fs)
        assert 3*cost==k*(2*k*k+1)
        cubic.append({'pairs':k,'neutral_segregated_field_energy':cost,'local_dipole_cost':k})
    return {'evaluations_not_simulations':rows,'neutrality_alone_counterexample':cubic,
            'inherited_S9_threshold':2**28,'receiver_modes_at_threshold':2**62,
            'formation_time_at_threshold':2**58,'original_trial_lower_target':'1/4',
            'independent_validation':False}


def run():
    start=time.perf_counter()
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
    root=HERE.parents[2]
    for path,digest in pins['sources'].items():
        assert hashlib.sha256((root/path).read_bytes()).hexdigest()==digest
    out={'arithmetic':'ordinary numpy/scipy double precision, not outward rounded',
         'tensor':tensor_control(),'pulses':pulse_control(),'fixed_path':fixed_control(),
         'complete_small_geometry':geometry_control(),'formulas':formulas()}
    for path,digest in pins['sources'].items():
        assert hashlib.sha256((root/path).read_bytes()).hexdigest()==digest
    out['source_pins_unchanged']=len(pins['sources'])
    out['resources']={'wall_seconds':time.perf_counter()-start,
                      'peak_rss_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024}
    assert out['resources']['wall_seconds']<60
    assert out['resources']['peak_rss_MiB']<1024
    return out


if __name__=='__main__':
    result=run()
    (HERE/'CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
