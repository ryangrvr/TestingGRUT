#!/usr/bin/env python3
"""Bounded exact finite-system calculation of the fixed GRUT-M1 proposal.

No random states, fitted Hamiltonian, time integrator, Lyapunov extraction,
or claim of gravity is used. Raw failed acceptance conditions are retained.
"""
import argparse
import hashlib
import json
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import eigh

HERE = Path(__file__).resolve().parent
BASE = dict(J=1., alpha=.5, hx=1.05, hz=.5, K=1., Delta=1., epsilon=.05)


def hamiltonian(L, p):
    n = 1 << (2*L)
    b = np.arange(n, dtype=np.int64)
    z = np.array([1.-2*((b >> j) & 1) for j in range(2*L)])
    diagonal = np.zeros(n)
    H = np.zeros((n, n))
    for i in range(L):
        j = (i+1) % L
        diagonal += p['J']*(1+p['alpha']*z[L+i])*z[i]*z[j]
        diagonal += p['hz']*z[i] - p['K']*z[L+i]*z[L+j]
        H[b ^ (1 << i), b] += p['hx']
        F = p['epsilon'] + (1-p['epsilon'])*(1-z[L+(i-1)%L]*z[L+j])/2
        H[b ^ (1 << (L+i)), b] -= p['Delta']*F
    H[b, b] += diagonal
    return H, z


def kron_hamiltonian(L, p):
    """Independent tensor-product construction; bit0 is the rightmost factor."""
    I = np.eye(2)
    X = np.array([[0., 1.], [1., 0.]])
    Z = np.diag([1., -1.])
    def word(ops):
        out = np.ones((1, 1))
        for k in reversed(range(2*L)):
            out = np.kron(out, ops.get(k, I))
        return out
    H = np.zeros((1 << (2*L),)*2)
    for i in range(L):
        j, prev = (i+1)%L, (i-1)%L
        H += p['J']*word({i:Z, j:Z})
        H += p['alpha']*p['J']*word({i:Z, j:Z, L+i:Z})
        H += p['hx']*word({i:X}) + p['hz']*word({i:Z})
        H -= p['K']*word({L+i:Z, L+j:Z})
        H -= p['Delta']*(1+p['epsilon'])/2*word({L+i:X})
        H += p['Delta']*(1-p['epsilon'])/2*word({L+prev:Z, L+i:X, L+j:Z})
    return H


def product(vectors):
    out = np.array([1.+0j])
    for v in reversed(vectors):
        out = np.kron(out, np.asarray(v, dtype=complex))
    return out


def state(L, mode, r=1):
    plusx = np.array([1., 1.])/np.sqrt(2)
    zbit = np.array([1., 0.]) if r == 1 else np.array([0., 1.])
    if mode == 'retention':
        return product([plusx]*L + [zbit]*L)
    if mode == 'readout':
        matter = [np.array([1., 0.]) for _ in range(L)]
        matter[1] = plusx
        return product(matter + [zbit]*L)
    if mode == 'formation':
        _, vec = eigh(np.array([[BASE['hz'], BASE['hx']], [BASE['hx'], -BASE['hz']]]))
        return product([vec[:, 0]]*L + [plusx]*L)
    raise ValueError(mode)


def evolve(E, V, psi, ts):
    coefficients = V.T.conj() @ psi
    return V @ (coefficients[:, None]*np.exp(-1j*E[:, None]*ts[None, :]))


def mean_diag(psi, diagonal):
    return np.einsum('nt,n,nt->t', psi.conj(), diagonal, psi).real


def mean_y(psi, bit, z):
    flip = np.arange(len(psi)) ^ (1 << bit)
    return np.sum(psi.conj()*(-1j*z[bit, :, None])*psi[flip], axis=0).real


def information(ts, means):
    success = .5 + np.abs(np.array(means[0])-np.array(means[1]))/4
    crossing = np.flatnonzero(success < .9)
    return dict(times=ts.tolist(), means=means, success=success.tolist(),
                first_sample_below_0_9=None if not len(crossing) else float(ts[crossing[0]]))


def prepared_runs(L, p, E, V, z):
    tr = np.linspace(0, 40, 201)
    tw = np.linspace(0, 2, 101)
    means, ymeans, read_mem = [], [], []
    preparation = []
    norm_error = energy_error = 0.
    for r in [1, -1]:
        for mode, ts in [('retention', tr), ('readout', tw)]:
            psi = state(L, mode, r)
            c = V.T.conj() @ psi
            initial_energy = float(np.sum(np.abs(c)**2*E))
            initial_energy_variance = float(np.sum(np.abs(c)**2*(E-initial_energy)**2))
            evolved = evolve(E, V, psi, ts)
            norm_error = max(norm_error, float(np.max(np.abs(np.sum(np.abs(evolved)**2, axis=0)-1))))
            # Evolve back into the energy basis, instead of just reusing phases.
            back = V.T.conj() @ evolved[:, ::20]
            energy_error = max(energy_error, float(np.max(np.abs(np.sum(np.abs(back)**2*E[:, None], axis=0)-initial_energy))))
            preparation.append(dict(mode=mode, bit=r, mean_energy=initial_energy,
                                    energy_variance=initial_energy_variance,
                                    excitation_above_ground=initial_energy-float(E[0])))
            memory = mean_diag(evolved, z[L]).tolist()
            if mode == 'retention':
                means.append(memory)
            else:
                ymeans.append(mean_y(evolved, 1, z).tolist())
                read_mem.append(memory)
    retention = information(tr, means)
    readout = information(tw, ymeans)
    read_memory = information(tw, read_mem)
    k = int(np.argmax(readout['success']))
    readout['best_registered_time'] = float(tw[k])
    readout['best_success'] = float(readout['success'][k])
    readout['memory_success_at_best'] = float(read_memory['success'][k])
    readout['memory_success'] = read_memory['success']
    formation = evolve(E, V, state(L, 'formation'), tr)
    order_op = sum(z[L+i]*z[L+(i+1)%L] for i in range(L))/L
    order = mean_diag(formation, order_op)
    return dict(retention=retention, readout=readout,
                formation=dict(times=tr.tolist(), nearest_neighbor_order=order.tolist(),
                               initial=float(order[0]), maximum=float(np.max(order)),
                               late_window_mean=float(np.mean(order[tr >= 20]))),
                preparation=preparation, norm_error=norm_error, energy_error=energy_error)


def equilibrium(L, E, V, z, beta=1.):
    ts = np.linspace(0, 4, 41)
    tm = np.linspace(0, 40, 201)
    p = np.exp(-beta*(E-E[0])); p /= p.sum()
    n = len(E)
    indices = np.arange(n)
    A = V.T @ V[indices ^ 1]
    dE = E[:, None]-E[None, :]
    y = p**.25
    sqrtp = np.sqrt(p)
    output = []
    for distance in [1, 2]:
        B = V.T @ V[indices ^ (1 << distance)]
        fd = float(np.sum(sqrtp[:, None]*sqrtp[None, :]*A*A) *
                   np.sum(sqrtp[:, None]*sqrtp[None, :]*B*B))
        C, F = [], []
        for t in ts:
            At = A*np.exp(1j*dE*t)
            AB = At @ B
            C.append(float(1-np.einsum('n,nm,mn->', p, AB, AB).real))
            D = (y[:, None]*At*y[None, :]) @ B
            F.append(float(np.einsum('nm,mn->', D, D).real))
        output.append(dict(distance=distance, times=ts.tolist(), squared_commutator=C,
                           regulated_F=F, factorized_reference=fd,
                           initial_relative_factorization_defect=abs(F[0]-fd)/max(abs(fd),1e-15)))
    T = V.T @ (z[L, :, None]*V)
    energy = float(p @ E)
    varianceE = float(p @ (E-energy)**2)
    meanT = float(p @ np.diag(T))
    overlapE = float(p @ (np.diag(T)*(E-energy)))/varianceE
    T[indices, indices] -= meanT+overlapE*(E-energy)
    weight = p[:, None]*np.abs(T)**2
    varianceT = float(weight.sum())
    Cmem = [float(np.sum(weight*np.cos(dE*t))/varianceT) for t in tm]
    plateau = float(np.sum(weight*(np.abs(dE) < 1e-10))/varianceT)
    fdt_error = float(np.max(np.abs(np.log(p[:, None]/p[None, :])+beta*dE)))
    return dict(beta=beta, energy=energy, variance_energy=varianceE,
                otoc=output, memory=dict(times=tm.tolist(), autocorrelation=Cmem,
                removed_mean=meanT, removed_energy_projection=overlapE,
                residual_variance=varianceT, finite_system_degenerate_plateau=plateau),
                detailed_balance_log_residual=fdt_error,
                lyapunov_exponent=None,
                exponent_reason="No large-size factorization or controlled exponential window is established; no MSS ratio inferred from these finite curves.")


def kinematic_checks(L, p, H, z):
    n = len(H)
    ids = np.arange(n)
    Z = np.diag(z[L])
    Y = np.zeros((n, n), dtype=complex)
    Y[ids, ids ^ (1 << 1)] = -1j*z[1]
    X0 = np.eye(n)[ids ^ 1]
    X1 = np.eye(n)[ids ^ 2]
    D = 1j*(H @ X0-X0 @ H)
    comm = D @ X1-X1 @ D
    out = []
    for r in [1, -1]:
        psi = state(L, 'retention', r)
        double = H @ (H @ Z-Z @ H)-(H @ Z-Z @ H) @ H
        second = float(-np.vdot(psi, double @ psi).real/r)
        phi = state(L, 'readout', r)
        first = float((1j*np.vdot(phi, (H @ Y-Y @ H) @ phi)).real)
        # Maximally mixed sigma, definite uniform tau pattern.
        selected = np.arange(1 << L)+(0 if r == 1 else ((1 << L)-1)*(1 << L))
        ct = .5*float(np.trace((comm.T.conj() @ comm)[np.ix_(selected, selected)]).real)/(1 << L)
        out.append(dict(bit=r, memory_second_derivative=second,
                        expected_memory_second=-4*p['Delta']**2*p['epsilon']**2,
                        readout_first_derivative=first,
                        expected_readout_first=4*p['J']*(1+p['alpha']*r)+2*p['hz'],
                        matter_commutator_t2_coefficient=ct,
                        expected_matter_coefficient=8*p['J']**2*(1+p['alpha']*r)**2))
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='results_initial.json')
    args = parser.parse_args()
    target = HERE/args.output
    if target.exists():
        raise SystemExit('Refusing to overwrite an execution record: '+str(target))
    start = time.monotonic()
    configs = [('candidate_L3',3,BASE), ('candidate_L4',4,BASE), ('candidate_L5',5,BASE),
               ('unfacilitated_L4',4,{**BASE,'epsilon':1.}),
               ('frozen_L4',4,{**BASE,'epsilon':0.}),
               ('no_feedback_L4',4,{**BASE,'alpha':0.}),
               ('half_softness_L4',4,{**BASE,'epsilon':.025})]
    result = dict(protocol_sha256=hashlib.sha256((HERE/'protocol.json').read_bytes()).hexdigest(),
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
                  runs={}, checks={})
    for name,L,p in configs:
        t0 = time.monotonic()
        H,z = hamiltonian(L,p)
        E,V = eigh(H, check_finite=True, driver='evr')
        run = dict(L=L, dimension=len(E), parameters=p, ground_energy=float(E[0]),
                   maximum_energy=float(E[-1]), hermiticity_error=float(np.max(np.abs(H-H.T))),
                   eigen_residual=float(np.max(np.abs(H @ V-V*E[None,:]))))
        run.update(prepared_runs(L,p,E,V,z))
        if name in ['candidate_L4','unfacilitated_L4','no_feedback_L4']:
            run['equilibrium'] = equilibrium(L,E,V,z)
        if name == 'candidate_L3':
            run['kronecker_max_error'] = float(np.max(np.abs(H-kron_hamiltonian(L,p))))
            run['kinematic_checks'] = kinematic_checks(L,p,H,z)
        run['wall_seconds'] = time.monotonic()-t0
        result['runs'][name] = run
        # Checkpoint each completed run; a crash preserves completed data.
        target.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(dict(run=name,seconds=round(run['wall_seconds'],3),
              min_memory_t10=min(v for t,v in zip(run['retention']['times'],run['retention']['success']) if t<=10),
              best_readout=run['readout']['best_success'],
              retained_at_readout=run['readout']['memory_success_at_best'])),flush=True)
    runs = result['runs']; candidate = runs['candidate_L4']
    result['checks'] = {
        'all_hamiltonians_hermitian':all(r['hermiticity_error']<1e-13 for r in runs.values()),
        'independent_kronecker_match':runs['candidate_L3']['kronecker_max_error']<1e-13,
        'all_eigen_residuals':all(r['eigen_residual']<1e-10 for r in runs.values()),
        'all_norms_conserved':all(r['norm_error']<1e-10 for r in runs.values()),
        'all_energies_conserved':all(r['energy_error']<1e-9 for r in runs.values()),
        'record_main_t10':min(v for t,v in zip(candidate['retention']['times'],candidate['retention']['success']) if t<=10)>=.9,
        'readout_main':candidate['readout']['best_success']>.6 and candidate['readout']['memory_success_at_best']>=.9,
        'active_matter_interactions':max(candidate['equilibrium']['otoc'][0]['squared_commutator'])>.1,
        'no_feedback_no_readout':max(runs['no_feedback_L4']['readout']['success'])<.500000001,
        'frozen_record':min(runs['frozen_L4']['retention']['success'])>1-1e-10,
        'short_time_derivatives':all(abs(c[k]-c[e])<1e-10 for c in runs['candidate_L3']['kinematic_checks'] for k,e in [('memory_second_derivative','expected_memory_second'),('readout_first_derivative','expected_readout_first'),('matter_commutator_t2_coefficient','expected_matter_coefficient')]),
        'zero_time_matter_commutators':all(abs(o['squared_commutator'][0])<1e-10 for r in runs.values() if 'equilibrium' in r for o in r['equilibrium']['otoc']),
        'finite_size_retention_t10':all(min(v for t,v in zip(runs[name]['retention']['times'],runs[name]['retention']['success']) if t<=10)>=.9 for name in ['candidate_L3','candidate_L5'])
    }
    result['all_registered_acceptance_pass'] = all(result['checks'].values())
    result['total_wall_seconds'] = time.monotonic()-start
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(checks=result['checks'],all_pass=result['all_registered_acceptance_pass'],seconds=result['total_wall_seconds'])),flush=True)


if __name__ == '__main__':
    main()
