"""Independent finite checks for the history-exact active receiver method."""
from itertools import combinations_with_replacement
from math import comb, factorial, sqrt
from pathlib import Path
import importlib.util

import numpy as np
from scipy.fft import dst
from scipy.sparse import diags

import active as A

_path = Path(__file__).resolve().parent.parent / "wm2_prepared_disorder_20261016/model.py"
_spec = importlib.util.spec_from_file_location("frozen_s4_model", _path)
M = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(M)


def field_isometry(frame, photons):
    physical_modes, rank = frame.shape
    physical = [list(combinations_with_replacement(range(physical_modes), k)) for k in range(photons+1)]
    lookup = [{o: i for i, o in enumerate(a)} for a in physical]
    source = list(combinations_with_replacement(range(rank), photons))
    out = np.zeros((len(physical[-1]), len(source)), complex)
    for column, occupation in enumerate(source):
        c = np.array([1+0j])
        seen = {}
        for step, mode in enumerate(occupation):
            seen[mode] = seen.get(mode, 0)+1
            y = np.zeros(len(physical[step+1]), complex)
            for j, old in enumerate(physical[step]):
                for r in range(physical_modes):
                    new = tuple(sorted(old+(r,)))
                    y[lookup[step+1][new]] += (c[j]*frame[r, mode]
                                               *sqrt((old.count(r)+1)/seen[mode]))
            c = y
        out[:, column] = c
    return out


def lift(evolution, exact):
    layout = evolution.layout
    frame = dst(evolution.w, type=1, norm="ortho", axis=0)
    out = np.zeros(len(exact["states"]), complex)
    frames = [field_isometry(frame, k) for k in range(layout.M+1)]
    for wi, word in enumerate(layout.words):
        k = len(word)
        amplitudes = frames[k] @ evolution.c[layout.offset[wi]:layout.offset[wi]+len(layout.occ[k])]
        physical_occ = list(combinations_with_replacement(range(evolution.modes), k))
        for j, o in enumerate(physical_occ):
            out[exact["index"][word, o]] = amplitudes[j]
    return out


def algebra_checks():
    rng = np.random.default_rng(20261016)
    rows = []
    for rank, photons in [(2, 1), (3, 2), (4, 3)]:
        raw = rng.normal(size=(rank, rank))+1j*rng.normal(size=(rank, rank))
        frame = np.linalg.qr(raw)[0]
        iso = field_isometry(frame, photons)
        assert np.linalg.norm(iso.conj().T@iso-np.eye(len(iso))) < 2e-12
        layout = A.Layout([photons], [0], rank)
        c = rng.normal(size=layout.size)+1j*rng.normal(size=layout.size)
        c /= np.linalg.norm(c)
        r = np.linalg.qr(rng.normal(size=(2, 2))+1j*rng.normal(size=(2, 2)))[0]
        pair = np.eye(rank, dtype=complex)
        pair[np.ix_([0, rank-1], [0, rank-1])] = r
        rotated = layout.beam(c.copy(), 0, rank-1, r)
        error = 0
        for k in range(photons+1):
            expected = field_isometry(pair, k)@layout.block(c, k)[0]
            error += np.linalg.norm(expected-layout.block(rotated, k)[0])**2
        error = sqrt(error)
        assert error < 2e-12
        rho = layout.density(c)
        assert abs(np.trace(rho).real-sum(k*np.linalg.norm(layout.block(c, k))**2 for k in range(photons+1))) < 2e-12
        rows.append(dict(rank=rank, maximum_photons=photons,
                         unitary_Fock_transform_error=error,
                         one_body_trace=float(np.trace(rho).real)))
    return rows


def commutator_checks():
    rows = []
    for n, delta, g in [(3, [.5, .25, .75], .9), (4, [0, -.25, .25, -.25], 1.1)]:
        m = M.vacuum(n, 4, delta, g, classes=True)
        b = m["H_emission"].toarray()
        a = m["H"].toarray()-b
        comm = a@b-b@a
        first, second = a@comm-comm@a, -b@comm+comm@b
        bounds = A.commutator_bounds(g)
        norms = [float(np.linalg.norm(x, 2)) for x in [b, first, second]]
        assert norms[0] <= g+1e-11
        assert norms[1] <= bounds["AAB"]+1e-11
        assert norms[2] <= bounds["BBA"]+1e-11
        rows.append(dict(N=n, receiver_modes=4, dimension=len(a), g=g,
                         emission_norm=norms[0], AAB_norm=norms[1], BBA_norm=norms[2], bounds=bounds))
    return rows


def complete_wavefunction_checks():
    rows = []
    for n, modes, delta in [(3, 5, [0, -.25, .25]), (4, 4, [0, -.25, .25, -.25])]:
        values = sorted(set(delta[1:]))
        counts = [delta[1:].count(x) for x in values]
        exact = M.vacuum(n, modes, delta, .9, classes=True)
        for time in [.5, 2, 8]:
            e = A.Evolution(counts, values, modes, .9, rank_cap=modes)
            r = e.run(time, strang_budget=.001, compression_budget=.01, max_seconds=60)
            assert r["status"] == "COMPLETE"
            y = M.evolve(exact["H"], M.unit(len(exact["states"]), exact["seed"]), time)
            approximation = lift(e, exact)
            error = float(np.linalg.norm(y-approximation))
            assert error <= r["analytic_state_error_account"]+2e-10
            rows.append(dict(N=n, modes=modes, time=time, full_dimension=len(y),
                             complete_wavefunction_error=error,
                             analytic_state_error_account=r["analytic_state_error_account"],
                             maximum_active_rank=r["maximum_active_rank"],
                             projections=r["projections"]))
    return rows


def quadratic_reference(evolution, time, delta):
    n, modes, g = evolution.N, evolution.modes, evolution.g
    off = np.r_[-g/n, -np.ones(modes-1)]
    h = diags([off, np.r_[0, np.full(modes, delta)], off], [-1, 0, 1], format="csr")
    z = M.evolve(h, M.unit(modes+1, 0), time)
    spectral = dst(z[1:], type=1, norm="ortho")
    mode = evolution.w.conj().T@spectral
    projection = np.zeros_like(evolution.c)
    m = n-1
    for k in range(n):
        for j, o in enumerate(evolution.layout.occ[k]):
            a = sqrt(factorial(m)/factorial(m-k)) * z[0]**(m-k)
            for r in o:
                a *= mode[r]
            for r in set(o):
                a /= sqrt(factorial(o.count(r)))
            evolution.layout.block(projection, k)[0, j] = a
    error = sqrt(max(0.0, 1-float(np.vdot(projection, projection).real))
                 + float(np.linalg.norm(evolution.c-projection)**2))
    return dict(state_error=error, exact_pump_survival=float(abs(z[0])**2),
                exact_macro_probability=float(1-sum(comb(m, k)*(1-abs(z[0])**2)**k*abs(z[0])**(2*(m-k))
                                                  for k in range(max(2, (n+1)//2)-1))))


def sign_symmetry_checks():
    rows = []
    for n, modes, delta in [(3, 4, [.5, .25, .75]), (4, 3, [0, -.25, 1/12, .25])]:
        positive = M.vacuum(n, modes, delta, .9, classes=False)
        negative = M.vacuum(n, modes, [-v for v in delta], .9, classes=False)
        for t in [.5, 2, 8]:
            x = M.evolve(positive["H"], M.unit(len(positive["states"]), positive["seed"]), t)
            y = M.evolve(negative["H"], M.unit(len(negative["states"]), negative["seed"]), t)
            error = float(np.max(abs(abs(x)**2-abs(y)**2)))
            assert error < 1e-11
            rows.append(dict(N=n, modes=modes, t=t, complete_basis_population_difference=error))
    return rows


def physical_reference_check():
    n, modes, duration = 3, 576, 36
    exact = M.vacuum(n, modes, [0, -.25, .25], .9, classes=True)
    y = M.evolve(exact["H"], M.unit(len(exact["states"]), exact["seed"]), duration)
    e = A.Evolution([1, 1], [-.25, .25], modes, .9)
    r = e.run(duration, max_seconds=120)
    assert r["status"] == "COMPLETE"
    frame = dst(e.w, type=1, norm="ortho", axis=0)
    projection = np.zeros_like(e.c)
    for wi, word in enumerate(e.layout.words):
        k = len(word)
        if k == 0:
            projection[e.layout.offset[wi]] = y[exact["index"][word, ()]]
        elif k == 1:
            field = np.asarray([y[exact["index"][word, (a,)]] for a in range(modes)])
            projection[e.layout.offset[wi]:e.layout.offset[wi]+e.layout.rank] = frame.conj().T@field
        else:
            tensor = np.zeros((modes, modes), complex)
            for a in range(modes):
                for b in range(a, modes):
                    value = y[exact["index"][word, (a,b)]]
                    tensor[a,b] = tensor[b,a] = value if a == b else value/sqrt(2)
            reduced = frame.conj().T@tensor@frame.conj()
            for j, (a,b) in enumerate(e.layout.occ[2]):
                projection[e.layout.offset[wi]+j] = reduced[a,b]*(1 if a == b else sqrt(2))
    missing = max(0.0, float(np.vdot(y,y).real-np.vdot(projection,projection).real))
    error = sqrt(float(np.linalg.norm(e.c-projection)**2)+missing)
    assert error <= r["analytic_state_error_account"]+2e-9
    return dict(N=n, physical_receiver_modes=modes, time=duration,
                complete_physical_matrix_dimension=len(y), complete_wavefunction_error=error,
                analytic_state_error_account=r["analytic_state_error_account"],
                maximum_active_rank=r["maximum_active_rank"],
                maximum_coefficient_dimension=r["maximum_coefficient_dimension"],
                complete_reference_norm_error=abs(float(np.vdot(y,y).real)-1),
                comparison_uses_full_two_photon_tensor=True)
