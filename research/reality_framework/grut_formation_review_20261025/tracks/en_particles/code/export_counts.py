"""en_particles part 1 (own code): recompute charge, energy and reference counts from SPECTRAL_EXPORT.npz.

Independent of Codex's analyze.py/export_covariance.py: reads only the exported arrays and the published
trace values (nQ, rho) for the denominators. Counts here are bin-level (50 bins); the per-mode exact
numbers are recomputed from a regenerated snapshot in snapshot_counts.py.

Usage: python3 -I export_counts.py <SPECTRAL_EXPORT.npz> <FIELD_RESULTS.json> <out.json>
"""
import sys, json
import numpy as np

LX, G2, LP = 1e-3, 1e-2, 3e-2


def fd4_symbol(N, L):
    a = 2*np.pi*np.fft.fftfreq(N); b = 2*np.pi*np.fft.rfftfreq(N); dx = L/N
    f = lambda x: 4*np.sin(x/2)**2*(1+np.sin(x/2)**2/3)/dx**2
    return f(a)[:, None, None]+f(a)[None, :, None]+f(b)[None, None, :]


def bin_mode_counts(N, L, edges):
    """Number of lattice modes (full k space, k=0 excluded) and mean FD4 k^2 per export bin."""
    f = np.fft.fftfreq(N)*2*np.pi*N/L; last = np.fft.rfftfreq(N)*2*np.pi*N/L
    kp = np.sqrt(f[:, None, None]**2+f[None, :, None]**2+last[None, None, :]**2)
    w = np.full(N//2+1, 2.); w[[0, -1]] = 1.; w = np.broadcast_to(w[None, None, :], kp.shape).copy(); w[0, 0, 0] = 0
    ids = np.clip(np.digitize(kp.ravel(), edges)-1, 0, len(edges)-2)
    cnt = np.bincount(ids, weights=w.ravel(), minlength=len(edges)-1)
    k2 = fd4_symbol(N, L)
    k2m = np.bincount(ids, weights=(w*k2).ravel(), minlength=len(edges)-1)/np.maximum(cnt, 1)
    kc = np.bincount(ids, weights=(w*kp).ravel(), minlength=len(edges)-1)/np.maximum(cnt, 1)
    return cnt, k2m, kc


def symplectic_eigs(C):
    """Symplectic eigenvalues of a real symmetric 2n x 2n matrix ordered (q1..qn, p1..pn)."""
    n = C.shape[0]//2
    J = np.block([[np.zeros((n, n)), np.eye(n)], [-np.eye(n), np.zeros((n, n))]])
    ev = np.linalg.eigvals(1j*J@C)
    return np.sort(np.abs(ev.real))[::2][::-1]  # each appears as +-nu


def main(npz, fr, out, tag='N48'):
    d = np.load(npz); res = json.load(open(fr))['cases'][tag]
    cov = d[tag+'_covariance']; ck = d[tag+'_k2_covariance']; edges = d[tag+'_k_edges']
    mq = d[tag+'_mean_q']; mp = d[tag+'_mean_p']; U = float(d[tag+'_potential_density']); L = float(d[tag+'_L'])
    N = res['N']; nQ = res['final']['nQ']; rho = res['final']['rho']
    # --- charge and energy reconstruction
    Q0 = mq[0]*mp[1]-mq[1]*mp[0]
    Qb = cov[:, 0, 3]-cov[:, 1, 2]
    Qnz = Qb.sum()
    Ekin = .5*(cov[:, 2, 2]+cov[:, 3, 3]+cov[:, 5, 5]).sum()+.5*np.sum(mp**2)
    Egrad = .5*(ck[:, 0, 0]+ck[:, 1, 1]+ck[:, 4, 4]).sum()
    Etot = Ekin+Egrad+U
    # --- Hartree mass from export
    r2 = mq[0]**2+mq[1]**2+(cov[:, 0, 0]+cov[:, 1, 1]).sum()
    ph2 = mq[2]**2+cov[:, 4, 4].sum()
    mH2 = 2*LX*r2+G2*ph2
    cnt, k2m, kc = bin_mode_counts(N, L, edges)
    Sqq = cov[:, 0, 0]+cov[:, 1, 1]; Spp = cov[:, 2, 2]+cov[:, 3, 3]; Skqq = ck[:, 0, 0]+ck[:, 1, 1]
    ok = Sqq > 0
    rows = {}
    for name, m2 in [('massless', 0.), ('isotropic_Hartree', mH2)]:
        for wname, w2 in [('omega2_from_k2cov', np.where(ok, Skqq/np.where(ok, Sqq, 1), 0)+m2),
                          ('omega2_from_bin_mean_k2', k2m+m2)]:
            om = np.sqrt(np.maximum(w2, 1e-300))
            ntot = np.where(ok, (Spp/om+om*Sqq)/2, 0)
            nm = .5*(ntot-Qb); npl = .5*(ntot+Qb)
            rows[name+'/'+wname] = dict(mass2=m2, nplus=npl.sum(), nminus=nm.sum(), nminus_over_net=nm.sum()/nQ)
    # --- lower envelope over a common omega per bin (Cauchy-Schwarz): ntot_b >= sqrt(Spp Sqq)
    nmin = np.sqrt(Spp*Sqq); nm_min = .5*(nmin-Qb)
    rows['per_bin_optimal_omega_lower_envelope'] = dict(nminus=nm_min.sum(), nminus_over_net=nm_min.sum()/nQ,
        note='min over a common omega per bin of the chosen-omega count; omega_b^2=Spp/Sqq')
    # --- basis-free Gaussian (symplectic) count of the charged 4x4 block per bin
    nu = np.zeros((50, 2)); qcheck = []
    for b in range(50):
        if not ok[b]: continue
        C = cov[b][np.ix_([0, 1, 2, 3], [0, 1, 2, 3])]
        nu[b] = symplectic_eigs(C)
    # If the per-bin block were U(1)-symmetric, nu = (n+, n-) with n+ - n- = Q_b.
    rows['symplectic_charged_block'] = dict(nu_large=nu[:, 0].sum(), nu_small=nu[:, 1].sum(),
        nu_small_over_net=nu[:, 1].sum()/nQ, sum_nu_over_net=nu.sum()/nQ,
        nu_large_minus_small_over_Qnz=(nu[:, 0]-nu[:, 1]).sum()/Qnz,
        note='basis-free lower bound for any linear canonical (Bogoliubov) single-bin basis; n- interpretation needs U(1) symmetry')
    # U(1) asymmetry diagnostic per bin: <chi1^2>-<chi2^2>, <chi1 chi2>
    asym = np.where(ok, np.hypot(cov[:, 0, 0]-cov[:, 1, 1], 2*cov[:, 0, 1])/np.where(ok, Sqq, 1), 0)
    # --- neutral references (count of phi quanta per net charge)
    Sff = cov[:, 4, 4]; Spf = cov[:, 5, 5]; Skff = ck[:, 4, 4]; okf = Sff > 0
    mphi2 = G2*r2+3*LP*ph2
    neut = {}
    for name, m2 in [('massless', 0.), ('Hartree', mphi2)]:
        om = np.sqrt(np.where(okf, Skff/np.where(okf, Sff, 1), 0)+m2)
        neut[name] = float(np.where(okf, (Spf/np.maximum(om, 1e-300)+om*Sff)/2, 0).sum()/nQ)
    neut['symplectic'] = float(np.sqrt(np.maximum(Sff*Spf-cov[:, 4, 5]**2, 0)).sum()/nQ)
    # --- seeded half-quanta (zero point) per net charge
    kmin, kmax = 1.8, 8.0
    zp_cont = .5*(kmax**3-kmin**3)/(6*np.pi**2)  # per complex species per unit volume, n- only
    out_d = dict(
        source='own code on SPECTRAL_EXPORT.npz (bin-level); denominators from FIELD_RESULTS.json final nQ, rho',
        tag=tag, nQ=nQ, rho=rho,
        charge=dict(Q0=Q0, Q_nonzero=Qnz, Q_total=Q0+Qnz, Q_total_over_nQ_minus1=(Q0+Qnz)/nQ-1,
                    Q0_over_nQ=Q0/nQ, Q_nonzero_over_nQ=Qnz/nQ),
        energy=dict(kinetic=Ekin, gradient=Egrad, potential=U, total=Etot, total_over_rho_minus1=Etot/rho-1,
                    zero_mode_kinetic=.5*np.sum(mp**2)),
        hartree=dict(mean_X2=r2, mean_phi2=ph2, mH2=mH2, mphi_Hartree2=mphi2),
        charged_counts=rows,
        U1_asymmetry_energy_weighted=float((asym*Sqq).sum()/Sqq.sum()),
        neutral_quanta_per_net_charge=neut,
        zero_point=dict(continuum_band_half_quanta_density_per_species=zp_cont,
                        half_quanta_nminus_per_net_charge=zp_cont/nQ,
                        all_lattice_modes_half_quanta_per_species=.5*N**3/L**3,
                        all_lattice_half_quanta_per_net_charge=.5*N**3/L**3/nQ),
        per_bin=dict(k_center=kc, modes=cnt, Q=Qb, Sqq=Sqq, Spp=Spp, Skqq=Skqq, nu=nu,
                     omega_eff_equal_time=np.sqrt(np.where(ok, Spp/np.where(ok, Sqq, 1), 0)),
                     omega_massless=np.sqrt(np.where(ok, Skqq/np.where(ok, Sqq, 1), 0))))
    def cv(v):
        if isinstance(v, np.ndarray): return v.tolist()
        if isinstance(v, np.generic): return v.item()
        raise TypeError(type(v))
    json.dump(out_d, open(out, 'w'), indent=1, default=cv)
    print(json.dumps({k: v for k, v in out_d.items() if k != 'per_bin'}, indent=1, default=cv))


if __name__ == '__main__':
    main(*sys.argv[1:4])
