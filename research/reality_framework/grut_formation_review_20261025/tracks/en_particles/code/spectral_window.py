"""en_particles part 2 analysis (own code): unequal-time power spectra and measured-frequency counts [E].

Input: window.npz from window_run.py (uniform conformal-time samples of psi_k=FFT(chi1+i chi2), pi_psi, phi_k).
Convention: psi_k(s) ~ exp(+i w s) with w>0 carries positive charge (Q=Im(psi* dpsi/ds)), i.e. particle (n+);
w<0 is the antiparticle branch (n-). Time FT X(w)=sum_j h_j psi(s_j) exp(-i w s_j), Hann taper h, normalised
P(w)=|X|^2/(M sum h^2) so that sum_w P = h^2-weighted time average of |psi_k|^2.
Measured-frequency (spectral) count per mode: n_pm = sum_{w>0 / w<0} |w| P(w). For a free field this equals the
standard occupation; for a broad spectrum it is the classical analogue of int dw/2pi |w| F(w,k) and assumes the
canonical sum rule int dw/2pi w rho = 1 (checked separately with the response twin run).
Usage: python3 -I spectral_window.py <window.npz> <out.json> [nQ=8113.477086664805] [pad=4] [i0=0] [i1=M]
"""
import sys, json
import numpy as np

LX, G2, LP = 1e-3, 1e-2, 3e-2


def fd4(N, L, full=True):
    a = 2*np.pi*np.fft.fftfreq(N); dx = L/N
    f = lambda x: 4*np.sin(x/2)**2*(1+np.sin(x/2)**2/3)/dx**2
    b = a if full else 2*np.pi*np.fft.rfftfreq(N)
    return f(a)[:, None, None]+f(a)[None, :, None]+f(b)[None, None, :]


def kcont(N, L, full=True):
    a = np.fft.fftfreq(N)*N*2*np.pi/L
    b = a if full else np.fft.rfftfreq(N)*N*2*np.pi/L
    return np.sqrt(a[:, None, None]**2+a[None, :, None]**2+b[None, None, :]**2)


def peak_stats(w, P, sign):
    sel = (w > 0) if sign > 0 else (w < 0)
    ws, Ps = w[sel], P[sel]
    if Ps.sum() <= 0: return dict(omega=None)
    i = int(np.argmax(Ps)); half = Ps[i]/2
    lo = i
    while lo > 0 and Ps[lo] > half: lo -= 1
    hi = i
    while hi < len(Ps)-1 and Ps[hi] > half: hi += 1
    fw = abs(ws[hi]-ws[lo])
    inside = np.abs(ws-ws[i]) <= max(fw, 1e-9)  # +-1 FWHM around the peak
    mean_abs = float((np.abs(ws)*Ps).sum()/Ps.sum())
    return dict(omega=float(ws[i]), fwhm=float(fw), power_fraction_within_1fwhm=float(Ps[inside].sum()/Ps.sum()),
                mean_abs_omega=mean_abs, rms_abs_omega=float(np.sqrt((ws**2*Ps).sum()/Ps.sum())))


def main(path, out, nQ=8113.477086664805, pad=4, i0=0, i1=None):
    d = np.load(path); sl = slice(i0, i1); s = d['s'][sl]; M = len(s); ds = s[1]-s[0]; N = int(d['N']); L = float(d['L'])
    psi = d['psi'][sl]; ppsi = d['ppsi'][sl]; phi = d['phi'][sl]; pphi = d['pphi'][sl]; zero = d['zero'][sl]
    h = np.hanning(M); h2 = (h**2).sum()
    wgt = h**2/h2  # weights for the matching equal-time averages
    # ---------------- equal-time references along the window (own code, full-FFT psi formulation)
    k2 = fd4(N, L); nz = k2 > 0
    k2r = fd4(N, L, full=False)
    wr = np.full(N//2+1, 2.); wr[[0, -1]] = 1.; wr = np.broadcast_to(wr[None, None, :], k2r.shape)
    eq = dict(massless=[], hartree=[], mH2=[], X2=[], phi2=[], Qnz=[], Q0=[], nplus_massless=[], phi_massless=[],
              phi_hartree=[], Q=[])
    for j in range(M):
        a = psi[j].astype(np.complex128); b = ppsi[j].astype(np.complex128)
        X2 = float((np.abs(a)**2).sum()); f2 = float((wr*np.abs(phi[j])**2).sum())
        mH2 = 2*LX*X2+G2*f2; mf2 = G2*X2+3*LP*f2
        Qk = np.imag(np.conj(a)*b)
        for name, m2 in [('massless', 0.), ('hartree', mH2)]:
            om = np.sqrt(np.where(nz, k2+m2, 1.))
            Nk = (np.abs(b)**2/om+om*np.abs(a)**2)/2
            nm = float((.5*(Nk-Qk))[nz].sum()); eq[name].append(nm/nQ)
            if name == 'massless': eq['nplus_massless'].append(float((.5*(Nk+Qk))[nz].sum())/nQ)
        for name, m2 in [('phi_massless', 0.), ('phi_hartree', mf2)]:
            om = np.sqrt(np.where(k2r > 0, k2r+m2, 1.))
            nf = (wr*(np.abs(pphi[j])**2/om+om*np.abs(phi[j])**2)/2)[k2r > 0].sum()
            eq[name].append(float(nf)/nQ)
        eq['mH2'].append(mH2); eq['X2'].append(X2); eq['phi2'].append(f2)
        eq['Qnz'].append(float(Qk[nz].sum())/nQ); eq['Q0'].append(float(Qk[~nz].sum())/nQ); eq['Q'].append(float(Qk.sum()))
    eqa = {k: np.array(v) for k, v in eq.items()}
    eq_summary = {k: dict(window_weighted_mean=float((wgt*v).sum()), min=float(v.min()), max=float(v.max()),
                          first=float(v[0]), last=float(v[-1])) for k, v in eqa.items()}
    # ---------------- zero mode spectrum
    Mp = pad*M
    w = 2*np.pi*np.fft.fftfreq(M, ds); wp = 2*np.pi*np.fft.fftfreq(Mp, ds)
    z = zero[:, 0]
    Z = np.fft.fft(h*z); P0 = np.abs(Z)**2/(M*h2)
    Zp = np.fft.fft(h*z, Mp); P0p = np.abs(Zp)**2/(M*h2)/pad
    i0 = int(np.argmax(P0p)); omega0 = float(wp[i0])
    zero_mode = dict(omega0_lab=omega0, n_plus=float((np.abs(w)*P0)[w > 0].sum()/nQ),
                     n_minus=float((np.abs(w)*P0)[w < 0].sum()/nQ),
                     charge_spectral=float((w*P0).sum()/nQ), charge_equal_time=float((wgt*eqa['Q0']).sum()),
                     neg=peak_stats(wp, P0p, -1), pos=peak_stats(wp, P0p, +1),
                     radius_mean=float(np.sqrt((wgt*np.abs(z)**2).sum())), radius_min=float(np.abs(z).min()),
                     radius_max=float(np.abs(z).max()))
    # ---------------- nonzero modes: shells of width 2pi/L in continuum |k|
    kc = kcont(N, L); dk = 2*np.pi/L
    sid = np.floor(kc/dk+.5).astype(int)  # shell n: |k| in [(n-.5)dk,(n+.5)dk)
    nsh = sid.max()+1
    flat_psi = psi.reshape(M, -1); sidf = sid.ravel(); k2f = k2.ravel()
    shells = []; tot = dict(n_minus=0., n_plus=0., charge_spec=0., charge_eq=0., n_minus_branchfreq=0.,
                            n_plus_branchfreq=0., n_minus_rms=0.)
    Pall = np.zeros(Mp)
    corot = np.zeros(Mp)
    for n in range(1, nsh):
        idx = np.flatnonzero(sidf == n)
        if len(idx) == 0: continue
        ts = flat_psi[:, idx].astype(np.complex128)*h[:, None]
        X = np.fft.fft(ts, axis=0); P = (np.abs(X)**2).sum(axis=1)/(M*h2)
        Xp = np.fft.fft(ts, Mp, axis=0); Pp = (np.abs(Xp)**2).sum(axis=1)/(M*h2)/pad
        Pall += Pp
        # co-rotating frame: psi exp(-i omega0 s) shifts w -> w - omega0
        npl = float((np.abs(w)*P)[w > 0].sum()); nmi = float((np.abs(w)*P)[w < 0].sum())
        csp = float((w*P).sum())
        # equal-time charge (h^2-weighted) for this shell
        a = flat_psi[:, idx].astype(np.complex128); b = ppsi.reshape(M, -1)[:, idx].astype(np.complex128)
        ceq = float((wgt[:, None]*np.imag(np.conj(a)*b)).sum())
        X2sh = float((wgt[:, None]*np.abs(a)**2).sum()); P2sh = float((wgt[:, None]*np.abs(b)**2).sum())
        pk_p = peak_stats(wp, Pp, +1); pk_m = peak_stats(wp, Pp, -1)
        pow_p = float(P[w > 0].sum()); pow_m = float(P[w < 0].sum())
        nbp = pk_p['omega']*pow_p if pk_p['omega'] else 0.; nbm = -pk_m['omega']*pow_m if pk_m['omega'] else 0.
        rms_m = pk_m.get('rms_abs_omega', 0.)*pow_m if pk_m['omega'] else 0.
        keff = float(np.sqrt(k2f[idx].mean()))
        shells.append(dict(shell=n, k_cont=n*dk, k_fd4_rms=keff, modes=len(idx), n_plus=npl/nQ, n_minus=nmi/nQ,
                           charge_spectral=csp/nQ, charge_equal_time=ceq/nQ, power_pos=pow_p, power_neg=pow_m,
                           omega_eq_time=float(np.sqrt(P2sh/X2sh)), pos=pk_p, neg=pk_m,
                           n_minus_branchfreq=nbm/nQ, n_plus_branchfreq=nbp/nQ))
        tot['n_minus'] += nmi/nQ; tot['n_plus'] += npl/nQ; tot['charge_spec'] += csp/nQ; tot['charge_eq'] += ceq/nQ
        tot['n_minus_branchfreq'] += nbm/nQ; tot['n_plus_branchfreq'] += nbp/nQ; tot['n_minus_rms'] += rms_m/nQ
    # ---------------- phi (real field): total quanta via sum_k sum_w |w| P_k(w)
    flat_phi = phi.reshape(M, -1); wrf = np.asarray(wr).ravel(); k2rf = k2r.ravel()
    sel = np.flatnonzero(k2rf > 0); nphi = 0.; phi_sh = []
    kcr = kcont(N, L, full=False).ravel(); sidr = np.floor(kcr/dk+.5).astype(int)
    for n in range(1, sidr.max()+1):
        idx = np.intersect1d(np.flatnonzero(sidr == n), sel)
        if len(idx) == 0: continue
        ts = flat_phi[:, idx].astype(np.complex128)*h[:, None]
        X = np.fft.fft(ts, axis=0); P = ((np.abs(X)**2)*wrf[idx][None, :]).sum(axis=1)/(M*h2)
        Xp = np.fft.fft(ts, Mp, axis=0); Pp = ((np.abs(Xp)**2)*wrf[idx][None, :]).sum(axis=1)/(M*h2)/pad
        # rfft half space with weight: both signs of w appear in -k partner -> symmetrise
        Ps = .5*(P+P[(-np.arange(M)) % M]); Pps = .5*(Pp+Pp[(-np.arange(Mp)) % Mp])
        nn = float((np.abs(w)*Ps).sum()); nphi += nn
        phi_sh.append(dict(shell=n, k_cont=n*dk, n=nn/nQ, peak=peak_stats(wp, Pps, +1)))
    res = dict(source='own analysis of a reproduction-type regenerated N48 window (Codex engine, unmodified)',
               window=dict(s_start=float(s[0]), s_end=float(s[-1]), samples=M, ds=float(ds),
                           frequency_resolution=float(2*np.pi/(M*ds)), hann_fwhm=float(1.44*2*np.pi/(M*ds))),
               nQ=nQ, equal_time=eq_summary, zero_mode=zero_mode,
               totals=dict(spectral_n_minus_per_net=tot['n_minus'], spectral_n_plus_per_net=tot['n_plus'],
                           spectral_charge_nonzero_per_net=tot['charge_spec'],
                           equal_time_charge_nonzero_per_net=tot['charge_eq'],
                           branchfreq_n_minus_per_net=tot['n_minus_branchfreq'],
                           branchfreq_n_plus_per_net=tot['n_plus_branchfreq'],
                           rms_freq_n_minus_per_net=tot['n_minus_rms'],
                           phi_quanta_spectral_per_net=nphi/nQ),
               shells=shells, phi_shells=phi_sh,
               lab_spectrum_all_nonzero=dict(omega=wp[np.argsort(wp)], P=Pall[np.argsort(wp)]))
    def cv(v):
        if isinstance(v, np.ndarray): return v.tolist()
        if isinstance(v, np.generic): return v.item()
        raise TypeError(type(v))
    json.dump(res, open(out, 'w'), default=cv)
    brief = {k: res[k] for k in ['window', 'equal_time', 'zero_mode', 'totals']}
    print(json.dumps(brief, indent=1, default=cv))
    print('shell k  modes  n-  n+  Q_spec Q_eq  w+peak fwhm  w-peak fwhm  w_eq  k_fd4')
    for r in shells[:24]:
        print(f"{r['shell']:3d} {r['k_cont']:5.2f} {r['modes']:5d} {r['n_minus']:.4f} {r['n_plus']:.4f} {r['charge_spectral']:+.4f} {r['charge_equal_time']:+.4f} "
              f"{r['pos'].get('omega') or 0:6.2f} {r['pos'].get('fwhm') or 0:5.2f} {r['neg'].get('omega') or 0:6.2f} {r['neg'].get('fwhm') or 0:5.2f} {r['omega_eq_time']:5.2f} {r['k_fd4_rms']:5.2f}")


if __name__ == '__main__':
    a = sys.argv[1:]
    kw = dict(x.split('=') for x in a[2:])
    main(a[0], a[1], nQ=float(kw.get('nQ', 8113.477086664805)), pad=int(kw.get('pad', 4)),
         i0=int(kw.get('i0', 0)), i1=int(kw['i1']) if 'i1' in kw else None)
