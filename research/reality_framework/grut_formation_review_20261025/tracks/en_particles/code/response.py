"""en_particles part 2b (own code): retarded response of psi from +-eps kicked twin runs [E].

G_k(D) = sum_shell dpsi_k(D) eta_k^* / sum_shell |eta_k|^2, dpsi = (psi(+eps)-psi(-eps))/2, kick in pi_psi=pi1+i pi2.
G(t,t') = -1/2 {psi(t), psi*(t')}, so rho(-D) = -rho(D)^*; rho_std(w) = -2 Im Gt(w), Gt = int_0^T G e^{-iwD} taper dD.
Branch weights Z_pm = int_{w>0 / w<0} dw/2pi w rho_std (free field: 1/2 each; sum rule Z_+ + Z_- = 1).
Classical FDT occupation per branch: nbar_pm = n_count_pm / (2 Z_pm) (n_count from spectral_window.py).
Usage: python3 -I response.py <K_plus/window.npz> <K_minus/window.npz> <kick.npy> <out.json> [stat_json]
"""
import sys, json
import numpy as np


def kcont(N, L):
    a = np.fft.fftfreq(N)*N*2*np.pi/L
    return np.sqrt(a[:, None, None]**2+a[None, :, None]**2+a[None, None, :]**2)


def main(fp, fm, fk, out, stat=None):
    P = np.load(fp); Mn = np.load(fm); eps = float(P['eps']); assert abs(float(Mn['eps'])+eps) < 1e-15
    s = P['s']; D = s-s[0]; ds = D[1]; N = int(P['N']); L = float(P['L'])
    eta = np.fft.fftn(np.load(fk), norm='forward')  # this is the +eps kick
    dpsi = (P['psi']-Mn['psi'])/2
    kc = kcont(N, L); dk = 2*np.pi/L; sid = np.floor(kc/dk+.5).astype(int).ravel()
    fl = dpsi.reshape(len(s), -1); ef = eta.ravel()
    T = D[-1]; taper = np.cos(np.pi*D/(2*T))**2  # 1 at D=0 with zero slope, 0 at D=T
    wgrid = np.linspace(-20, 20, 4001); dw = wgrid[1]-wgrid[0]
    E = np.exp(-1j*np.outer(wgrid, D))
    trap = np.full(len(D), ds); trap[[0, -1]] = ds/2
    st = None
    if stat: st = {r['shell']: r for r in json.load(open(stat))['shells']}
    rows = []
    for n in range(1, sid.max()+1):
        idx = np.flatnonzero(sid == n)
        if len(idx) == 0: continue
        G = (fl[:, idx]*np.conj(ef[idx])[None, :]).sum(axis=1)/(np.abs(ef[idx])**2).sum()
        # anomalous check: correlation with eta at -k (should be small after averaging)
        Gt = E@(G*taper*trap)
        rho = -2*Gt.imag
        Zp = float((wgrid*rho)[wgrid > 0].sum()*dw/(2*np.pi)); Zm = float((wgrid*rho)[wgrid < 0].sum()*dw/(2*np.pi))
        ip = np.argmax(np.where(wgrid > 0, rho, -np.inf)); im = np.argmax(np.where(wgrid < 0, -rho, -np.inf))
        row = dict(shell=n, k_cont=n*dk, modes=len(idx), Gprime0=float((G[1]/ds).real), G0=float(abs(G[0])),
                   Z_plus=Zp, Z_minus=Zm, Z_sum=Zp+Zm, rho_peak_pos=float(wgrid[ip]), rho_peak_neg=float(wgrid[im]),
                   rho_neg_weight_positive_part=float(((wgrid*rho)[(wgrid < 0)] > 0).mean()))
        if st and n in st:
            r = st[n]; row['n_count_minus'] = r['n_minus']; row['n_count_plus'] = r['n_plus']
            row['nbar_minus'] = r['n_minus']/(2*Zm) if Zm > 0 else None
            row['nbar_plus'] = r['n_plus']/(2*Zp) if Zp > 0 else None
            row['F_peak_pos'] = r['pos'].get('omega'); row['F_peak_neg'] = r['neg'].get('omega')
        rows.append(row)
    res = dict(source='own analysis of reproduction-type kicked twin runs (Codex engine unmodified)', eps=eps,
               s0=float(s[0]), T=float(T), shells=rows)
    if st:
        res['n_count_minus_total'] = sum(r['n_count_minus'] for r in rows if 'n_count_minus' in r)
        res['nbar_minus_total'] = sum(r['nbar_minus'] for r in rows if r.get('nbar_minus') is not None)
        res['nbar_plus_total'] = sum(r['nbar_plus'] for r in rows if r.get('nbar_plus') is not None)
    json.dump(res, open(out, 'w'), indent=1)
    print('eps', eps, 's0', s[0], 'T', T)
    print('sh  k     Gp0    Z+     Z-    sum   rho+    rho-   F+     F-    ncnt-  nbar-')
    for r in rows[:22]:
        print(f"{r['shell']:2d} {r['k_cont']:5.2f} {r['Gprime0']:6.3f} {r['Z_plus']:6.3f} {r['Z_minus']:6.3f} {r['Z_sum']:5.3f} "
              f"{r['rho_peak_pos']:6.2f} {r['rho_peak_neg']:6.2f} {r.get('F_peak_pos') or 0:6.2f} {r.get('F_peak_neg') or 0:6.2f} "
              f"{r.get('n_count_minus', 0):.4f} {r.get('nbar_minus') or 0:.4f}")
    if st: print('totals n_count- %.4f nbar- %.4f nbar+ %.4f' % (res['n_count_minus_total'], res['nbar_minus_total'], res['nbar_plus_total']))


if __name__ == '__main__':
    main(*sys.argv[1:])
