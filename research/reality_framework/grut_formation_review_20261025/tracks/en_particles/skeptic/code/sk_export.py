"""en_particles SKEPTIC (own code): export-only checks (charge/energy rebuild, bin-level references) and the tau=20
per-mode occupation in the seeded band. Usage: python3 -I sk_export.py <SPECTRAL_EXPORT.npz> <FIELD_RESULTS.json> <state_tau20.bin> <out.json>"""
import sys, json, struct
import numpy as np
LX, G2, LP = 1e-3, 1e-2, 3e-2

def main(npz, fr, st20, out):
    d = np.load(npz); fin = json.load(open(fr))['cases']['N48']['final']; nQ, rho = fin['nQ'], fin['rho']
    C = d['N48_covariance']; K = d['N48_k2_covariance']; mq = d['N48_mean_q']; mp = d['N48_mean_p']; U = float(d['N48_potential_density'])
    Qb = C[:, 0, 3]-C[:, 1, 2]; Q = mq[0]*mp[1]-mq[1]*mp[0]+Qb.sum()
    E = .5*(np.trace(C[:, 2:4, 2:4], axis1=1, axis2=2)+C[:, 5, 5]).sum()+.5*(mp**2).sum()+.5*(K[:, 0, 0]+K[:, 1, 1]+K[:, 4, 4]).sum()+U
    Sqq = C[:, 0, 0]+C[:, 1, 1]; Spp = C[:, 2, 2]+C[:, 3, 3]; Skq = K[:, 0, 0]+K[:, 1, 1]
    X2 = mq[0]**2+mq[1]**2+Sqq.sum(); f2 = mq[2]**2+C[:, 4, 4].sum(); mH2 = 2*LX*X2+G2*f2
    ok = Sqq > 0; k2b = np.where(ok, Skq/np.where(ok, Sqq, 1), 0)
    res = dict(charge_rebuild_rel=float(Q/nQ-1), energy_rebuild_rel=float(E/rho-1), mH2=float(mH2))
    for name, m2 in [('massless', 0.), ('hartree', mH2)]:
        w = np.sqrt(k2b+m2); Nb = np.where(ok, (Spp/np.where(ok, w, 1)+w*Sqq)/2, 0)
        res['binlevel_nminus_'+name] = float((.5*(Nb-Qb)).sum()/nQ)
    # tau=20 seeded-band occupation per mode (amplitude form, full FFT)
    with open(st20, 'rb') as f: N, L, s, bits = struct.unpack('<iddi', f.read(24))
    a = np.fromfile(st20, dtype='<f8', offset=24).reshape(6, N, N, N)
    psi = np.fft.fftn(a[0]+1j*a[1], norm='forward'); pi = np.fft.fftn(a[3]+1j*a[4], norm='forward')
    n = np.fft.fftfreq(N)*N; kc = 2*np.pi/L*np.sqrt(n[:, None, None]**2+n[None, :, None]**2+n[None, None, :]**2)
    th = 2*np.pi*np.fft.fftfreq(N); dx = L/N; one = ((4/3)*(2-2*np.cos(th))-(2-2*np.cos(2*th))/12)/dx**2
    k2 = one[:, None, None]+one[None, :, None]+one[None, None, :]
    mH20 = 2*LX*float((np.abs(psi)**2).sum())+G2*float(np.mean(a[2]**2))
    occ = {}
    for name, m2 in [('massless', 0.), ('hartree', mH20)]:
        w = np.sqrt(np.where(k2 > 0, k2+m2, 1.)); nm = np.abs(w*psi+1j*pi)**2/(4*w)*L**3  # quanta per mode
        rows = []
        for lo, hi in [(1.8, 2.6), (2.6, 3.5), (3.5, 5.5), (5.5, 7.5), (7.5, 8.0), (8.5, 12.)]:
            sel = (kc >= lo) & (kc < hi)
            rows.append([lo, hi, int(sel.sum()), float(nm[sel].mean())])
        occ[name] = rows
    res['tau20_mH2'] = mH20; res['tau20_nminus_per_mode_by_k'] = occ
    res['zero_point_band_per_net'] = .5*(8**3-1.8**3)/(6*np.pi**2)/nQ
    json.dump(res, open(out, 'w'), indent=1); print(json.dumps(res, indent=1))

if __name__ == '__main__': main(*sys.argv[1:])
