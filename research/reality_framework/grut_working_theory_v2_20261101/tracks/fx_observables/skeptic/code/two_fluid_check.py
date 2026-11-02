"""fx_observables SKEPTIC (Claude, 1 Nov 2026): is the single-fluid (tight-coupling) treatment of the equilibrium
vapour valid where it produces the claimed acoustic suppression?  [C] own code.

The track (code/vapour_transfer.py) treats objects + vapour as ONE fluid with c_eff^2 = gamma n_v T/rho_DM until
momentum decoupling a_m.  That requires the vapour mean free path to be short compared with the wavelength,
eps_TC = k v_bar / (a Gamma_m) << 1 (comoving k, physical Gamma_m), during the acoustic phase k r_s(a) >~ 1.  Here the
vapour moves at v_bar ~ 0.07-0.15 c while the loaded fluid's sound speed is c_eff ~ 0.01 c, so the diffusion length can
exceed the sound horizon.

This file (i) evaluates eps_TC at acoustic onset for each scenario, and (ii) re-integrates the pre-a_m era as TWO fluids
(pressureless objects + vapour fluid), coupled by momentum exchange Gamma_m (vapour) / Gamma_m rho_v/rho_o (objects)
and by chemical relaxation Gamma_c toward local two-phase equilibrium delta_v = r delta_o, with the vapour pressure
normalised so that the tight-coupling limit is exactly the track's c_eff^2.  Variant 'visc' adds a vapour viscosity
nu = c_v^2/(a Gamma_m + k c_v) (diffusive in the collisional limit, free-streaming rate in the collisionless one).
After a_m the track's own Gilbert marching is used unchanged in form (copied), but started from the separate object and
vapour states instead of a common one.  The vapour history is the track's history_B (copied loop, full arrays kept).
The CDM reference runs through the same two-fluid pipeline with f_v = 0.
"""
import importlib.util, json, math, os, sys, time
import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
TRACK = os.path.abspath(os.path.join(HERE, '..', '..'))
spec = importlib.util.spec_from_file_location('vt', os.path.join(TRACK, 'code', 'vapour_transfer.py'))
vt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(vt)          # read-only import of the track's module (main() not run)


def history_full(Q, SQ, stick, sigma_m=None, a_lo=1e-14, a_hi=1.0, n=8000):
    """copy of vt.history_B's loop, keeping the full arrays."""
    p = vt.obj_props(Q)
    Delta = 1 - p['Omega']
    R = p['R']
    sig_c = stick * math.pi * R * R
    sig_m = sigma_m if sigma_m is not None else math.pi * R * R
    lna = np.linspace(math.log(a_lo), math.log(a_hi), n)
    Nv = None; Tv = None; af = None; am = None; Tf = None
    rows = []; frozen = False
    for i, la in enumerate(lna):
        a = math.exp(la)
        nQ = vt.NQ0 / a**3
        Teq, nveq = vt.vapour_T(nQ, SQ, Delta)
        if Teq is None or (Nv is None and Teq > 0.05):
            continue
        Neq = min(nveq / nQ, 1.0)
        if Nv is None:
            Nv, Tv = Neq, Teq
        nobj = (1 - Nv) * nQ / Q
        Hh = vt.H_mX(a)
        if not frozen:
            Tv = Teq
            Gc = nobj * sig_c * math.sqrt(8 * Tv / math.pi) / Hh
            dl = lna[i] - lna[i - 1] if i > 0 else 0.0
            Nv = Neq + (Nv - Neq) * math.exp(-Gc * dl) if Gc * dl < 700 else Neq
            if Gc < 1.0:
                frozen = True; af = a; Tf = Tv
        else:
            Tv = Tf * (af / a)**2
            Gc = nobj * sig_c * math.sqrt(8 * Tv / math.pi) / Hh
        Gm = nobj * sig_m * math.sqrt(8 * Tv / math.pi) / Hh
        if am is None and frozen and Gm < 1.0:
            am = a
        rows.append((a, Nv, Tv, Gc, Gm))
        if am is not None and a > 10 * am:
            break
    r = np.array(rows)
    hist = dict(Q=Q, SQ=SQ, s=stick, Delta=Delta, EQ=p['EQ'], af=af, am=am, Tf=Tf,
                a=r[:, 0], Nv=r[:, 1], Tv=r[:, 2], Gc=r[:, 3], Gm=r[:, 4])
    # the track's own history (for f_v_mass, sigma_p_today, _ceff2, r_s) to make sure we use identical inputs
    hB = vt.history_B(Q, SQ, stick, sigma_m=sigma_m)
    hist['track'] = hB
    return hist


def make_funcs(H):
    la = np.log(H['a'])
    a0 = H['a'][0]
    EQ, Delta, af = H['EQ'], H['Delta'], H['af']

    def at(a):
        if a <= a0:
            i0 = 0
            Nv, Tv, Gc, Gm = H['Nv'][0], H['Tv'][0], H['Gc'][0] * a0 / a, H['Gm'][0] * a0 / a   # Gamma/H ~ 1/a in RD
        else:
            l = math.log(a)
            Nv = float(np.interp(l, la, H['Nv'])); Tv = float(np.exp(np.interp(l, la, np.log(H['Tv']))))
            Gc = float(np.exp(np.interp(l, la, np.log(np.maximum(H['Gc'], 1e-300)))))
            Gm = float(np.exp(np.interp(l, la, np.log(H['Gm']))))
        rho = (1 - Nv) * EQ + Nv * (1 + 1.5 * Tv)
        fv = Nv * (1 + 1.5 * Tv) / rho
        x = Delta / Tv
        if a < af:
            gam = (2.5 + x) / (1.5 + x - x / (2.5 + x))
            r = (1.5 + x) / (1.5 + x - x / (2.5 + x))
            cv2 = gam * Tv / r * (Nv / rho) / fv          # so that fv * cv2 * r = gam Nv Tv / rho (track's c_eff^2)
            Gch = Gc
        else:
            cv2 = 5 / 3 * Tv * (Nv / rho) / fv
            r = 1.0
            Gch = 0.0
        return fv, cv2, r, Gch, Gm, Tv
    return at


def growth2(k, H, reference=False, visc=False, dN=0.004, boost=1.0, rtol=1e-8):
    at = make_funcs(H)
    tr = H['track']
    a_i = vt.a_of_ktau(k, 10.0)
    a_m = tr['a_m']
    fb, fc = vt.Ob / vt.Om, vt.Oc / vt.Om

    def rhs(N, y):
        a = math.exp(N)
        Hc = a * vt.H_Mpc(a)
        do, Uo, dv, Uv = y
        if reference:
            fv, cv2, r, Gch, Gm = 0.0, 0.0, 1.0, 0.0, 0.0
        else:
            fv, cv2, r, Gch, Gm, Tv = at(a)
            Gch, Gm = Gch * boost, Gm * boost
        src = fc * ((1 - fv) * do + fv * dv)
        acc = 1.5 * Hc * Hc * vt.Om_a(a)
        Gmo = Gm * fv / (1 - fv)
        ex = Gch * (dv - r * do)
        ddo = Uo / (a * Hc) + ex * fv / (1 - fv)
        dUo = a / Hc * acc * src - Gmo * (Uo - Uv)
        ddv = Uv / (a * Hc) - ex
        dUv = a / Hc * (acc * src - cv2 * k * k * dv) - Gm * (Uv - Uo)
        if visc and not reference:
            cv = math.sqrt(cv2)
            nu = cv2 / (Gm * Hc + k * cv)            # conformal: a Gamma_m = Gm * a H = Gm * Hc
            dUv -= nu * k * k * Uv / Hc
        return [ddo, dUo, ddv, dUv]

    Ni = math.log(a_i)
    Hci = a_i * vt.H_Mpc(a_i)
    d0 = math.log(6.0)
    r0 = 1.0 if reference else at(a_i)[2]
    y0 = [d0, a_i * Hci, r0 * d0, a_i * Hci * r0]
    Nm = math.log(a_m)
    if Ni < Nm:
        sol = solve_ivp(rhs, (Ni, Nm), y0, method='Radau', rtol=rtol, atol=1e-12)
        ym = sol.y[:, -1]
    else:
        Nm = Ni                       # mode starts after a_m: no pre-a_m era (as in the track)
        ym = np.array(y0)
    # ---- after a_m: the track's marching (copied), objects cold, gas collisionless (Gilbert), baryons from a_rec
    fv = 0.0 if reference else tr['f_v_mass']
    sigp = 0.0 if reference else tr['sigma_p_today']
    Ngrid = np.linspace(Nm, 0.0, max(int(round((0.0 - Nm) / dN)), 50) + 1)
    dNs = Ngrid[1] - Ngrid[0]
    a = np.exp(Ngrid)
    Hc = a * np.array([vt.H_Mpc(z) for z in a])
    Oma = np.array([vt.Om_a(z) for z in a])
    C = np.concatenate([[0.0], np.cumsum(0.5 * (1 / (a[1:] * Hc[1:]) + 1 / (a[:-1] * Hc[:-1])) * dNs)])
    w_src = a * Hc * 1.5 * Oma
    A = a / Hc * 1.5 * Hc * Hc * Oma
    B = 1 / (a * Hc)
    n = len(Ngrid)
    do = np.zeros(n); dg = np.zeros(n); db = np.zeros(n); src = np.zeros(n)
    uo, ub = ym[1], 0.0
    do[0] = ym[0]; dg[0] = ym[2]; db[0] = 0.0
    fo, fg = fc * (1 - fv), fc * fv
    bar = a >= vt.A_REC
    src[0] = fo * do[0] + fg * dg[0]
    for i in range(1, n):
        S = C[i] - C[:i + 1]
        damp = np.exp(-0.5 * (k * sigp * S)**2)
        free = (ym[2] + ym[3] * S[0]) * damp[0]
        dg[i] = free + np.trapezoid(w_src[:i + 1] * src[:i + 1] * S * damp, dx=dNs)
        dop = do[i - 1] + dNs * B[i - 1] * uo
        uop = uo + dNs * A[i - 1] * src[i - 1]
        if bar[i - 1] or bar[i]:
            dbp = db[i - 1] + dNs * B[i - 1] * ub
            ubp = ub + dNs * A[i - 1] * src[i - 1]
        else:
            dbp, ubp = 0.0, 0.0
        srcp = fo * dop + fg * dg[i] + (fb * dbp if bar[i] else 0.0)
        do[i] = do[i - 1] + 0.5 * dNs * (B[i - 1] * uo + B[i] * uop)
        uo = uo + 0.5 * dNs * (A[i - 1] * src[i - 1] + A[i] * srcp)
        if bar[i - 1] or bar[i]:
            db[i] = db[i - 1] + 0.5 * dNs * (B[i - 1] * ub + B[i] * ubp)
            ub = ub + 0.5 * dNs * (A[i - 1] * src[i - 1] + A[i] * srcp)
        src[i] = fo * do[i] + fg * dg[i] + (fb * db[i] if bar[i] else 0.0)
    return (vt.Oc * (1 - fv) * do[-1] + vt.Oc * fv * dg[-1] + vt.Ob * db[-1]) / vt.Om, dict(do_am=float(ym[0]), dv_am=float(ym[2]))


def tc_diagnostic(H, ks):
    """eps_TC = k v_bar/(a Gamma_m) at acoustic onset k r_s(a) = 1 (or at a_m if onset is later)."""
    at = make_funcs(H)
    tr = H['track']
    am = tr['a_m']
    aa = np.geomspace(H['a'][0], am, 3000)
    c2 = np.array([at(z)[1] * at(z)[0] * at(z)[2] for z in aa])        # fv cv2 r = c_eff^2
    integ = np.sqrt(c2) / (aa * aa * np.array([vt.H_Mpc(z) for z in aa]))
    rs = np.concatenate([[0.0], np.cumsum(0.5 * (integ[1:] + integ[:-1]) * np.diff(aa))])
    out = {}
    for k in ks:
        j = np.searchsorted(k * rs, 1.0)
        onset = j < len(aa)
        a_on = aa[min(j, len(aa) - 1)]
        fv, cv2, r, Gch, Gm, Tv = at(a_on)
        vbar = math.sqrt(8 * Tv / math.pi)
        eps = k * vbar / (a_on * Gm * vt.H_Mpc(a_on))
        # vapour diffusion length up to a_on: lambda_D^2 = int v^2/(3 a Gamma_m) dtau
        sel = aa <= a_on
        vv = np.array([8 * at(z)[5] / math.pi for z in aa[sel]])
        Gmz = np.array([at(z)[4] for z in aa[sel]])
        Hc = aa[sel] * np.array([vt.H_Mpc(z) for z in aa[sel]])
        dtau = 1 / (aa[sel] * Hc)                                     # d tau / d a
        lamD = math.sqrt(np.trapezoid(vv / (3 * Gmz * Hc) * dtau, aa[sel])) if sel.sum() > 2 else 0.0
        out[str(k)] = dict(acoustic_onset_before_am=bool(onset), a_onset=float(a_on), eps_TC_at_onset=float(eps),
                           k_lambda_D_at_onset=float(k * lamD), r_s_full_Mpc=float(rs[-1]))
    return out


def control(fn):
    """forced tight coupling (rates x1e6) must reproduce the track's single-fluid P(k); also eps_TC at a_m."""
    trk = {}
    for r in json.load(open(os.path.join(TRACK, 'results', 'vapour_transfer.json')))['B']:
        if 'transfer' in r:
            trk[(r['Q'], r['S_over_Q'], r['sticking'])] = {str(t['k']): t['P_ratio'] for t in r['transfer']}
    out = []
    for Q, SQ, s in ((1e14, 1.235, 1.0), (1e12, 1.235, 1.0), (2e5, 1.235, 1.0)):
        H = history_full(Q, SQ, s, sigma_m=332.0 if Q == 2e5 else None)
        for k in (10.0, 50.0, 100.0):
            d0, _ = growth2(k, H, reference=True)
            d1, _ = growth2(k, H, boost=1e6)
            d2, _ = growth2(k, H, boost=1e6, dN=0.002)
            d0b, _ = growth2(k, H, reference=True, dN=0.002)
            row = dict(Q=Q, S_over_Q=SQ, k=k, two_fluid_forced_TC=(d1 / d0)**2, half_step=(d2 / d0b)**2,
                       track=trk[(Q, SQ, s)][str(k)], ref_vs_track_ref=d0 / vt.growth(k, H['track'], reference=True))
            print(row, flush=True)
            out.append(row)
    json.dump(out, open(fn, 'w'), indent=1, default=float)


def main():
    if sys.argv[1].endswith('control.json'):
        return control(sys.argv[1])
    t0 = time.time()
    cases = [(2e5, 1.235, 1.0), (2e5, 1.235, 1e-2), (1e6, 1.235, 1.0), (1e8, 1.235, 1.0), (1e10, 1.235, 1.0),
             (1e12, 1.235, 1.0), (1e14, 1.235, 1.0), (2e5, 2.07, 1.0), (1e8, 2.07, 1.0), (1e12, 2.07, 1.0)]
    if len(sys.argv) > 2:
        sel = [int(x) for x in sys.argv[2].split(',')]
        cases = [cases[i] for i in sel]
    ks = [1.0, 3.0, 10.0, 20.0, 50.0, 100.0, 300.0]
    trk = {}
    for r in json.load(open(os.path.join(TRACK, 'results', 'vapour_transfer.json')))['B']:
        if 'transfer' in r:
            trk[(r['Q'], r['S_over_Q'], r['sticking'])] = {str(t['k']): t['P_ratio'] for t in r['transfer']}
    out = []
    for Q, SQ, s in cases:
        H = history_full(Q, SQ, s, sigma_m=332.0 if Q == 2e5 else None)
        tc = tc_diagnostic(H, ks)
        row = dict(Q=Q, S_over_Q=SQ, sticking=s, a_f=H['af'], a_m=H['track']['a_m'], f_v_mass=H['track']['f_v_mass'],
                   r_s_track_Mpc=H['track']['r_s_Mpc'], tc=tc, P={})
        for k in ks:
            d0, _ = growth2(k, H, reference=True)
            d1, st1 = growth2(k, H)
            d2, st2 = growth2(k, H, visc=True)
            row['P'][str(k)] = dict(two_fluid=(d1 / d0)**2, two_fluid_visc=(d2 / d0)**2,
                                    track_single_fluid=trk.get((Q, SQ, s), {}).get(str(k)),
                                    delta_v_over_delta_o_at_am=st1['dv_am'] / st1['do_am'])
            print('Q=%.0e SQ=%.3f s=%g k=%g: track %.4f | 2-fluid %.4f | +visc %.4f | eps_TC(onset)=%.3g kLD=%.3g dv/do(am)=%.3f  %.0fs' % (
                Q, SQ, s, k, row['P'][str(k)]['track_single_fluid'] or -1, row['P'][str(k)]['two_fluid'],
                row['P'][str(k)]['two_fluid_visc'], tc[str(k)]['eps_TC_at_onset'], tc[str(k)]['k_lambda_D_at_onset'],
                row['P'][str(k)]['delta_v_over_delta_o_at_am'], time.time() - t0), flush=True)
        out.append(row)
    fn = sys.argv[1]
    json.dump(dict(scope=__doc__.strip().splitlines()[0], rows=out), open(fn, 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
