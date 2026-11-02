"""qm_transfer: linear Boltzmann code for Q-matter objects + a kinetic X vapour on LCDM (own code, Claude, 2 Nov 2026). [C]

Synchronous gauge, Ma & Bertschinger 1995 (astro-ph/9506072) conventions [R].  Species: photons (temperature F_l and
polarisation G_l hierarchies, Thomson coupling), baryons, massless neutrinos (hierarchy), and a dark-matter sector of
one of three kinds:
  'cdm'   collisionless cold DM (LCDM reference),
  'tca'   objects + vapour as ONE fluid with the two-phase adiabatic sound speed c_eff^2 (only while the vapour is
          tightly coupled and in chemical equilibrium),
  'kin'   cold object fluid (charge contrast d_oc, number contrast d_N, velocity th_o) + the vapour X as a
          momentum-binned Boltzmann hierarchy dN_l(q_j), l = 0..L, on a comoving-momentum grid (per-charge
          normalisation N(q) dq = vapour charge per total charge), with the collision operator
            elastic:     Gamma_el(q) [ <f>_angle - f ]                in the object rest frame,
            absorption:  -Gamma_a(q) f,  emission +Gamma_a(q) f_eq(q; T_o)  in the object rest frame,
          T_o set by the objects' energy balance (zero heat capacity), linearised exactly (incl. the boost terms
          and the perturbation of the rates delta_Gamma = d_N/3 + 2 d_oc/3).  Momentum and charge are exchanged
          with the object fluid; energy is conserved (objects take (1-Delta) m_X per absorbed charge).
Approximation schemes (CLASS-like, Blas, Lesgourgues & Tram 2011 [R]): photon-baryon tight coupling at zeroth order
early; radiation-streaming approximation (RSA) for photons and neutrinos once k tau > K_RSA (photons: before
recombination the slow photon-baryon solution with baryons tied, theta_b = theta_g = -(h'+6eta')/2; after release the
baryons are free with Thomson drag); after a_end(k) = max(k tau = K_END, end of vapour collisions) a sub-horizon
Newtonian phase: cold objects and baryons by predictor-corrector marching, the vapour by the exact collisionless
Volterra (Gilbert) integral with its q-resolved state at a_end in the free term.
Time variable N = ln a.  Units: Mpc, c = 1.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import spherical_jn
from scipy.interpolate import CubicSpline

K_RSA = 50.0
K_END = 100.0


class CDM:
    kind = 'cdm'
    a_coll_end = 0.0

    def __init__(self):
        pass


class Mode:
    """one k-mode through all phases."""

    def __init__(self, cos, dm, k, Lg=10, Ln=12, rtol=1e-5, verbose=False, a_out=(1.0,), K_end=K_END, K_rsa=K_RSA,
                 dN_late=0.004):
        self.c, self.dm, self.k = cos, dm, k
        self.Lg, self.Ln = Lg, Ln
        self.rtol = rtol
        self.verbose = verbose
        self.a_out = sorted(a_out)
        self.K_end, self.K_rsa, self.dN_late = K_end, K_rsa, dN_late
        self.nfev = 0

    # ------------------------------------------------------------------ layout
    def layout(self, ph, nu, dmk):
        idx = {}
        n = 0

        def add(name, m=1):
            nonlocal n
            idx[name] = slice(n, n + m); n += m
        add('eta')
        if dmk == 'cdm':
            add('dD')
        elif dmk == 'tca':
            add('dD'); add('thD')
        else:
            add('dD'); add('dN'); add('thD')
        add('db')
        if ph == 'tca':
            add('dg'); add('thg')                  # theta_b = theta_g
        elif ph == 'full':
            add('thb'); add('Fg', self.Lg + 1); add('Gg', self.Lg + 1)
        elif ph == 'slow_free':
            add('thb')
        # slow_tied: baryon velocity slaved
        if nu == 'full':
            add('Fn', self.Ln + 1)
        if dmk == 'kin':
            add('V', self.dm.Nq * (self.dm.Lv + 1))
        idx['_n'] = n
        return idx

    # ------------------------------------------------------------------ RHS
    def rhs_factory(self, ph, nu, dmk, idx):
        c, dm, k = self.c, self.dm, self.k
        k2 = k * k
        Lg, Ln = self.Lg, self.Ln
        if dmk == 'kin':
            Nq, Lv = dm.Nq, dm.Lv
            wq = dm.wq
            qg = dm.q
            ell = np.arange(Lv + 1, dtype=float)

        nstate = idx['_n']

        def f(N, Y):
            self.nfev += 1
            Y2 = Y.reshape(nstate, -1)
            m = Y2.shape[1]
            a = math.exp(N)
            cH = c.calH(a)
            Adm, Ab, Ag, An = c.A_species(a)
            tau = c.tau(a)
            out = np.zeros_like(Y2)
            eta = Y2[idx['eta']][0]
            # ---- DM density/momentum (per rho_dm)
            if dmk == 'cdm':
                dD = Y2[idx['dD']][0]; drho_dm = dD; mom_dm = 0.0
            elif dmk == 'tca':
                dD = Y2[idx['dD']][0]; thD = Y2[idx['thD']][0]
                drho_dm = dD; mom_dm = thD
            else:
                K = dm.kin_at(a)
                doc = Y2[idx['dD']][0]; dNo = Y2[idx['dN']][0]; tho = Y2[idx['thD']][0]
                V = Y2[idx['V']].reshape(Nq, Lv + 1, m) / wq[:, None, None]     # state = bin populations w_j dN_j
                drho_dm = (K['eo'] * doc + (wq * K['E']) @ V[:, 0, :]) / K['eDM']
                mom_dm = (K['eo'] * tho + k * ((wq * qg / a) @ V[:, 1, :])) / K['eDM']
            # ---- non-RSA sums
            Sd = Adm * drho_dm + Ab * Y2[idx['db']][0]
            St = Adm * mom_dm
            Bd = 0.0; Bt = 0.0
            if ph == 'tca':
                dg = Y2[idx['dg']][0]; thg = Y2[idx['thg']][0]
                Sd = Sd + Ag * dg
                St = St + (4 / 3) * Ag * thg + Ab * thg
            elif ph == 'full':
                Fg = Y2[idx['Fg']]; Gg = Y2[idx['Gg']]; thb = Y2[idx['thb']][0]
                thg = 0.75 * k * Fg[1]
                Sd = Sd + Ag * Fg[0]
                St = St + (4 / 3) * Ag * thg + Ab * thb
            elif ph == 'slow_free':
                thb = Y2[idx['thb']][0]
                St = St + Ab * thb
                Bd += Ag; Bt += (4 / 3) * Ag
            else:  # slow_tied
                Bd += Ag; Bt += (4 / 3) * Ag + Ab
            if nu == 'full':
                Fn = Y2[idx['Fn']]
                thn = 0.75 * k * Fn[1]
                Sd = Sd + An * Fn[0]
                St = St + (4 / 3) * An * thn
            else:
                Bd += An; Bt += (4 / 3) * An
            # h' from the 00 constraint (RSA species: delta_r = 4 (calH h' - k^2 eta)/k^2)
            hp = (k2 * eta * (1 - 4 * Bd / k2) + Sd) / (0.5 * cH - 4 * cH * Bd / k2)
            etap = (St - 0.5 * Bt * hp) / (k2 + 3 * Bt)
            th_slow = -0.5 * (hp + 6 * etap)
            out[idx['eta']][0] = etap
            if Y.ndim == 1:
                self._last = (float(np.ravel(hp)[0]), float(np.ravel(etap)[0]))
            # ---- baryons and photons
            dk = c.dkappa(a)
            R = c.Rgb(a)
            cs2 = c.cs2(a)
            db = Y2[idx['db']][0]
            if ph == 'tca':
                out[idx['dg']][0] = -(4 / 3) * thg - (2 / 3) * hp
                out[idx['db']][0] = -thg - 0.5 * hp
                out[idx['thg']][0] = (-cH * thg + cs2 * k2 * db + R * k2 * dg / 4) / (1 + R)
            elif ph == 'full':
                Fo = out[idx['Fg']]; Go = out[idx['Gg']]
                out[idx['db']][0] = -thb - 0.5 * hp
                out[idx['thb']][0] = -cH * thb + cs2 * k2 * db + R * dk * (thg - thb)
                FL1 = (2 * Lg + 1) / (k * tau) * Fg[Lg] - Fg[Lg - 1]
                GL1 = (2 * Lg + 1) / (k * tau) * Gg[Lg] - Gg[Lg - 1]
                Fo[0] = -k * Fg[1] - (2 / 3) * hp
                Fo[1] = (k / 3) * (Fg[0] - 2 * Fg[2]) + (4 / (3 * k)) * dk * (thb - thg)
                Pi = Fg[2] + Gg[0] + Gg[2]
                Fo[2] = (2 / 5) * k * Fg[1] - (3 / 5) * k * Fg[3] + (4 / 15) * hp + (8 / 5) * etap - (9 / 10) * dk * Fg[2] + 0.1 * dk * (Gg[0] + Gg[2])
                for l in range(3, Lg):
                    Fo[l] = k / (2 * l + 1) * (l * Fg[l - 1] - (l + 1) * Fg[l + 1]) - dk * Fg[l]
                Fo[Lg] = k / (2 * Lg + 1) * (Lg * Fg[Lg - 1] - (Lg + 1) * FL1) - dk * Fg[Lg]
                Go[0] = -k * Gg[1] + dk * (-Gg[0] + 0.5 * Pi)
                for l in range(1, Lg):
                    Go[l] = k / (2 * l + 1) * (l * Gg[l - 1] - (l + 1) * Gg[l + 1]) - dk * Gg[l]
                Go[2] += dk * 0.1 * Pi
                Go[Lg] = k / (2 * Lg + 1) * (Lg * Gg[Lg - 1] - (Lg + 1) * GL1) - dk * Gg[Lg]
            elif ph == 'slow_free':
                out[idx['db']][0] = -thb - 0.5 * hp
                out[idx['thb']][0] = -cH * thb + cs2 * k2 * db + R * dk * (th_slow - thb)
            else:
                out[idx['db']][0] = -th_slow - 0.5 * hp
            # ---- neutrinos
            if nu == 'full':
                Fo = out[idx['Fn']]
                FL1 = (2 * Ln + 1) / (k * tau) * Fn[Ln] - Fn[Ln - 1]
                Fo[0] = -k * Fn[1] - (2 / 3) * hp
                Fo[1] = (k / 3) * (Fn[0] - 2 * Fn[2])
                Fo[2] = (2 / 5) * k * Fn[1] - (3 / 5) * k * Fn[3] + (4 / 15) * hp + (8 / 5) * etap
                for l in range(3, Ln):
                    Fo[l] = k / (2 * l + 1) * (l * Fn[l - 1] - (l + 1) * Fn[l + 1])
                Fo[Ln] = k / (2 * Ln + 1) * (Ln * Fn[Ln - 1] - (Ln + 1) * FL1)
            # ---- dark matter
            if dmk == 'cdm':
                out[idx['dD']][0] = -0.5 * hp
            elif dmk == 'tca':
                out[idx['dD']][0] = -thD - 0.5 * hp
                out[idx['thD']][0] = -cH * thD + dm.ceff2(a) * k2 * dD
            else:
                v = K['v']                       # (Nq,)
                ga, ge = K['ga'], K['ge']        # conformal rates [1/Mpc]
                D, Deq, dTN, Nb, Neq, Ek = K['D'], K['Deq'], K['dTNeq'], K['N'], K['Neq'], K['Ek']
                Vo = np.zeros_like(V)
                kv = (k * v)[:, None]
                # free streaming
                Vo[:, 0, :] = -kv * V[:, 1, :]
                for l in range(1, Lv):
                    Vo[:, l, :] = kv / (2 * l + 1) * (l * V[:, l - 1, :] - (l + 1) * V[:, l + 1, :])
                VL1 = (2 * Lv + 1) / (kv * tau) * V[:, Lv, :] - V[:, Lv - 1, :]
                Vo[:, Lv, :] = kv / (2 * Lv + 1) * (Lv * V[:, Lv - 1, :] - (Lv + 1) * VL1)
                # metric
                Vo[:, 0, :] += (hp / 6) * D[:, None]
                if Lv >= 2:
                    Vo[:, 2, :] -= (hp / 15 + 0.4 * etap) * D[:, None]
                # collisions
                wE = wq * ga * (K['Delta'] + Ek)
                dTo = (wE @ V[:, 0, :]) / K['dTden']
                dGam = dNo / 3 + 2 * doc / 3
                c0 = ga[:, None] * (dTN[:, None] * dTo + (Neq - Nb)[:, None] * dGam - V[:, 0, :])
                gt = (ga + ge)[:, None]
                src1 = (ge * D + ga * Deq + ga * K['dlnGa'] * (Neq - Nb)) / (3 * k * v)
                c1 = -src1[:, None] * tho - gt * V[:, 1, :]
                Vo[:, 0, :] += c0
                Vo[:, 1, :] += c1
                Vo[:, 2:, :] -= (ga + ge)[:, None, None] * V[:, 2:, :]
                djc = -(wq @ c0)
                Smom = -k * ((wq * qg / a) @ c1)
                Co, eo_, jc = K['Co'], K['eo'], K['jc']
                out[idx['dD']][0] = -tho - 0.5 * hp + (djc - doc * jc) / Co
                out[idx['dN']][0] = -tho - 0.5 * hp
                out[idx['thD']][0] = -cH * tho + (Smom - K['eps_o'] * jc * tho) / eo_
                out[idx['V']] = (Vo * wq[:, None, None]).reshape(Nq * (Lv + 1), m)
            out /= cH
            if Y.ndim == 1:
                return out[:, 0]
            return out

        def jac(N, Y):
            return f(N, np.eye(nstate))

        return f, jac

    # ------------------------------------------------------------------ metric helpers (scalar state)
    def metric(self, ph, nu, dmk, idx, N, y):
        """h', eta' at the (scalar) state y"""
        f, _ = self.rhs_factory(ph, nu, dmk, idx)
        f(N, y)
        return self._last

    # ------------------------------------------------------------------ schedule
    def schedule(self):
        c, dm, k = self.c, self.dm, self.k
        lna = c.lna_tab
        taus = c.tau_tab
        a_i = c.a_of_tau(0.05 / k)
        a_i = min(a_i, 1e-6 * min(1.0, 1.0))
        t_rsa = self.K_rsa / k
        a_rsa = c.a_of_tau(t_rsa) if t_rsa < c.tau0 else 1.0
        # photon TCA end: k tau_c > 0.005 or tau_c/tau > 0.005
        lg = np.linspace(math.log(a_i), 0.0, 3000)
        aa = np.exp(lg)
        dk = np.array([c.dkappa(x) for x in aa])
        tt = np.array([c.tau(x) for x in aa])
        bad = (k / dk > 0.005) | (1.0 / (dk * tt) > 0.005)
        a_tcaend = float(aa[np.argmax(bad)]) if bad.any() else 1.0
        # baryon release: R * dkappa / calH < 300
        rel = np.array([c.Rgb(x) * c.dkappa(x) / c.calH(x) for x in aa]) < 300
        a_rel = float(aa[np.argmax(rel)]) if rel.any() else 1.0
        a_end = max(c.a_of_tau(self.K_end / k) if self.K_end / k < c.tau0 else 1.0, getattr(dm, 'a_coll_end', 0.0), a_rsa)
        a_end = min(a_end, 1.0)
        self.a_rel = a_rel
        # DM switch
        if dm.kind == 'cdm':
            a_kin = None
        else:
            a_kin = dm.a_kin_for_k(k, c)
            if a_kin is not None and a_kin <= a_i:
                a_kin = a_i
        events = {'a_i': a_i, 'a_rsa': a_rsa, 'a_tcaend': a_tcaend, 'a_rel': a_rel, 'a_end': a_end, 'a_kin': a_kin}
        return events

    def phase_at(self, a, ev):
        # photon mode
        if a < min(ev['a_tcaend'], ev['a_rsa']):
            ph = 'tca'
        elif a < ev['a_rsa']:
            ph = 'full'
        else:
            ph = 'slow_tied' if a < ev['a_rel'] else 'slow_free'
        nu = 'full' if a < ev['a_rsa'] else 'rsa'
        if self.dm.kind == 'cdm':
            dmk = 'cdm'
        elif self.dm.kind == 'frozen':
            dmk = 'kin'
        else:
            dmk = 'kin' if (ev['a_kin'] is not None and a >= ev['a_kin']) else 'tca'
        return ph, nu, dmk

    # ------------------------------------------------------------------ initial conditions and remaps
    def initial(self, a, idx, ph, nu, dmk):
        c, k, dm = self.c, self.k, self.dm
        x = k * c.tau(a)
        Rn = c.Onu / c.Or
        y = np.zeros(idx['_n'])
        C = 1.0
        y[idx['eta']] = 2 * C - (5 + 4 * Rn) / (6 * (15 + 4 * Rn)) * C * x * x
        dc = -0.5 * C * x * x
        y[idx['dD']] = dc
        y[idx['db']] = dc
        thg = -C * k * x**3 / 18.0
        if ph == 'tca':
            y[idx['dg']] = -(2 / 3) * C * x * x
            y[idx['thg']] = thg
        if nu == 'full':
            Fn = np.zeros(self.Ln + 1)
            Fn[0] = -(2 / 3) * C * x * x
            thn = -(23 + 4 * Rn) / (18 * (15 + 4 * Rn)) * C * k * x**3
            Fn[1] = 4 * thn / (3 * k)
            Fn[2] = 2 * 4 * C * x * x / (3 * (15 + 4 * Rn))
            y[idx['Fn']] = Fn
        if dmk == 'kin':
            K = dm.kin_at(a)
            dQ = dc
            V = np.zeros((dm.Nq, dm.Lv + 1))
            V[:, 0] = -(K['dNdlna'] / 3.0 + K['D'] / 3.0) * dQ
            if dm.kind == 'frozen':
                v = K['v']
                thn = -(23 + 4 * Rn) / (18 * (15 + 4 * Rn)) * C * k * x**3
                sig = 4 * C * x * x / (3 * (15 + 4 * Rn))
                V[:, 1] = -K['D'] / (3 * k * v) * thn * v**2
                if dm.Lv >= 2:
                    V[:, 2] = -0.5 * K['D'] * sig * v**2
            y[idx['V']] = (V * dm.wq[:, None]).reshape(-1)
            Co = K['Co']
            y[idx['dD']] = (dQ - dm.wq @ V[:, 0]) / Co if Co > 0 else dQ
            y[idx['dN']] = dQ
            y[idx['thD']] = 0.0
        return y

    def remap(self, y, idx_old, mode_old, idx_new, mode_new, N):
        c, k, dm = self.c, self.k, self.dm
        a = math.exp(N)
        ph0, nu0, dm0 = mode_old
        ph1, nu1, dm1 = mode_new
        z = np.zeros(idx_new['_n'])
        z[idx_new['eta']] = y[idx_old['eta']]
        z[idx_new['db']] = y[idx_old['db']]
        # photons
        if ph0 == 'tca' and ph1 == 'full':
            Fg = np.zeros(self.Lg + 1)
            Fg[0] = y[idx_old['dg']][0]
            thg = y[idx_old['thg']][0]
            Fg[1] = 4 * thg / (3 * k)
            z[idx_new['Fg']] = Fg
            z[idx_new['thb']] = thg
        elif ph0 == ph1 and ph0 == 'full':
            z[idx_new['Fg']] = y[idx_old['Fg']]; z[idx_new['Gg']] = y[idx_old['Gg']]; z[idx_new['thb']] = y[idx_old['thb']]
        elif ph0 == ph1 and ph0 == 'tca':
            z[idx_new['dg']] = y[idx_old['dg']]; z[idx_new['thg']] = y[idx_old['thg']]
        elif ph1 == 'slow_free':
            if ph0 == 'full':
                z[idx_new['thb']] = y[idx_old['thb']]
            elif ph0 == 'tca':
                z[idx_new['thb']] = y[idx_old['thg']]
            elif ph0 == 'slow_free':
                z[idx_new['thb']] = y[idx_old['thb']]
            else:   # from tied: theta_b = theta_slow
                hp, etap = self.metric(ph0, nu0, dm0, idx_old, N, y)
                z[idx_new['thb']] = -0.5 * (hp + 6 * etap)
        if nu0 == 'full' and nu1 == 'full':
            z[idx_new['Fn']] = y[idx_old['Fn']]
        # dark matter
        if dm0 == dm1:
            for key in ('dD', 'dN', 'thD', 'V'):
                if key in idx_old and key in idx_new:
                    z[idx_new[key]] = y[idx_old[key]]
        elif dm0 == 'tca' and dm1 == 'kin':
            dD = y[idx_old['dD']][0]; thD = y[idx_old['thD']][0]
            K = dm.kin_at(a)
            V = np.zeros((dm.Nq, dm.Lv + 1))
            # local two-phase equilibrium: dN = N [r dD + (q^2/(2 a^2 T) - 3/2) dT/T] (dT/T = dD / Dfac), comoving dipole
            r, Dfac, T = dm.tca_local(a)
            V[:, 0] = K['N'] * (r * dD + (K['Ek'] / T - 1.5) * dD / Dfac)
            V[:, 1] = -K['D'] / (3 * k * K['v']) * thD
            z[idx_new['V']] = (V * dm.wq[:, None]).reshape(-1)
            z[idx_new['dD']] = (dD - dm.wq @ V[:, 0]) / K['Co']
            z[idx_new['dN']] = dD
            z[idx_new['thD']] = thD
        else:
            raise ValueError('bad DM remap %s -> %s' % (dm0, dm1))
        return z

    # ------------------------------------------------------------------ driver
    def run(self):
        c, k, dm = self.c, self.k, self.dm
        ev = self.schedule()
        self.ev = ev
        pts = sorted(set([ev['a_i'], ev['a_end']] + [ev[x] for x in ('a_rsa', 'a_tcaend', 'a_rel', 'a_kin')
                                                         if ev[x] is not None and ev['a_i'] < ev[x] < ev['a_end']]))
        mode = self.phase_at(pts[0] * 1.0000001, ev)
        idx = self.layout(*mode)
        y = self.initial(pts[0], idx, *mode)
        self.log = []
        for j in range(len(pts) - 1):
            N0, N1 = math.log(pts[j]), math.log(pts[j + 1])
            newmode = self.phase_at(pts[j] * 1.0000001, ev)
            if newmode != mode:
                idx_new = self.layout(*newmode)
                y = self.remap(y, idx, mode, idx_new, newmode, N0)
                idx, mode = idx_new, newmode
            f, jac = self.rhs_factory(*mode, idx)
            scale = max(abs(y[idx['dD']][0]), 1e-30)
            sol = solve_ivp(f, (N0, N1), y, method='BDF', jac=jac, rtol=self.rtol, atol=1e-8 * scale)
            if not sol.success:
                raise RuntimeError('solve_ivp failed k=%g phase %s: %s' % (k, mode, sol.message))
            y = sol.y[:, -1]
            self.log.append(dict(a0=pts[j], a1=pts[j + 1], mode=mode, nsteps=int(sol.t.size), n=idx['_n'], njev=int(sol.njev), nlu=int(sol.nlu)))
        self.y_end, self.idx_end, self.mode_end = y, idx, mode
        return self.late(y, idx, mode, math.log(ev['a_end']))

    # ------------------------------------------------------------------ Newtonian late phase
    def late(self, y, idx, mode, Nend):
        c, k, dm = self.c, self.k, self.dm
        a_end = math.exp(Nend)
        f, _ = self.rhs_factory(*mode, idx)
        dy = f(Nend, y) * c.calH(a_end)          # d/dtau
        hp, etap = self.metric(*mode, idx, Nend, y)
        fb, fdm = c.Ob / c.Om, c.Odm / c.Om
        dmk = mode[2]
        if dmk == 'kin':
            K = dm.kin_at(a_end)
            Vend = y[idx['V']].reshape(dm.Nq, dm.Lv + 1) / dm.wq[:, None]
            eDM = K['eDM']
            fo = fdm * K['eo'] / eDM                     # object share of matter
            fg = fdm * float(dm.wq @ (K['N'] * K['E'])) / eDM
            Nv = float(dm.wq @ K['N'])
            do0, dpo0 = y[idx['dD']][0], dy[idx['dD']][0]
            # Newtonian frame shift of the vapour dipole: theta_N = theta_S + h'/2
            Vend[:, 1] += -K['D'] / (3 * k * K['v']) * (0.5 * hp)
        else:
            fo, fg = fdm, 0.0
            do0, dpo0 = y[idx['dD']][0], dy[idx['dD']][0]
        db0 = y[idx['db']][0]
        dpb0 = dy[idx['db']][0]
        if mode[0] == 'slow_tied':
            dpb0 = 0.0
        a_rel = self.a_rel
        dN = self.dN_late
        Ngrid = np.linspace(Nend, 0.0, max(int(round(-Nend / dN)), 60) + 1)
        dNs = Ngrid[1] - Ngrid[0]
        a = np.exp(Ngrid)
        Hc = c.calH_v(a)
        Oma = np.array([c.Omega_m_a(z) for z in a])
        cs2 = np.array([c.cs2(z) for z in a])
        Cc = np.concatenate([[0.0], np.cumsum(0.5 * (1 / (a[1:] * Hc[1:]) + 1 / (a[:-1] * Hc[:-1])) * dNs)])   # S = int dtau / a
        w_src = a * Hc * 1.5 * Oma
        A = a / Hc * 1.5 * Hc * Hc * Oma
        B = 1 / (a * Hc)
        P = a / Hc * cs2 * k * k
        n = len(Ngrid)
        do = np.zeros(n); dg = np.zeros(n); db = np.zeros(n); src = np.zeros(n)
        uo = a[0] * dpo0
        ub = a[0] * dpb0
        do[0], db[0] = do0, db0
        bar_free = a >= a_rel
        if fg > 0:
            # free term: exact free streaming of the q-resolved (discrete-velocity) state at a_end, on the SAME q grid
            # as the hierarchy (the consistent continuation of the discrete system; moments exact by construction)
            lq = np.log(dm.q)
            qf = np.exp(np.linspace(lq[0], lq[-1], 6 * dm.Nq))
            wf = qf * (np.log(qf[1]) - np.log(qf[0]))
            Nf = np.exp(CubicSpline(lq, np.log(np.maximum(K['N'], 1e-300)))(np.log(qf)))
            ns = min(n, 800)
            Ss = np.linspace(0.0, Cc[-1], ns)
            fs = np.zeros(ns)
            for l in range(dm.Lv + 1):
                J = spherical_jn(l, k * np.outer(Ss, dm.q))
                fs += (2 * l + 1) * (-1)**l * (J @ (dm.wq * Vend[:, l]))
            free = np.interp(Cc, Ss, fs) / Nv
            # kernel W(y) = int N j0(k q y) / Nv, tabulated in S
            Sgrid = np.linspace(0, Cc[-1] * 1.0001, 4000)
            W = (np.sinc(np.outer(Sgrid, k * qf) / math.pi) @ (wf * Nf)) / float(wf @ Nf)
            dg[0] = free[0]
        src[0] = fo * do[0] + fg * dg[0] + fb * db[0]
        outs = {}
        for i in range(1, n):
            if fg > 0:
                S = Cc[i] - Cc[:i + 1]
                Wv = np.interp(S, Sgrid, W)
                dg[i] = free[i] + np.trapezoid(w_src[:i + 1] * src[:i + 1] * S * Wv, dx=dNs)
            dop = do[i - 1] + dNs * B[i - 1] * uo
            uop = uo + dNs * A[i - 1] * src[i - 1]
            if bar_free[i - 1] or bar_free[i]:
                dbp = db[i - 1] + dNs * B[i - 1] * ub
                ubp = ub + dNs * (A[i - 1] * src[i - 1] - P[i - 1] * db[i - 1])
            else:
                dbp, ubp = db[i - 1], 0.0
            srcp = fo * dop + fg * dg[i] + fb * dbp
            do[i] = do[i - 1] + 0.5 * dNs * (B[i - 1] * uo + B[i] * uop)
            uo = uo + 0.5 * dNs * (A[i - 1] * src[i - 1] + A[i] * srcp)
            if bar_free[i - 1] or bar_free[i]:
                db[i] = db[i - 1] + 0.5 * dNs * (B[i - 1] * ub + B[i] * ubp)
                ub = ub + 0.5 * dNs * (A[i - 1] * src[i - 1] - P[i - 1] * db[i - 1] + A[i] * srcp - P[i] * dbp)
            else:
                db[i] = db[i - 1]
            src[i] = fo * do[i] + fg * dg[i] + fb * db[i]
        self.late_arrays = dict(a=a, do=do, dg=dg, db=db, fo=fo, fg=fg, fb=fb)
        res = {}
        for ao in self.a_out:
            j = int(np.argmin(abs(a - ao)))
            res[ao] = float(fo * do[j] + fg * dg[j] + fb * db[j])
        return res


def transfer_ratio(cos, dm, ks, ref_cache=None, **kw):
    """P_dm-model / P_LCDM at the a_out values (same initial amplitude)."""
    out = []
    for k in ks:
        m = Mode(cos, dm, k, **kw)
        r = m.run()
        if ref_cache is not None and k in ref_cache:
            r0 = ref_cache[k]
        else:
            m0 = Mode(cos, CDM(), k, **kw)
            r0 = m0.run()
            if ref_cache is not None:
                ref_cache[k] = r0
        out.append({str(ao): (r[ao] / r0[ao])**2 for ao in r})
    return out
