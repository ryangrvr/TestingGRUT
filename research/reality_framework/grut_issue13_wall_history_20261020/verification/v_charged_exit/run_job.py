"""Run one named verification job.  usage: python3 run_job.py <name>
Writes results/<name>.json (config + time series); field states go to the scratch directory."""
import json, os, sys, time
import numpy as np
import vlat
from vlat import Lattice

HERE = os.path.dirname(os.path.abspath(__file__))
SCR = os.environ.get('VCE_SCRATCH',
                     '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/vce_states')
os.makedirs(SCR, exist_ok=True)

DT = 0.05
JOBS = {}
def job(name, **kw):
    JOBS[name] = kw

# charged exit (A), kappa = 2, tau_i = 40, a: 1 -> 3, 1024^2 dx = 0.25 (L = 256)
for s in (101, 102, 103, 104, 105, 106):
    job('A_k2_t40_N1024_s%d' % s, kind='A', kappa=2.0, tau_i=40.0, N=1024, dx=0.25, seed=s, a_end=3.0, dtau_out=1.0)
# neutral quench (B)
for s in (201, 202):
    job('B_t40_N1024_s%d' % s, kind='B', tau_i=40.0, N=1024, dx=0.25, seed=s, a_end=3.0, dtau_out=1.0)
# sub-critical +v control, n = 0.9 n_c, tau_i = 150, a 1 -> 2.2
job('C1_k0.9_t150_N1024_s401', kind='C1', frac=0.9, tau_i=150.0, N=1024, dx=0.25, seed=401, a_end=2.2, dtau_out=2.5)
# decharged tests at the module's box (2048^2, dx = 0.25, L = 512)
job('A_k2_t150_N2048_s301_to1.6', kind='A', kappa=2.0, tau_i=150.0, N=2048, dx=0.25, seed=301, a_end=1.6, dtau_out=2.5,
    save_end=True)
job('CC_k2_t150_N2048_s301_from1.6', kind='CONT', src='A_k2_t150_N2048_s301_to1.6', decharge=False, a_end=2.5, dtau_out=2.5)
job('DC_k2_t150_N2048_s301_from1.6', kind='CONT', src='A_k2_t150_N2048_s301_to1.6', decharge=True, a_end=2.5, dtau_out=2.5)
job('A_k2_t40_N2048_s302_to1.8', kind='A', kappa=2.0, tau_i=40.0, N=2048, dx=0.25, seed=302, a_end=1.8, dtau_out=1.0,
    save_end=True)
job('CC_k2_t40_N2048_s302_from1.8', kind='CONT', src='A_k2_t40_N2048_s302_to1.8', decharge=False, a_end=3.0, dtau_out=1.0)
job('DC_k2_t40_N2048_s302_from1.8', kind='CONT', src='A_k2_t40_N2048_s302_to1.8', decharge=True, a_end=3.0, dtau_out=1.0)
# dt robustness (statistics) at 512^2
job('A_k2_t40_N512_s501_dt05', kind='A', kappa=2.0, tau_i=40.0, N=512, dx=0.25, seed=501, a_end=3.0, dtau_out=1.0)
job('A_k2_t40_N512_s501_dt025', kind='A', kappa=2.0, tau_i=40.0, N=512, dx=0.25, seed=501, a_end=3.0, dtau_out=1.0, dt=0.025)
# small 3D charged exit and neutral (caveat test of the 3D claim)
job('F_k2_t150_N2048_s301_freeze1.6', kind='CONT', src='A_k2_t150_N2048_s301_to1.6', decharge=False, freeze=True,
    tau_end=370.0, dtau_out=5.0)
job('F_dc_k2_t150_N2048_s301_freeze1.6', kind='CONT', src='A_k2_t150_N2048_s301_to1.6', decharge=True, freeze=True,
    tau_end=370.0, dtau_out=5.0)
job('A_k2_t150_N2048_s303_to1.6', kind='A', kappa=2.0, tau_i=150.0, N=2048, dx=0.25, seed=303, a_end=1.6, dtau_out=2.5,
    save_end=True)
job('F_k2_t150_N2048_s303_freeze1.6', kind='CONT', src='A_k2_t150_N2048_s303_to1.6', decharge=False, freeze=True,
    tau_end=370.0, dtau_out=5.0)
job('F_dc_k2_t150_N2048_s303_freeze1.6', kind='CONT', src='A_k2_t150_N2048_s303_to1.6', decharge=True, freeze=True,
    tau_end=370.0, dtau_out=5.0)
job('D3_k2_t40_N128_s601', kind='A', kappa=2.0, tau_i=40.0, N=128, dx=0.3, seed=601, a_end=2.2, dtau_out=2.0, ndim=3, dt=0.06)
job('D3_k2_t40_N192_s602', kind='A', kappa=2.0, tau_i=40.0, N=192, dx=0.3, seed=602, a_end=2.2, dtau_out=2.0, ndim=3, dt=0.06)


def run(name):
    cfg = dict(JOBS[name])
    t0 = time.time()
    dt = cfg.get('dt', DT)
    meta = {}
    if cfg['kind'] == 'CONT':
        lat = Lattice.load(os.path.join(SCR, cfg['src'] + '.npz'))
        if cfg['decharge']:
            lat.with_X = False
            for f in (lat.q1, lat.q2, lat.p1, lat.p2):
                f[...] = 0.0
        tau_i = lat.tau_i
        if cfg.get('freeze'):
            lat.a_freeze = lat.tau / lat.tau_i
    else:
        nd = cfg.get('ndim', 2)
        lat = Lattice(cfg['N'], cfg['dx'], cfg['tau_i'], ndim=nd, with_X=(cfg['kind'] != 'B'))
        tau_i = cfg['tau_i']
        if cfg['kind'] == 'A':
            meta = vlat.init_charged_exit(lat, cfg['kappa'], cfg['seed'])
        elif cfg['kind'] == 'B':
            vlat.init_neutral(lat, cfg['seed'])
        elif cfg['kind'] == 'C1':
            meta = vlat.init_subcritical_plus(lat, cfg['frac'], cfg['seed'])
    Q0 = lat.charge()
    tau_end = cfg['tau_end'] if 'tau_end' in cfg else tau_i * cfg['a_end']
    taus = np.arange(lat.tau, tau_end + 1e-9, cfg['dtau_out'])
    if taus[-1] < tau_end - 1e-9:
        taus = np.append(taus, tau_end)
    ts = []
    qmax = 0.0
    for tau in taus:
        lat.evolve(tau, dt)
        d = vlat.diagnostics(lat)
        if Q0 != 0:
            d['Q_rel_dev'] = lat.charge() / Q0 - 1
            qmax = max(qmax, abs(d['Q_rel_dev']))
        ts.append(d)
        print('%s tau=%.2f a=%.3f f+=%.3f f-=%.3f Lb=%.4g Lraw=%.4g wall=%.0fs' % (
            name, lat.tau, lat.a, d['f_plus'], d['f_minus'], d['L_basin'], d['L_raw'], time.time() - t0), flush=True)
    if cfg.get('save_end'):
        lat.save(os.path.join(SCR, name + '.npz'))
    out = dict(name=name, config=cfg, dt=dt, meta=meta, N=lat.N, dx=lat.dx, tau_i=tau_i, ndim=lat.ndim,
               L_box=lat.N * lat.dx, Q0=Q0, max_Q_rel_dev=qmax, wall_seconds=time.time() - t0, series=ts)
    with open(os.path.join(HERE, 'results', name + '.json'), 'w') as fh:
        json.dump(out, fh, indent=0)


if __name__ == '__main__':
    for nm in sys.argv[1:]:
        run(nm)
