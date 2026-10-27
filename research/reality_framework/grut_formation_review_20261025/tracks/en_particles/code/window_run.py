"""en_particles part 2 driver: dense conformal-time windows of the regenerated N48 history.

Uses Codex's UNMODIFIED early3d.cpp engine (reproduction-type use) from a COPY of the packet. Trick: the engine's
clock only labels output times, the Hamiltonian is autonomous in conformal time s. Restarting from a saved 64-bit
state with clock (tau0=s_start, a0=1, C=0) makes the engine's 'tau' equal to s, so samples/snapshots are uniform in
conformal time. Each chunk's snapshots are reduced immediately to Fourier data (own code) and deleted.

Usage: python3 -B window_run.py <packet_code_dir> <start_state.bin> <out_dir> <duration_s> [ds=0.05] [chunk=2.0]
       [threads=2] [hmax=0.004] [kick_eps=0] [kick_seed=7]
kick_eps!=0: perturbs the start state's charged momentum by a random-phase complex field of rms kick_eps
(linear-response twin run); the same kick is saved for the response analysis.
"""
import sys, os, struct, subprocess, json, glob, time
from pathlib import Path
import numpy as np


def read_state(path):
    with open(path, 'rb') as f:
        N, L, s, bits = struct.unpack('<iddi', f.read(24))
    dt = '<f8' if bits == 64 else '<f4'
    a = np.fromfile(path, dtype=dt, offset=24).reshape(6, N, N, N).astype(np.float64)
    return N, L, s, a


PREC = np.complex64


def reduce_snapshot(a):
    """Own-code Fourier reduction: full FFT of psi=chi1+i chi2 and pi_psi; rfft of phi and pi_phi (forward norm)."""
    psi = np.fft.fftn(a[0]+1j*a[1], norm='forward'); ppsi = np.fft.fftn(a[3]+1j*a[4], norm='forward')
    phi = np.fft.rfftn(a[2], norm='forward'); pphi = np.fft.rfftn(a[5], norm='forward')
    return psi.astype(PREC), ppsi.astype(PREC), phi.astype(np.complex64), pphi.astype(np.complex64), \
        complex(psi[0, 0, 0]), complex(ppsi[0, 0, 0]), float(phi[0, 0, 0].real), float(pphi[0, 0, 0].real)


def main():
    code = Path(sys.argv[1]); start = Path(sys.argv[2]); out = Path(sys.argv[3]); dur = float(sys.argv[4])
    ds = float(sys.argv[5]) if len(sys.argv) > 5 else .05
    chunk = float(sys.argv[6]) if len(sys.argv) > 6 else 2.
    threads = int(sys.argv[7]) if len(sys.argv) > 7 else 2
    hmax = float(sys.argv[8]) if len(sys.argv) > 8 else .004
    eps = float(sys.argv[9]) if len(sys.argv) > 9 else 0.
    kseed = int(sys.argv[10]) if len(sys.argv) > 10 else 7
    sys.path.insert(0, str(code)); import en_common as E
    global PREC
    if eps != 0: PREC = np.complex128  # response twin runs need the small differences
    exe = E.build(); out.mkdir(parents=True, exist_ok=True); work = out/'work'; work.mkdir(exist_ok=True)
    N, L, s0, a = read_state(start)
    if eps != 0:
        rng = np.random.default_rng(kseed)
        kick = (rng.normal(size=(N, N, N))+1j*rng.normal(size=(N, N, N)))/np.sqrt(2)*eps
        kick -= kick.mean()
        a[3] += kick.real; a[4] += kick.imag
        np.save(out/'kick.npy', kick.astype(np.complex128))
    nsteps = int(round(dur/ds)); nchunk = int(round(chunk/ds))
    store = dict(psi=[], ppsi=[], phi=[], pphi=[]); zm = []; s_list = []; diag = []
    cur = a; s_cur = s0; done = 0; t0 = time.monotonic()
    # first sample: the start state itself
    r = reduce_snapshot(cur)
    for key, v in zip(['psi', 'ppsi', 'phi', 'pphi'], r[:4]): store[key].append(v)
    zm.append(r[4:]); s_list.append(s_cur)
    while done < nsteps:
        m = min(nchunk, nsteps-done); D = m*ds
        init = work/'initial.bin'; E.write_initial(init, cur[:3], cur[3:], L)
        for f in glob.glob(str(work/'state_s*.bin')): os.unlink(f)
        args = [str(exe), str(init), str(work/'state'), str(hmax), str(D), str(ds), str(ds), '4', str(threads), '64',
                str(s_cur), '1', '0']
        subprocess.run(args, check=True, stdout=subprocess.DEVNULL)
        tr = np.genfromtxt(str(work/'state')+'.tsv', names=True, delimiter='\t')
        diag.append(np.array([tr['eta']+s_cur, tr['rho'], tr['nQ'], tr['X2'], tr['phi2']]).T[1:])
        files = sorted(glob.glob(str(work/'state_s*.bin')))
        files = files[1:]  # index 0 is the start state (already stored)
        seen = set()
        for f in files:
            Nn, Ll, ss, arr = read_state(f)
            key = round(ss, 9)
            if key in seen: continue  # final snapshot may duplicate the last sample
            seen.add(key)
            r = reduce_snapshot(arr)
            for kk, v in zip(['psi', 'ppsi', 'phi', 'pphi'], r[:4]): store[kk].append(v)
            zm.append(r[4:]); s_list.append(s_cur+ss); cur = arr
        s_cur = s_list[-1]; done += m
        print(f'chunk to s={s_cur:.4f} samples={len(s_list)} wall={time.monotonic()-t0:.1f}s', flush=True)
        for f in glob.glob(str(work/'state_s*.bin')): os.unlink(f)
    s_arr = np.array(s_list)
    assert np.allclose(np.diff(s_arr), ds, atol=1e-7), np.diff(s_arr)[:5]
    np.savez(out/'window.npz', s=s_arr, psi=np.array(store['psi']), ppsi=np.array(store['ppsi']),
             phi=np.array(store['phi']), pphi=np.array(store['pphi']),
             zero=np.array([[z[0], z[1], z[2], z[3]] for z in zm], dtype=complex), N=N, L=L,
             diag=np.concatenate(diag), start=str(start), eps=eps)
    np.save(out/'final_state.npy', cur)
    (out/'meta.json').write_text(json.dumps(dict(start=str(start), s0=s0, duration=dur, ds=ds, hmax=hmax, eps=eps,
        threads=threads, samples=len(s_list), wall=time.monotonic()-t0), indent=1))


if __name__ == '__main__':
    main()
