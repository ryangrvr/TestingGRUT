"""qm_micro item 1 table (Claude, 2 Nov 2026): s(T, theta) for X, and the Xbar channel fractions, on the planar M1
interface (scatter1d.solve, Numerov h = 0.01).  Bulk channels are labelled by the continuum dispersion branch."""
import json, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import scatter1d as sc
I = sc.Interface(h=0.01); Om = I.Om
Ts = (0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0); ths = (0, 30, 45, 60, 75, 85)
X, B = [], []
for T in Ts:
    E = 1 + T; p = math.sqrt(E * E - 1)
    row = dict(T=T, s={})
    for th in ths:
        r = sc.solve(I, E - Om, p * math.sin(math.radians(th)), 'X')
        row['s']['%d' % th] = 1 - r['reflected']
        if th == 0:
            row['channels_normal'] = dict(reflected_X=r['reflected'], into_bulk=r['trans_total'], phi=r['out'].get('out_phi', 0.0), Xbar=r['out'].get('out_Xbar', 0.0),
                                          bulk_modes=[(t['k'], t['frac']) for t in r['trans']], w=E - Om)
    X.append(row)
    for th in (0, 60):
        q = sc.solve(I, E + Om, p * math.sin(math.radians(th)), 'Xbar')
        B.append(dict(T=T, theta=th, w=E + Om, reflected_Xbar=q['reflected'], destroyed=1 - q['reflected'], converted_to_X=q['out'].get('out_X', 0.0),
                      phi_emission=q['out'].get('out_phi', 0.0), into_bulk=q['trans_total'], bulk_modes=[(t['k'], t['frac']) for t in q['trans']],
                      unitarity=q['unitarity']))
    print('T=%.2f X s(theta=%s) = %s | normal: R=%.5f bulk=%.5f' % (T, ths, ['%.4f' % row['s']['%d' % t] for t in ths], row['channels_normal']['reflected_X'], row['channels_normal']['into_bulk']))
for b in B:
    print('Xbar T=%.2f th=%d: destroyed %.5f (to X %.2e, phi %.2e, bulk %.5f: %s)' % (b['T'], b['theta'], b['destroyed'], b['converted_to_X'], b['phi_emission'], b['into_bulk'], [(round(k, 3), round(f, 4)) for k, f in b['bulk_modes']]))
json.dump(dict(X=X, Xbar=B, notes='bulk_modes: (k_x of the bulk channel, flux fraction); X at T <= 1 has only the phonon open inside'), open(sys.argv[1], 'w'), indent=1, default=float)
