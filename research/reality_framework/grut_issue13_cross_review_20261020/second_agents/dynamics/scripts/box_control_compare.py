"""Compare negative-seed pair in frozen box (L24,Lz48) vs enlarged box (L36,Lz60), coarse h=.75, plus published fine run.
Usage: python -I box_control_compare.py <small.csv> <big.csv> <published_fine.csv>"""
import sys, numpy as np
s, b, f = (np.genfromtxt(p, delimiter=',', names=True) for p in sys.argv[1:4])
Q0 = s['Q'][0]
def ret(a):
    return (a['c0_Q12'] + a['c1_Q12']) / a['Q'][0]
def gap(a, L, Lz):
    return np.min([np.min(np.stack([L-abs(a[f'c{j}_x']), L-abs(a[f'c{j}_y']), Lz-abs(a[f'c{j}_z'])], 1), 1)-12 for j in (0, 1)], axis=0)
m = s['t'] >= 120
mb = (b['t'] >= 120) & (b['t'] <= 160)
print('initial Q small/big/fine: %.9f %.9f %.9f ; initial E: %.6f %.6f' % (s['Q'][0], b['Q'][0], f['Q'][0], s['E'][0], b['E'][0]))
print('max |dQ|/Q0 big box run: %.2e, max |dE|/E0: %.2e' % (np.max(abs(b['Q']/b['Q'][0]-1)), np.max(abs(b['E']/b['E'][0]-1))))
print('late (120..160) min retention: small-box coarse %.5f%%  big-box coarse %.5f%%  published fine (small box) %.5f%%' % (100*ret(s)[m].min(), 100*ret(b)[mb].min(), 100*ret(f)[f['t'] >= 120].min()))
d = [np.max(abs(s[f'c{j}_Q12'][m]-b[f'c{j}_Q12'][mb]))/Q0 for j in (0, 1)]
print('max late |dQ12|/Q0 small vs big box (coarse): core0 %.5f%% core1 %.5f%%' % (100*d[0], 100*d[1]))
print('max late centre shift small vs big box: %.4f' % max(np.max(np.sqrt(sum((s[f'c{j}_{x}'][m]-b[f'c{j}_{x}'][mb])**2 for x in 'xyz'))) for j in (0, 1)))
print('small-box min R12 wall gap late: %.2f ; big-box min gap 120..240: %.2f' % (gap(s, 24, 48)[m].min(), gap(b, 36, 60)[b['t'] >= 120].min()))
print(' t   ret_small  ret_big   c0(x,y,z)_big           Q12_big(c0,c1)   sep_big')
for r in b:
    if r['t'] >= 100 and r['t'] % 10 == 0:
        k = np.where(s['t'] == r['t'])[0]
        rs = '%.5f' % (100*(s['c0_Q12'][k[0]]+s['c1_Q12'][k[0]])/Q0) if len(k) else '   -   '
        sep = np.sqrt(sum((r[f'c0_{x}']-r[f'c1_{x}'])**2 for x in 'xyz'))
        print('%3d  %9s  %8.5f  (%6.2f,%6.2f,%6.2f)  (%.2f,%.2f)  %.2f' % (r['t'], rs, 100*(r['c0_Q12']+r['c1_Q12'])/Q0, r['c0_x'], r['c0_y'], r['c0_z'], r['c0_Q12'], r['c1_Q12'], sep))
lt = b['t'] >= 160
print('big-box retention slope over 160..240: %.4f%% per 40 tau' % (100*40*np.polyfit(b['t'][lt], ret(b)[lt], 1)[0]))
