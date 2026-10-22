"""v_scenario reproduction driver: downloads the public inputs to the scratchpad (never into the repository), then
runs v1..v7. About 6 minutes on 4 cores (v2 dominates).  Run from anywhere:  python3 -I code/run_all.py
"""
import os
import subprocess
import sys
import tarfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DL = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/qtheory/v_scenario_dl'
URLS = {
    'desi_mean.txt': 'https://raw.githubusercontent.com/CobayaSampler/bao_data/master/desi_bao_dr2/'
                     'desi_gaussian_bao_ALL_GCcomb_mean.txt',
    'desi_cov.txt': 'https://raw.githubusercontent.com/CobayaSampler/bao_data/master/desi_bao_dr2/'
                    'desi_gaussian_bao_ALL_GCcomb_cov.txt',
    'pplus.dat': 'https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/'
                 '4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat',
    'pplus_cov.cov': 'https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/'
                     '4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES_STAT%2BSYS.cov',
    'ss/src.tar': 'https://arxiv.org/src/1803.01038',
    'ss/eos_w.dat': 'https://member.ipmu.jp/satoshi.shirai/eos_w.dat',
}


def fetch():
    for rel, url in URLS.items():
        p = os.path.join(DL, rel)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        if not os.path.exists(p):
            urllib.request.urlretrieve(url, p)
    x = os.path.join(DL, 'ss', 'x')
    if not os.path.exists(os.path.join(x, 'anc', 'standardmodel2018.dat')):
        os.makedirs(x, exist_ok=True)
        with tarfile.open(os.path.join(DL, 'ss', 'src.tar')) as t:
            t.extractall(x, filter='data')


if __name__ == '__main__':
    fetch()
    for s in ('v1_symbolic.py', 'v3_realistic.py', 'v4_light_q.py', 'v5_early.py', 'v6_data.py',
              'v2_toy_full.py', 'v7_offsets.py'):
        print('==', s, flush=True)
        subprocess.run([sys.executable, '-I', os.path.join(HERE, s)], check=True)
