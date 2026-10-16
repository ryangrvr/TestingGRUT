"""Verify saved evidence, or regenerate the complete bounded packet."""
from pathlib import Path
import argparse
import subprocess
import sys
P=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--full',action='store_true'); args=ap.parse_args()
    scripts=['controls.py','report.py','plot.py'] if args.full else []
    for script in scripts:subprocess.run([sys.executable,str(P/script)],check=True)
    subprocess.run([sys.executable,str(P/'verify.py')]+(['--refresh'] if args.full else []),check=True)


if __name__=='__main__':main()
