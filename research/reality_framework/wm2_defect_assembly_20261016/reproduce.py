"""Check saved evidence or regenerate exact configuration analysis first."""
import argparse
from pathlib import Path
import subprocess
import sys
P=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--full',action='store_true');args=ap.parse_args()
    if args.full:
        for script in ['controls.py','report.py','plot.py']:
            subprocess.run([sys.executable,str(P/script)],check=True)
    subprocess.run([sys.executable,str(P/'verify.py')]+(['--refresh'] if args.full else []),check=True)


if __name__=='__main__':main()
