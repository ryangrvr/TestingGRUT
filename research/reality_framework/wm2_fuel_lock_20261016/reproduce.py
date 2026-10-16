"""Check saved evidence, or regenerate the complete bounded calculation."""
import argparse
from pathlib import Path
import subprocess
import sys

P=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--recompute',action='store_true')
    args=parser.parse_args()
    if args.recompute:
        for name in ['controls.py','report.py','plot.py']:
            subprocess.run([sys.executable,str(P/name)],check=True)
    subprocess.run([sys.executable,str(P/'verify.py')]+(['--refresh'] if args.recompute else []),check=True)


if __name__=='__main__':main()
