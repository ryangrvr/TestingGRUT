"""Import the archived, unmodified B2 implementation."""
import sys,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
B2=HERE.parent/'grut_detector_reconstruction_20261108'
sys.path.insert(1,str(B2))
from model import reconstruction,distances,key,tv
from point import write,run as run_b2_point
from precheck import precheck

def card():return json.loads((HERE/'protocol.json').read_text())
