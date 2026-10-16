"""Frozen dependencies; all newly written deliverables stay in this packet."""
from pathlib import Path
import sys
P=Path(__file__).resolve().parent;PREVIOUS=P.parent/'cp12_source_resolved_propagation_20261015'
sys.path.insert(0,str(PREVIOUS))
from shared import *
from source_basis import extend,save_cache
from response_run import response,ratio
from neutral_krylov import reduced
P=Path(__file__).resolve().parent


def neutral_path(n,q,b,t):return cache_root()/f'HOST_N{n}_Q{q}_B{b}_W{t:g}.npz'
