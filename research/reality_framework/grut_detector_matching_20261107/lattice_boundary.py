"""M0-02 wrapper: analytic ghost-endpoint hypothesis, unchanged lattice solver."""
import argparse
import copy
import hashlib
import json
from pathlib import Path

import numpy as np

import lattice

HERE=Path(__file__).resolve().parent


def run(card,n,smoothed=False,refined=False):
    geometry=copy.deepcopy(card)
    a=card["length"]/(n+1)
    geometry["length"]=card["length"]*(n-1)/(n+1)
    relative=(np.arange(n-1)-(n-2)/2)*a
    weights=np.exp(-0.5*(relative/card["sigma"])**2)
    weights/=weights.sum()
    delta=float(np.dot(weights,(-1.)**np.arange(n-1))/4)
    if smoothed:
        geometry["phase"]=card["phase"]-card["coupling"]*delta
    result=lattice.run(geometry,n,refined=refined)
    result.update(target_length=card["length"],occupied_site_span=geometry["length"],
                  ghost_endpoint_span=(n+1)*a,detector_smoothed=smoothed,
                  alternating_detector_constant=delta,
                  continuum_observable_mean=result["mean"]-(delta if smoothed else 0))
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--N",type=int,required=True)
    parser.add_argument("--smoothed",action="store_true")
    parser.add_argument("--refined",action="store_true")
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    path=HERE/"protocol_m0_02.json"
    card=json.loads(path.read_text())
    result=run(card,args.N,args.smoothed,args.refined)
    result["protocol_sha256"]=hashlib.sha256(path.read_bytes()).hexdigest()
    with args.output.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")
    print(json.dumps({k:result[k] for k in ["N","mean","gap","detector_smoothed","seconds"]}))
