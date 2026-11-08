"""Check P1 preservation and arithmetic, not continuum truth."""
import argparse
import copy
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import tempfile

import numpy as np
from analyze import analyze

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate(out):
    out = Path(out)
    provenance = json.loads((out/"provenance.json").read_text())
    require(provenance["clean_source"], "Dirty production source")
    hashes = provenance["source_and_input_sha256"]
    for relative, expected in hashes.items():
        require(hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()==expected,
                "Changed frozen source/input: "+relative)
    recomputed = analyze(out)
    stored = json.loads((out/"analysis.json").read_text())
    require(recomputed==stored, "Stored analysis disagrees with frozen analysis")
    count, maxnorm, minprob = 0, 0., 1.
    for point in provenance["card"]["points"]:
        path = out/(point["label"]+".json")
        if not path.exists():
            continue
        data = json.loads(path.read_text())
        event = json.loads((out/("event_"+point["label"]+".json")).read_text())
        require(event["status"]=="complete", "Point without completed event")
        require(data["label"]==point["label"] and data["cutoff"]==point["cutoff"]
                and data["quadrature"]==point["quadrature"], "Wrong point dictionary")
        require(data["dimension"]<=provenance["card"]["dimension_cap"], "Dimension cap exceeded")
        checkpoint = np.load(path.with_suffix(".npz"), allow_pickle=False)
        states, energy, ground = [checkpoint[k] for k in ("basis","energies","ground")]
        require(states.shape==(data["dimension"],data["occupied_modes"]), "Wrong basis shape")
        require(len({tuple(row) for row in states})==len(states), "Duplicate occupation tuple")
        require(np.all(states>=0) and np.all(states==np.floor(states)), "Invalid occupation tuple")
        omega = np.sqrt(1/np.pi+(np.arange(1,states.shape[1]+1)*np.pi/5)**2)
        require(float(np.max(states@omega))<=point["cutoff"]+1e-10, "State above declared cutoff")
        require(np.all(np.isfinite(energy)) and np.all(np.diff(energy)>=-1e-12), "Bad spectrum")
        require(np.max(abs(energy[:6]-data["energies_low"]))<1e-13, "Checkpoint spectrum disagreement")
        require(abs(np.vdot(ground,ground).real-1)<1e-10, "Unnormalized preparation")
        require(abs(data["gap"]-(energy[1]-energy[0]))<1e-13, "Wrong gap")
        for key, p in data["histories"].items():
            n = len(key.split(","))
            expected = {"".join(s) for s in itertools.product("+-",repeat=n)}|{"overflow"}
            require(set(p)==expected, "Missing, conditioned or extra outcomes")
            values = np.asarray(list(p.values()))
            require(np.all(np.isfinite(values)) and values.min()>=-1e-10, "Invalid probability")
            error = abs(float(values.sum())-1)
            require(error<1e-9, "Non-normalized complete history")
            losses = data["overflow_by_step"][key]
            require(len(losses)==n and min(losses)>=-1e-10, "Bad absorbing overflow ledger")
            require(abs(sum(losses)-p["overflow"])<1e-12, "Overflow ledger does not close")
            maxnorm, minprob = max(maxnorm,error), min(minprob,float(values.min()))
            count += 1
    require(count==12 if stored["complete"] else count<12, "Wrong number of completed histories")
    return {"verification":"saved_source_checkpoint_probability_and_frozen_analysis_checks_pass",
            "source_and_input_hashes":len(hashes), "complete_history_distributions":count,
            "max_probability_norm_error":maxnorm,"minimum_probability":minprob,
            "physical_verdict":stored["verdict"],"scientific_gates":stored.get("gates"),
            "independent_physical_implementation":False,
            "scope":"Preservation/arithmetic checks; not a continuum certificate or interval enclosure"}


def negative_controls(out):
    answers = {}
    with tempfile.TemporaryDirectory(prefix="grut_p1_verify_") as directory:
        trial = Path(directory)/"raw"
        shutil.copytree(out, trial)
        original = (trial/"analysis.json").read_text()
        altered = json.loads(original)
        altered["verdict"] = "continuum_certified"
        (trial/"analysis.json").write_text(json.dumps(altered))
        try:
            validate(trial)
            raise AssertionError("Corrupt verdict was accepted")
        except ValueError:
            answers["corrupt_verdict_rejected"] = True
        (trial/"analysis.json").write_text(original)
        point = trial/"E12_q128.json"
        data = json.loads(point.read_text())
        key = next(iter(data["histories"]))
        data["histories"][key].pop("overflow")
        point.write_text(json.dumps(data))
        # Update derived arithmetic so rejection must also inspect normalization/outcomes.
        (trial/"analysis.json").write_text(json.dumps(analyze(trial)))
        try:
            validate(trial)
            raise AssertionError("Missing overflow was accepted")
        except ValueError:
            answers["missing_overflow_rejected"] = True
    return answers


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, default=HERE/"raw")
    parser.add_argument("--negative-controls", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = validate(args.raw)
    if args.negative_controls:
        result["negative_controls"] = negative_controls(args.raw)
    if args.output:
        with args.output.open("x") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")
    print(json.dumps(result,indent=2))
