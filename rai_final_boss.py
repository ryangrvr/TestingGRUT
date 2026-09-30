#!/usr/bin/env python3
"""RAI-FINAL-BOSS -- gates the committed record. Read-only. Banks nothing. W-0.

Journal-read AND de-pinned (the twice-caught lesson): no verdict value is a pass condition;
gates assert well-formedness and record/journal consistency via dynamic needles; the blocked
stage must be DISCLOSED, not hidden; free-text quality is not mechanically certifiable.
"""
import json, hashlib, subprocess, sys, os, re

ROOT = os.path.dirname(os.path.abspath(__file__))
# Journal path resolution order: $RAI_JOURNAL env var > first CLI arg > repo-relative default.
def _default_journal_path():
    for wf in sorted(os.listdir(os.path.join(ROOT, "subagents", "workflows"))) \
            if os.path.isdir(os.path.join(ROOT, "subagents", "workflows")) else []:
        p = os.path.join(ROOT, "subagents", "workflows", wf, "journal.jsonl")
        if os.path.exists(p): return p
    return None

def _resolve_journal_path():
    if os.environ.get("RAI_JOURNAL"): return os.environ["RAI_JOURNAL"]
    if len(sys.argv) > 1 and not sys.argv[1].startswith("-"): return sys.argv[1]
    return _default_journal_path()

JPATH = _resolve_journal_path()

def read(p):
    with open(os.path.join(ROOT, p), encoding="utf-8") as f: return f.read()

FAILS, N = [], 0
def gate(cond, label, kind="JOURNAL"):
    global N; N += 1
    ok = bool(cond)
    print(f"  [{'PASS' if ok else 'FAIL'}] {kind:9s} {label}")
    if not ok: FAILS.append(label)
    return ok

print("\n== A. GOVERNANCE ==")
reg = json.loads(read("provenance/claims.json"))
claims = reg["claims"] if isinstance(reg, dict) else reg
gate(len(claims) == 74, "register 74 nodes", "GRAPH")
gate(hashlib.sha256(read("provenance/claims.json").encode()).hexdigest()
     .startswith("beaeb84e8a6f8468"), "register sha256 unchanged", "GRAPH")
_rh = subprocess.run(["git","-C",ROOT,"rev-parse","HEAD"],capture_output=True,text=True)
_rv = subprocess.run(["git","-C",ROOT,"rev-parse","origin/v4"],capture_output=True,text=True)
if _rh.returncode != 0 or _rv.returncode != 0:
    gate(False, f"HEAD == origin/v4 by ref identity "
                f"(git error: {(_rh.stderr or _rv.stderr).strip()})", "GRAPH")
else:
    gate(_rh.stdout.strip() == _rv.stdout.strip(), "HEAD == origin/v4 by ref identity", "GRAPH")
_rst = subprocess.run(["git","-C",ROOT,"status","--porcelain"],capture_output=True,text=True)
if _rst.returncode != 0:
    gate(False, f"no prior result modified (git status error: {_rst.stderr.strip()})", "GRAPH")
else:
    porc = _rst.stdout
    gate(not any(l[:2] in (" M","M ","MM","D ") for l in porc.splitlines()),
         "no prior result modified", "GRAPH")
BY = {c["id"]: c for c in claims}
for nid in ("u3_split_origin","u4_constitutive_origin"):
    gate(BY[nid].get("tier") == "to-derive", f"{nid} fence unmoved", "GRAPH")

print("\n== B. JOURNAL (well-formedness; no verdict value is a pass condition) ==")
res, failed = [], []
try:
    with open(JPATH, encoding="utf-8") as f:
        for line in f:
            if line.strip():
                d = json.loads(line)
                if d.get("type") == "result" and isinstance(d.get("result"), dict): res.append(d["result"])
                if d.get("type") == "failed": failed.append(d.get("key",""))
except (OSError, json.JSONDecodeError) as e:
    gate(False, f"journal readable ({e})")
    print(f"\nBATTERY: {N-len(FAILS)}/{N}  FAILURES: {FAILS}")
    sys.exit(1)
def last(pred):
    xs=[r for r in res if pred(r)]; return xs[-1] if xs else None
recon_c = last(lambda r: "deletion_table" in r)
recon_w = last(lambda r: "arrow_test" in r)
att_c   = last(lambda r: "discharge_class" in r and "verification" in r)
syn     = last(lambda r: "clpw_final" in r)
fin     = last(lambda r: "o1_clpw" in r)
sents   = [r for r in res if "sentence" in r and "grades" in r]
hosts   = [r for r in res if "hidden_input" in r and "attack" in r]
resd    = last(lambda r: "residue_verdict" in r)
gate(recon_c is not None and recon_w is not None, "both reconstructions present")
gate(att_c is not None, "clpw attack present")
gate(len(failed) >= 3, f"the wiesbrock:attack refusals are ON THE RECORD ({len(failed)} failed events)")
gate(len(hosts) >= 2, f"double hostile ran ({len(hosts)} hostiles)")
gate(syn is not None and syn.get("clpw_final") and syn.get("wiesbrock_final") and syn.get("final_state"),
     f"synthesis well-formed (journal: {syn.get('clpw_final')} / {syn.get('wiesbrock_final')} / {syn.get('final_state')})")
gate(len(sents) == 2, "blind sentence pair present")
gate(fin is not None and all(fin.get(k) for k in
     ("o1_clpw","o8_1space","o14_empirical_discriminator","o15_mathematical_discriminator",
      "o17_primary_sentence","o18_hostile_sentence","o19_sentence_comparison","final_classification")),
     f"final output complete (classification: {fin.get('final_classification') if fin else None})")
gate(resd is not None and resd.get("residue_verdict") in ("PHYSICAL","GAUGE","INPUT","MIXED","UNDETERMINED"),
     f"residue verdict well-formed (journal: {resd.get('residue_verdict') if resd else None})")

print("\n== C. MAIN-LOOP DIFF (the blocked stage's replacement, computed here) ==")
rw = json.dumps(recon_w) if recon_w else ""
hb = " ".join(json.dumps(x) for x in hosts)
gate("Remark (19)" in rw or "(19)" in rw, "recon: both hsm orientations located in sources", "DIFF")
gate("inclusion order alone" in rw or "orientation-free" in rw,
     "recon: positivity from inclusion order alone", "DIFF")
gate("epsilon" in hb or "sign" in hb.lower(),
     "hostile B independently names the signed nesting datum", "DIFF")
gate(("hsm" in rw.lower() or "half-sided" in rw.lower()) and ("hsm" in hb.lower() or "half-sided" in hb.lower()),
     "both independent passes attribute the orientation to the half-sided clause", "DIFF")

print("\n== D. RECORD/JOURNAL CONSISTENCY (dynamic needles) ==")
rec = "RAI_FINAL_BOSS.md"
if not (syn is not None and fin is not None and resd is not None):
    gate(False, "syn/fin/resd records present before consistency check", "GRAPH")
    print("  [ -- ] consistency loop skipped: journal missing a required record")
elif os.path.exists(os.path.join(ROOT, rec)):
    MD = re.sub(r"\s+"," ", read(rec).replace(">"," ").replace("*",""))
    gate("[[" not in MD, "no template token", "GRAPH")
    for val, lbl in ((syn.get("clpw_final"), "record states the journal's CLPW class"),
                     (syn.get("wiesbrock_final"), "record states the journal's Wiesbrock class"),
                     (fin.get("final_classification"), "record states the final classification"),
                     (resd.get("residue_verdict"), "record states the residue verdict")):
        gate(val is not None and (val in MD or val.replace("-", " ") in MD),
             f"{lbl} ({val})", "GRAPH")
    gate("classifier" in MD.lower() and ("refused" in MD.lower() or "blocked" in MD.lower()),
         "record DISCLOSES the blocked stage", "GRAPH")
    gate("own derivation" in MD.lower() or "campaign's own derivations" in MD.lower()
         or "OWN-DERIVATION" in MD,
         "record flags the campaign-own derivations as unverified-by-literature", "GRAPH")
else:
    print("  [ -- ] record not yet written")

print(f"\nBATTERY: {N-len(FAILS)}/{N}" + (f"  FAILURES: {FAILS}" if FAILS else ""))
if syn is not None and fin is not None:
    print(f"FROM JOURNAL: CLPW={syn.get('clpw_final')} WIES={syn.get('wiesbrock_final')} "
          f"STATE={syn.get('final_state')} FINAL={fin.get('final_classification')}")
print("W-0 -- reported, NOT banked. Register unmodified. Fences unmoved.")
sys.exit(1 if FAILS else 0)
