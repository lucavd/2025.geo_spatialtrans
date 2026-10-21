# tools/S1.5_review_mutants_analyse.py - lettura dei log dei mutanti del revisore rieseguiti (copia di review_code/analyse_mutants.py)

import re, glob, os, json
RC="/mnt/micron/geo_spatialtrans/S1.5/postreview"
def parse(f):
    L=open(f).read().splitlines(); P=set(); F=set()
    for l in L:
        m=re.match(r'^(PASS|FAIL) (\S+)', l)
        if m: (P if m.group(1)=="PASS" else F).add(m.group(2))
    res=[l for l in L if l.startswith("RESULT")]; ex=[l for l in L if l.startswith("exit")]; tm=[l for l in L if re.match(r'^[0-9.]+ s [0-9]+ kB',l)]
    halted = any("Execution halted" in l for l in L)
    return P,F,(res[-1] if res else "no RESULT"),(ex[-1] if ex else ""),(tm[-1] if tm else ""),halted
bP,bF,br,be,bt,bh=parse(RC+"/logs/base.log")
print("BASE", br, be, bt, "FAIL:", sorted(bF), "n assert", len(bP)+len(bF))
desc=json.load(open(RC+"/mutants.json"))
out=[]
for t in sorted(os.path.basename(x)[:-4] for x in glob.glob(RC+"/logs/*.log")):
    if t=="base": continue
    P,F,r,e,tm,h=parse(RC+f"/logs/{t}.log")
    new=sorted((F - bF) & bP)
    missing=sorted(bP - P - F)
    out.append(dict(id=t, desc=desc.get(t, "mutante dell'autore"), detected=bool(new), by=new, result=r, halted=h, n_not_run=len(missing), time=tm))
    print(t, "DET" if new else "---", r, "halted" if h else "", "not_run=%d"%len(missing), new[:8])
import csv
w=csv.DictWriter(open("results/S1.5/S1.5_review_mutants_after.csv","w",newline=""), fieldnames=["id","desc","detected","by","result","halted","n_not_run"]); w.writeheader()
for o in out: w.writerow({k:(" ".join(o[k]) if k=="by" else o[k]) for k in w.fieldnames})
json.dump(out, open(RC+"/mutants_result.json","w"), ensure_ascii=False, indent=1)
