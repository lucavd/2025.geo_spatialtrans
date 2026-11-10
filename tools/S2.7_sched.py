#!/usr/bin/env python3
"""tools/S2.7_sched.py - S2.7: esegue una lista di compiti con un tetto di processi e una guardia della temperatura (come S2.6).
Uso (root del repo): .venv/bin/python tools/S2.7_sched.py <tasks.json> <log> [slots] [memoria GB]
tasks.json: lista di {"name", "cmd", "cost", "mem"} (cost = processi del compito, mem = memoria prevista in GB). Il tetto si rilegge da /mnt/micron/geo_spatialtrans/S2.7/SLOTS a
ogni giro se esiste; con Tctl > 88 C non parte nessun compito nuovo. Un compito il cui file "done" (campo facoltativo) esiste si salta."""
import json, os, re, subprocess, sys, time
T = json.load(open(sys.argv[1])); LOG = sys.argv[2]; SL = int(sys.argv[3]) if len(sys.argv) > 3 else 16; MEM = float(sys.argv[4]) if len(sys.argv) > 4 else 150.0
SF = "/mnt/micron/geo_spatialtrans/S2.7/SLOTS"
def log(m):
    with open(LOG, "a") as f: f.write(time.strftime("%Y-%m-%d %H:%M:%S ") + m + "\n")
def tctl():
    try:
        o = subprocess.run(["sensors"], capture_output=True, text=True).stdout; m = re.search(r"Tctl:\s*\+([0-9.]+)", o); return float(m.group(1)) if m else 0.0
    except Exception: return 0.0
todo = [t for t in T if not (t.get("done") and os.path.exists(t["done"]))]
log(f"compiti {len(T)}, da fare {len(todo)}")
run = []; fails = 0
while todo or run:
    for p, t, f in list(run):
        if p.poll() is not None:
            run.remove((p, t, f)); f.close(); log(f"fine {t['name']} exit {p.returncode}"); fails += p.returncode != 0
    slots = SL
    if os.path.exists(SF):
        try: slots = int(open(SF).read().strip())
        except Exception: pass
    used = sum(t["cost"] for _, t, _ in run); mem = sum(t.get("mem", 4) for _, t, _ in run); temp = tctl()
    while todo and used + todo[0]["cost"] <= slots and mem + todo[0].get("mem", 4) <= MEM and temp <= 88.0:
        t = todo.pop(0); f = open(t.get("log", os.devnull), "a")
        p = subprocess.Popen(t["cmd"], shell=True, stdout=f, stderr=subprocess.STDOUT, executable="/bin/bash")
        run.append((p, t, f)); used += t["cost"]; mem += t.get("mem", 4); log(f"avvio {t['name']} (costo {t['cost']}, in uso {used}/{slots}, Tctl {temp:.1f})")
    time.sleep(20)
log(f"finito; compiti falliti {fails}")
sys.exit(1 if fails else 0)
