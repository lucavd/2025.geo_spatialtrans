"""tools/S2.6_run.py - S2.6: esecuzione del pannello (domande 1, 2, 3) con limiti di processi e di memoria. Definizioni: docs/steps/S2.6.md.
Uso (root del repo): .venv/bin/python tools/S2.6_run.py [plan|run] [q1,q2,q3]
Ogni compito scrive un file di uscita stabile in /mnt/micron/geo_spatialtrans/S2.6/; un compito fatto si salta (rilanciabile). Processi
attivi <= contenuto del file WORKERS (guardia della temperatura); memoria prevista dei compiti attivi <= MEM_GB (180). I compiti grandi in
un pezzo solo delle sezioni girano uno alla volta; un controllo ferma un processo R oltre 200 GB di memoria residente."""
import os, sys, json, time, subprocess, random
M = "/mnt/micron/geo_spatialtrans/S2.6"; WK = f"{M}/WORKERS"; MEM_GB = 180; KILL_GB = 200
RS = "Rscript --vanilla"
ENV = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
A6 = ["A1", "A2", "A3", "A4", "A5", "A6"]; SIDES = [1000, 2000, 3000]
def draws(seed, n):  # stessi numeri di R: RNGkind Mersenne-Twister; set.seed(seed); sample.int(1e8, n)
    out = subprocess.run(["Rscript", "--vanilla", "-e", f'RNGkind("Mersenne-Twister","Inversion","Rejection"); set.seed({seed}); cat(sample.int(1e8,{n}))'],
                         capture_output=True, text=True, check=True).stdout.split()
    return [int(x) for x in out]
# memoria prevista (GB) dalle esplorazioni E-2, E-5: picco a 1 mm2 per tessuto, parte fissa 2.15 GB, crescita lineare con l'area
P1 = {"A1": 5.84, "A2": 8.85, "A3": 5.31, "A4": 3.16, "A5": 7.09, "A6": 3.52}
GEN = {"A1": 38.9, "A2": 69.3, "A3": 35.0, "A4": 13.1, "A5": 56.0, "A6": 2.9}   # milioni di trascritti generati per mm2 (E-5; stroma stimato)
TISSUE = {"A1": 22.52, "A2": 33.08, "A4": 22.47, "A5": 25.58, "A6": 42.63}
def mem_one(A, mm2, stripes=True): s = P1[A] - 2.15; return (2.15 + s * mm2) * (1.0 if stripes else 0.5) + 1
def tasks(which):
    T = []
    sL = draws(20261209, 36); sO = draws(20261210, 36); sD = draws(20261212, 24); sQ2 = draws(20261211, 72); sS = draws(20261213, 5); sS2 = draws(20261214, 10)
    lay = lambda A, side, k: f"{M}/layers/pc_{A}_{side}um_px2_s{sL[(A6.index(A) * 3 + SIDES.index(side)) * 2 + k]}.rds"
    for A in A6:
        for side in SIDES:
            for k in (0, 1):
                i = (A6.index(A) * 3 + SIDES.index(side)) * 2 + k; L = lay(A, side, k)
                T.append(dict(id=f"L_{A}_{side}_{k}", q="q1", cmd=f"{RS} tools/S2.6_explore.R pc {A} {side} 2 {sL[i]} {M}/layers", out=L, mem=3, cpu=1, dep=[]))
    for A in A6:
        for side in SIDES:
            for k in (0, 1):
                i = (A6.index(A) * 3 + SIDES.index(side)) * 2 + k; mm2 = (side / 1000) ** 2
                T.append(dict(id=f"O1_{A}_{side}_{k}", q="q1", cmd=f"{RS} tools/S2.6_pieces.R one {lay(A, side, k)} {sO[i]} {M}/q1 one_{A}_{side}_k{k} 1 1",
                              out=f"{M}/q1/one_{A}_{side}_k{k}/run.json", mem=mem_one(A, mm2), cpu=1, dep=[f"L_{A}_{side}_{k}"]))
        # quattro pezzi interni del disegno (504 um, bordo 280 um) del quadrato di 9 mm2, k = 0: livelli dei pezzi e simulazione in processi separati
        T.append(dict(id=f"DC_{A}", q="q1", cmd=f"{RS} tools/S2.6_pieces.R cut {lay(A, 3000, 0)} 504 280 15,16,21,22 {M}/q1/design_{A}",
                      out=f"{M}/q1/design_{A}/cut.json", mem=4, cpu=1, dep=[f"L_{A}_3000_0"]))
        for j, p in enumerate((15, 16, 21, 22)):
            T.append(dict(id=f"D_{A}_{p}", q="q1", cmd=f"{RS} tools/S2.6_pieces.R one {M}/q1/design_{A}/piece_p{p:03d}.rds {sD[A6.index(A) * 4 + j]} {M}/q1/design_{A} p{p:03d} 1 1",
                          out=f"{M}/q1/design_{A}/p{p:03d}/run.json", mem=mem_one(A, 1.132), cpu=1, dep=[f"DC_{A}"]))
    for A in A6:
        L = lay(A, 2000, 0)
        for d, (Tu, H) in enumerate([(None, None), (1000, 280), (504, 280), (1000, 0)]):
            arm = ["ONE", "P1000", "P504", "P1000H0"][d]
            for k in (0, 1, 2):
                s = sQ2[(A6.index(A) * 4 + d) * 3 + k]; tag = f"{arm}_{A}_k{k}"; out = f"{M}/q2/{tag}/run.json"
                if arm == "ONE":
                    T.append(dict(id=f"Q2_{tag}", q="q2", cmd=f"{RS} tools/S2.6_pieces.R one {L} {s} {M}/q2 {tag} 0 1", out=out, mem=mem_one(A, 4, False), cpu=1, dep=[f"L_{A}_2000_0"]))
                else:
                    W = 4
                    T.append(dict(id=f"Q2_{tag}", q="q2", cmd=f"{RS} tools/S2.6_pieces.R pieces {L} {Tu} {H} {s} {M}/q2 {tag} 0 {W}", out=out,
                                  mem=W * mem_one(A, ((Tu + 2 * H) / 1000) ** 2, False) + 2, cpu=W, dep=[f"L_{A}_2000_0"]))
                Ts = "1000 504" if arm == "ONE" else str(Tu)
                T.append(dict(id=f"M2_{tag}", q="q2", cmd=f"{RS} tools/S2.6_measures.R {M}/q2/{tag} {Ts}", out=f"{M}/q2/{tag}/measures_T{Ts.split()[-1]}.csv", mem=6, cpu=1, dep=[f"Q2_{tag}"]))
    SEC = ["A1", "A2", "A4", "A5", "A6"]
    for i, A in enumerate(SEC):
        Ls = f"{M}/sections/sec_{A}_s{sS[i]}.rds"
        T.append(dict(id=f"SL_{A}", q="q3", cmd=f"{RS} tools/S2.6_section_layer.R {A} {sS[i]} {M}/sections", out=Ls, mem=14, cpu=1, dep=[]))
        if A in ("A1", "A4", "A6"):
            gb = 0.092 * GEN[A] * TISSUE[A] + 0.87 * TISSUE[A] + 3
            T.append(dict(id=f"S1_{A}", q="q3", cmd=f"{RS} tools/S2.6_pieces.R one {Ls} {sS2[2 * i]} {M}/q3 one_{A} 1 0", out=f"{M}/q3/one_{A}/run.json",
                          mem=gb, cpu=1, dep=[f"SL_{A}"], exclusive_big=True))
        pp = mem_one(A, 2.434) + 4; W = max(2, min(16, int(120 // pp)))
        T.append(dict(id=f"SP_{A}", q="q3", cmd=f"{RS} tools/S2.6_pieces.R pieces {Ls} 1000 280 {sS2[2 * i + 1]} {M}/q3 pieces_{A} 1 {W}", out=f"{M}/q3/pieces_{A}/run.json",
                      mem=W * pp + 8, cpu=W, dep=[f"SL_{A}"]))
    return [t for t in T if t["q"] in which]
def workers():
    try: return max(1, min(28, int(open(WK).read().split()[0])))
    except Exception: return 16
def rss_guard(procs):
    for t, p in procs:
        try:
            out = subprocess.run(["ps", "-o", "rss=", "--ppid", str(p.pid)], capture_output=True, text=True).stdout.split()
            out += subprocess.run(["ps", "-o", "rss=", "-p", str(p.pid)], capture_output=True, text=True).stdout.split()
            gb = max([int(x) for x in out] + [0]) / 1048576
            if gb > KILL_GB: p.kill(); print(f"[guardia] {t['id']} fermato: {gb:.0f} GB", flush=True)
        except Exception: pass
def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "plan"; which = (sys.argv[2] if len(sys.argv) > 2 else "q1,q2,q3").split(",")
    os.makedirs(M, exist_ok=True)
    if not os.path.exists(WK): open(WK, "w").write("16\n")
    ALL = tasks(["q1", "q2", "q3"]); T = [t for t in ALL if t["q"] in which]
    done = {t["id"] for t in ALL if os.path.exists(t["out"])}   # le dipendenze possono stare in un'altra domanda (correzione dopo il primo lancio)
    print(f"{len(T)} compiti, {len(done)} gia' fatti; memoria prevista massima per compito {max(t['mem'] for t in T):.0f} GB", flush=True)
    if mode == "plan":
        for t in T: print(t["id"], f"{t['mem']:.1f} GB", t["cpu"], "fatto" if t["id"] in done else "", t["cmd"])
        return
    os.makedirs(f"{M}/logs", exist_ok=True)
    todo = [t for t in T if t["id"] not in done]; run = []; failed = set()
    while todo or run:
        for t, p in list(run):
            if p.poll() is not None:
                run.remove((t, p))
                if p.returncode == 0 and os.path.exists(t["out"]): done.add(t["id"]); print(time.strftime("%H:%M:%S"), "fatto", t["id"], flush=True)
                else: failed.add(t["id"]); print(time.strftime("%H:%M:%S"), "FALLITO", t["id"], p.returncode, flush=True)
        rss_guard(run)
        cpu = sum(t["cpu"] for t, _ in run); mem = sum(t["mem"] for t, _ in run); big = any(t.get("exclusive_big") for t, _ in run)
        for t in list(todo):
            if any(d in failed for d in t["dep"]): todo.remove(t); failed.add(t["id"]); print("SALTATO", t["id"], "(dipendenza fallita)", flush=True); continue
            if not all(d in done for d in t["dep"]): continue
            if t.get("exclusive_big") and big: continue
            if cpu + t["cpu"] > workers() and run: continue
            if mem + t["mem"] > MEM_GB and run: continue
            lg = open(f"{M}/logs/{t['id']}.log", "a")
            p = subprocess.Popen(t["cmd"], shell=True, stdout=lg, stderr=subprocess.STDOUT, env=ENV)
            run.append((t, p)); todo.remove(t); cpu += t["cpu"]; mem += t["mem"]; big = big or t.get("exclusive_big", False)
            print(time.strftime("%H:%M:%S"), "avvio", t["id"], f"({t['mem']:.0f} GB, {t['cpu']} proc; attivi {len(run)}, {cpu} proc, {mem:.0f} GB)", flush=True)
        time.sleep(10)
    print("FINE", len(done), "fatti,", len(failed), "falliti:", sorted(failed), flush=True)
    sys.exit(1 if failed else 0)
if __name__ == "__main__": main()
