"""tools/S1.3d_arbiter_batch.py — S1.3d: arbitro esatto in parallelo su un elenco di chiavi (salta le verità già scritte).
Uso: .venv/bin/python tools/S1.3d_arbiter_batch.py <outdir> <file_chiavi> [--o2-keys file] [--cores 40] [--no-o1]"""
import sys, os, json, argparse, time, traceback
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib.util
spec = importlib.util.spec_from_file_location("arb", os.path.join(os.path.dirname(os.path.abspath(__file__)), "S1.3d_arbiter.py")); arb = importlib.util.module_from_spec(spec); spec.loader.exec_module(arb)

def one(args):
    outdir, key, do_o1, do_o2 = args
    geom, truth = os.path.join(outdir, "geom"), os.path.join(outdir, "truth")
    need = []
    if do_o1 and not os.path.exists(os.path.join(truth, key + "_o1.json")): need.append("o1")
    if do_o2 and not os.path.exists(os.path.join(truth, key + "_o2.json")): need.append("o2")
    if not need: return key, "skip", 0
    t0 = time.time()
    try:
        old = sys.argv
        for what in need:
            sys.argv = ["x", geom, truth, key] + (["--o2", "--no-o1"] if what == "o2" else [])
            import io, contextlib
            with contextlib.redirect_stdout(io.StringIO()): arb.main()
        sys.argv = old
        return key, "ok", time.time() - t0
    except Exception as e:
        return key, "ERRORE " + repr(e) + " " + traceback.format_exc()[-500:], time.time() - t0

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("outdir"); ap.add_argument("keys"); ap.add_argument("--o2-keys", default=None)
    ap.add_argument("--cores", type=int, default=40); ap.add_argument("--no-o1", action="store_true"); a = ap.parse_args()
    keys = [k.strip() for k in open(a.keys) if k.strip()]
    o2 = set(k.strip() for k in open(a.o2_keys)) if a.o2_keys else set()
    jobs = [(a.outdir, k, not a.no_o1, k in o2) for k in keys]
    n_err = 0
    with Pool(a.cores) as p:
        for key, st, dt in p.imap_unordered(one, jobs):
            if st.startswith("ERRORE"): n_err += 1; print(key, st, flush=True)
    print(f"arbitro: {len(jobs)} chiavi, errori {n_err}", flush=True)
    sys.exit(1 if n_err else 0)
