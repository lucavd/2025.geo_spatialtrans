#!/usr/bin/env python3
"""tools/F3.2_download.py - F3.2: manifest e scaricamento riprendibile dei dati di verifica e del ponte di versione.

Uso:  python3 tools/F3.2_download.py manifest   -> data/holdout/F3.2_download_manifest.csv (dall'inventario di R1 + CDN per P1/P5)
      python3 tools/F3.2_download.py download   -> scarica, verifica, scrive data/holdout/F3.2_download_status.csv
      python3 tools/F3.2_download.py extract    -> estrae i tar accanto all'originale (<nome>_x/, marcatore .extracted)
Integrita' (cancello G-2 a): md5 uguale al pubblicato (inventario di R1) quando c'e'; altrimenti ETag multiparte del CDN
ricalcolato (md5 delle md5 delle parti, parti da 8 MiB; la dimensione della parte si ricava dal numero di parti dell'ETag).
Dati di verifica in /mnt/micron/geo_spatialtrans/data_holdout/<etichetta>/, ponte di versione in .../data_bridge/<etichetta>/.
"""
import csv, hashlib, math, os, subprocess, sys, tarfile, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
INV = ROOT / "data/real/inventory/files_10x_pages_2026-09-18.csv"
OUT = ROOT / "data/holdout"; OUT.mkdir(parents=True, exist_ok=True)
MAN = OUT / "F3.2_download_manifest.csv"; STATUS = OUT / "F3.2_download_status.csv"
MICRON = Path("/mnt/micron/geo_spatialtrans")
CDN = "https://cf.10xgenomics.com/samples/spatial-exp"
# etichetta, ruolo, dataset_id 10x, versione Space Ranger, tessuto, livello
SETS = [
    ("CRC_P1_SR300", "holdout", "Visium_HD_Human_Colon_Cancer_P1", "3.0.0", "tumore e stroma del colon (P1)", "L2"),
    ("CRC_P5_SR300", "holdout", "Visium_HD_Human_Colon_Cancer_P5", "3.0.0", "tumore e stroma del colon (P5)", "L2"),
    ("CRC_6p5mm_SR410", "holdout", "Visium_HD_6p5mm_Human_Colon_Cancer", "4.1.0", "tumore e stroma del colon (6,5 mm)", "L2"),
    ("MouseBrain_FFPE_SR401", "holdout", "Visium_HD_Mouse_Brain", "4.0.1", "cervello di topo, seconda sezione", "L2"),
    ("HumanKidney_FFPE_SR401", "holdout", "Visium_HD_Human_Kidney_FFPE", "4.0.1", "rene umano", "L3"),
    ("RatLiver_6p5mm_SR410", "holdout", "Visium_HD_6p5mm_Rat_Liver", "4.1.0", "fegato di ratto", "L3"),
    ("CRC_P2_SR300", "bridge", "Visium_HD_Human_Colon_Cancer", "3.0.0", "tumore e stroma del colon (P2, taratura)", "ponte"),
]
KEEP = ("_tissue_image.", "_binned_outputs.tar.gz", "_spatial.tar.gz", "_segmented_outputs.tar.gz", "_barcode_mappings.parquet",
        "_metrics_summary.csv", "_web_summary.html", "_image.tif")
BRIDGE_SKIP = ("_tissue_image.", "_image.tif", "_segmented_outputs.tar.gz")   # ponte: l'immagine si confronta con l'md5 pubblicato

def head(url):
    r = subprocess.run(["curl", "-sI", "--max-time", "20", url], capture_output=True, text=True).stdout
    code, h = None, {}
    for ln in r.replace("\r", "").split("\n"):
        if ln.startswith("HTTP/"): code = int(ln.split()[1])
        elif ":" in ln: k, v = ln.split(":", 1); h[k.strip().lower()] = v.strip()
    return code, h

def manifest():
    inv = list(csv.DictReader(open(INV)))
    p2_names = [r["filename"] for r in inv if r["dataset_id"] == "Visium_HD_Human_Colon_Cancer" and "/3.0.0/" in r["url_download"]]
    rows = []
    for lab, role, ds, ver, tis, lev in SETS:
        cand = [r for r in inv if r["dataset_id"] == ds and f"/{ver}/" in r["url_download"]]
        if cand:
            items = [(r["filename"], r["md5"], r["url_download"], "inventario R1") for r in cand]
        else:   # P1, P5: stessi nomi dei file di P2 in 3.0.0, sul CDN
            items = [(fn.replace("Visium_HD_Human_Colon_Cancer", ds), "", f"{CDN}/{ver}/{ds}/" + fn.replace("Visium_HD_Human_Colon_Cancer", ds), "CDN, nomi di P2 3.0.0") for fn in p2_names]
        for fn, md5, url, src in items:
            if not any(k in fn for k in KEEP): continue
            if role == "bridge" and any(k in fn for k in BRIDGE_SKIP): continue
            code, h = head(url)
            if code != 200: print("ASSENTE", code, url); continue
            dest = MICRON / ("data_bridge" if role == "bridge" else "data_holdout") / lab / fn
            rows.append(dict(label=lab, role=role, dataset_id=ds, spaceranger=ver, tissue=tis, level=lev, filename=fn,
                             bytes=h.get("content-length", ""), md5_published=md5, etag=h.get("etag", "").strip('"'),
                             url=url, source=src, dest=str(dest)))
    with open(MAN, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    tot = sum(int(r["bytes"] or 0) for r in rows)
    print(f"manifest: {len(rows)} file, {tot/1e9:.1f} GB -> {MAN}")

def digests(path, etag):
    n_parts = int(etag.split("-")[1]) if "-" in etag else 1
    size = path.stat().st_size
    part = 8 * 1024 * 1024
    if "-" in etag and math.ceil(size / part) != n_parts:      # dimensione della parte diversa da 8 MiB: si ricava dal numero di parti
        for mb in (5, 10, 16, 25, 32, 50, 64, 100, 128):
            if math.ceil(size / (mb * 1024 * 1024)) == n_parts: part = mb * 1024 * 1024; break
    whole, parts = hashlib.md5(), []
    with open(path, "rb") as f:
        while True:
            b = f.read(part)
            if not b: break
            whole.update(b); parts.append(hashlib.md5(b).digest())
    et = hashlib.md5(b"".join(parts)).hexdigest() + f"-{len(parts)}" if "-" in etag else whole.hexdigest()
    return whole.hexdigest(), et, part

def download():
    rows = list(csv.DictReader(open(MAN)))
    done = {}
    if STATUS.exists():
        for r in csv.DictReader(open(STATUS)): done[r["dest"]] = r
    out = list(done.values())
    for r in rows:
        dest = Path(r["dest"]); dest.parent.mkdir(parents=True, exist_ok=True)
        if r["dest"] in done and done[r["dest"]]["status"] == "OK" and dest.exists(): print("SKIP", dest.name, flush=True); continue
        ok, info = False, {}
        for att in (1, 2, 3):
            print(f"[{time.strftime('%FT%T')}] GET {dest.name} ({att})", flush=True)
            rc = subprocess.run(["curl", "-sS", "-L", "--retry", "5", "--retry-delay", "10", "-C", "-", "-o", str(dest), r["url"]]).returncode
            if rc not in (0, 33): time.sleep(20); continue
            code, h = head(r["url"]); et_now = h.get("etag", "").strip('"')
            md5, et, part = digests(dest, r["etag"])
            if r["md5_published"]:
                ok = md5 == r["md5_published"]; how = "md5 pubblicato"
            else:
                ok = (et == r["etag"]) and (et_now == r["etag"]); how = "ETag multiparte"
            info = dict(md5_observed=md5, etag_recomputed=et, part_bytes=part, etag_at_end=et_now, check=how)
            if ok: break
            print("MISMATCH", dest.name, info, flush=True); dest.unlink(missing_ok=True)
        rec = dict(r, **info, bytes_observed=dest.stat().st_size if dest.exists() else 0, status="OK" if ok else "FAIL",
                   attempts=att, finished=time.strftime("%FT%T%z"))
        out = [o for o in out if o["dest"] != r["dest"]] + [rec]
        with open(STATUS, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rec)); w.writeheader(); w.writerows(out)
        print(f"[{time.strftime('%FT%T')}] {rec['status']} {dest.name} {rec['bytes_observed']} B ({info.get('check')})", flush=True)
    n_fail = sum(o["status"] != "OK" for o in out)
    print(f"RIEPILOGO: OK={len(out)-n_fail} FAIL={n_fail}"); sys.exit(1 if n_fail else 0)

def extract():
    for r in csv.DictReader(open(STATUS)):
        p = Path(r["dest"])
        if not p.name.endswith(".tar.gz") or r["status"] != "OK": continue
        x = p.parent / (p.name.split(r["dataset_id"] + "_")[1][:-len(".tar.gz")] + "_x")
        if (x / ".extracted").exists(): print("SKIP", x); continue
        x.mkdir(exist_ok=True); print("EXTRACT", p.name, "->", x, flush=True)
        subprocess.run(["tar", "-xzf", str(p), "-C", str(x)], check=True); (x / ".extracted").write_text(time.strftime("%FT%T%z"))

if __name__ == "__main__":
    {"manifest": manifest, "download": download, "extract": extract}[sys.argv[1]]()
