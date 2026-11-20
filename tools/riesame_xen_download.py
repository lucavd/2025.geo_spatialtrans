"""tools/riesame_xen_download.py - Riesame Parte B, braccio Xenium: scarica lo zip Xenium di P2 (GEO GSE280314) e la registrazione
pubblicata da 10x (GitHub), verifica l'integrita' (dimensione HTTP, CRC di tutti i membri dello zip, sha1 git del file GitHub), estrae
solo i file che servono e scrive l'inventario (md5) in <dest>/inventory.tsv. Rilanciabile (curl -C -). Uso: python3 tools/riesame_xen_download.py"""
import os, sys, json, hashlib, subprocess, zipfile, time
DEST = "/mnt/micron/geo_spatialtrans/data_xenium/P2"
ZIPN = "GSE280314_Xenium_V1_Human_Colon_Cancer_P2_CRC_Add_on_FFPE_outs.zip"
URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE280nnn/GSE280314/suppl/" + ZIPN
SIZE = 15874246655                       # Content-Length di GEO (HEAD del 2026-11-20)
GH = "https://raw.githubusercontent.com/10XGenomics/HumanColonCancer_VisiumHD/main/Figures/Xenium_Visium_Alignment/Xenium_P2_cell_info.csv.gz"
GH_SHA = "156bb0f3d3b33db64ed7f59a0b683891103eef15"; GH_SIZE = 48860122   # API GitHub contents (blob sha1)
WANT = ["transcripts.parquet", "nucleus_boundaries.parquet", "cell_boundaries.parquet", "cells.parquet", "experiment.xenium",
        "gene_panel.json", "metrics_summary.csv"]
def md5(p, bs=1 << 24):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(bs), b""): h.update(b)
    return h.hexdigest()
def log(*a): print(time.strftime("%H:%M:%S"), *a, flush=True)
os.makedirs(f"{DEST}/outs", exist_ok=True); os.makedirs(f"{DEST}/tenx_alignment", exist_ok=True)
zp = f"{DEST}/{ZIPN}"
if not (os.path.exists(zp) and os.path.getsize(zp) == SIZE):
    for k in range(5):
        r = subprocess.run(["curl", "-sS", "-L", "-C", "-", "--retry", "5", "-o", zp, URL])
        if os.path.exists(zp) and os.path.getsize(zp) == SIZE: break
        log("curl attempt", k, "rc", r.returncode, "size", os.path.getsize(zp) if os.path.exists(zp) else None)
assert os.path.getsize(zp) == SIZE, "size mismatch"
log("zip size OK", SIZE)
gp = f"{DEST}/tenx_alignment/Xenium_P2_cell_info.csv.gz"
if not os.path.exists(gp): subprocess.run(["curl", "-sS", "-L", "-o", gp, GH], check=True)
data = open(gp, "rb").read(); gsha = hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()
assert len(data) == GH_SIZE and gsha == GH_SHA, ("github file mismatch", len(data), gsha)
log("github file OK, git blob sha1", gsha)
inv = [("zip", ZIPN, SIZE, None), ("github", "Xenium_P2_cell_info.csv.gz", GH_SIZE, md5(gp))]
z = zipfile.ZipFile(zp)
names = z.namelist(); json.dump([(i.filename, i.file_size, i.compress_size, "%08x" % i.CRC) for i in z.infolist()], open(f"{DEST}/zip_listing.json", "w"))
log("members", len(names))
bad = z.testzip()                      # legge ogni membro e controlla il CRC32
assert bad is None, f"CRC error in {bad}"
log("testzip OK: CRC32 of all", len(names), "members")
for w in WANT:
    cand = [n for n in names if n.split("/")[-1] == w]
    cand = sorted(cand, key=lambda n: n.count("/"))[:1]        # il membro meno annidato (cartella outs)
    assert cand, f"missing {w}"
    n = cand[0]; out = f"{DEST}/outs/{w}"
    with z.open(n) as src, open(out, "wb") as dst:
        while True:
            b = src.read(1 << 24)
            if not b: break
            dst.write(b)
    inv.append(("member", n, os.path.getsize(out), md5(out)))
    log("extracted", n, os.path.getsize(out))
inv[0] = ("zip", ZIPN, SIZE, md5(zp))
with open(f"{DEST}/inventory.tsv", "w") as f:
    f.write("kind\tname\tbytes\tmd5\n")
    for r in inv: f.write("\t".join(map(str, r)) + "\n")
log("inventory written"); print(open(f"{DEST}/inventory.tsv").read())
