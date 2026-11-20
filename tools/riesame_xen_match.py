"""tools/riesame_xen_match.py - Riesame Parte B: criterio XD-25 (docs/riesame/DECISIONI_B_XENIUM.md). Quale file di registrazione di 10x
(GitHub, Xenium_P{1,2,5}_cell_info.csv.gz) corrisponde a quale pacchetto Xenium di GEO GSE280314 (P1, P2, P5)? Per ogni coppia il criterio
XD-24. GitHub: download con verifica dello sha1 git (API contents). GEO P1 e P5: lettura parziale dello zip (indice + cells.parquet) con
richieste HTTP Range, CRC32 del membro verificato da zipfile. Uso: .venv/bin/python -B tools/riesame_xen_match.py"""
import os, io, json, hashlib, zipfile, urllib.request, numpy as np, pandas as pd
M = "/mnt/micron/geo_spatialtrans/data_xenium"; GHAPI = "https://api.github.com/repos/10XGenomics/HumanColonCancer_VisiumHD/contents/Figures/Xenium_Visium_Alignment"
GHRAW = "https://raw.githubusercontent.com/10XGenomics/HumanColonCancer_VisiumHD/main/Figures/Xenium_Visium_Alignment/"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE280nnn/GSE280314/suppl/GSE280314_Xenium_V1_Human_Colon_Cancer_{p}_CRC_Add_on_FFPE_outs.zip"
PS = ["P1", "P2", "P5"]; info = {"bytes_fetched": {}}
def md5(p):
    h = hashlib.md5(); f = open(p, "rb")
    for b in iter(lambda: f.read(1 << 24), b""): h.update(b)
    return h.hexdigest()
class RangeFile(io.RawIOBase):
    def __init__(s, url):
        s.url = url; s.pos = 0; s.n = 0
        r = urllib.request.urlopen(urllib.request.Request(url, method="HEAD")); s.size = int(r.headers["Content-Length"])
    def seekable(s): return True
    def readable(s): return True
    def tell(s): return s.pos
    def seek(s, off, wh=0):
        s.pos = off if wh == 0 else (s.pos + off if wh == 1 else s.size + off); return s.pos
    def read(s, k=-1):
        if k is None or k < 0: k = s.size - s.pos
        if k == 0 or s.pos >= s.size: return b""
        end = min(s.size, s.pos + k) - 1
        d = urllib.request.urlopen(urllib.request.Request(s.url, headers={"Range": f"bytes={s.pos}-{end}"})).read()
        s.pos += len(d); s.n += len(d); return d
    def readinto(s, b):
        d = s.read(len(b)); b[:len(d)] = d; return len(d)
# 1) GitHub
meta = {x["name"]: x for x in json.load(urllib.request.urlopen(GHAPI))}
gh = {}
for p in PS:
    name = f"Xenium_{p}_cell_info.csv.gz"; d = f"{M}/{p}/tenx_alignment"; os.makedirs(d, exist_ok=True); fp = f"{d}/{name}"
    if not os.path.exists(fp): urllib.request.urlretrieve(GHRAW + name, fp); info["bytes_fetched"][name] = os.path.getsize(fp)
    data = open(fp, "rb").read(); sha = hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()
    assert sha == meta[name]["sha"] and len(data) == meta[name]["size"], (name, sha)
    gh[p] = pd.read_csv(fp); info[f"github_{p}"] = {"rows": int(len(gh[p])), "md5": md5(fp), "git_sha1": sha}
# 2) GEO: P2 dal pacchetto gia' estratto, P1 e P5 per lettura parziale
geo = {"P2": pd.read_parquet(f"{M}/P2/outs/cells.parquet")}; info["geo_P2"] = {"rows": int(len(geo["P2"])), "source": "extracted zip (inventory.tsv)"}
for p in ("P1", "P5"):
    out = f"{M}/{p}/outs/cells.parquet"; os.makedirs(os.path.dirname(out), exist_ok=True)
    rf = RangeFile(GEO.format(p=p)); z = zipfile.ZipFile(io.BufferedReader(rf, buffer_size=1 << 20))
    lst = [(i.filename, i.file_size, i.compress_size, "%08x" % i.CRC) for i in z.infolist()]
    json.dump(lst, open(f"{M}/{p}/zip_listing.json", "w"))
    cand = sorted([n for n, *_ in lst if n.split("/")[-1] == "cells.parquet"], key=lambda n: n.count("/"))[0]
    with z.open(cand) as src: data = src.read()          # zipfile verifica il CRC32 alla fine della lettura
    open(out, "wb").write(data); geo[p] = pd.read_parquet(out)
    info[f"geo_{p}"] = {"rows": int(len(geo[p])), "zip_bytes": rf.size, "member": cand, "member_bytes": len(data), "md5": md5(out), "n_members": len(lst)}
    info["bytes_fetched"][f"GEO_{p}_partial"] = rf.n
# 3) XD-24 per ogni coppia
def norm(s): return s.map(lambda v: v.decode() if isinstance(v, (bytes, bytearray)) else str(v)).str.strip().str.lower()
rows = []
for gp in PS:
    G = geo[gp].assign(key=norm(geo[gp].cell_id)).drop_duplicates("key").set_index("key")
    for hp in PS:
        H = gh[hp].assign(key=norm(gh[hp].cell_id)).drop_duplicates("key").set_index("key")
        com = G.index.intersection(H.index); r = {"geo": gp, "registration": hp, "n_geo": len(geo[gp]), "n_registration": len(gh[hp]),
                                                  "n_common": int(len(com)), "share_geo_in_registration": len(com) / len(geo[gp])}
        if len(com) > 1:
            d = H.loc[com, ["x_centroid", "y_centroid"]].to_numpy(float) - G.loc[com, ["x_centroid", "y_centroid"]].to_numpy(float)
            r.update({"dx_mean": d[:, 0].mean(), "dy_mean": d[:, 1].mean(), "dx_sd": d[:, 0].std(ddof=1), "dy_sd": d[:, 1].std(ddof=1)})
            for c in ("transcript_counts", "total_counts"): r[f"share_equal_{c}"] = float((G.loc[com, c].to_numpy() == H.loc[com, c].to_numpy()).mean())
        r["c1_ids"] = r["share_geo_in_registration"] >= 0.99; r["c2_constant"] = bool(len(com) > 1 and r["dx_sd"] < 1 and r["dy_sd"] < 1)
        r["c3_counts"] = bool(len(com) > 1 and r["share_equal_transcript_counts"] >= 0.99 and r["share_equal_total_counts"] >= 0.99)
        r["match"] = bool(r["c1_ids"] and r["c2_constant"] and r["c3_counts"]); rows.append(r)
T = pd.DataFrame(rows); os.makedirs("results/riesame", exist_ok=True); T.to_csv("results/riesame/xen_match.csv", index=False)
info["matches_per_geo"] = {gp: T[(T.geo == gp) & T.match].registration.tolist() for gp in PS}
json.dump(info, open("results/riesame/xen_match.json", "w"), indent=1, default=str)
print(json.dumps(info, indent=1, default=str)); pd.set_option("display.width", 250); print(T.round(4).to_string(index=False))
