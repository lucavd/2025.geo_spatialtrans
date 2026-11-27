import json, h5py, collections
from g1_common import XEN, BIN
t = json.load(open(f"{XEN}_gene_panel.json"))["payload"]["targets"]
desc = collections.Counter(x.get("type", {}).get("descriptor") for x in t); print("descriptors", desc)
ens = [(x["type"]["data"].get("name"), x["type"]["data"].get("id"), x["type"].get("descriptor")) for x in t if str(x.get("type", {}).get("data", {}).get("id", "")).startswith("ENSG")]
print("targets with ENSG", len(ens), "unique ids", len(set(e[1] for e in ens)), "unique names", len(set(e[0] for e in ens)))
with h5py.File(f"{BIN}/square_002um/filtered_feature_bc_matrix.h5", "r") as f:
    fid = set(x.decode() for x in f["matrix/features/id"][:])
with h5py.File(f"{BIN}/square_002um/raw_feature_bc_matrix.h5", "r") as f:
    fid_raw = set(x.decode() for x in f["matrix/features/id"][:])
ids = set(e[1] for e in ens)
print("shared filtered", len(ids & fid), "shared raw", len(ids & fid_raw), "n features", len(fid))
print("not in visium:", sorted((e[0], e[1]) for e in ens if e[1] not in fid))
dup = [k for k, v in collections.Counter(e[0] for e in ens).items() if v > 1]; print("dup names", dup)
