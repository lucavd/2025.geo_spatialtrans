# 09_cell2location.py — cell2location (GPU) on the files written by 09a_export_c2l.R (ADEMP §4, §6;
# DECISIONS 2026-10-30, 2026-10-31, 2026-11-04). First use: timing of one replicate -> n_c2l amendment.
#
# Settings fixed in ADEMP §4 (not changed here):
#   reference gene filter  filter_genes(cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
#   N_cells_per_location   true cells per spot of the condition (oracle information)
#   detection_alpha        200 (primary), 20 (sensitivity, --alpha 20)
#   reference of arm P     retrained in every replicate; all other arms use the real reference, trained once
# Settings NOT in ADEMP, taken from the official tutorial and to be written in the n_c2l amendment before
# any result: RegressionModel batch_key = donor_id, max_epochs 250, posterior 1000 samples;
# Cell2location max_epochs 30000, batch_size None, train_size 1, posterior 1000 samples.
# Seed: scvi.settings.seed = 1000 * rep + 100 * condition + dataset slot (slots of 05_deconvolve.R);
# real reference seed 0, reference of P seed 1000 * rep.
#
# Output per dataset: <c2l dir>/r<rep>/<cond>/<dataset>/props_a<alpha>.csv (spots x types, cell
# proportions = posterior mean abundance / its sum over types). Proportions are cleaned and scored in R
# with the functions of 05_deconvolve.R (clean_props, rmse, jsd), in a later script.
# Timing: <c2l dir>/c2l_timing.csv, one row per step (seconds, peak RAM of the process, peak GPU memory).
#
# Usage (working directory = thesis/):
#   .venv-c2l/bin/python scripts/09_cell2location.py <tissue> [main|pilot] [rep_from] [rep_to]
#       [--alpha 200] [--conds c02_d0174,...] [--datasets anchor,...]

import argparse, os, resource, sys, time
import numpy as np, pandas as pd, scipy.io, anndata
import torch, scvi, cell2location
from cell2location.models import RegressionModel, Cell2location
from cell2location.utils.filtering import filter_genes

ap = argparse.ArgumentParser()
ap.add_argument("tissue", choices=["colon", "lymphoid", "brain"])
ap.add_argument("mode", nargs="?", default="main", choices=["main", "pilot"])
ap.add_argument("rep_from", nargs="?", type=int, default=1)
ap.add_argument("rep_to", nargs="?", type=int, default=None)
ap.add_argument("--alpha", type=float, default=200)
ap.add_argument("--conds", default=None, help="comma-separated subset of conditions (test runs)")
ap.add_argument("--datasets", default=None, help="comma-separated subset of datasets (test runs)")
a = ap.parse_args()
rep_to = a.rep_to if a.rep_to is not None else a.rep_from

data_dir = os.environ.get("TESI_DATA", "")
if not data_dir or not os.path.isdir(data_dir):
    sys.exit("TESI_DATA non impostata")
base = os.path.join(data_dir, "c2l", a.mode, a.tissue)
types = open(os.path.join(base, "types.txt")).read().split()
hd_depth = {"colon": 174, "lymphoid": 62, "brain": 121}[a.tissue]
cond_keys = [f"c{n:02d}_d{hd_depth:04d}" for n in (2, 10)]          # as in 05_deconvolve.R
slot = {"anchor": 1, "E": 2, "Cminus": 3, "Eprime": 4, "O": 5, "P": 6, "G": 7}
conds = a.conds.split(",") if a.conds else cond_keys
only = set(a.datasets.split(",")) if a.datasets else None
timing_file = os.path.join(base, "c2l_timing.csv")
print(f"cell2location {cell2location.__version__}, scvi-tools {scvi.__version__}, torch {torch.__version__}, "
      f"GPU {torch.cuda.is_available()}", flush=True)
if not torch.cuda.is_available():
    sys.exit("No GPU visible to torch")


def read_counts(d, obs_file):
    m = scipy.io.mmread(os.path.join(d, "counts.mtx")).T.tocsr().astype(np.float32)   # obs x genes
    genes = open(os.path.join(d, "genes.txt")).read().split("\n")[: m.shape[1]]
    if obs_file.endswith(".csv"):
        obs = pd.read_csv(os.path.join(d, obs_file)).set_index("cell")
    else:
        obs = pd.DataFrame(index=open(os.path.join(d, obs_file)).read().split("\n")[: m.shape[0]])
    return anndata.AnnData(X=m, obs=obs, var=pd.DataFrame(index=genes))


def log_time(step, rep, cond, dataset, seconds):
    peak_ram_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2   # Linux: KB
    peak_gpu_gb = torch.cuda.max_memory_allocated() / 1024 ** 3
    row = pd.DataFrame([dict(tissue=a.tissue, rep=rep, cond=cond, dataset=dataset, step=step, alpha=a.alpha,
                             seconds=round(seconds, 1), peak_ram_gb=round(peak_ram_gb, 2),
                             peak_gpu_gb=round(peak_gpu_gb, 2), time=time.strftime("%Y-%m-%d %H:%M:%S"))])
    row.to_csv(timing_file, mode="a", header=not os.path.exists(timing_file), index=False)
    print(f"{row.time[0]} {step} rep {rep} {cond} {dataset}: {seconds / 60:.1f} min, "
          f"RAM {peak_ram_gb:.1f} GB, GPU {peak_gpu_gb:.1f} GB", flush=True)


def signatures(d, seed, rep):
    """Reference signatures (genes x types); trained once per reference and cached in d."""
    f = os.path.join(d, "signatures.csv")
    if os.path.exists(f):
        return pd.read_csv(f, index_col=0)
    t0 = time.time()
    scvi.settings.seed = seed
    ad = read_counts(d, "cells.csv")
    ad = ad[:, filter_genes(ad, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)].copy()
    RegressionModel.setup_anndata(adata=ad, batch_key="donor_id", labels_key="cell_type")
    mod = RegressionModel(ad)
    mod.train(max_epochs=250)
    ad = mod.export_posterior(ad, sample_kwargs={"num_samples": 1000, "batch_size": 2500})
    fn = ad.uns["mod"]["factor_names"]
    sig = ad.varm["means_per_cluster_mu_fg"][[f"means_per_cluster_mu_fg_{k}" for k in fn]].copy()
    sig.columns = fn
    sig.to_csv(f + ".tmp"); os.replace(f + ".tmp", f)
    log_time("reference", rep, "", os.path.basename(d), time.time() - t0)
    return sig


sig_real = signatures(os.path.join(base, "ref_real"), 0, 0)
for rep in range(a.rep_from, rep_to + 1):
    rd = os.path.join(base, f"r{rep:03d}")
    sig_p = None
    for ci, cond in enumerate(cond_keys, start=1):
        if cond not in conds:
            continue
        for ds in sorted(os.listdir(os.path.join(rd, cond)), key=lambda s: slot.get(s, 99)):
            if only and ds not in only:
                continue
            dd = os.path.join(rd, cond, ds)
            f_out = os.path.join(dd, f"props_a{a.alpha:g}.csv")
            if os.path.exists(f_out):
                continue
            if ds == "P" and sig_p is None:
                sig_p = signatures(os.path.join(rd, "ref_P"), 1000 * rep, rep)
            sig = sig_p if ds == "P" else sig_real
            t0 = time.time()
            scvi.settings.seed = 1000 * rep + 100 * ci + slot[ds]
            torch.cuda.reset_peak_memory_stats()
            sp = read_counts(dd, "spots.txt")
            n_cells = int(open(os.path.join(dd, "n_cells.txt")).read())
            genes = np.intersect1d(sp.var_names, sig.index)
            sp = sp[:, genes].copy()
            Cell2location.setup_anndata(adata=sp)
            mod = Cell2location(sp, cell_state_df=sig.loc[genes], N_cells_per_location=n_cells,
                                detection_alpha=a.alpha)
            mod.train(max_epochs=30000, batch_size=None, train_size=1)
            sp = mod.export_posterior(sp, sample_kwargs={"num_samples": 1000, "batch_size": mod.adata.n_obs})
            ab = pd.DataFrame(sp.obsm["means_cell_abundance_w_sf"].values, index=sp.obs_names,
                              columns=sp.uns["mod"]["factor_names"])
            props = ab.div(ab.sum(axis=1), axis=0).reindex(columns=types, fill_value=0.0)
            props.to_csv(f_out + ".tmp"); os.replace(f_out + ".tmp", f_out)
            log_time("mapping", rep, cond, ds, time.time() - t0)
print("Done", flush=True)
