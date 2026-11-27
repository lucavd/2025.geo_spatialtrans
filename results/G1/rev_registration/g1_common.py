"""Independent re-computation of the S1 registration gate (X1). Common helpers.
Written from the rule text (docs/sonda/REGOLA_S1.md); no project code imported."""
import json, numpy as np
D = "/mnt/micron/geo_spatialtrans/data_lung/X1"
XEN = f"{D}/xenium/Xenium_V1_Human_Lung_Cancer_FFPE"
BIN = f"{D}/visium/binned/binned_outputs"
OUT = "/mnt/micron/geo_spatialtrans/G1_review/rev_registration"
XPIX_UM = 0.2125          # Xenium morphology pixel size (rule, chain step 1)
MARGIN_BINS = 75          # 150 um margin around each 1 mm region (2 um bins)

def to_h(p):
    return np.c_[p, np.ones(len(p))]

def proj_apply(H, p):
    q = to_h(np.asarray(p, float)) @ H.T
    return q[:, :2] / q[:, 2:3]

def fit_projective_lsq(src, dst):
    """Linear least squares for an 8-parameter projective map dst ~ H src,
    in centred/scaled coordinates (h33 fixed to 1)."""
    ms, ss = src.mean(0), src.std()
    md, sd = dst.mean(0), dst.std()
    s = (src - ms) / ss; d = (dst - md) / sd
    x, y = s[:, 0], s[:, 1]; u, v = d[:, 0], d[:, 1]
    n = len(s); z = np.zeros(n); o = np.ones(n)
    A = np.vstack([np.c_[x, y, o, z, z, z, -u * x, -u * y],
                   np.c_[z, z, z, x, y, o, -v * x, -v * y]])
    b = np.r_[u, v]
    h = np.linalg.lstsq(A, b, rcond=None)[0]
    Hn = np.array([[h[0], h[1], h[2]], [h[3], h[4], h[5]], [h[6], h[7], 1.0]])
    Ts = np.array([[1 / ss, 0, -ms[0] / ss], [0, 1 / ss, -ms[1] / ss], [0, 0, 1]])
    Td_inv = np.array([[sd, 0, md[0]], [0, sd, md[1]], [0, 0, 1]])
    H = Td_inv @ Hn @ Ts
    return H / H[2, 2]

def load_chain():
    return json.load(open(f"{OUT}/chain.json"))

def xenium_to_array(xy_um, chain):
    """Xenium um -> H&E full-res pixel (inverse of he_imagealignment, Xenium px 0.2125 um)
    -> continuous 2 um array coordinates (col, row), bin centres at integers."""
    Minv = np.array(chain["Minv"]); Hp2a = np.array(chain["H_px2arr"])
    px = proj_apply(Minv, np.asarray(xy_um, float) / XPIX_UM)
    return proj_apply(Hp2a, px)
