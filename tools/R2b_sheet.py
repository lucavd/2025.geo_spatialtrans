
import pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image
df = pd.read_csv("results/R2b/R2b_windows.csv")
fig, axes = plt.subplots(6, 4, figsize=(16, 24), facecolor="white")
for ax, (_, r) in zip(axes.ravel(), df.iterrows()):
    im = np.asarray(Image.open(r.png)); ax.imshow(im); ax.set_axis_off()
    ax.set_title(f"{r.win_id} ({r.roi_id}) {r.side_um:.0f} um  CP={r.n_cellpose_rgb:.0f} SD={r.n_stardist_he:.0f}", fontsize=10)
fig.tight_layout(); fig.savefig("results/R2b/R2b_windows_sheet.png", dpi=60)

