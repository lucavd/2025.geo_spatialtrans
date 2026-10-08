import pandas as pd, numpy as np
nuc = pd.read_parquet("results/R2/R2_nuclei_all.parquet", columns=["archetype","roi_id","method","scale","keep","x_um","y_um","um_per_px","area_um2","orientation","eccentricity"])
nuc = nuc[(nuc.scale==1)&nuc.keep]
rois = pd.read_csv("results/R2/R2_rois_checked.csv", float_precision="round_trip")
rows=[]
for (A,roi,me),d in nuc.groupby(["archetype","roi_id","method"]):
    if me not in ("cellpose_rgb","spaceranger","stardist_he"): continue
    upp=d.um_per_px.iloc[0]; fx=(d.x_um/upp)%1; fy=(d.y_um/upp)%1
    r=rois[(rois.archetype==A)&(rois.roi_id==roi)].iloc[0]; S=r.side_px*r.um_per_px
    Wm = 912 if A=="A5" else 913; px=S/Wm
    on_edge = (np.abs(d.y_um/px-np.round(d.y_um/px))<1e-9)|(np.abs(d.x_um/px-np.round(d.x_um/px))<1e-9)
    # quantization step of centroid in native px: fraction lying on multiples of 1/k
    q2=np.mean((np.abs(fx*2-np.round(fx*2))<1e-6)&(np.abs(fy*2-np.round(fy*2))<1e-6))
    ar=np.round(d.area_um2/upp**2)  # n pixels
    th=np.degrees(np.pi/2-d.orientation)%180; near_ax=np.mean((np.minimum(th%90,90-th%90))<5)
    rows.append(dict(A=A,roi=roi,method=me,n=len(d),frac_halfpx_grid=q2,n_unique_area=d.area_um2.nunique(),on_mask_edge=int(on_edge.sum()),frac_theta_pq_within5_of_axes=near_ax))
o=pd.DataFrame(rows); pd.set_option("display.width",200)
print(o.groupby(["A","method"]).agg(n=("n","sum"),halfpx=("frac_halfpx_grid","mean"),uniq_area=("n_unique_area","median"),on_edge=("on_mask_edge","sum"),near_ax=("frac_theta_pq_within5_of_axes","mean")).to_string())
o.to_csv("results/R3/review/code/rc_quantization.csv",index=False)
