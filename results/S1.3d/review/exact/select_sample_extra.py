
import random, glob, os
D="/mnt/micron/geo_spatialtrans/S1.3d"
H=open(f"{D}/keys_H.txt").read().split()
h2=[k for k in H if os.path.exists(f"{D}/truth/{k}_o2_s2.parquet")]
rng=random.Random(20261016)
o2=rng.sample(sorted(h2),3); o1=rng.sample(sorted(H),4)
with open(f"{D}/review_exact/keys_sample_extra.txt","w") as f:
    for i,k in enumerate(o2+["S_cperf"]): f.write(f"o2 {k} {3000+i}\n")
    for i,k in enumerate(o1+["S_cperf"]): f.write(f"o1 {k} {4000+i}\n")
print(len(h2), o2, o1)
