
import random, collections
D="/mnt/micron/geo_spatialtrans/S1.3d"
rng=random.Random(20261015)
O2=[k for k in open(f"{D}/keys_O2.txt").read().split()]
N=open(f"{D}/keys_N.txt").read().split()
by=collections.defaultdict(list)
for k in O2: by[k[:2]].append(k)
quota={"A1":2,"A2":2,"A3":2,"A4":2,"A5":1,"A6":1}
o2=[k for a in sorted(quota) for k in rng.sample(sorted(by[a]),quota[a])]
byN=collections.defaultdict(list)
for k in N: byN[k[:2]].append(k)
q1={"A1":4,"A2":4,"A3":3,"A4":3,"A5":3,"A6":3}
o1=[k for a in sorted(q1) for k in rng.sample(sorted(byN[a]),q1[a])]
with open(f"{D}/review_exact/keys_sample.txt","w") as f:
    for i,k in enumerate(o2): f.write(f"o2 {k} {1000+i}\n")
    for i,k in enumerate(o1): f.write(f"o1 {k} {2000+i}\n")
print(o2); print(o1)
