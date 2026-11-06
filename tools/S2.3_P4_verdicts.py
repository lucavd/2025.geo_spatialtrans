"""tools/S2.3_P4_verdicts.py — S2.3: esiti di P-4 (intervallo [0.8, 1.25]) e innesco di R-2, come scritto e post hoc PH-1."""
import pandas as pd
out = []
for suf, lab in (("", "come scritto"), ("_PH1_nocorr", "PH-1 senza correzioni")):
    T = pd.read_csv(f"results/S2.3/S2.3_P4_glass{suf}.csv")
    T["verdict"] = ["confermata" if (lo >= 0.8 and hi <= 1.25) else ("smentita per eccesso" if lo > 1.25 else ("smentita per difetto" if hi < 0.8 else "indeterminata")) for lo, hi in zip(T.lo, T.hi)]
    T["variant"] = lab; out.append(T)
    far = T[T.band.isin(["32-64", "64-128", "128-256"]) & (T.verdict == "smentita per eccesso")].groupby("key").size()
    print(lab, "R-2: sezioni con >= 2 fasce smentite per eccesso:", dict(far[far >= 2]), "-> innesco", int((far >= 2).sum()) >= 2)
    F = pd.read_csv(f"results/S2.3/S2.3_R2_fits{suf}.csv"); w = F.loc[F.groupby("key").aic.idxmin(), ["key", "candidate"]]
    print("  vincente per sezione:", dict(zip(w.key, w.candidate)))
P = pd.concat(out); P.to_csv("results/S2.3/S2.3_P4_verdicts.csv", index=False)
print(P[["variant", "key", "band", "ratio", "lo", "hi", "verdict"]].round(3).to_string())
