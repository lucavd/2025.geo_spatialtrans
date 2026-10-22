import sys, time, pandas as pd
out = "/mnt/micron/geo_spatialtrans/R6/"
try:
    import decoupler as dc
    for org in ["human", "mouse"]:
        df = dc.op.resource("PanglaoDB", organism=org)
        df.to_csv(out + f"PanglaoDB_omnipath_{org}.tsv", sep="\t", index=False)
        print(org, df.shape, list(df.columns)[:12], flush=True)
except Exception as e:
    print("decoupler failed:", repr(e)[:500])
