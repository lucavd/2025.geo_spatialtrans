
import pandas as pd
pd.set_option('display.width',250); pd.set_option('display.max_columns',30); pd.set_option('display.max_rows',300)
d=pd.read_csv('results/S1.5b/S1.5b_verdicts.csv')
print(d.arm.value_counts().to_string()); print(d.id.value_counts().to_string())
print(d.to_string())
a=pd.read_csv('results/S1.5b/S1.5b_arm_rule.csv'); print(a.to_string())
r=pd.read_csv('results/S1.5b/S1.5b_real_refs_roi_type.csv'); print(r.shape); print(r.columns.tolist()); print(r.head(3).to_string())
