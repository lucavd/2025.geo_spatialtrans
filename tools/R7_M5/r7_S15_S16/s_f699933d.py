
import pandas as pd
pd.set_option('display.width',250); pd.set_option('display.max_columns',30); pd.set_option('display.max_rows',400); pd.set_option('display.max_colwidth',80)
d=pd.read_csv('results/S1.5c/S1.5c_verdicts.csv')
print(d.groupby(['arm']).size().to_string()); print(d.groupby('id').size().to_string())
print(d.source.value_counts().to_string())
print(d.verdict.value_counts().to_string())
print(d[d.arm==d.arm.unique()[0]].head(5).to_string())
print(d.arm.unique())
w=d[d.arm=='N-a'] if 'N-a' in d.arm.unique() else d[d.arm==d.arm.unique()[0]]
print(w.to_string())
print(pd.read_csv('results/S1.5c/S1.5c_rules.csv').to_string())
cp=pd.read_csv('results/S1.5c/S1.5c_cp1.csv'); print(cp.shape); print(cp.head(20).to_string())
