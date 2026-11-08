
import pandas as pd
pd.set_option('display.width',250); pd.set_option('display.max_columns',30); pd.set_option('display.max_rows',400); pd.set_option('display.max_colwidth',90)
d=pd.read_csv('results/S1.5d/S1.5d_verdicts.csv')
print(d.groupby('arm').size().to_string()); print(d.groupby('id').size().to_string())
print(d.groupby(['id','source']).size().to_string())
print(d.groupby('unit').size().shape)
print(d[(d.arm=='RV')&(d.archetype=='A2')].to_string())
for f in ['ruleA','ruleB','ruleC','ruleC_sensitivity_U3c','counts','cp']:
    x=pd.read_csv(f'results/S1.5d/S1.5d_{f}.csv'); print('==',f, x.shape); print(x.to_string())
