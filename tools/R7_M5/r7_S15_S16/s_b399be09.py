
import pandas as pd
pd.set_option('display.width',250); pd.set_option('display.max_columns',30); pd.set_option('display.max_rows',400); pd.set_option('display.max_colwidth',90)
for st in ['S1.5e','S1.6']:
    d=pd.read_csv(f'results/{st}/{st}_verdicts.csv')
    print('=====',st, d.shape); print(d.groupby('arm').size().to_string()); print(d.groupby(['id','source']).size().to_string())
for f in ['results/S1.5e/S1.5e_ruleA.csv','results/S1.5e/S1.5e_ruleB.csv','results/S1.5e/S1.5e_ruleQ.csv','results/S1.5e/S1.5e_counts.csv','results/S1.5e/S1.5e_cp.csv','results/S1.6/S1.6_counts.csv','results/S1.6/S1.6_CP2.csv','results/S1.5/S1.5_CP1_verdicts.csv']:
    x=pd.read_csv(f); print('==',f,x.shape); print(x.to_string())
for f in ['results/S1.6/S1.6_P1.csv','results/S1.6/S1.6_CP1.csv']:
    x=pd.read_csv(f); print('==',f,x.shape); print(x.columns.tolist()); print(x.head(12).to_string())
