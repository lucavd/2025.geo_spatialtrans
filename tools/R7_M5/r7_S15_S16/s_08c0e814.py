
import pandas as pd, numpy as np
pd.set_option('display.width',250); pd.set_option('display.max_columns',30); pd.set_option('display.max_colwidth',40)
T=pd.read_pickle('/tmp/r7_S15_S16/rows.pkl')
T['arm']=T.notes.str.extract(r'arm=([^ ;]+)')[0]
c=['step','prediction_id','archetype','arm','quantity','est','lo','hi','a','b','lo_new','hi_new','a_new','b_new','kappa_est','kappa_margin','verdict_recorded','verdict_new','verdict_new_stress']
SZC=['P-8a','P-8b','P-8c','P-10','P-14']; SP=['P-1','P-2','P-3','P-4','P-5','P-6','P-7']
x=T[(T.step=='S1.5b')&(T.prediction_id.isin(['P-1','P-2']))&(T.arm.isin(['group','recentred']))&(T.verdict_new!=T.verdict_recorded)]
print(x[c].to_string())
for arm in ['RV','W3','U']:
    x=T[(T.step=='S1.5d')&(T.arm==arm)&(T.prediction_id.isin(SZC+SP))&((T.verdict_new=='confermata')!=(T.verdict_recorded=='confermata'))]
    print('== S1.5d',arm); print(x[c].to_string())
for arm in ['Q','RV']:
    x=T[(T.step=='S1.5e')&(T.arm==arm)&(T.prediction_id.isin(SP))&((T.verdict_new=='confermata')!=(T.verdict_recorded=='confermata'))]
    print('== S1.5e',arm); print(x[c].to_string())
