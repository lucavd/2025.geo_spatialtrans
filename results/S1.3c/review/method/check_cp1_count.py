
import numpy as np, pandas as pd
cp=pd.read_csv('results/S1.3c/S1.3c_cp_roi.csv'); c1=cp[cp.cp=='CP-1']
t=c1.assign(p=c1['pass'].astype(str)=='True').groupby(['archetype','k']).p.agg(['sum','size','mean'])
print(t.to_string())
# riproduzione del meccanismo: colonna object con numpy.bool_ e NaN
s=pd.Series([np.bool_(True),np.bool_(True),False,np.bool_(True),np.nan],dtype=object)
print('sum object con np.bool_:', s.dropna().sum(), int(s.dropna().sum()), '| mean:', s.dropna().mean())
ba_r, ba_p = np.float64(0.6), np.float64(0.5); ok=(0.9<=0.95<=1.1) and (ba_p<ba_r); print(type(ok))
