
import pandas as pd
T=pd.read_pickle('/tmp/r7_S15_S16/rows.pkl')
T['clip']=T.notes.str.contains('clipped')
print(T.groupby('step').clip.sum().to_string())
print(T.groupby(['step','transition']).size().to_string())
print((T.verdict_new_stress!=T.verdict_new).groupby(T.step).sum().to_string())
print(T.groupby(['est_unit','margin_unit']).size().to_string())
