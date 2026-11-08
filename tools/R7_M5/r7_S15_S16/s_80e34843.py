
import pandas as pd
pd.set_option('display.width',250); pd.set_option('display.max_columns',30); pd.set_option('display.max_rows',200)
d=pd.read_csv('results/S1.5/S1.5_verdicts.csv')
print(d.drop(columns=['tissue']).to_string())
