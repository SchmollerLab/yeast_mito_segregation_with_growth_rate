import os

import pandas as pd

cwd_path = os.path.dirname(os.path.abspath(__file__))

df_rel_path = r'\tables_out\start_cell_1010---_len32_final.csv'

df_path = f'{cwd_path}{os.sep}{df_rel_path}'

df = pd.read_csv(df_path, index_col=['strain', 'mtdna_ratio', 'simulation_index'])

time_col = df.loc[('WT and WT',1.0,0), 'time'].values

idxs = df.index.unique()

for idx in idxs:
    df.loc[idx, 'time'] = time_col
    
df.reset_index().to_csv(df_path, index=False)