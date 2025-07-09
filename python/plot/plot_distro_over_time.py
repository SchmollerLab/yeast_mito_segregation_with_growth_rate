import os

import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

cwd_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
tables_path = os.path.join(cwd_path, 'segregation_simulator', 'tables_out')

nuc = np.zeros(32, dtype=int) # 1, 0, 1, 0, etc. 
nuc[::2] = 1
nuc = tuple(nuc)
nuc_format = ''.join([str(n) for n in nuc])

table_basename = f'start_cell_{nuc_format}.csv'

df_filepath = os.path.join(tables_path, table_basename)

df = pd.read_csv(df_filepath)

ncols = 4
nrows = 4

for growth_rate_ratio in df['growth_rate_ratio'].unique():
    df_grr = df[df['growth_rate_ratio'] == growth_rate_ratio].copy()
    fig, ax = plt.subplots(nrows, ncols, figsize=(12, 10))
    ax = ax.flatten()
    for i, time_i in enumerate(df_grr['time'].unique()):
        axes = ax[i]
        df_i = df_grr[df_grr['time'] == time_i]
        sns.histplot(
            data=df_i, 
            ax=axes,
            x='h',
            stat='probability',
            # hue='growth_rate_ratio'
        )
    fig.suptitle(f'Growth rate ratio = {growth_rate_ratio}')
    
plt.show()