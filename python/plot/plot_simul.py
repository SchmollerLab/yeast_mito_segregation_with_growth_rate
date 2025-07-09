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

data = df.groupby(
    ['growth_rate_ratio' , 'simulation_index', 'time']
).agg(mean_h=('h', 'mean')).reset_index()

sns.boxplot(
    data=data, 
    x='time',
    y='mean_h',
    hue='growth_rate_ratio'
)
plt.show()