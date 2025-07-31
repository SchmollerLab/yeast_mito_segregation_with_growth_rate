import os

import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import ks_2samp

from python.segregation_simulator import (
    wt_doubling_time,
    startbud,
    ngen,
    ndau,
    nspl,
    number_simulations,
    number_of_cells,
    tables_path
)

from python.segregation_simulator.utils import (
    get_cell_inital_state, 
    get_table_filenames, 
    get_single_cells_filename
)

start_cell_type = '00...11' # '11...00', '1010...', '00...11

nuc, nuc_format = get_cell_inital_state(start_cell_type)
table_basename, table_filename, single_cells_filename = get_table_filenames(
    nuc_format, number_of_cells, number_simulations=number_simulations
)

print(f'Loading table {table_filename}')

df_filepath = os.path.join(tables_path, table_filename)

df = pd.read_csv(df_filepath)

data = df.groupby(
    ['growth_rate_ratio' , 'simulation_index', 'time']
).agg(mean_h=('mean_h', 'mean')).reset_index()

data['time'] = data['time'].round(2)

last_timepoint = data['time'].max()
data_last_timepoint = data[data['time'] == last_timepoint]

print('-'*100)
print('Groups statistics:')
print(data.groupby(['time', 'growth_rate_ratio']).describe())
print('-'*100)

fig, ax = plt.subplots(figsize=(10, 6))

sns.boxplot(
    data=data, 
    x='time',
    y='mean_h',
    hue='growth_rate_ratio'
)
fig.suptitle(f'Start cell = {start_cell_type}')

plt.show()