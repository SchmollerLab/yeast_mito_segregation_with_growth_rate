import os

import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import ks_2samp

from yeast_mito_segregation.segregation_simulator import (
    wt_doubling_time,
    startbud,
    ngen,
    ndau,
    nspl,
    number_simulations,
    number_of_cells,
    tables_path,
    table_endname
)

from yeast_mito_segregation.segregation_simulator.utils import (
    get_cell_inital_state, 
    get_table_filenames, 
    get_single_cells_filename
)

start_cell_type = '1010...' # '11...00', '1010...', '00...11
num_rows = 3

nuc, nuc_format = get_cell_inital_state(start_cell_type)
table_basename, table_filename, single_cells_filename = get_table_filenames(
    nuc_format, appended_text=table_endname
)

print(f'Loading table {table_filename}')

df_filepath = os.path.join(tables_path, table_filename)

df = pd.read_csv(df_filepath)

df = df.groupby(
    ['growth_rate_ratio' , 'mtdna_ratio', 'simulation_index', 'time']
).agg(mean_h=('mean_h', 'mean')).reset_index()

df['time'] = df['time'].round(1)

df_WT = df[df['growth_rate_ratio'] == 1.0]
df_WT_ratio_one = df_WT[df_WT['mtdna_ratio'] == 1.0]

growth_rate_ratios = df['growth_rate_ratio'].unique()

for growth_rate_ratio in growth_rate_ratios:
    if growth_rate_ratio == 1.0:
        continue

    df_growth_rate = df[df['growth_rate_ratio'] == growth_rate_ratio]
    mtDNA_ratios = df_growth_rate['mtdna_ratio'].unique()
    
    num_cols = (
        len(mtDNA_ratios) // num_rows + int(len(mtDNA_ratios) % num_rows > 0)
    )
    
    if num_cols == 0:
        import pdb; pdb.set_trace()
        continue
    elif num_cols == 1:
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        ax = [ax]
    else:
        fig, ax = plt.subplots(num_rows, num_cols, figsize=(19, 10.5))
        ax = ax.flatten()
        fig.subplots_adjust(
            top=0.92, 
            bottom=0.05, 
            left=0.05, 
            right=0.95, 
            hspace=0.4, 
            wspace=0.3
        )
    
    for a, mtDNA_ratio in enumerate(mtDNA_ratios):
        df_growth_rate_ratio = df_growth_rate[
            df_growth_rate['mtdna_ratio'] == mtDNA_ratio
        ]
        data = pd.concat([df_WT_ratio_one, df_growth_rate_ratio])
        
        # print('-'*100)
        # print('Groups statistics:')
        # print(data.groupby(['time', 'growth_rate_ratio']).describe())
        # print('-'*100)

        axes = ax[a]
        sns.boxplot(
            data=data, 
            x='time',
            y='mean_h',
            hue='growth_rate_ratio',
            ax=axes
        )
        axes.set_title(
            f'mtDNA ratio = {mtDNA_ratio}'
        )
        
    fig.suptitle(
        f'Start cell = {start_cell_type}'
    )

plt.show()