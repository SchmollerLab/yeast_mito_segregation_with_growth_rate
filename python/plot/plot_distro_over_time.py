import os

import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

from python.segregation_simulator import (
    other_strain_growth_rate,
    wt_doubling_time,
    startbud,
    ngen,
    ndau,
    nspl,
    number_simulations,
    number_of_cells,
)

from python.segregation_simulator.utils import (
    get_cell_inital_state, 
    get_table_filenames, 
    get_single_cells_filename
)

start_cell_type = '11...00' # '11...00', '1010...'

cwd_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
tables_path = os.path.join(cwd_path, 'segregation_simulator', 'tables_out')

nuc, nuc_format = get_cell_inital_state(start_cell_type)
table_basename, table_filename, single_cells_filename = get_table_filenames(
    nuc_format, number_of_cells
)

ncols = 4
nrows = 4

for growth_rate_ratio in (1.0, other_strain_growth_rate):
    df_cell_filename = get_single_cells_filename(
        single_cells_filename, growth_rate_ratio, 0, 0
    )

    df_filepath = os.path.join(tables_path, df_cell_filename)

    df_grr = pd.read_csv(df_filepath)

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