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
 
CELL_KEYS = ('s0_c0_WT_and_WT', 's0_c0_WT_ic_and_Datp6')

start_cell_type = '1010...' # '11...00', '1010...'

cwd_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
tables_path = os.path.join(cwd_path, 'segregation_simulator', 'tables_out')

nuc, nuc_format = get_cell_inital_state(start_cell_type)
table_basename, table_filename, single_cells_filename = get_table_filenames(
    nuc_format, number_of_cells
)

single_cells_filepath = os.path.join(tables_path, single_cells_filename)

ncols = 4
nrows = 4

for cell_key in CELL_KEYS:
    df_cell = pd.read_hdf(single_cells_filepath, key=cell_key)

    fig, ax = plt.subplots(nrows, ncols, figsize=(12, 10))
    ax = ax.flatten()
    for i, (time, df_time) in enumerate(df_cell.groupby(level=0)):
        axes = ax[i]
        sns.histplot(
            data=df_time, 
            ax=axes,
            x='h',
            stat='probability',
            # hue='growth_rate_ratio'
        )
    fig.suptitle(f'Cell key = {cell_key}, start cell = {start_cell_type}')
    
plt.show()