import os

import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

from python.segregation_simulator import (
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

simul_index = 9
cell_index = 1

CELL_KEYS = (
    f's{simul_index}_c{cell_index}_WT_and_WT', 
    f's{simul_index}_c{cell_index}_WT_ic_and_Datp6'
)

start_cell_type = '00...11' # '11...00', '1010...', '00...11

cwd_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
tables_path = os.path.join(cwd_path, 'segregation_simulator', 'tables_out')

nuc, nuc_format = get_cell_inital_state(start_cell_type)
table_basename, table_filename, single_cells_filename = get_table_filenames(
    nuc_format, number_of_cells, number_simulations=number_simulations
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