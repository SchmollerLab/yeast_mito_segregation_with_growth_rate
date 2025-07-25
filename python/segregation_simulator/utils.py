from typing import Literal

import numpy as np
import pandas as pd

from python.segregation_simulator import (
    growth_rate_wt_atp6_neongreen
)

def get_cell_inital_state(start_cell_type: Literal['1010...', '11...00']):
    if start_cell_type == '1010...':
        nuc = np.zeros(32, dtype=int) # 1, 0, 1, 0, etc. 
        nuc[::2] = 1
    elif start_cell_type == '11...00':
        nuc = np.ones(32, dtype=int) # 1, 1, etc., 0, 0 
        nuc[16:] = 0
    
    nuc_format = ''.join([str(n) for n in nuc])
    return tuple(nuc), nuc_format

def get_table_filenames(nuc_format, number_of_cells):
    table_basename = f'start_cell_{nuc_format}'
    single_cells_filename = f'{table_basename}_single_cell_data.h5'
    table_filename = f'{table_basename}_num_cells_per_colony_{number_of_cells}.csv'
    return table_basename, table_filename, single_cells_filename

def get_single_cells_filename(single_cells_filename, growth_rate_ratio, s, c):
    filename = f'{single_cells_filename}_s{s}_c{c}_gr{growth_rate_ratio}.csv'
    return filename

def calc_growth_rate_ratios(df_post_growth_mating_filepath):
    df_pgm = pd.read_csv(df_post_growth_mating_filepath)
    df_pgm['WT_growth_rate_hours'] = growth_rate_wt_atp6_neongreen
    
    hours_exp = 20
    
    df_pgm['growth_rate_hours'] = (
        (np.log(df_pgm['Ratio']/(100-df_pgm['Ratio'])) 
        + hours_exp*df_pgm['WT_growth_rate_hours'])
        / hours_exp
    )
    
    df_pgm['growth_rate_ratio'] = (
        df_pgm['growth_rate_hours'] / df_pgm['WT_growth_rate_hours']
    )
    
    growth_rate_ratios_mean = (
        df_pgm.groupby('Strain')['growth_rate_ratio'].mean().to_dict()
    )

    return growth_rate_ratios_mean