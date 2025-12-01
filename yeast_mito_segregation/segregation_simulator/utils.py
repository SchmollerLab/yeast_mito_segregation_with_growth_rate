from typing import Literal

import numpy as np
import pandas as pd

from yeast_mito_segregation.segregation_simulator import (
    growth_rate_wt_atp6_neongreen
)

def get_cell_inital_state(
        start_cell_type: Literal['1010...', '11...00', '00...11'],
        len_nuc: int = 32
    ):
    if start_cell_type == '1010...':
        nuc = np.zeros(len_nuc, dtype=int) 
        nuc[::2] = 1 # 1, 0, 1, 0, etc. 
    elif start_cell_type == '11...00':
        nuc = np.ones(len_nuc, dtype=int) 
        nuc[16:] = 0 # 1, 1, etc., 0, 0 
    elif start_cell_type == '00...11':
        nuc = np.ones(len_nuc, dtype=int) 
        nuc[:16] = 0 # 0, 0, etc., 1, 1
    
    nuc_format = start_cell_type.replace('.', '-')
    nuc_format = f'{nuc_format}_len{len_nuc}'
    return tuple([int(val) for val in nuc]), nuc_format

def get_table_filenames(nuc_format, appended_text=''):
    table_basename = f'start_cell_{nuc_format}'
    single_cells_filename = (
        f'{table_basename}_{appended_text}.h5'
    )
    table_filename = single_cells_filename.replace('.h5', '.csv')
    return table_basename, table_filename, single_cells_filename

def get_single_cells_filename(single_cells_filename, growth_rate_ratio, s, c):
    filename = f'{single_cells_filename}_s{s}_c{c}_gr{growth_rate_ratio}.csv'
    return filename

def calc_growth_rate_ratios(df_post_growth_mating_filepath):
    df_pgm = pd.read_csv(df_post_growth_mating_filepath)
    # df_pgm['WT_growth_rate_hours'] = growth_rate_wt_atp6_neongreen
    
    hours_exp = 20
    
    # For 'WT_growth_rate' -->
    # see "yeast_mito_segregation/experimental_data/add_WT_growth_rate_hours.py"
    df_pgm['growth_rate_hours'] = (
        (
            np.log(df_pgm['Ratio']/(100-df_pgm['Ratio'])) 
            + hours_exp*df_pgm['WT_growth_rate']
        )
        / hours_exp
    )
    
    df_pgm['growth_rate_ratio'] = (
        df_pgm['growth_rate_hours'] / df_pgm['WT_growth_rate']
    )

    growth_rate_ratios_mean = (
        df_pgm.groupby('Strain')['growth_rate_ratio'].mean().to_dict()
    )

    return growth_rate_ratios_mean

def set_mtdna_amount(nuc, mtdna_ratio):
    # Here we adjust the number of 1s in `nuc` to simulate different amounts 
    # of mtDNA in the strains (see experimental_data/qPCR.png) 
    # where WT is 0s.
    # For example, if mtdna_ratio in strain is 1.3, we need more 1s to 
    # make the distribution of 1s = 1.3 higher than 0s
    # e.g., number of 0s = 32/2.3 = 13.91 = 14 --> 18 1s
    # The number of 0s to change to 1 is selected randomly
    num_zeros_required = int(round(len(nuc) / (mtdna_ratio + 1)))
    current_num_zeros = np.sum(np.array(nuc) == 0)
    
    if current_num_zeros == num_zeros_required:
        return nuc
    
    nuc = np.array(nuc)
    num_ones_required = len(nuc) - num_zeros_required
    current_num_ones = len(nuc) - current_num_zeros
    
    if num_ones_required > num_zeros_required:
        # need to change some 0s to 1s
        num_zeros_to_change = current_num_zeros - num_zeros_required
        zeros_indices = np.where(nuc == 0)[0]
        indices_to_change = np.random.choice(
            zeros_indices, size=num_zeros_to_change, replace=False
        )
        nuc[indices_to_change] = 1        
    else:
        # need to change some 1s to 0s
        num_ones_to_change = current_num_ones - num_ones_required
        ones_indices = np.where(nuc == 1)[0]
        try:
            indices_to_change = np.random.choice(
                ones_indices, size=num_ones_to_change, replace=False
            )
        except ValueError:
            import pdb; pdb.set_trace()
            
        nuc[indices_to_change] = 0
    
    return tuple(nuc)
