import os

from tqdm import tqdm

import numpy as np
import pandas as pd

from yeast_mito_sim import (
    cell, 
    family_simulator, 
    print_family_table
)

cwd_path = os.path.dirname(os.path.abspath(__file__))
tables_path = os.path.join(cwd_path, 'tables_out')

number_simulations = 100

nuc = np.zeros(32, dtype=int) # 1, 0, 1, 0, etc. 
nuc[::2] = 1
nuc = tuple(nuc)
nuc_format = ''.join([str(n) for n in nuc])

other_strain_growth_rate = 0.9314
wt_doubling_time = 1.48369631
startbud = 32
ngen = 14
ndau = 14
nspl = 4

table_basename = f'start_cell_{nuc_format}.csv'

dfs = {}

for growth_rate_ratio in (1.0, other_strain_growth_rate):
    pbar = tqdm(total=number_simulations, ncols=100)
    for n in range(number_simulations):
        start_cell = cell(
            nuc, 
            -1, 
            startbud=startbud, 
            nspl=nspl, 
            ndau=ndau, 
            maxnaddexp=0
        )
        f, g = family_simulator(
            start_cell, ngen, 
            other_strain_growth_rate=growth_rate_ratio
        )
        
        table_filename = f'{n}_{table_basename}'
        table_out_filepath = os.path.join(
            tables_path, table_filename
        )
        
        print_family_table(
            f, g, table_out_filepath, add=False, 
            doubling_time=wt_doubling_time
        )
        
        df = pd.read_csv(table_out_filepath, index_col='time')

        df['nspl'] = nspl
        df['ndau'] = ndau
        df['ngen'] = ngen
        df['startbud'] = startbud
        df['growth_rate_ratio'] = growth_rate_ratio

        dfs[(growth_rate_ratio, n)] = df
        
        os.remove(table_out_filepath)
        
        pbar.update()

    pbar.close()
    
final_df = pd.concat(dfs, names=['growth_rate_ratio', 'simulation_index'])
final_table_out_filepath = os.path.join(tables_path, table_basename)
final_df.to_csv(final_table_out_filepath)
    
    