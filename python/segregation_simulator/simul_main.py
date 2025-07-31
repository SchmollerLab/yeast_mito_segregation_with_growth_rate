import os

from tqdm import tqdm

import numpy as np
import pandas as pd

from yeast_mito_sim import (
    cell, 
    family_simulator, 
    print_family_table, 
    get_family_table
)
from python.segregation_simulator import (
    # other_strain_growth_rate,
    wt_doubling_time,
    startbud,
    ngen,
    ndau,
    nspl,
    number_simulations,
    number_of_cells,
    tables_path,
    df_post_growth_mating_filepath
)
from python.segregation_simulator.utils import (
    get_cell_inital_state, 
    get_table_filenames, 
    calc_growth_rate_ratios
)

STRAINS_TO_SIMUL = (
    'atp6'
)
while True:
    answer = input(
        'Do you want to save the data? [y/n] '
    )
    if answer.lower() == 'q':
        exit()
        
    if answer.lower() == 'y':
        SAVE = True
        break
    elif answer.lower() == 'n':
        SAVE = False
        break
    else:
        print(f'{answer} is not a valid answer. Please enter "y" or "n".')

start_cell_types = ('00...11',) 
# start_cell_types = ('00...11', '11...00', '1010..')

# start_cell_type = '1010...' # '11...00', '1010...'
growth_rate_ratios_mapper = calc_growth_rate_ratios(df_post_growth_mating_filepath)

growth_rate_ratios_to_simul_mapper = {}
for strain_to_simul in STRAINS_TO_SIMUL:
    growth_rate_ratio_strain_to_simul_mapper = {
        strain: value for strain, value in growth_rate_ratios_mapper.items()
        if strain.endswith(strain_to_simul)
    }
    growth_rate_ratios_to_simul_mapper = {
        **growth_rate_ratios_to_simul_mapper, 
        **growth_rate_ratio_strain_to_simul_mapper
    }

growth_rate_ratios_to_simul_mapper = {
    'WT and WT': 1.0,
    **growth_rate_ratios_to_simul_mapper
}

pbar_cell_type = tqdm(
    total=len(start_cell_types), 
    ncols=100, 
    desc='Start cell', 
    position=0
)
for start_cell_type in start_cell_types:
    nuc, nuc_format = get_cell_inital_state(start_cell_type)
    table_basename, table_filename, single_cells_filename = get_table_filenames(
        nuc_format, number_of_cells, number_simulations=number_simulations
    )

    example_cell_filepath = os.path.join(tables_path, single_cells_filename)
    
    if SAVE:
        hdf_store_cells = pd.HDFStore(example_cell_filepath, mode='w')
    
    dfs = {}

    pbar_strains = tqdm(
        total=len(growth_rate_ratios_to_simul_mapper), 
        ncols=100, 
        leave=False, 
        desc='Strain',
        position=1
    )
    for strain, growth_rate_ratio in growth_rate_ratios_to_simul_mapper.items():
        pbar_simul = tqdm(
            total=number_simulations, 
            ncols=100, 
            leave=False, 
            desc='Simul ',
            position=2
        )
        for s in range(number_simulations):
            pbar_cells = tqdm(
                total=number_of_cells, 
                ncols=100, 
                leave=False, 
                desc='Cells ',
                position=3
            )
            df_cells = {}
            for c in range(number_of_cells):
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
                    other_strain_growth_rate=growth_rate_ratio,
                    start_cell_id=0
                )

                df_cell = get_family_table(
                    f, g, doubling_time=wt_doubling_time, 
                    index_col='time'
                )
                
                df_cells[c] = df_cell
                
                cell_key = (
                    f's{s}_c{c}_{strain}'
                    .replace(' ', '_')
                    .replace('∆', 'D')
                    .replace('(ic)', 'ic')
                )
                
                if SAVE:
                    hdf_store_cells[cell_key] = df_cell

                pbar_cells.update()

            pbar_cells.close()
            df_cells = pd.concat(df_cells, names=['cell_idx']).reset_index()
            
            df_colony = df_cells.groupby('time').agg(mean_h=('h', 'mean'))
            
            df_colony['nspl'] = nspl
            df_colony['ndau'] = ndau
            df_colony['ngen'] = ngen
            df_colony['startbud'] = startbud
            df_colony['growth_rate_ratio'] = growth_rate_ratio
            
            dfs[(growth_rate_ratio, s)] = df_colony
            
            pbar_simul.update()
        
        pbar_simul.close()
        pbar_strains.update()

    pbar_strains.close()

    final_df = pd.concat(dfs, names=['growth_rate_ratio', 'simulation_index'])
    final_table_out_filepath = os.path.join(
        tables_path, table_filename
    )
    if SAVE:
        final_df.to_csv(final_table_out_filepath)

    if SAVE:
        hdf_store_cells.close()
    pbar_cell_type.update()

pbar_cell_type.close()
    