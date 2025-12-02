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
from yeast_mito_segregation.segregation_simulator import (
    # other_strain_growth_rate,
    wt_doubling_time,
    startbud,
    ngen,
    ndau,
    nspl,
    number_simulations,
    number_of_cells,
    tables_path,
    df_post_growth_mating_filepath,
    df_qpcr_data_filepath,
    df_mtdna_ratio_test_data_filepath,
    force_number_of_runs,
    table_endname
)
from yeast_mito_segregation.segregation_simulator.utils import (
    get_cell_inital_state, 
    get_table_filenames, 
    calc_growth_rate_ratios,
    set_mtdna_amount
)

DF_MTDNA_RATIO_FILEPATH = (
    df_mtdna_ratio_test_data_filepath
    # df_qpcr_data_filepath
    # df_mtdna_ratio_test_data_filepath
)

STRAINS_TO_SIMUL = (
    # 'WT',
    # 'atp6',
    # '(ic)',
    # '(il)',
    # 'cob',
    # 'cox2',
    '∆cox4∆atp6',
    '∆cox4∆cob',
    '∆cox4∆cox2',
    '∆rip1∆atp6',
    '∆rip1∆cob',
    '∆rip1∆cox2'
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

while True:
    answer = input(
        'Do you want to concatenate to existing table? [y/n] '
    )
    if answer.lower() == 'q':
        exit()
        
    if answer.lower() == 'y':
        CONCAT = True
        break
    elif answer.lower() == 'n':
        CONCAT = False
        break
    else:
        print(f'{answer} is not a valid answer. Please enter "y" or "n".')

start_cell_types = ('1010...',) 
# start_cell_types = ('00...11', '11...00', '1010..')

# start_cell_type = '1010...' # '11...00', '1010...'
growth_rate_ratios_mapper = calc_growth_rate_ratios(df_post_growth_mating_filepath)

mtDNA_amounts_df = pd.read_csv(DF_MTDNA_RATIO_FILEPATH, index_col='strain')

growth_rate_ratios_to_simul_mapper = {}
for strain_to_simul in STRAINS_TO_SIMUL:
    if strain_to_simul == 'WT':
        growth_rate_ratio_strain_to_simul_mapper = {'WT and WT': 1.0}
    else:
        growth_rate_ratio_strain_to_simul_mapper = {
            strain: value for strain, value in growth_rate_ratios_mapper.items()
            if strain.endswith(strain_to_simul)
        }
    growth_rate_ratios_to_simul_mapper = {
        **growth_rate_ratios_to_simul_mapper, 
        **growth_rate_ratio_strain_to_simul_mapper
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
        nuc_format, appended_text=table_endname
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
        if strain not in mtDNA_amounts_df.index:
            pbar_strains.update()
            continue 
        
        mtdna_ratios = mtDNA_amounts_df.loc[[strain]]
        pbar_mtdna_ratio = tqdm(
            total=len(mtdna_ratios), 
            ncols=100, 
            leave=False, 
            desc='mtDNA ratio',
            position=2
        )
        for mtdna_ratio_row in mtdna_ratios.itertuples():
            if force_number_of_runs:
                num_cells = number_of_cells
                num_simul = number_simulations
            else:
                try:
                    num_simul = int(mtdna_ratio_row.number_simulations)
                except Exception as err:
                    num_simul = number_simulations
                try:
                    num_cells = int(mtdna_ratio_row.number_of_cells)
                except Exception as err:
                    num_cells = number_of_cells

            mtdna_ratio = mtdna_ratio_row.value
            
            pbar_simul = tqdm(
                total=num_simul, 
                ncols=100, 
                leave=False, 
                desc='Simul ',
                position=3
            )
            for s in range(num_simul):
                pbar_cells = tqdm(
                    total=num_cells, 
                    ncols=100, 
                    leave=False, 
                    desc='Cells ',
                    position=4
                )
                df_cells = {}
                for c in range(num_cells):                
                    nuc_distr = set_mtdna_amount(nuc, mtdna_ratio)
                    start_cell = cell(
                        nuc_distr, 
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
                        .replace('(il)', 'il')
                    )
                    
                    if SAVE:
                        hdf_store_cells[cell_key] = df_cell

                    pbar_cells.update()

                pbar_cells.close()
                df_cells = pd.concat(df_cells, names=['cell_idx']).reset_index()
                
                df_colony = df_cells.groupby('time').agg(mean_h=('h', 'mean'))
                
                import pdb; pdb.set_trace()
                
                df_colony['nspl'] = nspl
                df_colony['ndau'] = ndau
                df_colony['ngen'] = ngen
                df_colony['startbud'] = startbud
                
                dfs[(strain, growth_rate_ratio, mtdna_ratio, s)] = df_colony
                
                pbar_simul.update()
            
            pbar_simul.close()
            pbar_mtdna_ratio.update()
        
        pbar_mtdna_ratio.close()
        pbar_strains.update()

    pbar_strains.close()

    names = ['strain', 'growth_rate_ratio', 'mtdna_ratio', 'simulation_index']
    final_df = pd.concat(dfs, names=names)
    final_table_out_filepath = os.path.join(
        tables_path, table_filename
    )
    
    if CONCAT and os.path.exists(final_table_out_filepath):
        saved_df = pd.read_csv(final_table_out_filepath, index_col=names)
        final_df = pd.concat([saved_df, final_df])

    if SAVE:
        final_df.to_csv(final_table_out_filepath)

    if SAVE:
        hdf_store_cells.close()
    pbar_cell_type.update()

pbar_cell_type.close()
    