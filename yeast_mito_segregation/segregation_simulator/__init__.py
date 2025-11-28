import os

import numpy as np

cwd_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
tables_path = os.path.join(cwd_path, 'segregation_simulator', 'tables_out')
experimental_data_path = os.path.join(cwd_path, 'experimental_data')
df_post_growth_mating_filepath = os.path.join(
    experimental_data_path, 'PostGrowthMatingAssayWT.csv'
)
df_fast_filepath = os.path.join(
    experimental_data_path, 'FAST_WT.csv'
)
df_qpcr_data_filepath = os.path.join( experimental_data_path, 'qPCR_values.csv')
df_mtdna_ratio_test_data_filepath  = os.path.join(
    experimental_data_path, 'mtDNA_ratio_test_data.csv'
)

os.makedirs(tables_path, exist_ok=True)

# other_strain_growth_rate = 0.9314
wt_doubling_time = 1.48369631 # Measured, see manuscript (search "1.48")
startbud = 32
ngen = 14
ndau = 11
nspl = 5
number_simulations = 2
number_of_cells = 5

# If False, `number_simulations` and `number_of_cells` are ignored and the
# values from the experimental data file are used 
# (see "yeast_mito_segregation\experimental_data\mtDNA_ratio_test_data.csv").
force_number_of_runs = False

# Append this text to the end of the table filenames
table_endname = 'final'

growth_rate_wt_atp6_neongreen = np.log(2)/wt_doubling_time # 0.46718