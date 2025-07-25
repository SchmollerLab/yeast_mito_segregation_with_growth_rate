import os

cwd_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
tables_path = os.path.join(cwd_path, 'segregation_simulator', 'tables_out')
experimental_data_path = os.path.join(cwd_path, 'experimental_data')
df_post_growth_mating_filepath = os.path.join(
    experimental_data_path, 'PostGrowthMatingAssayWT.csv'
)

other_strain_growth_rate = 0.9314
wt_doubling_time = 1.48369631
startbud = 32
ngen = 14
ndau = 11
nspl = 5
number_simulations = 100
number_of_cells = 30

growth_rate_wt_atp6_neongreen = 0.46718