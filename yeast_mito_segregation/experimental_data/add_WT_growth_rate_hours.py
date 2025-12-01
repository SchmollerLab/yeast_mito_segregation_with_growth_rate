import numpy as np
import pandas as pd

from yeast_mito_segregation.segregation_simulator import (
    df_post_growth_mating_filepath
)

# Measured, see manuscript (search "1.48")
ATP6_NG_growth_rate_hours = 1.48369631

# Measured, see these tables: 
# - "yeast_mito_segregation/experimental_data/cox4-growthOnPlate.csv"
# - "yeast_mito_segregation/experimental_data/rip1-growthOnPlate.csv"
rip1_del_ATP6_NG_growth_rate_hours = 1.629396
cox4_del_ATP6_NG_growth_rate_hours = 1.632216

df = pd.read_csv(df_post_growth_mating_filepath)
df['WT_doubling_time'] = ATP6_NG_growth_rate_hours
df['WT_growth_rate'] = np.log(2)/ATP6_NG_growth_rate_hours

cox4_del_mask = df['Strain'].str.startswith('∆cox4')
df.loc[cox4_del_mask, 'WT_doubling_time'] = cox4_del_ATP6_NG_growth_rate_hours
df.loc[cox4_del_mask, 'WT_growth_rate'] = (
    np.log(2)/cox4_del_ATP6_NG_growth_rate_hours
)

rip1_del_mask = df['Strain'].str.startswith('∆rip1')
df.loc[rip1_del_mask, 'WT_doubling_time'] = rip1_del_ATP6_NG_growth_rate_hours
df.loc[rip1_del_mask, 'WT_growth_rate'] = (
    np.log(2)/rip1_del_ATP6_NG_growth_rate_hours
)

df.to_csv(df_post_growth_mating_filepath, index=False)