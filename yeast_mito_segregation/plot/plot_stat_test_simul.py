import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import beta, combine_pvalues

from yeast_mito_segregation.segregation_simulator import (
    table_endname, tables_path, df_post_growth_mating_filepath
)
    

from yeast_mito_segregation.segregation_simulator.utils import (
    get_cell_inital_state, 
    get_table_filenames, 
    get_single_cells_filename
)

def monte_carlo_one_sided_p(obs, sims, direction='less', add_one_correction=True):
    """
    Compute one-sided Monte Carlo p-value comparing observed value to simulated null distribution.

    Args:
        obs : float
            Observed statistic (percent dark, e.g. 0.42 for 42%).
        sims : array-like
            Simulated statistics under null (length N).
        direction : 'less' or 'greater'
            'less' => test obs is smaller than null (left-tailed).
            'greater' => test obs is larger than null (right-tailed).
        add_one_correction : bool
            Use (k+1)/(N+1) correction to avoid p=0.

    Returns:
        p : float
            Monte Carlo p-value.
        k : int
            Number of simulated stats as or more extreme than obs (depends on direction).
        N : int
            Number of simulations.
    """
    sims = np.asarray(sims)
    N = sims.size
    if direction == 'less':
        k = np.sum(sims <= obs)
    elif direction == 'greater':
        k = np.sum(sims >= obs)
    else:
        raise ValueError("direction must be 'less' or 'greater'")

    if add_one_correction:
        p = (k + 1) / (N + 1)
    else:
        p = k / N
    return p, k, N

def clopper_pearson_ci(k, N, alpha=0.05):
    """
    Clopper-Pearson (exact binomial) CI for a proportion p = k/N.
    Here used to give CI for the Monte Carlo p estimate (k successes in N sims).
    Returns (lower, upper).
    """
    lower = 0.0 if k == 0 else beta.ppf(alpha/2, k, N - k + 1)
    upper = 1.0 if k == N else beta.ppf(1 - alpha/2, k + 1, N - k)
    return lower, upper

STRAIN = 'WT (ic) and ∆atp6'
mtDNA_RATIO = 0.71
add_correction = True

start_cell_type = '1010...' # '11...00', '1010...', '00...11

nuc, nuc_format = get_cell_inital_state(start_cell_type)
table_basename, table_filename, single_cells_filename = get_table_filenames(
    nuc_format, appended_text=table_endname
)

print(f'Loading table {table_filename}')

df_filepath = os.path.join(tables_path, table_filename)

df_simul = pd.read_csv(df_filepath, index_col='strain')
df_simul_strain = df_simul.loc[STRAIN].reset_index()
df_simul_strain['time'] = df_simul_strain['time'].round(1)
last_time = df_simul_strain['time'].max()
df_simul_strain_last_timepoint = (
    df_simul_strain[df_simul_strain['time'] == last_time]
)
df_simul_strain_last_timepoint = (
    df_simul_strain_last_timepoint
        [df_simul_strain_last_timepoint['mtdna_ratio'] == mtDNA_RATIO]
)

df_exp = pd.read_csv(df_post_growth_mating_filepath, index_col='Strain')

sim_percents = df_simul_strain_last_timepoint['mean_h'].values
experimental_percents = df_exp.loc[STRAIN]['Ratio'].values/100
p_values = []

print(sim_percents.min(), sim_percents.max())

for obs_percent in experimental_percents:
    # Compute one-sided p-value for 'less' (obs smaller than simulated -> selection against)
    p, k, N = monte_carlo_one_sided_p(
        obs=obs_percent, 
        sims=sim_percents, 
        direction='less', 
        add_one_correction=add_correction
    )
    # note: use k+1,N+1 if correction used
    ci_lower, ci_upper = clopper_pearson_ci(
        k+1 if add_correction else k, 
        N+1 if add_correction else N,
    )  

    print(f"Observed percent dark: {obs_percent:.3f}")
    print(f"Simulations: N = {N}, #sim <= obs = {k}")
    print(f"One-sided Monte Carlo p-value (left-tailed, +1 correction): {p:.4f}")
    print(
        f"Approx 95% CI for p (Clopper-Pearson): "
        f"[{ci_lower:.4f}, {ci_upper:.4f}]"
    )
    
    p_values.append(p)

    # Visualization
    plt.hist(sim_percents, bins=25, alpha=0.8)
    plt.axvline(obs_percent, color='red', linestyle='--', linewidth=2, label=f'Observed = {obs_percent:.3f}')
    plt.xlabel('Percent dark (simulated null)')
    plt.ylabel('Count')
    plt.title('Simulated null distribution vs observed')
    plt.legend()
    plt.show()

fisher_result = combine_pvalues(p_values)
fisher_pval = fisher_result.pvalue

print(f'Fisher p-value = {fisher_pval:.3f}')