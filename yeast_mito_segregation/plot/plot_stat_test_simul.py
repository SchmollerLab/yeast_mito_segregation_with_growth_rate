import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import beta, combine_pvalues

from collections import defaultdict

from yeast_mito_segregation.segregation_simulator import (
    table_endname, tables_path, df_post_growth_mating_filepath,
    df_fast_filepath, df_qpcr_data_filepath
)
    

from yeast_mito_segregation.segregation_simulator.utils import (
    get_cell_inital_state, 
    get_table_filenames, 
    get_single_cells_filename
)

def monte_carlo_p_value(obs, sims, one_sided=True, direction='less', add_one_correction=True):
    """
    Compute Monte Carlo p-value comparing observed value to simulated null distribution.

    Args:
        obs : float
            Observed statistic (percent dark, e.g. 0.42 for 42%).
        one_sided : bool
            If True, compute one-sided p-value.
            If False, compute two-sided p-value.
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
    if not one_sided:
        k = np.sum(np.abs(sims - np.mean(sims)) >= np.abs(obs - np.mean(sims)))
    elif direction == 'less':
        k = np.sum(sims <= obs)
    elif direction == 'greater':
        k = np.sum(sims >= obs)
    else:
        raise ValueError("direction must be 'less' or 'greater' or one_sided=False")

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

df_qcr = pd.read_csv(df_qpcr_data_filepath, index_col='strain')

SHOW_PLOT = True
add_correction = True

start_cell_type = '1010...' # '11...00', '1010...', '00...11

nuc, nuc_format = get_cell_inital_state(start_cell_type)
table_basename, table_filename, single_cells_filename = get_table_filenames(
    nuc_format, appended_text=table_endname
)

print(f'Loading table {table_filename}')

df_filepath = os.path.join(tables_path, table_filename)

df_simul = pd.read_csv(df_filepath)
df_pvalues_data = defaultdict(list)
for strain, df_simul_strain in df_simul.groupby('strain'):
    df_simul_strain = df_simul_strain.reset_index()
    df_simul_strain['time'] = df_simul_strain['time'].round(1)
    last_time = df_simul_strain['time'].max()
    df_simul_strain_last_timepoint = (
        df_simul_strain[df_simul_strain['time'] == last_time]
    )
    mtDNA_RATIO = df_qcr.at[strain, 'value']
    
    df_simul_strain_last_timepoint = (
        df_simul_strain_last_timepoint
            [df_simul_strain_last_timepoint['mtdna_ratio'] == mtDNA_RATIO]
    )

    df_exp = pd.read_csv(df_fast_filepath, index_col='Mating')

    sim_percents = df_simul_strain_last_timepoint['mean_h'].values
    
    if len(sim_percents) == 0:
        import pdb; pdb.set_trace()

    strain_x = strain.replace('and', 'x')
    if strain_x not in df_exp.index:
        continue

    df_exp_strain = df_exp.loc[strain.replace('and', 'x')]
    p_values_ts = []
    p_values_os = []

    # print(sim_percents.min(), sim_percents.max())

    # experimental_percents = [experimental_percents.mean()]
    
    fig, ax = plt.subplots(1, len(df_exp_strain), figsize=(15, 6))
    ax = ax.flatten()
    
    fig.suptitle(f'Simulated null distribution vs observed - {strain}')
    
    print('*'*100)
    print(f'STRAIN: {strain}\n')
    for a, row in enumerate(df_exp_strain.itertuples()):
        
        axes: plt.Axes = ax[a]
        
        obs_percent = row.Ratio/100
        
        # Compute one-sided p-value for 'less' (obs smaller than simulated -> selection against)
        p_os, k_os, N = monte_carlo_p_value(
            obs=obs_percent, 
            sims=sim_percents, 
            one_sided=True,
            direction='less', 
            add_one_correction=add_correction
        )
        # note: use k+1,N+1 if correction used
        ci_lower_os, ci_upper_os = clopper_pearson_ci(
            k_os+1 if add_correction else k_os, 
            N+1 if add_correction else N,
        )  
        
        # Compute two-sided p-value
        p_ts, k_ts, N = monte_carlo_p_value(
            obs=obs_percent, 
            sims=sim_percents, 
            one_sided=False,
            add_one_correction=add_correction
        )
        # note: use k+1,N+1 if correction used
        ci_lower_ts, ci_upper_ts = clopper_pearson_ci(
            k_ts+1 if add_correction else k_ts, 
            N+1 if add_correction else N,
        )

        print(f"Observed percent dark: {obs_percent:.3f}")
        print(f"Simulations: N = {N}, #sim <= obs = {k_os}")
        print(f"One-sided Monte Carlo p-value (left-tailed, +1 correction): {p_os:.4f}")
        print(
            f"Approx 95% CI for p (Clopper-Pearson): "
            f"[{ci_lower_os:.4f}, {ci_upper_os:.4f}]"
        )
        print(f"Two-sided Monte Carlo p-value (left-tailed, +1 correction): {p_ts:.4f}")
        print(
            f"Approx 95% CI for p (Clopper-Pearson): "
            f"[{ci_lower_ts:.4f}, {ci_upper_ts:.4f}]"
        )
        print('-'*100)

        p_values_ts.append(p_ts)
        p_values_os.append(p_os)
        df_pvalues_data['strain'].append(strain)
        df_pvalues_data['Replicate'].append(row.Replicate)
        df_pvalues_data['two_sided_monte_carlo_p_value'].append(p_ts)
        df_pvalues_data['one_sided_monte_carlo_p_value'].append(p_os)

        # Visualization
        axes.hist(sim_percents, bins=25, alpha=0.8)
        axes.axvline(obs_percent, color='red', linestyle='--', linewidth=2, label=f'Observed = {obs_percent:.3f}')
        axes.set_xlabel('Percent dark (simulated null)')
        axes.set_ylabel('Count')
        axes.set_title(f'Replicate {row.Replicate}')
        axes.legend()
    
    fisher_result_os = combine_pvalues(p_values_os)
    fisher_pval_os = fisher_result_os.pvalue

    df_pvalues_data['combined_one_sided_p_value_fisher'].extend([fisher_pval_os]*len(p_values_os))
    
    fisher_result_ts = combine_pvalues(p_values_ts)
    fisher_pval_ts = fisher_result_ts.pvalue
    df_pvalues_data['combined_two_sided_p_value_fisher'].extend([fisher_pval_ts]*len(p_values_ts))
    
    print('='*100)

    print(f'Fisher one-sided p-value = {fisher_pval_os:.3f}')
    print(f'Fisher two-sided p-value = {fisher_pval_ts:.3f}')

    print('*'*100)
    
    if SHOW_PLOT:
        plt.show()
        import pdb; pdb.set_trace()
    else:
        plt.close()

    

df_pvalues = pd.DataFrame(df_pvalues_data)

print(df_pvalues)

pvalues_tablename = table_filename.replace(
    '.csv', '_monte_carlo_statistical_analysis.csv'
)
pvalues_table_filepath = os.path.join(tables_path, pvalues_tablename)

answer = input(f'Save p-values table to {pvalues_table_filepath}? (y/n): ')
if answer.lower() != 'y':
    print('Saving cancelled')
    exit(0)
    
df_pvalues.to_csv(pvalues_table_filepath, index=False)

print(f'Table saved at {pvalues_table_filepath}')