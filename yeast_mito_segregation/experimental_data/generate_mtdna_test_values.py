import numpy as np
import pandas as pd

from yeast_mito_segregation.segregation_simulator import (
    df_mtdna_ratio_test_data_filepath
)

data = {
    'WT (ic) and ∆atp6': np.linspace(1.0, 0.1, 10),
    'WT (ic) and ∆cob': np.linspace(1.0, 0.1, 10),
    'WT (ic) and ∆cox2': np.linspace(1.0, 0.1, 10),
}

df = (
    pd.concat(
        {k: pd.DataFrame(v, columns=['value']) for k, v in data.items()},
        names=['strain'] 
    )
    .reset_index()
    .drop(columns='level_1')
    .set_index('strain')
)

df.to_csv(df_mtdna_ratio_test_data_filepath)