"""Load csv files into one pandas dataframe."""
import pandas as pd

IMPORTANT_COLUMNS = [
    "frame_i",
    "Cell_ID",
    "relationship",
    "relative_ID",
    "generation_num",
    "emerg_frame_i",
    "division_frame_i",
    "gui_NG_amount_autoBkgr_maxProj",
    "gui_mKate_amount_autoBkgr_maxProj",
]


def load_dataframe(neon_green_csv_file: str, m_kate_csv_file: str) -> pd.DataFrame:
    """Load both csv files in one pandas dataframe.

    Only the important columns are kept.

    Args:
        neon_green_csv_file (str): path to the neon green csv file
        m_kate_csv_file (str): path to the mKate csv file

    Returns:
        pd.DataFrame: combined dataframe
    """

    neon_green_df = pd.read_csv(neon_green_csv_file)
    mkate_df = pd.read_csv(m_kate_csv_file)

    df = pd.concat([neon_green_df, mkate_df], axis=1).T.drop_duplicates().T

    # TODO: maybe make sure that dataframe is sorted by frame_i and cell_id,
    # such that the adding of cells to the tree works correctly every time

    return df[IMPORTANT_COLUMNS]
