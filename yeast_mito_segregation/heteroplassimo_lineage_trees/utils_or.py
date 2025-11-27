import pandas as pd
from heteroplassimo.cell_data import CellState


def add_signal_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add signal columns (boolean) for mKate and NG."""

    ng_signal = []
    mKate_signal = []
    cell_state = []

    # create new column to threshold upon: (visit every cell at every frame and divide the normalised_to_bkgr protein amount with its volume to get the concentration)
    df["mKate_conc_auto_vol_pxl"] = (
        df["mKate_amount_autoBkgr"] / df["cell_area_pxl"]
    )
    df["NG_conc_auto_vol_pxl"] = df["NG_amount_autoBkgr"] / df["cell_area_pxl"]

    for _, row in df.iterrows():

        zygotes = df[
            (df.frame_i == 0)
            & (df.relationship == "mother")
            & (df.experiment_foldername == row.experiment_foldername)
            & (df.Position_n == row.Position_n)
        ]

        if row.frame_i <= 20:
            df_ = df.loc[
                (df.Cell_ID.isin(list(zygotes.Cell_ID.values)))
                & (df.frame_i == row.frame_i)
            ]
        else:
            df_ = df.loc[
                (df.Cell_ID.isin(list(zygotes.Cell_ID.values)))
                & (df.frame_i.isin([19, 20]))
            ]

        zygote_ng_signal_in_same_frame = df_.NG_conc_auto_vol_pxl.median()
        zygote_ng_signal_in_same_frame_std = df_.NG_conc_auto_vol_pxl.std()
        zygote_mKate_signal_in_same_frame = df_.mKate_conc_auto_vol_pxl.median()
        zygote_mKate_signal_in_same_frame_std = df_.mKate_conc_auto_vol_pxl.std()
        # print("NORM. NG SIGNAL ZYGOTES:", zygote_ng_signal_in_same_frame, "NORM. NG SIGNAL CELL:", row.NG_conc_auto_vol_pxl)
        # print("NORM. MKATE SIGNAL ZYGOTES:", zygote_mKate_signal_in_same_frame, "NORM. MKATE SIGNAL CELL:", row.mKate_conc_auto_vol_pxl)

        ng_sig = (
            row.NG_conc_auto_vol_pxl > zygote_ng_signal_in_same_frame
        )  # - zygote_ng_signal_in_same_frame_std / 2
        ng_signal.append(ng_sig)

        mKate_sig = (
            row.mKate_conc_auto_vol_pxl > zygote_mKate_signal_in_same_frame
        )  # - zygote_mKate_signal_in_same_frame_std / 2
        mKate_signal.append(mKate_sig)

        if ng_sig and mKate_sig:

            cell_state.append(CellState.HETEROPLASMIC)
        elif ng_sig:

            cell_state.append(CellState.HOMOPLASMIC_NG)
        elif mKate_sig:

            cell_state.append(CellState.HOMOPLASMIC_MKATE)
        else:
            cell_state.append(CellState.NO_SIGNAL)

            # print(zygote_signal, ng_signal_std, row.NG_concentration_autoBkgr_from_vol_vox, row.NG_concentration_autoBkgr_from_vol_vox > zygote_signal - ng_signal_std)

    df["ng_signal"] = ng_signal
    df["mKate_signal"] = mKate_signal
    df["cell_state"] = cell_state

    return df
