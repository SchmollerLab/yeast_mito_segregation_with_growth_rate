import pandas as pd
import sklearn
import numpy as np
from sklearn.preprocessing import minmax_scale
from heteroplassimo.cell_data import CellState



# TODO: normalize
def add_signal_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add signal columns (boolean) for mKate and NG."""

    ng_signal = []
    mKate_signal = []
    cell_state = []

    # create new column to threshold upon:
    # (visit every cell at every frame and divide the normalised_to_bkgr protein amount with its volume to get the concentration)

    # df["mKate_conc_auto_vol_pxl"] =  df["mKate_amount_autoBkgr"]  / df["cell_area_pxl"]
    # df["NG_conc_auto_vol_pxl"] =  df["NG_amount_autoBkgr"]  / df["cell_area_pxl"]

    ### normalize the NG_conc and mKate_conc to the zygote values per position!
    ### the min of zygotes should be the 0 (smaller -, higher than the min +)

    ### Apply the signal threshold per position!
    ### Threshold = the minimum of the values from the first two frames of the zygotes (per position), i.e. 90 min after timelapse starts.

    # normalize the signal concentration columns to be between 0 and 1.
    df[
        [
            "NG_concentration_autoBkgr_from_vol_fl",
            "mKate_concentration_autoBkgr_from_vol_fl",
        ]
    ] = minmax_scale(
        df[
            [
                "NG_concentration_autoBkgr_from_vol_fl",
                "mKate_concentration_autoBkgr_from_vol_fl",
            ]
        ]
    )

    for _, row in df.iterrows():

        zygotes = df[
            (df.frame_i <= 1)
            & (df.relationship == "mother")
            & (df.experiment_foldername == row.experiment_foldername)
            & (df.Position_n == row.Position_n)
        ]

        zygote_ng_signal_in_first2frames = (
            zygotes.NG_concentration_autoBkgr_from_vol_fl.min()
        )
        zygote_mKate_signal_in_first2frames = (
            zygotes.mKate_concentration_autoBkgr_from_vol_fl.min()
        )

        #  if row.frame_i <= 20:
        #     df_ = df.loc[(df.Cell_ID.isin(list(zygotes.Cell_ID.values))) & (df.frame_i == row.frame_i)] #change this two only the first 2 frames
        #   else:
        #      df_ = df.loc[(df.Cell_ID.isin(list(zygotes.Cell_ID.values))) & (df.frame_i.isin([19, 20]))]

        ng_sig = (
            row.NG_concentration_autoBkgr_from_vol_fl
            >= zygote_ng_signal_in_first2frames
        )  # - zygote_ng_signal_in_same_frame_std / 2
        ng_signal.append(ng_sig)

        mKate_sig = (
            row.mKate_concentration_autoBkgr_from_vol_fl
            >= zygote_mKate_signal_in_first2frames
        )  # - zygote_mKate_signal_in_same_frame_std / 2
        mKate_signal.append(mKate_sig)

        if ng_sig == True and mKate_sig == True:
            cell_state.append(CellState.HETEROPLASMIC)
        elif ng_sig == True and mKate_sig == False:
            cell_state.append(CellState.HOMOPLASMIC_NG)
        elif mKate_sig == True and ng_sig == False:
            cell_state.append(CellState.HOMOPLASMIC_MKATE)
        else:
            cell_state.append(CellState.NO_SIGNAL)
            # if no signal -> get signal from prior cell in the lineage

            # print(zygote_signal, ng_signal_std, row.NG_concentration_autoBkgr_from_vol_vox, row.NG_concentration_autoBkgr_from_vol_vox > zygote_signal - ng_signal_std)

    df["ng_signal"] = ng_signal
    df["mKate_signal"] = mKate_signal
    df["cell_state"] = cell_state

    return df

def get_eucl_distance(df: pd.DataFrame) -> pd.DataFrame:

      # MAYBE APPLY THE Eucl Dist ON THE log(signal) columns  
      # normalize the signal concentration columns to be between 0 and 1.
    df[
        [
            "NG_concentration_autoBkgr_from_vol_fl",
            "mKate_concentration_autoBkgr_from_vol_fl",
        ]
    ] = minmax_scale(
        df[
            [
                "NG_concentration_autoBkgr_from_vol_fl",
                "mKate_concentration_autoBkgr_from_vol_fl",
            ]
        ]
    )
    

    df = pd.merge(df, df, left_on = ['frame_i','Cell_ID'], right_on =['frame_i','relative_ID'], how = 'left', suffixes =[ '', '_rel'],) \
                    .drop(['relative_ID_rel', 'relationship_rel'], axis=1)

    g = (df.NG_concentration_autoBkgr_from_vol_fl - df.NG_concentration_autoBkgr_from_vol_fl_rel)
    r = (df.mKate_concentration_autoBkgr_from_vol_fl - df.mKate_concentration_autoBkgr_from_vol_fl_rel)

    dist_per_frame = np.sqrt((np.square(g) + np.square(r)))
        
    df['eucl_dist'] = dist_per_frame

    return df

def frames_to_time(df: pd.DataFrame) -> pd.DataFrame:

    time = (15 * df["frame_i"]).div(60)

    df["time"] = time

    return df 

def get_signal_ratio(df: pd.DataFrame) -> pd.DataFrame: 

    #for _, row in df.iterrows():

    cell_signal_ratio = (  
        (df.mKate_concentration_autoBkgr_from_vol_fl) # adding +1 to each float so there is no ZeroDivisionError 
        / (df.NG_concentration_autoBkgr_from_vol_fl)
     ) 

    log_ratio = np.log((df.mKate_concentration_autoBkgr_from_vol_fl)/(df.NG_concentration_autoBkgr_from_vol_fl))
                           
    df["cell_signal_ratio"] = cell_signal_ratio
    df['log_ratio'] = log_ratio

    return df 

def get_log_signal_and_ratio(df: pd.DataFrame) -> pd.DataFrame:

    log_mk = np.log(df['mKate_concentration_autoBkgr_from_vol_fl'])
    log_ng = np.log(df['NG_concentration_autoBkgr_from_vol_fl'])
    ratio_log =  ((np.log(df['mKate_concentration_autoBkgr_from_vol_fl'])) / (np.log(df['NG_concentration_autoBkgr_from_vol_fl']))) 


    df["log_mk"] = log_mk
    df["log_ng"] = log_ng
    df["ratio_log"] = ratio_log

    return df 
