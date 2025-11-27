import pandas as pd
import matplotlib as plt
import matplotlib.pyplot as plt
from heteroplassimo.utils import CellState
from heteroplassimo.utils import add_signal_columns
from heteroplassimo.utils import get_eucl_distance
from heteroplassimo.utils import get_signal_ratio
from heteroplassimo.utils import get_log_signal_and_ratio
from heteroplassimo.utils import frames_to_time
import os


def combine_all_dataframes():

    files = os.listdir("./csv_files")

    df = pd.read_csv(f"./csv_files/{files[0]}")

    for file_name in files[1:]:
        if file_name.startswith("."):
            continue
        df = df.append(pd.read_csv(f"./csv_files/{file_name}"))

    df.to_csv(f"csv_files/all_stats.csv")

    stats_df = pd.DataFrame(
        columns=["depth", "hetero", "no_signal", "homo_ng", "homo_mkate"]
    )

    for depth in df.depth.unique():

        depth_df = df[df.depth == depth].drop_duplicates(
            subset=["exp", "pos", "tree", "cell_id"]
        )

        cell_state_percentage = (
            depth_df.cell_state.value_counts() / len(depth_df)
        ).to_frame()

        stats_df = stats_df.append(
            {
                "depth": depth,
                "hetero": cell_state_percentage[
                    cell_state_percentage.index == "CellState.HETEROPLASMIC"
                ].cell_state.values[0],
                "no_signal": cell_state_percentage[
                    cell_state_percentage.index == "CellState.NO_SIGNAL"
                ].cell_state.values[0]
                if len(
                    cell_state_percentage[
                        cell_state_percentage.index == "CellState.NO_SIGNAL"
                    ].cell_state
                )
                > 0
                else 0,
                "homo_ng": cell_state_percentage[
                    cell_state_percentage.index == "CellState.HOMOPLASMIC_NG"
                ].cell_state.values[0]
                if len(
                    cell_state_percentage[
                        cell_state_percentage.index == "CellState.HOMOPLASMIC_NG"
                    ].cell_state
                )
                > 0
                else 0,
                "homo_mkate": cell_state_percentage[
                    cell_state_percentage.index == "CellState.HOMOPLASMIC_MKATE"
                ].cell_state.values[0]
                if len(
                    cell_state_percentage[
                        cell_state_percentage.index == "CellState.HOMOPLASMIC_MKATE"
                    ].cell_state
                )
                > 0
                else 0,
            },
            ignore_index=True,
        )

    plt.plot(stats_df.depth, stats_df.hetero, label="HETEROPLASMIC")
    plt.plot(stats_df.depth, stats_df.no_signal, label="NO_SIGNAL")
    plt.plot(stats_df.depth, stats_df.homo_ng, label="HOMOPLASMIC_NG")
    plt.plot(stats_df.depth, stats_df.homo_mkate, label="HOMOPLASMIC_MKATE")

    plt.legend()

    plt.xlabel("lineage_depth")
    plt.ylabel("Percentage of cells")
    #plt.show()

    plt.savefig('cell_states_per_lineage_depth.svg')
    plt.savefig('cell_states_per_lineage_depth.pdf')
    plt.close()



# combine_all_dataframes()   # un-comment if you want the function to be executed!


def change_of_signal_in_time():


    df = pd.read_csv("multiExp_acdc_output__.csv")

    df = df[df.frame_i <= 32]

    ng_median = df.groupby(["frame_i"]).apply(lambda x: x['NG_concentration_autoBkgr_from_vol_fl'].median()/len(x))
    ng_mean = df.groupby(["frame_i"]).apply(lambda x: x['NG_concentration_autoBkgr_from_vol_fl'].mean()/len(x))

    mkate_median = df.groupby(["frame_i"]).apply(lambda x: x['mKate_concentration_autoBkgr_from_vol_fl'].median()/len(x))
    mkate_mean = df.groupby(["frame_i"]).apply(lambda x: x['mKate_concentration_autoBkgr_from_vol_fl'].mean()/len(x))


    #ng_median = df.groupby(["frame_i"]).apply(lambda x: x['NG_amount_autoBkgr'].median()/len(x))
    #ng_mean = df.groupby(["frame_i"]).apply(lambda x: x['NG_amount_autoBkgr'].mean()/len(x))

    #mkate_median = df.groupby(["frame_i"]).apply(lambda x: x['mKate_amount_autoBkgr'].median()/len(x))
    #mkate_mean = df.groupby(["frame_i"]).apply(lambda x: x['mKate_amount_autoBkgr'].mean()/len(x))

    # ng_median.plot()
    # ng_mean.plot()

    # mkate_median.plot()
    # mkate_mean.plot()

    #plots show the signal of mKate and ng in time, normalised to the amount of cells in each frame.

    #plt.show()


    plt.plot(ng_median, label="ng_median")
    plt.plot(ng_mean, label="ng_mean")
    plt.plot(mkate_median, label="mkate_median")
    plt.plot(mkate_mean, label="mkate_mean")

    plt.legend()

    #plt.savefig('NG_mKate_signal_intensity_per_frame_norm_means_medians.svg')
    #plt.savefig('NG_mKate_signal_intensity_per_frame_norm_means_medians.pdf')
    plt.close()

    # for exp in sorted(df.experiment_foldername.unique()):

    #     exp_df = df[df["experiment_foldername"] == exp]

    #     for pos in exp_df.Position_n.unique():

    #         df_ = exp_df[
    #             (exp_df["experiment_foldername"] == exp) & (exp_df.Position_n == pos)
    #         ]
    #         df_ = add_signal_columns(df_)
    #         df_ = df_.reset_index()

    #         for cell_id in df_.Cell_ID.unique():

    #             cell_df = df_[df_.Cell_ID == cell_id]


    #             plt.plot(cell_df.frame_i, cell_df.NG_concentration_autoBkgr_from_vol_fl, linewidth=0.2, alpha=0.5)


    # plt.show()

change_of_signal_in_time()

def combine_all_stats_v2():

         files = os.listdir("./csv_files")
         
         df_all = pd.DataFrame()
    
         for file_name in files[1:]:
            if not "v2" in file_name:
                continue

            print(df_all)
            df_all = df_all.append(pd.read_csv(f"./csv_files/{file_name}"))
            # df_all = df_all.reset_index()

         df_all.to_csv(f"csv_files/all_stats_v2.csv")
         print(df_all)

combine_all_stats_v2()


def change_of_state_from_parent():

    files = os.listdir("./csv_files")
    df = pd.read_csv(f"./csv_files/{files[0]}")
    
    for file_name in files[1:]:
        df = df.append(pd.read_csv(f"./csv_files/{file_name}"))
 
    perc_ng = []
    perc_mkate = []
    depths = []


    for depth in df.depth.unique():

        depth_df = df[df.depth == depth].drop_duplicates(subset=["exp", "pos", "tree", "cell_id"])

        depths.append(depth)
        perc_mkate.append(len(depth_df[(depth_df.cell_state == "CellState.HOMOPLASMIC_MKATE") & (depth_df.parent_cell_state == "CellState.HETEROPLASMIC")]) / len(depth_df))
        perc_ng.append(len(depth_df[(depth_df.cell_state == "CellState.HOMOPLASMIC_NG") & (depth_df.parent_cell_state == "CellState.HETEROPLASMIC")]) / len(depth_df))

    plt.bar(pd.Series(depths) - 0.2, perc_ng, 0.4, label="change_to_ng")
    plt.bar(pd.Series(depths) + 0.2, perc_mkate, 0.4, label="change_to_mK")


    plt.legend()

    plt.xlabel("Lineage_depth")
    plt.ylabel("Percentage of cell_state different to parent_cell_state")
    #plt.show()

    plt.savefig('change_of_cell_state_from parent_per_lineage_depth.svg')
    plt.savefig('change_of_cell_state_from parent_per_lineage_depth.pdf')
    plt.close()


change_of_state_from_parent()