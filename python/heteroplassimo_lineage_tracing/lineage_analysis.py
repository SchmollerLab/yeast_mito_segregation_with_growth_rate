from heteroplassimo.lineages_tree import LineagesTree
from heteroplassimo.file_loader import load_dataframe
from heteroplassimo.utils import add_signal_columns
from heteroplassimo.utils import get_eucl_distance
from heteroplassimo.utils import get_signal_ratio
from heteroplassimo.utils import get_log_signal_and_ratio
from heteroplassimo.utils import frames_to_time
import pandas as pd
import pydot


if __name__ == "__main__":

    #  ng_path = "NGAtp6_su9mKate001_crop_s08_acdc_output_NG.csv"
    #  mkate_path = "NGAtp6_su9mKate001_crop_s08_acdc_output_mKate.csv"
    df = pd.read_csv("multiExp_acdc_output__145.csv")

    # df = add_signal_columns(df)
    # df = load_dataframe(neon_green_csv_file=ng_path, m_kate_csv_file=mkate_path)

    for exp in sorted(df.experiment_foldername.unique()):

        exp_df = df[df["experiment_foldername"] == exp]

        for pos in exp_df.Position_n.unique():

            print(exp, pos)

            df_ = exp_df[
                (exp_df["experiment_foldername"] == exp) & (exp_df.Position_n == pos)
            ]
            df_ = frames_to_time(df_)
            df_ = add_signal_columns(df_)
            df_ = get_eucl_distance(df_)
            df_ = get_signal_ratio(df_)
            df_ = get_log_signal_and_ratio(df_)
            df_ = df_.reset_index()

            lineage_trees = LineagesTree(df=df_)
            lineage_trees.create_graphviz(
                f"dot_files/lineage_tree_{exp}_{pos}", show=True
            )
            lineage_trees.get_stats(exp=exp, pos=pos)
            # lineage_trees.loss_of_signal_stats(plot=True)
            # lineages_trees._get_bud_signal_state_counts(plot=True)
            lineage_trees.get_stats_v2(exp=exp, pos=pos)

            for i, tree in enumerate(lineage_trees.trees):
                graphs = pydot.graph_from_dot_file(
                    f"dot_files/lineage_tree_{exp}_{pos}_{i+1}.dot"
                )
                graph = graphs[0]
                graph.write_svg(f"svg_files/lineage_tree_{exp}_{pos}_{i+1}.svg")
      