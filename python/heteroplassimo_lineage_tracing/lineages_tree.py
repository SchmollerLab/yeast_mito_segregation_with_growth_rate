import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from treelib import Tree, Node
from io import StringIO
import codecs
from heteroplassimo.cell_data import CellData, CellState


def rename_children(tree, parent):
    for i, child in enumerate(tree.children(parent.identifier)):
        new_name = f"{parent.tag}_{i+1}"
        tree.update_node(parent.identifier, identifier=parent.identifier)
        rename_children(tree, child)


class LineagesTree:
    def __init__(self, df: pd.DataFrame) -> None:
        # df["signal"] = df.NG_h5_amount_autoBkgr_maxProj > 25_000  # TODO make more general
        self.df = df
        self.trees = self._construct_tree_from_df(df=self.df)

    @property
    def n_zygotes(self) -> int:
        """Get number of zygotes."""
        return len(
            self.df.loc[
                (self.df.frame_i == 0)
                & (self.df.emerg_frame_i == -1)
                & (self.df.relationship == "mother")
            ]
        )

    def _construct_tree_from_df(self, df: pd.DataFrame):

        trees = [Tree() for _ in range(self.n_zygotes)]

        print(f"{len(trees)} Trees", trees)

        for index, tree in enumerate(trees):
            print(f"Creating Zygote{index+1} Node")
            tree.create_node(f"Zygote{index+1}", f"Zygote{index+1}")

        for row in self.df.itertuples():

            if row.Index < self.n_zygotes:
                print("Create first cell", row.Index + 1)

                last_mum_mk_signal = self.df[(self.df.Cell_ID == row.Cell_ID) & (self.df.relationship == 'mother')].sort_values(by='frame_i', ascending=True)["mKate_concentration_autoBkgr_from_vol_fl"].values[-1]
                last_mum_ng_signal = self.df[(self.df.Cell_ID == row.Cell_ID) & (self.df.relationship == 'mother')].sort_values(by='frame_i', ascending=True)["NG_concentration_autoBkgr_from_vol_fl"].values[-1]

                trees[row.Index].create_node(
                    f"Cell_ID_{row.Index + 1}",
                    f"Cell_ID_{row.Index + 1}",
                    parent=f"Zygote{row.Index + 1}",
                    data=CellData(
                        cell_id=row.Cell_ID,
                        relative_id=row.relative_ID,
                        frame_id= row.frame_i,
                        time = row.time,
                        cell_state=row.cell_state,
                        ng_signal=row.ng_signal,
                        mKate_signal=row.mKate_signal,
                        mKate_concentration_autoBkgr_from_vol_fl=row.mKate_concentration_autoBkgr_from_vol_fl,
                        NG_concentration_autoBkgr_from_vol_fl=row.NG_concentration_autoBkgr_from_vol_fl,
                        generation_num=row.generation_num,
                        mk_int_last_as_bud= last_mum_mk_signal,
                        ng_int_last_as_bud= last_mum_ng_signal,
                        false_ng_signals_until_now=0,
                        true_ng_signals_until_now=0,
                        false_mKate_signals_until_now=0,
                        true_mKate_signals_until_now=0,
                        eucl_dist= row.eucl_dist,
                        cell_signal_ratio=row.cell_signal_ratio,
                        log_ratio = row.log_ratio,
                        log_mk = row.log_mk,
                        log_ng = row.log_ng,
                        ratio_log = row.ratio_log,
                    ),
                )
                continue

            if row.relationship == "bud" and row.relative_ID < row.Cell_ID:

                last_bud_mk_signal = self.df[(self.df.Cell_ID == row.Cell_ID) & (self.df.relationship == 'bud')].sort_values(by='frame_i', ascending=True)["mKate_concentration_autoBkgr_from_vol_fl"].values[-1]
                last_bud_ng_signal = self.df[(self.df.Cell_ID == row.Cell_ID) & (self.df.relationship == 'bud')].sort_values(by='frame_i', ascending=True)["NG_concentration_autoBkgr_from_vol_fl"].values[-1]    

                ng_signal = self._get_ng_signal(df=df, row=row)
                mkate_signal = self._get_mkate_signal(df=df, row=row)
                ng_signal_counts = self._get_bud_signal_state_counts(
                    df=df, row=row, column_name="ng_signal"
                )
                mKate_signal_counts = self._get_bud_signal_state_counts(
                    df=df, row=row, column_name="mKate_signal"
                )

                try:
                    n_ng_false = ng_signal_counts[False]
                except:
                    n_ng_false = 0

                try:
                    n_ng_true = ng_signal_counts[True]
                except:
                    n_ng_true = 0

                try:
                    n_mKate_false = mKate_signal_counts[False]
                except:
                    n_mKate_false = 0

                try:
                    n_mKate_true = mKate_signal_counts[True]
                except:
                    n_mKate_true = 0

                data = CellData(
                    cell_id=row.Cell_ID,
                    relative_id = row.relative_ID,
                    frame_id=row.frame_i,
                    time = row.time, 
                    cell_state=row.cell_state,
                    ng_signal=row.ng_signal,
                    mKate_signal=row.mKate_signal,
                    mKate_concentration_autoBkgr_from_vol_fl=row.mKate_concentration_autoBkgr_from_vol_fl,
                    NG_concentration_autoBkgr_from_vol_fl=row.NG_concentration_autoBkgr_from_vol_fl,
                    generation_num=row.generation_num,
                    mk_int_last_as_bud= last_bud_mk_signal,
                    ng_int_last_as_bud= last_bud_ng_signal,
                    false_ng_signals_until_now=n_ng_false,
                    true_ng_signals_until_now=n_ng_true,
                    false_mKate_signals_until_now=n_mKate_false,
                    true_mKate_signals_until_now=n_mKate_true,
                    eucl_dist= row.eucl_dist,
                    cell_signal_ratio = row.cell_signal_ratio,
                    log_ratio = row.log_ratio,
                    log_mk = row.log_mk,
                    log_ng = row.log_ng,
                    ratio_log = row.ratio_log,

                )

                fails = 0

                for tree in trees:
                    try:
                        tree.create_node(
                            f"Cell_ID_{row.Cell_ID}",
                            f"Cell_ID_{row.Cell_ID}",
                            parent=f"Cell_ID_{int(row.relative_ID)}",
                            data=data,
                        )
                    except:
                        continue

                if fails == len(trees):
                    print(
                        f"Could not add Cell with id {row.Cell_ID} to any Tree; row Index: {row.Index}"
                    )
                    fails += 1

                fails = 0

        return trees

    @staticmethod
    def _get_ng_signal(df: pd.DataFrame, row) -> bool:
        """Check whether the signal is True or False.

        The row is always a row with relationship = "bud"

        Args:
            df (pd.DataFrame): the initial Dataframe
            row (...): row with relationshiop = "bud"

        Returns:
            bool: Wheather signal of given bud row is True or False
        """
        cell = df.loc[(df.Cell_ID == row.Cell_ID)]

        # signal of last bud cell state
        return (
            row.ng_signal
            if len(cell) == 1
            else cell.loc[cell.relationship == "bud"].ng_signal.values[-1]
        )

    @staticmethod
    def _get_mkate_signal(df: pd.DataFrame, row) -> bool:
        """Check whether the signal is True or False.

        The row is always a row with relationship = "bud"

        Args:
            df (pd.DataFrame): the initial Dataframe
            row (...): row with relationshiop = "bud"

        Returns:
            bool: Wheather signal of given bud row is True or False
        """
        cell = df.loc[(df.Cell_ID == row.Cell_ID)]

        # signal of last bud cell state
        return (
            row.mKate_signal
            if len(cell) == 1
            else cell.loc[cell.relationship == "bud"].mKate_signal.values[-1]
        )

    @staticmethod
    def _get_bud_signal_state_counts(
        df: pd.DataFrame, row, column_name: str
    ) -> pd.DataFrame:

        cell = df.loc[(df.Cell_ID == row.Cell_ID)]  # & (df.frame_i <= row.frame_i)]

        bud_cells = cell.loc[cell.relationship == "bud"]

        return bud_cells[column_name].value_counts()

    def get_stats_v2(self, exp: str, pos: str):

        df = pd.DataFrame(
            columns=[
                "tree",
                "lineage",
                "depth",
                'time',
                "cell_id",
                'relative_id',
                'mKate_intensity',
                'ng_intensity',
                "ng_signal",
                "mKate_signal",
                'mk_int_last_as_bud', 
                'ng_int_last_as_bud',
                "cell_state",
                "parent_cell_state",
                "eucl_dist",
                'cell_signal_ratio',
                'log_ratio',
                'log_mk',
                'log_ng',
                'ratio_log',
            ]
        )

        for i_tree, tree in enumerate(self.trees):

            lineages = tree.paths_to_leaves()

            for i_lineage, lineage in enumerate(lineages):

                for depth, cell in enumerate(lineage[1:]):

                    cell_node = tree.get_node(cell)
                    cell_data: CellData = cell_node.data

                    parent_node = tree.parent(cell_node.identifier)
                    parent_cell_data: CellData = parent_node.data

                    df = df.append(
                        {
                            "exp": exp,
                            "pos": pos,
                            "tree": i_tree,
                            "lineage": i_lineage,
                            "depth": depth,
                            'time': cell_data.time, 
                            'cell_id': cell_data.cell_id,
                            'mKate_intensity':cell_data.mKate_concentration_autoBkgr_from_vol_fl,
                            'ng_intensity':cell_data.NG_concentration_autoBkgr_from_vol_fl,
                            'cell_signal_ratio': cell_data.cell_signal_ratio,
                            'log_ratio': cell_data.log_ratio,
                            'log_mk': cell_data.log_mk, 
                            'log_ng': cell_data.log_ng,
                            'ratio_log': cell_data.ratio_log,
                            "cell_state": cell_data.cell_state,
                            'mk_int_last_as_bud': cell_data.mk_int_last_as_bud,
                            'ng_int_last_as_bud': cell_data.ng_int_last_as_bud,
                            "relative_id": cell_data.relative_id,
                            "parent_cell_state": parent_cell_data.cell_state if parent_cell_data else None,
                            "parent_cell_ratio_log": parent_cell_data.ratio_log if parent_cell_data else None, 
                            "eucl_dist": cell_data.eucl_dist
                        },
                        ignore_index=True,
                    )

        df.to_csv(f"csv_files/stats_v2_{exp}_{pos}.csv")
        print(df)

    # def loss_of_signal_stats(self, plot: bool = True):

    #     max_tree_depth = max(tree.depth() for tree in self.trees)
    #     columns = [f"Depth {i + 1}" for i in range(max_tree_depth)]

    #     df_depth = pd.DataFrame(columns=columns)

    #     for tree in self.trees:
    #         loss_of_ng_signal = dict(zip(columns, [0 for _ in range(max_tree_depth)]))
    #         loss_of_mKate_signal = dict(zip(columns, [0 for _ in range(max_tree_depth)]))

    #         lineages = tree.paths_to_leaves()

    #         for lineage in lineages:

    #             ng_signals = [tree.get_node(l).data.ng_signal for l in lineage[1:]]
    #             mKate_signals = [tree.get_node(l).data.mKate_signal for l in lineage[1:]]

    #             if False in ng_signals:
    #                 # find index/deptch of first False cell signal in lineage
    #                 depth = ng_signals.index(False) + 1  # to account for zygote
    #                 loss_of_ng_signal[f"Depth {depth}"] = loss_of_ng_signal[f"Depth {depth}"] + 1

    #             if False in mKate_signals:
    #                 # find index/deptch of first False cell signal in lineage
    #                 depth = mKate_signals.index(False) + 1  # to account for zygote
    #                 loss_of_mKate_signal[f"Depth {depth}"] = loss_of_mKate_signal[f"Depth {depth}"] + 1

    #         df_depth = df_depth.append({loss_of_ng_signal, loss_of_mKate_signal}, ignore_index=True)

    #     ax = df_depth.T.plot.bar()
    #     plt.ylabel("# of Cells")
    #     plt.tight_layout()
    #     plt.legend(title="Zygote")
    #     plt.savefig(f"loss_of_signal.pdf")

    def get_stats(self, exp: str, pos: str):

        df_ = pd.DataFrame(
            columns=[
                "zygote_tree",
                "number_of_lineages",
                "last_ng_signal_true",
                "ng_percentage",
                "last_mKate_signal_true",
                "mKate_percentage",
            ]
        )

        for index, tree in enumerate(self.trees):

            last_ng_signal_true = 0
            last_mKate_signal_true = 0
            lineages = tree.paths_to_leaves()

            for lineage in lineages:
                ng_signal_last = tree.get_node(lineage[-1]).data.ng_signal # edit here instead of ng_signal to have the log_ratio
                mKate_signal_last = tree.get_node(lineage[-1]).data.mKate_signal

                # if the last signal of the lineage is True
                # we count it
                if ng_signal_last:
                    last_ng_signal_true += 1
                if mKate_signal_last:
                    last_mKate_signal_true += 1

            df_ = df_.append(
                {
                    "zygote_tree": index + 1,
                    "number_of_lineages": len(lineages),
                    "last_ng_signal_true": last_ng_signal_true,
                    "ng_percentage": round(last_ng_signal_true / len(lineages), 2)
                    * 100,
                    "last_mKate_signal_true": last_mKate_signal_true,
                    "mKate_percentage": round(last_mKate_signal_true / len(lineages), 2)
                    * 100,
                },
                ignore_index=True,
            )

        df_.zygote_tree = df_.zygote_tree.astype(int)

        df_.to_csv(f"csv_files/stats_{exp}_{pos}.csv")
        print(df_)
 


    def create_graphviz(
        self, filename=None, shape="doublecircle", graph="digraph", show: bool = False
    ):
        """
        Export the tree(s) in the dot format of the graphviz software.

        Can be rendered for example here: https://dreampuf.github.io/GraphvizOnline/

        Args:
            filename ([type], optional): [description]. Defaults to None.
            shape (str, optional): [description]. Defaults to "doublecircle".
            graph (str, optional): [description]. Defaults to "digraph".
            show (bool, optional): [description]. Defaults to False.
        """
        trees = self.trees

        for index, tree in enumerate(trees):
            nodes, connections = [], []

            if show:
                tree.show()

            for n in tree.expand_tree(mode=tree.WIDTH):
                nid = tree[n].identifier
                cell_data = tree[n].data
                if not cell_data:
                    state = f'"{nid}" [label="{tree[n].tag}", shape={shape}]'
                else:
                    if cell_data.cell_state == CellState.HETEROPLASMIC:
                        color = "#feb24c"
                    elif cell_data.cell_state == CellState.HOMOPLASMIC_NG:
                        color = "#238b45"
                    elif cell_data.cell_state == CellState.HOMOPLASMIC_MKATE:
                        color = "#bd0026"
                    elif cell_data.cell_state == CellState.NO_SIGNAL:
                        color = "#bdbdbd"
                    state = f'"{nid}" [label="{tree[n].tag}\nFrame: {cell_data.frame_id}\nCell State: {cell_data.cell_state.value}",shape={shape}, color="{color}", style=filled]'
                nodes.append(state)

                for c in tree.children(nid):
                    cid = c.identifier
                    connections.append(f'"{nid}" -> "{cid}"')

            # write nodes and connections to dot format
            is_plain_file = filename is not None
            if is_plain_file:
                f = codecs.open(f"{filename}_{index+1}.dot", "w", "utf-8")
            else:
                f = StringIO()

            f.write(graph + " tree {\n")
            for n in nodes:
                f.write("\t" + n + "\n")

            if connections:
                f.write("\n")

            for c in connections:
                f.write("\t" + c + "\n")

            f.write("}")

            if not is_plain_file:
                print(f.getvalue())

            f.close()

    # for every node(cell), get parent and check cell_state
    # barplot of percentage of cells that have a different id than the parent_id
    # stats how many are cell_state=hetero and have children homo_ng or homo_mKate
    # treelib.node.Node(tag=None, identifier=None, expanded=True, data=None)
    #    bpointer
    # The parent ID of a node. This attribute can be accessed and modified with . and = operator respectively.

    # plot the amount of cells per cell_state in relation to the tree.level
    # class treelib.tree.Tree(tree=None, deep=False, node_class=None)
    #    level(nid, filter=None)
    # Get the node level in this tree. The level is an integer starting with ‘0’ at the root.
    # In other words, the root lives at level ‘0’;

    # def bla(self):
    #    trees = self.trees
    #    tree = trees[0]
    #    root_node = tree.get_node(tree.root)
    #    first_node = tree.children(root_node.identifier)[0]
    #    rename_children(tree, first_node)
