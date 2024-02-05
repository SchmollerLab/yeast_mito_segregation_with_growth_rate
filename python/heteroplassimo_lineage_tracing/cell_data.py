from dataclasses import dataclass
from enum import Enum


class CellState(Enum):
    HETEROPLASMIC = "heteroplasmic"
    HOMOPLASMIC_NG = "homoplasmic_ng"
    HOMOPLASMIC_MKATE = "homoplasmic_mKate"
    NO_SIGNAL = "no_signal"


@dataclass
class CellData:
    cell_id: int
    frame_id: int
    time: float
    cell_state: CellState
    ng_signal: bool
    mKate_signal: bool
    relative_id: int
    mKate_concentration_autoBkgr_from_vol_fl: float
    NG_concentration_autoBkgr_from_vol_fl: float
    mk_int_last_as_bud: float
    ng_int_last_as_bud: float
    log_mk: float
    log_ng: float 
    eucl_dist: float
    log_ratio: float
    ratio_log: float 
    cell_signal_ratio: float
    generation_num: float
    false_ng_signals_until_now: int
    true_ng_signals_until_now: int
    false_mKate_signals_until_now: int
    true_mKate_signals_until_now: int
