"""F1b: extract the HRM validation series from the archived LJ loop (no MATLAB needed).

Reads models/LJ_test_sampled_vel_ACS_2_loop_1180.mat (which also holds Simulink
operating-point objects that only MATLAB can read) and writes the plain arrays the
Julia HRM replica needs to results/f1b_hrm_replica/archive_1180.mat:

  control_series        (2361×3) [t, u_ACS current [A], active_mpc], as fed to the HRM
                         root inports (row 1 is the initial zero row of the loop script)
  time_snapshots        (61)     MPC query times
  tw_max_V2P_trajectory (2362×2) [t, max_i T_w,i] logged by the HRM
  X_1_trajectory        (2362)   HRM leak flow (mdot_c_leak_signal) [kg/s]
  X_2_trajectory        (2362)   HRM Te_c of node 12 (Te_c_signal(12)) [K]
  sampled_A, sampled_B, sampled_mdot_c_drift, sampled_vel_ACS, Lambda_t, Mu_t

Layout of the logged trajectories (see matlab/DET_scheduler_v1_perc_online_LJ.m): row 1 is
a zero row; rows 2..2322 are the full first HRM run 0..1160 s; each later row is the END
value of the HRM re-simulation from t = 0 for the next query (1160.5..1180 s).
"""
from pathlib import Path
import warnings

import numpy as np
import scipy.io as sio

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "models" / "LJ_test_sampled_vel_ACS_2_loop_1180.mat"
DST = ROOT / "results" / "f1b_hrm_replica" / "archive_1180.mat"

KEYS = ["control_series", "time_snapshots", "tw_max_V2P_trajectory", "X_1_trajectory",
        "X_2_trajectory", "sampled_A", "sampled_B", "sampled_mdot_c_drift",
        "sampled_vel_ACS", "Lambda_t", "Mu_t", "u_V2P_trajectory"]


def main():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        d = sio.loadmat(SRC, variable_names=KEYS, squeeze_me=True)
    out = {k: np.asarray(d[k], dtype=float) for k in KEYS}
    DST.parent.mkdir(parents=True, exist_ok=True)
    sio.savemat(DST, out, do_compression=True)
    for k, v in out.items():
        print(f"{k:24s} {v.shape}")
    print("written", DST.relative_to(ROOT))


if __name__ == "__main__":
    main()
