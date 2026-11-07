"""F4: is the MC-dropout leak forecast calibrated?

The archived loop (LJ_test_sampled_vel_ACS_2_loop_1180.mat) stores the HRM leak and
Te_c sampled every 0.5 s: X_1_trajectory = [0; full root trajectory 0..1160 s
(2321 samples); one sample per later snapshot]. For every snapshot whose 25-step
horizon is fully inside the record (t_q <= 1167.5 s) we rebuild the BiLSTM window
ending at t_q, sample the ensemble and compare it with the realized leak.

Reported: coverage of the central 90 % band, PIT histogram, RMSE of the median vs the
spread, and the error of the integrated leak over the horizon (what drives the PCT,
since d(PCT)/d(inventory) ~ -1.6e-2 K/kg).

Usage: .venv\\Scripts\\python python/f4_forecast_calibration.py [n_traj] [seed]
"""
import os
import sys

import numpy as np
import scipy.io
import torch
import torch.nn as nn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from f4_leak_ensemble import ROOT, MODEL, Seq2SeqLSTM  # noqa: E402

N_TRAJ = int(sys.argv[1]) if len(sys.argv) > 1 else 500
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 1
H = 25            # forecast steps used by the MPC (t_q + 0.5 ... t_q + 12.5 s)
DT = 0.5
DPCT_DM = -1.6e-2  # K/kg, from F3


def main():
    sys.stdout.reconfigure(encoding="utf-8")   # Windows console (cp1252)
    a = scipy.io.loadmat(os.path.join(ROOT, "models", "LJ_test_sampled_vel_ACS_2_loop_1180.mat"),
                         variable_names=["X_1_trajectory", "X_2_trajectory", "X_minmax", "Y_minmax",
                                         "len_seq", "time_snapshots"])
    X1, X2 = a["X_1_trajectory"].ravel(), a["X_2_trajectory"].ravel()
    Xmm, Ymm, L = a["X_minmax"], a["Y_minmax"], int(a["len_seq"].item())
    ts = a["time_snapshots"].ravel()
    t_rec = np.concatenate([[np.nan], np.arange(2321) * DT, ts[1:len(X1) - 2322 + 1]])
    assert len(t_rec) == len(X1) and t_rec[-1] == 1180.0, (len(t_rec), t_rec[-1])

    model = Seq2SeqLSTM()
    model.load_state_dict(torch.load(MODEL, map_location="cpu"))
    for m in model.modules():
        (m.train() if isinstance(m, nn.LSTM) else m.eval())
    torch.manual_seed(SEED)

    rows, pits, inside = [], [], []
    print(" t_q     sd(1)  sd(25) | RMSE med | 90% cover | ΔM real-med [kg]  sd_ΔM [kg]  z     ΔPCT equiv [K]")
    for i in range(len(X1)):
        t_q = t_rec[i]
        if not (t_q >= 1160.0 and t_q + H * DT <= 1180.0):
            continue
        X = np.column_stack([X1[i - L + 1:i + 1], X2[i - L + 1:i + 1]])
        xt = torch.tensor((X - Xmm[:, 0]) / (Xmm[:, 1] - Xmm[:, 0]), dtype=torch.float32)
        with torch.no_grad():
            Y = np.stack([model(xt).numpy()[:H, 0] for _ in range(N_TRAJ)])
        Y = Y * (Ymm[0, 1] - Ymm[0, 0]) + Ymm[0, 0]          # (n, H) leak, kg/s
        real = X1[i + 1:i + 1 + H]
        lo, med, hi = np.percentile(Y, [5, 50, 95], axis=0)
        cov = np.mean((real >= lo) & (real <= hi))
        pit = (Y < real).mean(axis=0)
        dM = Y.sum(axis=1) * DT                               # integrated leak over the horizon
        dM_real = real.sum() * DT
        z = (dM_real - np.median(dM)) / dM.std()
        pits.append(pit); inside.append(cov)
        rows.append((t_q, Y.std(0)[0], Y.std(0)[-1], np.sqrt(np.mean((med - real) ** 2)), cov,
                     dM_real - np.median(dM), dM.std(), z))
        print(f"{t_q:7.1f}  {Y.std(0)[0]:.3f}  {Y.std(0)[-1]:.3f} |  {rows[-1][3]:.4f}  |   {cov:4.2f}    |"
              f"   {rows[-1][5]:+7.3f}        {dM.std():.3f}    {z:+5.2f}   {DPCT_DM * rows[-1][5]:+.4f}")

    pits = np.concatenate(pits)
    hist = np.histogram(pits, bins=10, range=(0, 1))[0] / len(pits)
    r = np.array(rows)
    print(f"\n{len(rows)} queries x {H} steps. Mean 90 % coverage: {np.mean(inside):.2f} (nominal 0.90)")
    print("PIT histogram (10 bins, uniform = 0.10 each):", np.round(hist, 3))
    print(f"Mean RMSE of the median: {r[:,3].mean():.4f} kg/s; mean ensemble sd: {np.mean(r[:,1:3]):.4f} kg/s")
    print(f"Integrated-leak error: |z| mean {np.abs(r[:,7]).mean():.2f}, max {np.abs(r[:,7]).max():.2f}; "
          f"max |ΔPCT| equivalent {np.abs(DPCT_DM * r[:,5]).max():.4f} K")
    os.makedirs(os.path.join(ROOT, "results", "f4_risk"), exist_ok=True)
    np.savetxt(os.path.join(ROOT, "results", "f4_risk", "forecast_calibration.csv"), r, delimiter=",",
               header="t_q,sd_step1,sd_step25,rmse_median,coverage90,dM_real_minus_median,sd_dM,z",
               comments="", fmt="%.6f")


if __name__ == "__main__":
    main()
