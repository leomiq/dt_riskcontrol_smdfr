"""F4: regenerate the MC-dropout leak-forecast ensemble of an archived query.

Rebuilds the BiLSTM input window exactly as DET_scheduler_v1_perc_online_LJ.m
(run_mc_dropout_inference): last len_seq samples of X_1_trajectory (HRM leak) and
X_2_trajectory (Te_c, node 12), min-max normalized with X_minmax; the output is
denormalized with Y_minmax. Dropout is sampled as in mc_inference.py (LSTM in train
mode, rest in eval).

Checks: (i) the archived forecast Y_denorm (one dropout sample, N_par = 1) must be a
plausible member of the ensemble; (ii) the deterministic forecast (dropout off).

Usage, from the repository root:
    .venv\\Scripts\\python python/f4_leak_ensemble.py [archive.mat] [n_traj] [seed]
Writes results/f4_risk/leak_ensemble_<t_q>.mat.
"""
import os
import sys

import numpy as np
import scipy.io
import torch
import torch.nn as nn

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ARCHIVE = os.path.join(ROOT, "models", "LJ_test_sampled_vel_ACS_2_loop_1180.mat")
TQ = float(sys.argv[1]) if len(sys.argv) > 1 else None
N_TRAJ = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 1
MODEL = os.path.join(ROOT, "models", "model_200x200_drop1_1.pth")


class Seq2SeqLSTM(nn.Module):  # same architecture as mc_inference.py
    def __init__(self, input_dim=2, output_dim=2, hidden_dim=32, num_layers=2, dropout=0.05):
        super().__init__()
        self.lstm = nn.LSTM(input_dim, hidden_dim, num_layers=num_layers,
                            dropout=dropout, batch_first=True, bidirectional=True)
        self.fc = nn.Sequential(nn.Linear(64, 10), nn.ReLU(), nn.Linear(10, output_dim))

    def forward(self, x):
        out, _ = self.lstm(x)
        return self.fc(out)


def main():
    sys.stdout.reconfigure(encoding="utf-8")   # Windows console (cp1252)
    a = scipy.io.loadmat(ARCHIVE, variable_names=[
        "X_1_trajectory", "X_2_trajectory", "X_minmax", "Y_minmax", "len_seq",
        "Y_denorm", "N_par", "time_snapshots", "snapshot_id", "t_sam"])
    L = int(a["len_seq"].item())
    Xmm, Ymm = a["X_minmax"], a["Y_minmax"]
    X1, X2 = a["X_1_trajectory"].ravel(), a["X_2_trajectory"].ravel()
    # record layout: [0; root trajectory 0..1160 s every 0.5 s (2321); one sample per later snapshot]
    ts = a["time_snapshots"].ravel()
    t_rec = np.concatenate([[np.nan], np.arange(2321) * 0.5, ts[1:len(X1) - 2321]])
    t_q = float(ts[int(a["snapshot_id"].item()) - 1]) if TQ is None else TQ
    i = int(np.flatnonzero(t_rec == t_q)[0])
    X = np.column_stack([X1[i - L + 1:i + 1], X2[i - L + 1:i + 1]])
    real = X1[i + 1:]                  # realized leak after t_q (may be shorter than the horizon)
    Xn = (X - Xmm[:, 0]) / (Xmm[:, 1] - Xmm[:, 0])
    denorm = lambda Yn: Yn * (Ymm[:, 1] - Ymm[:, 0]) + Ymm[:, 0]

    model = Seq2SeqLSTM()
    model.load_state_dict(torch.load(MODEL, map_location="cpu"))
    xt = torch.tensor(Xn, dtype=torch.float32)

    model.eval()
    with torch.no_grad():
        y_det = denorm(model(xt).numpy())

    torch.manual_seed(SEED)
    for m in model.modules():          # enable_mc_dropout
        (m.train() if isinstance(m, nn.LSTM) else m.eval())
    with torch.no_grad():
        Y = np.stack([denorm(model(xt).numpy()) for _ in range(N_TRAJ)])   # (n, L, 2)

    # the archived single sample exists only for the archived snapshot
    y_arch = a["Y_denorm"][0] if i == len(X1) - 1 else y_det
    leak = Y[:, :, 0]
    q = np.percentile(leak, [5, 50, 95], axis=0)
    # rank of the archived sample within the ensemble, per step (0..1); uniform if it is a member
    rank = (leak < y_arch[:, 0]).mean(axis=0)
    rms = lambda d: np.sqrt(np.mean(d ** 2))
    dist_arch = rms(y_arch[:26, 0] - q[1, :26])
    dist_members = np.array([rms(leak[i, :26] - q[1, :26]) for i in range(N_TRAJ)])

    print(f"archive {os.path.basename(ARCHIVE)}: t_q = {t_q} s, N_par archived = {int(a['N_par'].item())}, len_seq = {L}")
    print(f"window: leak [{X[:, 0].min():.4f}, {X[:, 0].max():.4f}] kg/s, Te_c [{X[:, 1].min():.2f}, {X[:, 1].max():.2f}]")
    print(f"ensemble: {N_TRAJ} trajectories, seed {SEED}")
    print(" k   leak p05    p50     p95   sd     | det.    archived  rank")
    for k in range(26):
        print(f"{k+1:2d}  {q[0,k]:7.4f} {q[1,k]:7.4f} {q[2,k]:7.4f} {leak[:,k].std():.4f} | {y_det[k,0]:7.4f} {y_arch[k,0]:7.4f}  {rank[k]:.2f}")
    print(f"RMS distance to the ensemble median over the horizon: archived {dist_arch:.4f} kg/s; "
          f"members: median {np.median(dist_members):.4f}, p95 {np.percentile(dist_members, 95):.4f}; "
          f"archived sample at member percentile {100 * (dist_members < dist_arch).mean():.0f}")

    out = os.path.join(ROOT, "results", "f4_risk")
    os.makedirs(out, exist_ok=True)
    scipy.io.savemat(os.path.join(out, f"leak_ensemble_{int(t_q)}.mat"), {
        "Y": Y, "y_det": y_det, "y_arch": y_arch, "t_q": t_q, "t_sam": float(a["t_sam"].item()),
        "seed": SEED, "model": os.path.basename(MODEL), "leak_now": X[-1, 0], "real": real})


if __name__ == "__main__":
    main()
