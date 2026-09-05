#%%
# mc_inference.py

import torch
import torch.nn as nn
import numpy as np
import scipy.io

import os

# ensure we load the .pth from this script’s folder:
#_HERE = os.path.dirname(os.path.abspath(__file__))
#MODEL_PATH = os.path.join(_HERE, "model_200x200_drop1_1.pth")

# Absolute path to your .pth
MODEL_PATH = "/g100/home/userexternal/lmiquele/master/test_python_3_8_12/gbm_rmm_cineca/model_200x200_drop1_1.pth"

# Debug print so you see it in MATLAB’s Command Window
print(f"DEBUG: loading model from {MODEL_PATH}; exists? {os.path.exists(MODEL_PATH)}; size={os.path.getsize(MODEL_PATH)}")

# 1) Re-declare the SAME model architecture
class Seq2SeqLSTM(nn.Module):
    def __init__(self, input_dim, output_dim, hidden_dim=32, num_layers=2, dropout=0.05):
        super().__init__()
        self.lstm = nn.LSTM(input_dim, hidden_dim, num_layers=num_layers,
                            dropout=dropout, batch_first=True, bidirectional=True)
        self.fc = nn.Sequential(
                nn.Linear(64, 10),  # 10 neurons
                nn.ReLU(),          # ReLU activation
                nn.Linear(10, output_dim)
            )

    def forward(self, x):
        out, _ = self.lstm(x)   # out: [batch_size, seq_len, 64]
        out = self.fc(out)          # apply FC(64->10->output_dim) + ReLU
        return out
    
# 2) Helper to enable dropout layers in train mode, but keep others in eval
def enable_mc_dropout(net):
    """Force dropout layers to be in train mode (to sample at inference),
       while keeping other layers (like LSTM, Linear) in eval mode if you wish."""
    for module in net.modules():
        if isinstance(module, nn.LSTM):
            #print(f"Enabling dropout for {module}")
            module.train()
        else:
            # If you prefer to keep everything else in eval mode:
            module.eval()

# 3) Load the trained model once at module-level
#    So we don't reload it every time we call 'mc_dropout_infer'.
model = Seq2SeqLSTM(input_dim=2, output_dim=2, hidden_dim=32, dropout=0.05)

#%% Loading weights from the trained model
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
#state_dict = torch.load("model_200x200_drop1_1.pth", map_location=device)
state_dict = torch.load(MODEL_PATH, map_location=device)

# Loading weights (dictionary) into the model architecture
model.load_state_dict(state_dict)

#%%
model.to(device)

#%%
model.eval()  # initial mode is eval()

#%% MC dropout inference

def mc_dropout_single_rollout(model, init_window, device="cuda"):
    window_tensor = torch.tensor(init_window, dtype=torch.float32).to(device)
    enable_mc_dropout(model)  # Turn on dropout, keep rest in eval
    
    with torch.no_grad():
        # forward pass
        out = model(window_tensor)
        # unsqueeze to add batch dimension
        out = out.unsqueeze(0)
        
    # move to cpu and convert to numpy
    preds_t = out.cpu().numpy().reshape(-1, init_window.shape[1])

    return preds_t

def mc_dropout_trajectories(model, init_window, n_trajectories=59, device="cuda"):
    #model.to(device)
    trajs_list = []
    for _ in range(n_trajectories):
        traj = mc_dropout_single_rollout(model, init_window, device)
        trajs_list.append(traj)  # shape => [steps, output_dim]

    # shape => [n_trajectories, steps, output_dim]
    result_np = np.stack(trajs_list, axis=0)

    #scipy.io.savemat("mc_results.mat", {"trajectories": result_np})

    return result_np

# %% Example usage

# init_window = np.random.rand(1,50,2)
# single_roll = mc_dropout_single_rollout(model, init_window)
# trajs = mc_dropout_trajectories(model, init_window, n_trajectories=59)

# %% plot for each output dimension
# import matplotlib.pyplot as plt

# for i in range(2):
#     plt.figure()
#     #plt.plot(init_window[0,:,i], label="Initial Window")
#     for j in range(trajs.shape[0]):
#         #plt.plot(single_roll[:,i], label="Single Rollout")
#         plt.plot(trajs[j,:,i])

#     plt.legend()
#     plt.title(f"Output Dimension {i}")
#     plt.show()

# %%
