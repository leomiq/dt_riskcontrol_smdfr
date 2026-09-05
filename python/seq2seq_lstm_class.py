import torch
import torch.nn as nn

class Seq2SeqLSTM(nn.Module):
    def __init__(self, input_dim, output_dim, hidden_dim=32, num_layers=2, dropout=0.1):
        super(Seq2SeqLSTM, self).__init__()
        self.hidden_dim = hidden_dim
        self.num_layers = num_layers
        
        self.lstm = nn.LSTM(
            input_dim,
            hidden_dim,
            num_layers=num_layers,
            dropout=dropout,
            batch_first=True,
            bidirectional=True
        )
        
        self.fc = nn.Sequential(
            nn.Linear(hidden_dim * 2, 10),
            nn.ReLU(),
            nn.Linear(10, output_dim)
        )
    
    def forward(self, x):
        out, _ = self.lstm(x)
        out = self.fc(out)
        return out
