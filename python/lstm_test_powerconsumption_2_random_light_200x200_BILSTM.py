#%% Packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import os
import random as rd
from numpy import array
from sklearn.metrics import mean_squared_error
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader
import copy

#%% Parameters
# variable failure time
# mission time (t_M), in seconds
t_M = 1000
# fixed sample time (t_sam), in seconds
t_sam = 0.5
# number of fixed observation points (length of input sequence)
n_input = 100
n_out = 100

#%% import dataset object from dataset_gbm.py

from dataset_gbm import dataset

#%% OPTION 1. Splitting and formatting a single sequence from dataset (sequential-based splitting)

# split each univariate dataset into train/test sets, sequentially
def split_dataset(dataset_list, n_input = 100, n_out = 100, percent_train = 0.8):
    train_list = list()
    test_list = list()
    for data in dataset_list:
        # split into standard weeks
        total_aval_samples = len(data) - n_input - n_out #from index 0 to len(data)-n_input-n_out
        if total_aval_samples >= 0:
            train_part = int(total_aval_samples*percent_train); test_part = total_aval_samples - train_part
            #mod_timesteps = perc_part%n_input

            # randomly taking the upper or lower sequence
            if rd.randint(0,1) == 1:
                # test_part = len(data) - (perc_part + mod_timesteps)
                train_samples_index = data[:train_part]
                test_samples_index = data[train_part:train_part + test_part]
            else:
                test_samples_index = data[:test_part]
                train_samples_index = data[test_part:test_part + train_part]

            train_samples = [data[j:j+n_input] for j in range(train_samples_index.shape[0])]
            train_list.extend(train_samples)

            test_samples = [data[j:j+n_input] for j in range(test_samples_index.shape[0])]
            test_list.extend(test_samples)

    #creating the arrays
    train_array = array(train_list); test_array = array(test_list)

    #shuffle the instances?
    np.random.shuffle(train_array); np.random.shuffle(test_array)

    #train_array = train_array.shuffle()
    return array(train_array), array(test_array)

#%% OPTION2 . Taking random samples for training and testing from the dataset (random-based splitting)

# split a univariate dataset into train/test sets
def select_random_dataset(dataset_list, n_input, n_out, total_scenarios=10, percent_train = 0.85, seed=42):
    #rd.Random(seed) for reproducibility
    to_train_list = list()
    #init_sample_ids = list()
    to_test_list = list()

    # taking random samples from dataset without replacement
    selected_scenarios = rd.Random(seed).sample(range(len(dataset_list)), total_scenarios)
    
    for ID_sim in selected_scenarios:
        dataset = dataset_list[ID_sim]
        n_dataset = dataset.shape[0]
        # discarding the last n_input + n_out samples
        n_obs = list(range(n_dataset - n_input - n_out))
        rd.Random(seed).shuffle(n_obs)
        n_samples_train = int(len(n_obs)*percent_train)
        n_obs_train = n_obs[:n_samples_train]
        n_obs_test = n_obs[n_samples_train:]
        train_samples = [dataset[j:j+n_input] for j in n_obs_train]
        to_train_list.extend(train_samples)
        test_samples = [dataset[j:j+n_input] for j in n_obs_test]
        to_test_list.extend(test_samples)

    return array(to_train_list), array(to_test_list)

#%% OPTION3. Splitting the simulations into training, validation and testing sets (simulation-based splitting)

# separate the simulation ids of dataset into training, validation and testing (preserving sequential data)
def select_random_sims(dataset_list, n_input, n_out, num_train=0.6, num_validation=0.3, num_test=0.1, seed=42):
    to_train_list = list()
    to_validation_list = list()
    to_test_list = list()
    # list of ID_SIM
    ID_sim_list = [i for i in range(len(dataset_list))]
    # split randomly the list of ID_SIM into training, validation and testing disjoint lists
    rd.Random(seed).shuffle(ID_sim_list)
    n_train = int(len(ID_sim_list)*num_train)
    n_validation = int(len(ID_sim_list)*num_validation)
    n_test = int(len(ID_sim_list)*num_test)
    id_train = ID_sim_list[0:n_train]
    id_validation = ID_sim_list[n_train:n_train+n_validation]
    id_test = ID_sim_list[n_train+n_validation:n_train+n_validation+n_test]
    # select the samples from the dataset for to_train_list
    for i in id_train:
        dataset = dataset_list[i]
        train_samples = [dataset[j:j+n_input] for j in range(dataset.shape[0]-n_input-n_out)]
        to_train_list.extend(train_samples)
    # select the samples from the dataset for to_validation_list
    for i in id_validation:
        dataset = dataset_list[i]
        validation_samples = [dataset[j:j+n_input] for j in range(dataset.shape[0]-n_input-n_out)]
        to_validation_list.extend(validation_samples)
    # select the samples from the dataset for to_test_list
    for i in id_test:
        dataset = dataset_list[i]
        test_samples = [dataset[j:j+n_input] for j in range(dataset.shape[0]-n_input-n_out)]
        to_test_list.extend(test_samples)

    train_array = array(to_train_list); validation_array = array(to_validation_list); test_array = array(to_test_list)
    #np.random.seed(seed)
    #np.random.shuffle(train_array); np.random.shuffle(validation_array); np.random.shuffle(test_array)

    return train_array, validation_array, test_array

#%% To training, to validation, to test
# convert history into inputs and outputs
def to_training(dataset, to_train, n_input, n_out):
    #x_1, x_2, x_3, x_4 = list(), list(), list(), list()
    X = list()
    Y = list()
    #y_1, y_2 = list(), list()
    #in_start = 0
    # step over the entire history one time step at a time
    for i in range(len(to_train)):
        # get the first time value of the sample
        ID_sim = int(to_train[i][0,0])
        tin_start = to_train[i][0,1]
        index = np.where(dataset[ID_sim][:,1] == tin_start)
        in_start = int(index[0])
        in_end = in_start + n_input
        out_end = in_end + n_out
        # ensure we have enough data for this instance
        if out_end <= dataset[ID_sim].shape[0]:
            #x_1_input = dataset[ID_sim][in_start:in_end, 2]
            #x_1_input = x_1_input.reshape((len(x_1_input), 1))
            #x_1.append(x_1_input)
            #x_2_input = dataset[ID_sim][in_start:in_end, 3]
            #x_2_input = x_2_input.reshape((len(x_2_input), 1))
            #x_2.append(x_2_input)
            #x_3_input = dataset[ID_sim][in_start:in_end, 4]
            #x_3_input = x_3_input.reshape((len(x_3_input), 1))
            #x_3.append(x_3_input)
            #x_4_input = dataset[ID_sim][in_start:in_end, 4]
            #x_4_input = x_4_input.reshape((len(x_4_input), 1))
            #x_4.append(x_4_input)
            X.append(dataset[ID_sim][in_start:in_end, 2:4])
            #concat_y = np.concatenate((data[in_end:out_end, 1], data[in_end:out_end, 3]))
            #y.append(concat_y)
            #y_1.append(dataset[ID_sim][in_end:out_end, 2])
            #y_2.append(dataset[ID_sim][in_end:out_end, -1])
            Y.append(dataset[ID_sim][in_end:out_end, [2,-1]])
    #return array(x_1), array(x_2), array(x_3), array(x_4), array(y_1), array(y_2)
    #return array(X), array(y_1), array(y_2)
    return array(X), array(Y)

# convert history into inputs and outputs
def to_validation(dataset, to_validation, n_input, n_out):
    #x_1, x_2, x_3, x_4 = list(), list(), list(), list()
    X = list()
    Y = list()
    #in_start = 0
    # step over the entire history one time step at a time
    for i in range(len(to_validation)):
        # get the first time value of the sample
        ID_sim = int(to_validation[i][0,0])
        tin_start = to_validation[i][0,1]
        index = np.where(dataset[ID_sim][:,1] == tin_start)
        in_start = int(index[0])
        in_end = in_start + n_input
        out_end = in_end + n_out
        # ensure we have enough data for this instance
        if out_end <= dataset[ID_sim].shape[0]:
            x_input = dataset[ID_sim][in_start:in_end, 2:4]
            #x_input = x_input.reshape((n_input, 2))
            X.append(x_input)
            #concat_y = np.concatenate((data[in_end:out_end, 1], data[in_end:out_end, 3]))
            #y.append(concat_y)
            y_output = dataset[ID_sim][in_end:out_end, [2,-1]]
            Y.append(y_output)
            #y_2.append(dataset[ID_sim][in_end:out_end, -1])
    #return array(x_1), array(x_2), array(x_3), array(x_4), array(y_1), array(y_2)
    return array(X), array(Y)

# convert history into inputs and outputs
def to_testing(dataset, to_test, n_input, n_out):
    #...
    history = list()
    test = list()
    for i in range(len(to_test)):
        # get the first time value of the sample
        ID_sim = int(to_test[i][0,0])
        tin_start = to_test[i][0,1]
        index = np.where(dataset[ID_sim][:,1] == tin_start)
        in_start = int(index[0])
        in_end = in_start + n_input
        out_end = in_end + n_out
        # ensure we have enough data for this instance
        if out_end <= dataset[ID_sim].shape[0]:
            x_input = dataset[ID_sim][in_start:in_end, 2:4]
            #x_input = x_input.reshape((n_input, 2))
            history.append(x_input)
            #concat_y = np.concatenate((dataset[i][in_end:out_end, 1], dataset[i][in_end:out_end, 3]))
            #y.append(concat_y)
            #y_1.append(dataset[i][in_end:out_end, 2])
            #y_2.append(dataset[i][in_end:out_end, 3])
            y_mat = dataset[ID_sim][in_end:out_end, [2,-1]]
            test.append(y_mat)
    return array(history), array(test)

#%% Seq2Seq LSTM Model (Pytorch)

class Seq2SeqLSTM(nn.Module):
    def __init__(self, input_dim, output_dim, hidden_dim=32, num_layers=2, dropout=0.2):
        super(Seq2SeqLSTM, self).__init__()
        self.hidden_dim = hidden_dim
        self.num_layers = num_layers
        
        # Bidirectional LSTM
        self.lstm = nn.LSTM(input_dim, hidden_dim, 
                            num_layers=num_layers, 
                            dropout=dropout, 
                            batch_first=True, 
                            bidirectional=True)
        
        # Fully connected layers
        #self.fc1 = nn.Linear(hidden_dim * 2, 10)  # 2x hidden_dim for bidirectional
        #self.relu = nn.ReLU()
        #self.fc2 = nn.Linear(10, output_dim)

    # def forward(self, x):
    #     # LSTM outputs a sequence of hidden states
    #     lstm_out, _ = self.lstm(x)  # shape: [batch_size, seq_len, 2*hidden_dim]
        
    #     # If you only need the last time step:
    #     last_output = lstm_out[:, -1, :]  # shape: [batch_size, 2*hidden_dim]
        
    #     # Pass the last hidden state through the fully connected layers
    #     x = self.relu(self.fc1(last_output))  # shape: [batch_size, 10]
    #     x = self.fc2(x)                       # shape: [batch_size, output_dim]

    #     # If you prefer [batch_size, 1, output_dim], do x = x.unsqueeze(1)
    #     x = x.unsqueeze(1)
    #     return x

        self.fc = nn.Sequential(
                nn.Linear(hidden_dim*2, 10),  # 10 neurons
                nn.ReLU(),          # ReLU activation
                nn.Linear(10, output_dim)
            )

    def forward(self, x):
        out, _ = self.lstm(x)   # out: [batch_size, seq_len, 64]
        out = self.fc(out)          # apply FC(64->10->output_dim) + ReLU
        return out


def build_model_seq2seq_pytorch(train_X, train_Y, val_X, val_Y,
                                n_epochs=100, batch_size=64, 
                                learning_rate=1e-3, device='cuda', patience=30):
    """
    Build and train a Seq2Seq LSTM model using DataLoader for mini-batching.
    
    Parameters:
        train_X: numpy array or torch tensor, shape [N, seq_len, n_features]
        train_Y: numpy array or torch tensor, shape [N, (1), output_dim]
        val_X:   validation features, same shape structure as train_X
        val_Y:   validation targets, same shape structure as train_Y
        n_epochs: number of epochs to train
        batch_size: how many samples per training step
        learning_rate: learning rate for optimizer
        device: 'cpu' or 'cuda' for GPU
    
    Returns:
        Trained PyTorch model
    """
    # 1) Create TensorDatasets from your train and validation data
    train_dataset = TensorDataset(torch.tensor(train_X, dtype=torch.float32),
                                  torch.tensor(train_Y, dtype=torch.float32))
    val_dataset   = TensorDataset(torch.tensor(val_X,   dtype=torch.float32),
                                  torch.tensor(val_Y,   dtype=torch.float32))

    # 2) Create DataLoaders
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader   = DataLoader(val_dataset,   batch_size=batch_size, shuffle=False)

    # 3) Model initialization
    n_features = train_X.shape[-1]      # For example, if train_X is [N, 50, 2], n_features=2
    output_dim = train_Y.shape[-1]      # If train_Y is [N, 1, 2], then output_dim=2

    model = Seq2SeqLSTM(input_dim=n_features, output_dim=output_dim).to(device)

    # 4) Define loss and optimizer
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)

    # Lists to store losses
    train_losses = []
    val_losses   = []
    best_val_loss = np.inf
    best_model_state = None

    epochs_no_improve = 0

    # 5) Training loop
    for epoch in range(n_epochs):
        model.train()  # put model in training mode
        running_loss = 0.0
        running_val_loss = 0.0

        # Iterate over batches from train_loader
        for x_batch, y_batch in train_loader:
            # Move data to the specified device (CPU or GPU)
            x_batch = x_batch.to(device)
            y_batch = y_batch.to(device)

            # Forward pass
            outputs = model(x_batch)
            loss = criterion(outputs, y_batch)

            # Backward pass and optimization
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            running_loss += loss.item()

        # Average training loss for this epoch
        train_loss = running_loss / len(train_loader)
        train_losses.append(train_loss)

        # 6) Validation loop
        model.eval()  # put model in evaluation mode
        with torch.no_grad():
            for x_val, y_val in val_loader:
                x_val = x_val.to(device)
                y_val = y_val.to(device)
                
                val_outputs = model(x_val)
                v_loss = criterion(val_outputs, y_val)
                running_val_loss += v_loss.item()

        # Average validation loss
        epoch_val_loss = running_val_loss / len(val_loader)
        val_losses.append(epoch_val_loss)

        print(f"Epoch [{epoch+1}/{n_epochs}] | "
              f"Train Loss: {train_loss:.4f} | Val Loss: {epoch_val_loss:.4f}")

        # Check if this is the best model so far
        if epoch_val_loss < best_val_loss:
            best_val_loss = epoch_val_loss
            # Save the best model state
            best_model_state = copy.deepcopy(model.state_dict())
            epochs_no_improve = 0  # ...and reset counter
        else:
            epochs_no_improve += 1

        # Early stopping check
        if epochs_no_improve >= patience:
            print(f"Early stopping triggered after {epoch+1} epochs.")
            break
    
    # 7) Restore the best model
    if best_model_state is not None:
        model.load_state_dict(best_model_state)

    return model, train_losses, val_losses

#%% Extra functions

# make a forecast
def forecast(model, history_sample, n_features):
    # history_sample is an array vector
    input_x = history_sample.reshape((1, len(history_sample), n_features))
    # forecast the two timeseries with the input
    yhat = model.predict(input_x, verbose=0)
    # we only want the vector forecast
    #yhat = yhat[0]
    return yhat
 
# evaluate a single model
def evaluate_model(model, history, test, n_input, n_out):
    # fit model
    # model, backup = build_model(train, n_input, n_out)
    # history is an array of (sample,timesteps_input1)
    # history = [x for x in train]
    # walk-forward validation over each sample
    predictions = list()
    for i in range(len(history)):
        # predict the timeseries with the input
        history_sample = history[i,:]
        yhat_sequence = forecast(model, history_sample, n_input)
        # store the predictions
        predictions.append(yhat_sequence)
    # evaluate predictions days for each week
    predictions = array(predictions)
    score, scores = evaluate_forecasts(test, predictions)
    return score, scores

# evaluate one or more weekly forecasts against expected values
def evaluate_forecasts(test, predictions):
    scores = list()
    # calculate an RMSE score for each day
    for i in range(test.shape[2]):
        scores[i] = mean_squared_error(test[:,:,i], predictions[:,:,i])
    score = mean_squared_error(test, predictions)
    return score, scores

# summarize scores
def summarize_scores(name, score, scores):
    s_scores = ', '.join(['%.1f' % s for s in scores])
    print('%s: [%.3f] %s' % (name, score, s_scores))

#%% Evaluate the model

# load the new file
#dataset = read_csv('household_power_consumption_days.csv', header=0, infer_datetime_format=True, parse_dates=['datetime'], index_col=['datetime'])

#train, test = split_dataset(dataset, n_input, n_out)
#train, test = select_random_dataset(dataset, n_input, n_out, total_scenarios = len(dataset), percent_train = 0.85, seed = 42)
train, validation, test = select_random_sims(dataset, n_input, n_out, num_train=0.6, num_validation=0.3, num_test=0.1, seed = 42)

to_train = to_training(dataset, train, n_input, n_out)

to_valid = to_validation(dataset, validation, n_input, n_out)

to_test = to_testing(dataset, test, n_input, n_out)

#%% Build the model

#model, backup = build_model(dataset, train, n_input, n_out)

start = time.time()
model, train_losses, val_losses = build_model_seq2seq_pytorch(to_train[0], to_train[1], to_valid[0], to_valid[1], n_epochs=300, batch_size=64, learning_rate=1e-3)
end = time.time()

elapsed_time = (end - start)/3600

#%% Plotting the training and validation losses

plt.plot(train_losses[1:], label='train')
plt.plot(val_losses[1:], label='validation')
#plt.title('model train vs validation loss')
plt.ylabel('loss (MSE) [-]')
plt.xlabel('epoch [-]')
plt.legend()
plt.show()

#%% Evaluate the model

x_input = np.random.rand(1,100,2)
y_pred = model(torch.tensor(x_input, dtype=torch.float32, device='cuda'))

#%% Save best model

#torch.save(model.state_dict(), 'model_200x200_drop1.pth')
torch.save(model.state_dict(), 'model_100x100_drop2_BILSTM.pth')
