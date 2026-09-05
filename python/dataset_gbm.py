#%% Packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import torch
import random as rd
from numpy import array

#%% Parameters
# variable failure time
# mission time (t_M), in seconds
t_M = 1000
# fixed sample time (t_sam), in seconds
t_sam = 0.5
# !!! number of fixed observation points (length of input sequence) (MAYBE TO BE CHANGED TO 50)
n_input = 100
n_out = 100

#%% import dataset object from dataset_gbm.py
#from dataset_gbm import dataset

#%% Concatenating available data from the simulations

available_sheets = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
df_list_rom = []; df_list_hrm = []; df_list_fault_times = []

# Training databases are failure-injection sweeps (large .xls batches) and are
# not redistributed with this repository; set TRAINING_DB_DIR to their location
# to re-run the training pipeline. The trained weights are shipped in models/.
abs_path_datasets = os.path.join(
    os.environ.get(
        'TRAINING_DB_DIR',
        os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                     'data', 'failure_injection_databases')),
    '')

for i in available_sheets:
    rom_str = abs_path_datasets + 'data_drift_rom_v2_batch_' + str(i) + '.xls'
    hrm_str = abs_path_datasets + 'data_drift_hrm_batch_' + str(i) + '.xls'
    fault_times_str = abs_path_datasets + 'data_fault_times_batch_' + str(i) + '.xls'
    xls_rom = pd.read_excel(rom_str, sheet_name=None)
    xls_hrm = pd.read_excel(hrm_str, sheet_name=None)
    xls_fault_times = pd.read_excel(fault_times_str, sheet_name=None)
    df_list_rom.extend([xls_rom[key] for key in xls_rom.keys()])
    df_list_hrm.extend([xls_hrm[key] for key in xls_hrm.keys()])
    df_list_fault_times.extend([xls_fault_times[key] for key in xls_fault_times.keys()])


#%% Concatenating the fault times from the simulations

df_fault_times = pd.concat(df_list_fault_times, axis=0)
df_fault_times = df_fault_times.reset_index(drop=True)
df_fault_times['Simulation'] = [i for i in range(df_fault_times.shape[0])]

#%% Plotting a single simulation

# fig, axs = plt.subplots(2,3, figsize=(15,8))
# #fault_time = 500

# id_sample = 0
# fault_time = df_fault_times.iloc[id_sample,1]

# xl_hrm_sample = df_list_hrm[id_sample][df_list_hrm[id_sample]['t[s]'] > fault_time - n_input*t_sam]
# xl_rom_sample = df_list_rom[id_sample][df_list_rom[id_sample]['t[s]'] > fault_time - n_input*t_sam]

# # Set font size
# title_fontsize = 16
# label_fontsize = 14

# axs[0, 0].plot(xl_hrm_sample['t[s]'], xl_hrm_sample['mdot_c_leak'])
# axs[0, 0].set_title('Mass flow leakage rate', fontsize=title_fontsize)
# axs[0, 0].set_xlabel('Time [s]', fontsize=label_fontsize)
# axs[0, 0].set_ylabel(r'$\dot{m}_{\mathrm{c,leak}}$ [kg\s]', fontsize=label_fontsize)

# axs[0, 1].plot(xl_hrm_sample['t[s]'], xl_hrm_sample['M_c_1'])
# axs[0, 1].set_title('Coolant mass in node 1', fontsize=title_fontsize)
# axs[0, 1].set_xlabel('Time [s]', fontsize=label_fontsize)
# axs[0, 1].set_ylabel(r'$M_{c}^{1}$ [kg]', fontsize=label_fontsize)

# for i in range(12):
#     axs[0, 2].plot(xl_hrm_sample['t[s]'], xl_hrm_sample['tw_' + str(i + 1)], label='Node ' + str(i + 1))
# axs[0, 2].set_title('Piping wall temperature per node', fontsize=title_fontsize)
# axs[0, 2].set_xlabel('Time [s]', fontsize=label_fontsize)
# axs[0, 2].set_ylabel(r'$T_{w}^{i}$ [K]', fontsize=label_fontsize)
# axs[0, 2].legend(loc='upper right')

# for i in range(12):
#     axs[1, 0].plot(xl_hrm_sample['t[s]'], xl_hrm_sample['tc_' + str(i + 1)], label='Node ' + str(i + 1))
# axs[1, 0].set_title('Coolant temperature per node', fontsize=title_fontsize)
# axs[1, 0].set_xlabel('Time [s]', fontsize=label_fontsize)
# axs[1, 0].set_ylabel(r'$T_{c}^{i}$ [K]', fontsize=label_fontsize)
# axs[1, 0].legend(loc='upper right')

# axs[1, 1].plot(xl_hrm_sample['t[s]'], xl_hrm_sample['Tp'])
# axs[1, 1].set_title('Peak cladding temperature', fontsize=title_fontsize)
# axs[1, 1].set_xlabel('Time [s]', fontsize=label_fontsize)
# axs[1, 1].set_ylabel(r'$T_p$ [K]', fontsize=label_fontsize)

# axs[1, 2].plot(xl_hrm_sample['t[s]'], xl_hrm_sample['P_tot'] / 1E8)
# axs[1, 2].set_title('Normalized thermal power', fontsize=title_fontsize)
# axs[1, 2].set_xlabel('Time [s]', fontsize=label_fontsize)
# axs[1, 2].set_ylabel(r'$P_N$ [-]', fontsize=label_fontsize)

# plt.tight_layout()
# plt.show()

#%% Plotting a single simulation. difference wrt HRM and ROM

# id_sample = 0
# fault_time = df_fault_times.iloc[id_sample,1]

# xl_hrm_sample = df_list_hrm[id_sample][df_list_hrm[id_sample]['t[s]'] > fault_time - n_input*t_sam]
# xl_rom_sample = df_list_rom[id_sample][df_list_rom[id_sample]['t[s]'] > fault_time - n_input*t_sam]

# fig, axs = plt.subplots(1,2, figsize=(15,8))
# #fault_time = 500

# axs[0].plot(xl_rom_sample['t[s]'], xl_rom_sample['mdot_c_leak'], label=r'${\dot{m}}_{\mathrm{leak}}$ [kg/s], ROM')
# axs[0].plot(xl_hrm_sample['t[s]'], xl_hrm_sample['mdot_c_leak'], label=r'${\dot{m}}_{\mathrm{leak}}$ [kg/s], HRM')
# axs[0].set_xlabel('Time [s]')
# axs[0].legend()

# axs[1].plot(xl_rom_sample['t[s]'], xl_rom_sample['Tp'], label=r'$T_{\mathrm{p}}$ [K], ROM')
# axs[1].plot(xl_hrm_sample['t[s]'], xl_hrm_sample['Tp'], label=r'$T_{\mathrm{p}}$ [K], HRM')
# axs[1].axvline(x=fault_time, color='k', linestyle='--', label=r'fault time [s]')
# axs[1].set_xlabel('Time [s]')
# axs[1].legend()

#%% Selecting the timeframe of interest (failure time as random event)

df_fault_times = pd.concat(df_list_fault_times, axis=0)
# reset indexes
df_fault_times = df_fault_times.reset_index(drop=True)
#reset simulation ids from 0 to n
df_fault_times['Simulation'] = [i for i in range(df_fault_times.shape[0])]

for j in range(df_fault_times.shape[0]):
    fault_time = df_fault_times.iloc[j,1]
    df_list_rom[j] = df_list_rom[j][df_list_rom[j]["t[s]"] > fault_time - n_input*t_sam] 
    df_list_hrm[j] = df_list_hrm[j][df_list_hrm[j]["t[s]"] > fault_time - n_input*t_sam]

#%% Adding the delta_Tp column to the HRM dataframes

# computing delta_Tp
for i in range(len(df_list_rom)):
    xl_rom = df_list_rom[i]
    xl_hrm = df_list_hrm[i]

    xl_hrm['delta_Tp'] = xl_hrm['Tp'] - xl_rom['Tp']

#%% Calculating maximum and minimum for each feature, both rom and hrm
# maximum and minimum for each feature
n_trials = len(df_list_hrm)

max_rom = pd.concat([df_list_rom[i].max() for i in range(n_trials)], axis=1).max(axis=1)
max_hrm = pd.concat([df_list_hrm[i].max() for i in range(n_trials)], axis=1).max(axis=1)
#max_vec = pd.concat([max_rom, max_hrm], axis=1).max(axis=1)

min_rom = pd.concat([df_list_rom[i].min() for i in range(n_trials)], axis=1).min(axis=1)
min_hrm = pd.concat([df_list_hrm[i].min() for i in range(n_trials)], axis=1).min(axis=1)
#min_vec = pd.concat([min_rom, min_hrm], axis=1).min(axis=1)

# create excel file with the maximum and minimum values
max_rom = pd.DataFrame(max_rom).T
max_hrm = pd.DataFrame(max_hrm).T
min_rom = pd.DataFrame(min_rom).T
min_hrm = pd.DataFrame(min_hrm).T

# Create an Excel writer object
with pd.ExcelWriter('max_min_valuess.xlsx') as writer:
    # Write each DataFrame to a different sheet
    max_rom.to_excel(writer, sheet_name='max_rom')
    max_hrm.to_excel(writer, sheet_name='max_hrm')
    min_rom.to_excel(writer, sheet_name='min_rom')
    min_hrm.to_excel(writer, sheet_name='min_hrm')

#%% GLOBALLY normalize the inputs and outputs

def normalize_df(df, max_vec, min_vec):
    df_norm = df.copy()
    for i in df.columns:
        max_val = max_vec[i][0]    
        min_val = min_vec[i][0]
        df_norm[i] = (df[i] - min_val) / (max_val - min_val)
    return df_norm

df_list_rom_norm = [normalize_df(df_list_rom[i], max_rom, min_rom) for i in range(n_trials)]
df_list_hrm_norm = [normalize_df(df_list_hrm[i], max_hrm, min_hrm) for i in range(n_trials)]

#%% denormalize the time columns
for i in range(n_trials):
    df_list_rom_norm[i]['t[s]'] = df_list_rom[i]['t[s]'] 
    df_list_hrm_norm[i]['t[s]'] = df_list_hrm[i]['t[s]']

#%% Plotting simulations
# fig, axs = plt.subplots(5,2, figsize=(6,15))
# #fault_time = 500

# for i in range(5):
#     for j in range(2):
#         #plot the new added columns with respect to time
#         id_list = 5*j + i
#         xl_hrm = df_list_hrm[-id_list]
#         xl_rom = df_list_rom[-id_list]
#         #plt.plot(xl_hrm['t[s]'], xl_hrm['tw_max'], label='High-resolution model (Physical Object)')
#         axs[i, j].plot(xl_hrm['t[s]'], xl_hrm['mdot_c_leak'], label=r'$\Delta T_{\mathrm{max}}$ [K]')
#         axs[i, j].plot(xl_rom['t[s]'], xl_rom['mdot_c_leak'], label=r'$\Delta T_{\mathrm{max}}$ [K]')
#         #axs[i,j].plot(xl_hrm['t[s]'], xl_hrm['delta_Tp'], label=r'$T_{\mathrm{max}}$ [K]')
#         #axs[i,j].plot(xl_hrm['t[s]'], xl_hrm['Tp'], label=r'$ T_{\mathrm{max}}$ [K]')
#         #axs[i,j].plot(xl_hrm['t[s]'], xl_hrm['mdot_c_1'], label=r'${\dot{m}}_{c}^{1}$ [kg/s]')
#         #axs[i,j].plot(xl_hrm['t[s]'], xl_hrm['M_c_1'], label=r'${M}_{c}^{1}$ [kg]')
#         axs[i, j].legend()
#         axs[i, j].set_xlim([900,930])  # Set the timespan for the plot
#         #axs[i,j].axvline(x=500, color='k', linestyle='--')

# plt.tight_layout()
# plt.show()

# # # Obs: STRONG CORRELATION BETWEEN MDOT_C_1 AND DELTA_TW_MAX

#%% CONSTRUCTING THE TRAINING DATASET (IMPORTANT!)

# time, input_X, output_Y
# input_X = mdot_c_leak, Te_c_12_hrm, T_e_f_12_hrm, Tin_c
# output_Y = mdot_c_leak, delta_Tp
dataset = []

for i in range(len(df_list_hrm_norm)):
    data_list = pd.DataFrame()
    data_list['ID_sim'] = [i for _ in range(len(df_list_hrm_norm[i]))]
    data_list['t[s]'] = df_list_hrm[i]['t[s]'].values
    #data_list['mdot_f_1'] = df_list_hrm_norm[i]['mdot_f_1'].values
    data_list['mdot_c_leak'] = df_list_hrm_norm[i]['mdot_c_leak'].values
    data_list['Te_c_12_hrm'] = df_list_hrm_norm[i]['Te_c_12'].values
    data_list['Te_f_12_hrm'] = df_list_hrm_norm[i]['Te_f_12'].values
    data_list['Tin_c'] = df_list_hrm_norm[i]['Tin_c'].values
    data_list['delta_Tp'] = df_list_hrm_norm[i]['delta_Tp'].values
    data_list['Tp'] = df_list_hrm_norm[i]['Tp'].values

    data_list = data_list.to_numpy()
    dataset.append(data_list)

#%% separate the simulation ids of dataset into training, validation and testing (preserving sequential data)
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

#%% Selecting the training, validation and test sets
train, validation, test = select_random_sims(dataset, n_input, n_out, num_train=0.6, num_validation=0.3, num_test=0.1, seed = 42)
