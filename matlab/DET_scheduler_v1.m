%% Useful files to save the data
%simulation_core_dualfluid_FI3_sim_rom_extract;
%save('simcore_DF_sim_ROM_extract');

%% Parallel trees
%Ncomps = 5;
%parpool;
N_par = 59; M_workers = 10;
% A set of (independent) trees is...
forest_data = {};
forest_roots = {};
forest_roots_hrm = {};
% tree_par_index: indexes for trees to be run in parallel workers

%% importing maximum and minimum values for BiLSTM model ((de)normalization)

%Importing table from current directorz
%max_min_values = readtable('max_min_valuess.xlsx');
X_sheets_str = {'min_hrm','max_hrm'};
Y_sheets_str = {'min_hrm','max_hrm'};
X_labels = {'mdot_c_leak', 'Te_c_12'};
Y_labels = {'mdot_c_leak', 'delta_Tp'};

X_minmax = zeros(2,2);
Y_minmax = zeros(2,2);

% input_X = mdot_c_leak, Te_c_12_hrm, T_e_f_12_hrm, Tin_c
% output_Y = mdot_c_leak, delta_Tp

% dim(X) = [n_features, min and max]
for i = 1:length(X_sheets_str)
    sheet = X_sheets_str{i};
    table = readtable('max_min_valuess.xlsx','Sheet',sheet);
    X_minmax(:,i) = table{1,X_labels}';
end

for i = 1:length(Y_sheets_str)
    sheet = Y_sheets_str{i};
    table = readtable('max_min_valuess.xlsx','Sheet',sheet);
    Y_minmax(:,i) = table{1,Y_labels}';
end

%% Running the "ground truth" simulation
% At each time step, we will compute this MCDET
% Current time snapshot
time_snapshot_root = 650;
len_seq = 100;
% sampling time
t_sam = 0.5;
%mission time
tsim = time_snapshot_root + len_seq*t_sam;

%Saving the branching and tree data in a cell array
snapshots_data = {};

%% Resampling: parameters
% Number of branches for each resampling
N_resampling = 2;
N_branches = 2;

%% Loading LOCA initial conditions (HRM)

%fixing the seed for reproducibility
rng(1,"twister");

% loading variables of interest: coolant mass flow rate (leakage)
Ncomps = 5;
name_root = 'SMDFR_hrm_MCDET';
name_Wks = 'simulation_core_dualfluid_FI3_sim_hrm_extract_MCDET';

load_system(name_root);
simRoot = Simulink.SimulationInput(name_root);

% Getting the variables from the model workspace
mdlWks = get_param(name_root,'ModelWorkspace');
unif_dist_drift = getVariable(mdlWks,'unif_dist_drift');

% Getting sampling parameters for B (LOCA nonlinear dynamics)
B_l = getVariable(mdlWks,'B_l');
B_s = getVariable(mdlWks, 'B_s');
sampled_B = B_l + (B_s - B_l)*rand();

% Getting sampling parameters for B (LOCA nonlinear dynamics)
A_l = getVariable(mdlWks,'A_l');
A_s = getVariable(mdlWks, 'A_s');
sampled_A = A_l + (A_s - A_l)*rand();

% Setting a new sample for mdot_c_drift
sampled_mdot_c_drift = unif_dist_drift(1,2) + (unif_dist_drift(2,2) - unif_dist_drift(1,2))*rand();
simRoot = setVariable(simRoot, 'B', sampled_B, 'Workspace', name_root);
simRoot = setVariable(simRoot, 'A', sampled_A, 'Workspace', name_root);
simRoot = setVariable(simRoot, 'mdot_c_drift', sampled_mdot_c_drift, 'Workspace', name_root);
%simRoot = setVariable(simRoot, 'A', 0.18*sampled_mdot_c_drift, 'Workspace', name_root);
%disp(string(sampled_mdot_c_drift));

% Getting the variables from the model workspace
vel_ACS_ub = getVariable(mdlWks,'vel_ACS_ub');
vel_ACS_lb = getVariable(mdlWks,'vel_ACS_lb');

% Setting a new sample for mdot_c_drift
sampled_vel_ACS = vel_ACS_lb + (vel_ACS_ub - vel_ACS_lb)*rand();
%disp(string(sampled_vel_ACS))
simRoot = setVariable(simRoot,'vel_ACS', sampled_vel_ACS, 'Workspace', name_root);

% Lambda_t / Mu_t: TARGET failure and repair rates
Lambda_t = getVariable(mdlWks,'Lambda_t');
Mu_t = getVariable(mdlWks,'Mu_t');

% Including the initial parent last state to start the simulation
%parent_last_state = [0 1];

% Saving information of root simulation (initial state could be useful...)
% Setting the stopping time to finish the simulation. Start time at
% t = 0
simRoot = setModelParameter(simRoot, "StopTime", string(time_snapshot_root));

% Setting intial state of SimulationInput object
simRoot = setModelParameter(simRoot, "SaveFinalState", "on", "SaveOperatingPoint", "on", "FinalStateName", "finalOPRoot");

% Running the model
outRoot = sim(simRoot);

%Saving the data
finalOP_Raiz = outRoot.finalOPRoot;

% data of the root. struct:{'id': time snapshot integer, 'data':out type, 'final_OP': OP dtype, 'init_state': list}
init_state = [0,1];
forest_roots_hrm{1} = struct('id', time_snapshot_root, 'data', outRoot, 'final_OPRoot', finalOP_Raiz, 'init_state', init_state, 'Lambda_t', Lambda_t, 'Mu_t', Mu_t);


%% MC dropout inference

% Define input data for prediction in MATLAB (n_features = 2)
X_1 = outRoot.mdot_c_leak_signal.Data(end-len_seq+1:end);
X_2 = outRoot.Te_c_signal.Data(12,:,end-len_seq+1:end);

X_input = zeros(len_seq,2);  
X_input(:,1) = X_1; X_input(:,2) = X_2;

X_norm = norm_function(X_input, X_minmax);

% importing the module
mod = py.importlib.import_module('mc_inference');
%py.importlib.reload(mod);

% Convert to Python
X_py = py.numpy.array(X_norm);
num_trajs = py.int(59); % 3 trajectories

% measure computational time 
dropout_time = tic;
result_py = mod.mc_dropout_trajectories(mod.model, X_py, num_trajs);%, py.str('cuda'));
dropout_finish = toc(dropout_time);

% Convert python list -> MATLAB double
Y_norm = double(result_py);

Y_denorm = denorm_function(Y_norm, Y_minmax);

%% IS-MCDET branching, for a single time snapshot
init_loop = tic; 
%parfor (tree_par_index = 1:N_par, M_workers)
for tree_par_index = 1:N_par
    % loading variables of interest: coolant mass flow rate (leakage)
    Ncomps = 5;
    name_root = 'SMDFR_rom_MCDET_v2';
    name_Wks = 'simulation_core_dualfluid_FI3_sim_rom_extract_MCDET_v2';
    load_system('SMDFR_rom_MCDET_v2.slx');
    simromRoot = Simulink.SimulationInput(name_root);
    simromRoot = setModelParameter(simromRoot,"SaveFinalState", "on", "SaveOperatingPoint", "on", "FinalStateName", "finalOP");

    %tsim = time_snapshot_root + length_seq*t_sam;
    %predicted window:seconds
    simromRoot = setModelParameter(simromRoot,"StopTime", string(time_snapshot_root));

    % Combine into a single input matrix (for use in Inport). 
    %1: mdot_c_leak; 2:tw_max
    inport_data = [transpose(0:t_sam:time_snapshot_root), outRoot.mdot_c_leak_signal.Data(:), zeros(1,length(0:t_sam:time_snapshot_root))'];
    simromRoot = simromRoot.setExternalInput(inport_data);

    % Get target distributions parameters
    mdl_Wks = get_param(name_root,'ModelWorkspace');
    Lambda_t = getVariable(mdl_Wks,'Lambda_t');
    Mu_t = getVariable(mdl_Wks,'Mu_t');
    %! set sampled ACS velocity
    simromRoot = setVariable(simromRoot,'vel_ACS', sampled_vel_ACS, 'Workspace', name_Wks);

    outrom = sim(simromRoot);
    close_system('SMDFR_rom_MCDET_v2.slx');

    romOP = outrom.finalOP;
    
    queue_idle = {};
    % structure of the queue: {snapshot_time,id_str,states}, where states = [demand of ps, availability of ps]
    % including first node as current snapshot of the system (i.e., the root)
    root_state = romOP.loggedStates{14}.Values.Data;
    queue_idle{1} = {0,0,root_state};
    
    tree = {};
    % including first element of the tree
    % structure: left (parent node), right (child nodes), 'id' (integer value)
    tree{1} = struct('left', 'R', 'right', [],'id', 0);
    
    % data of the search tree. struct:{'id': integer type, 'data':out type, 'final_OP': OP dtype, 'init_state': list}
    tree_data = {};
    
    % Setting up the simulation loop
    % We set a final snapshot of the operating point at the end of the
    % simulation (branching, and not necessarily the mission time)
    % time_snapshot = 0;
    %tsim = time_snapshot_root + length_seq*t_sam;
    % counter k_nodes is defined separately to add up new nodes
    k_nodes = 0;

    % copying forest roots for each tree generated on each predicted signal
    init_state_rom = [0,1];

    forest_roots{tree_par_index} = struct('id', time_snapshot_root, 'data', outrom, 'final_OPRoot', romOP, 'init_state', init_state_rom, 'Lambda_t', Lambda_t, 'Mu_t', Mu_t);

    while ~ isempty(queue_idle)
        branch = queue_idle{1};
        current_node = branch{2};
        %display(current_node)
        if current_node == 0
            load_system('SMDFR_rom_MCDET_v2.slx');
            simIn = Simulink.SimulationInput('SMDFR_rom_MCDET_v2');

            % feeding WBM with PO "design parameters?"
            for index = [1:2,8:23]
                vec_values = finalOP_Raiz.loggedStates{index}.Values.Data;
                romOP.loggedStates{index}.Values.Data = vec_values;
            end 
            
            % feeding WBM with PO "time-dependent measurable quantities"
            for index = [3,4,5]
                vec_values = finalOP_Raiz.loggedStates{index+1}.Values.Data;
                romOP.loggedStates{index}.Values.Data = mean(vec_values(1:12));
            end 
            
            % feeding WBM with PO "time-dependent measurable quantities?"
            romOP.loggedStates{6}.Values.Data = finalOP_Raiz.loggedStates{7}.Values.Data;
            romOP.loggedStates{7}.Values.Data = finalOP_Raiz.loggedStates{3}.Values.Data;

            % resetting new_ttf to force resampling times according to IS
            % in ROM
            romOP.loggedStates{21}.Values.Data = 0;

            %! 
            simIn = setInitialState(simIn, romOP);
            simIn = setModelParameter(simIn, "StopTime", string(tsim));
            simIn = setModelParameter(simIn, "SaveFinalState", "on", "SaveOperatingPoint", "on", "FinalStateName", "finalOP");

            % set predicted ts in matlab, IMPORT

            %concatenating initial value + predicted ts of BB
            pred_mdot_c_leak_arr = cat(2, outRoot.mdot_c_leak_signal.Data(end,1), Y_denorm(tree_par_index,:,1));

            %loading BB signal into WB (input)
            %pred_mdot_c_leak = timeseries(pred_mdot_c_leak_arr, time_snapshot_root:0.5:tsim);
            %pred_mdot_c_leak = cat(2, transpose(time_snapshot_root:0.5:tsim), transpose(pred_mdot_c_leak_arr));
            %simIn = setVariable(simIn, 'pred_mdot_c_leak', pred_mdot_c_leak, 'Workspace', name_root);

            %concatenating initial value + predicted ts of BB
            %rec_root_delta_tw_max = max(max([outRoot.tw_signal.Data(1,1,:), mean(outRoot.tw_signal.Data(2:11,1,:)), outRoot.tw_signal.Data(12,1,:)]));
            rec_rom_tw_max = max(outrom.tw_signal.Data(1,1,:));

            delta_tw_max_0 = outRoot.tw_max_signal.Data(end) - rec_rom_tw_max;
            %delta_tw_max_0 = Y_denorm(tree_par_index,1,2);
            pred_delta_tw_max_arr = cat(2, delta_tw_max_0, Y_denorm(tree_par_index,:,2));

            %loading BB signal into WB (output) for error correction
            %pred_delta_tw_max = timeseries(pred_delta_tw_max_arr, time_snapshot_root:0.5:tsim);
            %pred_delta_tw_max = cat(2, transpose(time_snapshot_root:0.5:tsim), transpose(pred_delta_tw_max_arr));
            %simIn = setVariable(simIn, 'pred_delta_tw_max', pred_delta_tw_max, 'Workspace', name_root);

            inport_data = [transpose(time_snapshot_root:t_sam:tsim), pred_mdot_c_leak_arr', pred_delta_tw_max_arr'];
            simIn = simIn.setExternalInput(inport_data);

            %! set sampled ACS velocity
            simIn = setVariable(simIn, 'vel_ACS', sampled_vel_ACS, 'Workspace', name_Wks);

            % RESAMPLING MTTF AND MTTR OF PS1 COMPONENTS 
            % Getting the variables from the model workspace
            mdl_Wks = get_param(name_root,'ModelWorkspace');

            Ncomps = getVariable(mdl_Wks,'Ncomps');
            Lambda_s_resample_lb = getVariable(mdl_Wks,'Lambda_s_resample_lb');
            Lambda_s_resample_ub = getVariable(mdl_Wks,'Lambda_s_resample_ub');
            Mu_s_resample_lb = getVariable(mdl_Wks,'Mu_s_resample_lb');
            Mu_s_resample_ub = getVariable(mdl_Wks,'Mu_s_resample_ub');

            % Get target distributions parameters
            Lambda_t = getVariable(mdl_Wks,'Lambda_t');
            Mu_t = getVariable(mdl_Wks,'Mu_t');

            %sampled_vel_ACS = getVariable(mdlWks,'vel_ACS');
            %disp([string(sampled_vel_ACS), string(tree_par_index)]);

            cell_trans = {};
            Lambda_s = [];
            Mu_s = [];

            % Resampling failure and repair rates
            for i = 1:Ncomps
                %Lambda_s(i) = Lambda_s_resample_lb(i) + (Lambda_s_resample_ub(i) - Lambda_s_resample_lb(i))*rand();
                %Mu_s(i) = Mu_s_resample_lb(i) + (Mu_s_resample_ub(i) - Mu_s_resample_lb(i))*rand();

                %!!! v(10): Multiple importance sampling. three
                %candidates for the importance distribution

                Lambda_s_ensemble = [2*Lambda_t(i) 1E3*Lambda_t(i) 5E5*Lambda_t(i)];
                Mu_s_ensemble = [2*Mu_t(i) 1E1*Mu_t(i) 1E3*Mu_t(i)];

                Lambda_s = Lambda_s_ensemble(randi(3));
                Mu_s = Mu_s_ensemble(randi(3));

                cell_trans{i} = [0 Lambda_s; Mu_s 0];
            end

            % Setting the resampled variables again into the model workspace
            mat_trans = cell2mat(cell_trans);
            mat_trans = reshape(mat_trans,2,2,5);
            
            simIn = setVariable(simIn,'Lambda_s', Lambda_s, 'Workspace', name_root);
            simIn = setVariable(simIn,'Mu_s', Mu_s, 'Workspace', name_root);
            simIn = setVariable(simIn,'mat_trans', mat_trans, 'Workspace', name_root);

            %%
            % start the simulation
            out = sim(simIn);
            nombre = get_param(bdroot, 'Name');

            %close_system('SMDFR_rom_MCDET.slx')

            % Saving information of current branch
            % oper_point = get_param(nombre,'CurrentOperatingPoint')
            init_state = branch{3};
            finalOP = out.finalOP;
            tree_data{end+1} = struct('id',current_node,'data',out,'final_OP',finalOP,'init_state',init_state,'Lambda_t',Lambda_t,'Mu_t', Mu_t);
            
            % Preparing new branch(es) to be stored in the queue
            time_snapshot = finalOP.snapshotTime;
            
            % End condition
            parent_last_state = finalOP_Raiz.loggedStates{14}.Values.Data;
            new_queue_state = finalOP.loggedStates{14}.Values.Data;
            state_var_ps = finalOP.loggedStates{11}.Values.Data;
            
            if time_snapshot < tsim
                % if the system is not demanding anymore the ps:
                if parent_last_state(1) == 1 && new_queue_state(1) == 0
                    % branch is reaching  a terminal condition
                    % do not add new branches
                elseif state_var_ps == 2
                    % branch is reaching a terminal condition
                    % do not add new branches
                else
                    for i = 1:N_resampling
                        % Generate new resamples based on new failure/repair rates
                        k_nodes = k_nodes + 1;
                        new_queue_event = {time_snapshot, k_nodes, new_queue_state};
                        % Adding new branch(es) in the queue list (last in, first out)
                        queue_idle{end+1} = new_queue_event;
                        %compl_branch_event_value = 0;
                        % Adding new branch(es) to the search tree
                        %
                        if mod(i,2) == 1
                            % proceed with the transition of preceding branch
                            compl_branch_event_value = 0;
                        else
                            % mantain the same system state of preceding
                            % branch (complement branch)
                            compl_branch_event_value = 1;
                        end
                        tree = add_node(tree, current_node, k_nodes, compl_branch_event_value);
                        
                        % Updating the internal state of the variables (discarded)
                        %if all(new_queue_state == [1 1])
                        % activate protection system by modifying the internal states of the
                        % last operating point
                        % secondary pump
                        %end
                    end
                end
            end
            
            % Discard the current branch from the queue list (last in, first out)
            queue_idle(1) = [];
        else
            % Detect the parent node
            parent_node = search_parent_node(tree,current_node);
            
            % Obtaining finalOP of the parent node
            parent_finalOP = search_parent_finalOP(tree_data,parent_node);
            
            % Obtaining initial state of the parent node
            parent_init_state = search_parent_init_state(tree_data,parent_node);
            
            % Obtaining final state of the parent node
            %parent_final_state = parent_finalOP.loggedStates{15}.Values.Data;
            parent_final_state = search_final_state(tree_data,parent_node,Ncomps);

            % Setting complementary branch event value, ON parent_finalOP
            compl_branch_event_value = search_node(tree,current_node).compl_branch_event;
            
            %parent_finalOP.loggedStates{21}.Values.Data = compl_branch_event_value;

            load_system('SMDFR_rom_MCDET_v2.slx');
            
            % SET SIMULATION PARAMETERS
            nombre = get_param(bdroot, 'Name');
            % Setting up the next branch to be simulated
            simIn2 = Simulink.SimulationInput(nombre);
            
            % Including state of the parent node in the environment
            % for DET scheduler purposes
            parent_last_state = parent_final_state;
            simIn2 = setVariable(simIn2,'parent_last_state', parent_last_state, 'Workspace', nombre);

            % Including complementary branches
            if compl_branch_event_value == 0
                %display(compl_branch_event_value);
                simIn2 = setVariable(simIn2,'initial_branch_state', parent_last_state, 'Workspace', nombre);
                simIn2 = setVariable(simIn2,'complementary_branch_flag', compl_branch_event_value, 'Workspace', nombre);
                parent_last_time = parent_finalOP.snapshotTime;
                simIn2 = setVariable(simIn2,'parent_last_time', parent_last_time, 'Workspace', nombre);
            else
                %display(compl_branch_event_value);
                if parent_last_state == [1 0]
                    initial_branch_state = [1 1];
                    simIn2 = setVariable(simIn2,'initial_branch_state', initial_branch_state, 'Workspace', nombre);
                    simIn2 = setVariable(simIn2,'complementary_branch_flag', compl_branch_event_value, 'Workspace', nombre);
                elseif parent_last_state == [1 1]
                    initial_branch_state = [1 0];
                    simIn2 = setVariable(simIn2,'initial_branch_state', initial_branch_state, 'Workspace', nombre);
                    simIn2 = setVariable(simIn2,'complementary_branch_flag', compl_branch_event_value, 'Workspace', nombre);
                else % [0 1] and [0 0]: terminal states. do not do anything
                    simIn2 = setVariable(simIn2,'initial_branch_state', parent_last_state, 'Workspace', nombre);
                    simIn2 = setVariable(simIn2,'complementary_branch_flag', 0, 'Workspace', nombre);
                end
                parent_last_time = parent_finalOP.snapshotTime;
                %display(parent_last_time);
                simIn2 = setVariable(simIn2,'parent_last_time', parent_last_time, 'Workspace', nombre);
            end

            % Setting the mission time
            simIn2 = setModelParameter(simIn2, "StopTime", string(tsim));
            
            %!!! UPDATING STATES OF THE CURRENT BRANCH INTO parent_finalOP
            % modifying parent_finalOP
            % branch_state = branch{3};
            % parent_finalOP.loggedStates{10}.Values.Data = [0,1];
            
            % Setting intial state of SimulationInput object
            % Setting step_counters for branching
            parent_finalOP.loggedStates{15}.Values.Data = 0;
            parent_finalOP.loggedStates{23}.Values.Data = 0;
            simIn2 = setInitialState(simIn2,parent_finalOP);

            % Resampling MTTF and MTTR of PS1 components (can be done rdly)
            if current_node >= 1
                % Getting the variables from the model workspace
                mdlWks = get_param(nombre,'ModelWorkspace');

                Ncomps = getVariable(mdlWks,'Ncomps');
                Lambda_s_resample_lb = getVariable(mdlWks,'Lambda_s_resample_lb');
                Lambda_s_resample_ub = getVariable(mdlWks,'Lambda_s_resample_ub');
                Mu_s_resample_lb = getVariable(mdlWks,'Mu_s_resample_lb');
                Mu_s_resample_ub = getVariable(mdlWks,'Mu_s_resample_ub');

                % Get target distributions parameters
                Lambda_t = getVariable(mdlWks,'Lambda_t');
                Mu_t = getVariable(mdlWks,'Mu_t');

                %sampled_vel_ACS = getVariable(mdlWks,'vel_ACS');
                %disp([string(sampled_vel_ACS), string(tree_par_index)]);

                cell_trans = {};
                Lambda_s = [];
                Mu_s = [];

                % Resampling failure and repair rates
                for i = 1:Ncomps
                    %Lambda_s(i) = Lambda_s_resample_lb(i) + (Lambda_s_resample_ub(i) - Lambda_s_resample_lb(i))*rand();
                    %Mu_s(i) = Mu_s_resample_lb(i) + (Mu_s_resample_ub(i) - Mu_s_resample_lb(i))*rand();

                    %!!! v(10): Multiple importance sampling. three
                    %candidates for the importance distribution

                    Lambda_s_ensemble = [2*Lambda_t(i) 1E3*Lambda_t(i) 5E5*Lambda_t(i)];
                    Mu_s_ensemble = [2*Mu_t(i) 1E1*Mu_t(i) 1E3*Mu_t(i)];

                    Lambda_s = Lambda_s_ensemble(randi(3));
                    Mu_s = Mu_s_ensemble(randi(3));

                    cell_trans{i} = [0 Lambda_s; Mu_s 0];
                end

                % Setting the resampled variables again into the model workspace
                mat_trans = cell2mat(cell_trans);
                mat_trans = reshape(mat_trans,2,2,5);

                simIn2 = setVariable(simIn2,'Lambda_s', Lambda_s, 'Workspace', nombre);
                simIn2 = setVariable(simIn2,'Mu_s', Mu_s, 'Workspace', nombre);
                simIn2 = setVariable(simIn2,'mat_trans', mat_trans, 'Workspace', nombre);
                
            end
            
            %! set predicted ts in matlab, INPORT

            %concatenating initial value + predicted ts of BB
            pred_mdot_c_leak_arr = cat(2, outRoot.mdot_c_leak_signal.Data(end,1), Y_denorm(tree_par_index,:,1));

            %loading BB signal into WB (input)
            %pred_mdot_c_leak = timeseries(pred_mdot_c_leak_arr, time_snapshot_root:0.5:tsim);
            %pred_mdot_c_leak = cat(2, transpose(time_snapshot_root:0.5:tsim), transpose(pred_mdot_c_leak_arr));
            %simIn = setVariable(simIn, 'pred_mdot_c_leak', pred_mdot_c_leak, 'Workspace', name_Wks);

            %concatenating initial value + predicted ts of BB
            rec_root_delta_tw_max = max(max([outRoot.tw_signal.Data(1,1,:), mean(outRoot.tw_signal.Data(2:11,1,:)), outRoot.tw_signal.Data(12,1,:)]));
            rec_rom_tw_max = max(outrom.tw_signal.Data(1,1,:));

            delta_tw_max_0 = outRoot.tw_max_signal.Data(end) - rec_delta_tw_max;
            %delta_tw_max_0 = Y_denorm(tree_par_index,1,2);
            pred_delta_tw_max_arr = cat(2, delta_tw_max_0, Y_denorm(tree_par_index,:,2));

            %loading BB signal into WB (output) for error correction
            %pred_delta_tw_max = timeseries(pred_delta_tw_max_arr, time_snapshot_root:0.5:tsim);
            %pred_delta_tw_max = cat(2, transpose(time_snapshot_root:0.5:tsim), transpose(pred_delta_tw_max_arr));
            %simIn = setVariable(simIn, 'pred_delta_tw_max', pred_delta_tw_max, 'Workspace', name_Wks);

            inport_data = [transpose(time_snapshot_root:t_sam:tsim), pred_mdot_c_leak_arr', pred_delta_tw_max_arr'];
            simIn2 = simIn2.setExternalInput(inport_data);

            %! set sampled ACS velocity
            simIn2 = setVariable(simIn2,'vel_ACS', sampled_vel_ACS, 'Workspace', name_Wks);
            
            %sampled_mdot_c_drift = leakage_rate_data{tree_par_index};
            %simIn2 = setVariable(simIn2,'mdot_c_drift', sampled_mdot_c_drift, 'Workspace', name_root);
            %
            simIn2 = setModelParameter(simIn2,"SaveFinalState", "on", "SaveOperatingPoint", "on", "FinalStateName", "finalOP");
            
            %SIMULATE THE CURRENT BRANCH
            out2 = sim(simIn2);

            %close_system('SMDFR_rom_MCDET.slx');
            
            % Saving information of current branch
            % oper_point = get_param(nombre,'CurrentOperatingPoint')
            init_state = branch{3};
            finalOP = out2.finalOP;
            tree_data{end+1} = struct('id',current_node,'data',out2,'final_OP',finalOP,'init_state',init_state,'Lambda_t',Lambda_t,'Mu_t',Mu_t);
            
            % Preparing new branch(es) to be stored in the queue
            time_snapshot = finalOP.snapshotTime;
            
            % End condition
            new_queue_state = finalOP.loggedStates{14}.Values.Data;
            state_var_ps = finalOP.loggedStates{11}.Values.Data;
            
            if time_snapshot < tsim
                % if the system is not demanding anymore the ps:
                if parent_last_state(1) == 1 && new_queue_state(1) == 0
                    % branch is reaching a terminal condition
                    % do not add new branches
                elseif state_var_ps == 2
                    % branch is reaching a terminal condition
                    % do not add new branches
                else
                    for i = 1:N_branches
                    % NOTE: IF THE RESAMPLING OCURR EVERY BRANCHING TIME:
                    % MODIFY ALSO THE TRANSITION STATE TO MODIFY THE TIME TO TRANSITION
                        k_nodes = k_nodes + 1;
                        new_queue_state = finalOP.loggedStates{14}.Values.Data;
                        new_queue_event = {time_snapshot, k_nodes, new_queue_state};
                        % Adding new branch(es) in the queue list (last in, first out)
                        queue_idle{end+1} = new_queue_event;
                        %
                        if mod(i,2) == 1
                            % proceed with the transition of preceding branch
                            compl_branch_event_value = 0;
                        else
                            % mantain the same system state of preceding
                            % branch (complement branch)
                            compl_branch_event_value = 1;
                        end
                        % Adding new branch(es) to the search tree
                        tree = add_node(tree, current_node, k_nodes, compl_branch_event_value);
                        
                        % Updating the internal state of the variables (discarded)
                        % if all(new_queue_state == [1 1])
                        % activate protection system by modifying the internal states of the
                        % last operating point
                        % secondary pump
                        % end
                    end
                end
            end
            
            % Discard the current branch from the queue list (last in, first out)
            queue_idle(1) = [];
        end
    end
    forest_data{tree_par_index} = struct('tree', tree, 'tree_data', tree_data);
    close_system('SMDFR_rom_MCDET_v2.slx');
    %time_snapshot_root = time_snapshot_root + dt_root;
    % END OF OUTER LOOP
end
end_loop = toc(init_loop); 
%snapshots_data{t_index} = struct('tree', tree, 'tree_data', tree_data);
%time_snapshot_root = time_snapshot_root + dt_root;
%end

%%
%save()
%%
% Functions
function tree = add_node(tree, parent_node, value, compl_value)
    new_node = struct('left', parent_node, 'right', [], 'id', value,'compl_branch_event',compl_value);
    % updating parent node
    for i = 1:length(tree)
        current_node = tree{i};
        if current_node.id == parent_node
            child_nodes_parent = current_node.right;
            child_nodes_parent(end+1) = value;
            current_node.right = child_nodes_parent;
            tree{i} = current_node;
        end
        % including new node
    end
    tree{end+1} = new_node;
end

function node = search_node(tree,value)
    for i = 1:length(tree)
        current_node = tree{i};
        if current_node.id == value
            node = current_node;
        end
    end
end

function final_state = search_final_state(tree_data,value,Ncomps)
    node = 0;
    for i = 1:length(tree_data)
        current_node = tree_data{i};
        if current_node.id == value
            node = current_node;
        end
    end
    current_sys_state = node.data.current_sys_state.Data(end,:);
    if all(current_sys_state == ones(Ncomps))
        availability_state = 1;
    else
        availability_state = 0;
    end
    protection_systems = node.data.protection_systems.Data(end);
    final_state = [protection_systems, availability_state];
end

function parent_node = search_parent_node(tree,value)
    for i = 1:length(tree)
        current_node = tree{i};
        if current_node.id == value
            parent_node = current_node.left;
        end
    end
end

function parent_finalOP = search_parent_finalOP(tree_data,parent_node)
    for i = 1:length(tree_data)
        branch_id = tree_data{i}.id;
        if branch_id == parent_node
            parent_finalOP = tree_data{i}.final_OP;
            break
        end
    end
end

function parent_init_state = search_parent_init_state(tree_data,parent_node)
    for i = 1:length(tree_data)
        branch_id = tree_data{i}.id;
        if branch_id == parent_node
            parent_init_state = tree_data{i}.init_state;
            break
        end
    end
end

function X_norm = norm_function(X_input, X_minmax)
    %copy X_input
    X_norm = zeros(size(X_input));
    %dim(X) = [n_samples,sequence_lenght, n_features]
    for i = 1:size(X_input,2)
        X_norm(:,i) = (X_input(:,i) - X_minmax(i,1))/(X_minmax(i,2) - X_minmax(i,1));
    end
end

function Y_denorm = denorm_function(Y_input, Y_minmax)
    %copy Y_input
    Y_denorm = zeros(size(Y_input));
    %dim(Y) = [n_samples,sequence_lenght, n_features]
    for i = 1:size(Y_input,3)
        Y_denorm(:,:,i) = Y_input(:,:,i)*(Y_minmax(i,2) - Y_minmax(i,1)) + Y_minmax(i,1);
    end
end
