%DET_2.m
function [forest_data, forest_roots] = DET_2(det_params, u_trial, d_mpc)  

    persistent isModelLoaded  % Ensure model stays loaded
    
    if isempty(isModelLoaded)
        load_system('rom_MCDET_v2_2022b_control_CMEX.slx');
        set_param('rom_MCDET_v2_2022b_control_CMEX', ...
            'SimulationMode', 'accelerator', ...
            "SaveFinalState", "on",...
            "SaveOperatingPoint", "on",...
            "FinalStateName", "finalOP",...
            'FastRestart', 'on');  % Faster successive runs
        isModelLoaded = true;
    end

    % Unpack necessary parameters
    N_par                   = det_params.N_par;
    M_workers               = det_params.M_workers;
    finalOP_Raiz            = det_params.finalOP_Raiz;
    finalOP2_Raiz           = det_params.finalOP2_Raiz;
    romOP                   = det_params.romOP;
    tsim                    = det_params.tsim;
    time_snapshot_root      = det_params.time_snapshot_root;
    t_sam                   = det_params.t_sam;
    outRoot                 = det_params.outRoot;
    outrom                  = det_params.outrom;
    Y_denorm                = det_params.Y_denorm;
    Y_denorm_50             = det_params.Y_denorm_50;
    sampled_vel_ACS         = det_params.sampled_vel_ACS;
    Lambda_t                = det_params.Lambda_t;
    Mu_t                    = det_params.Mu_t;
    Ncomps                  = det_params.Ncomps;
    time_snapshots          = det_params.time_snapshots;
    snapshot_id             = det_params.snapshot_id;
    N_resampling            = det_params.N_resampling;
    N_branches              = det_params.N_branches;
    name_root               = det_params.name_root;
    name_Wks                = det_params.name_Wks;
    delayed_time_rom        = det_params.delayed_time_rom;
    control_series          = det_params.control_series;

    % Optional DSM safety-margin penalty settings, forwarded to fitness_cost via
    % the root node of tree_data (see the note there): a global cannot be used
    % because the cost is evaluated inside a parfor in luus_jaakola_det.
    if isfield(det_params, 'dsm_thr') && isfield(det_params, 'lambda_dsm')
        dsm_thr_cfg    = det_params.dsm_thr;
        lambda_dsm_cfg = det_params.lambda_dsm;
    else
        dsm_thr_cfg    = [];
        lambda_dsm_cfg = [];
    end

    Raiz_to_Rom             = [1,2,4,5,6,7,3,8,9,10,22,23,24];

    %rng(42,"twister");
    %parfor (tree_par_index = 1:N_par, M_workers)
    for tree_par_index = 1:N_par
        queue_idle = {};
        % structure of the queue: {snapshot_time,id_str,states}, where states = [demand of ps, availability of ps]
        % including first node as current snapshot of the system (i.e., the root)
        %! rom or hrm state?
        %root_state = romOP.loggedStates{13}.Values.Data;
        root_state = finalOP2_Raiz.loggedStates{23}.Values.Data;
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
    
        %! rom or hrm state?
        % copying forest roots for each tree generated on each predicted signal
        %init_state_rom = root_state;
        init_state_rom = finalOP2_Raiz.loggedStates{23}.Values.Data;
    
        forest_roots{tree_par_index} = struct('id', time_snapshot_root - delayed_time_rom, 'data', outrom, 'final_OPRoot', romOP, 'init_state', init_state_rom, 'Lambda_t', Lambda_t, 'Mu_t', Mu_t);
    
        romOP_temp = romOP;
    
        while ~ isempty(queue_idle)
            branch = queue_idle{1};
            current_node = branch{2};
            %display(current_node)
            if current_node == 0
                load_system('rom_MCDET_v2_2022b_control_CMEX.slx');
                simIn = Simulink.SimulationInput('rom_MCDET_v2_2022b_control_CMEX');
    
                % feeding WBM with PO "design parameters?"
                for index = [1:2,6:13]
                    raiz2rom_index = Raiz_to_Rom(index);
                    vec_values = finalOP2_Raiz.loggedStates{raiz2rom_index}.Values.Data;
                    romOP_temp.loggedStates{index}.Values.Data = vec_values;
                end 
                
                % feeding WBM with PO "time-dependent measurable quantities"
                for index = [3,4,5]
                    raiz2rom_index = Raiz_to_Rom(index);
                    vec_values = finalOP2_Raiz.loggedStates{raiz2rom_index}.Values.Data;
                    romOP_temp.loggedStates{index}.Values.Data = mean(vec_values(1:12));
                end 
    
                % resetting new_ttf 'xS_5' to force resampling times according to IS
                % in ROM
                % romOP_temp.loggedStates{20}.Values.Data = 0;
    
                % setting the initial state of the simulation
                %simIn.set_param(name_root,'LoadInitialState','on','InitialState','romOP');
    
                %! 
                simIn = setInitialState(simIn, romOP_temp);
                simIn = setModelParameter(simIn, "StopTime", string(tsim));

                %DEACTIVATED UNDER PERSISTENT SIMULATION
                % simIn = setModelParameter(simIn, 'FastRestart', 'on',...
                %     "SaveFinalState", "on",...
                %     "SaveOperatingPoint", "on",...
                %     "FinalStateName", "finalOP");
                % set parameter 'parent_last_time' in the model workspace
    
                % set predicted ts in matlab, IMPORT
    
                %concatenating initial value + predicted ts of BB
                pred_mdot_c_leak_arr = cat(2, outRoot.mdot_c_leak_signal.Data(end,1), Y_denorm_50(:,1)');

                pred_delta_tw_max_arr = Y_denorm_50(:,2)';

                %! compensating effect of data assimilation / reduction
                %pred_delta_tw_max_comp = [17,11,5,2,0.8,1];
                % for i = 1:5
                %     pred_delta_tw_max_arr(i) = pred_delta_tw_max_arr(i); %- pred_delta_tw_max_comp(i);
                % end

                %TESTING ANTICIPATED ROM INITALIZATION, AVOIDING TRANSIENT
                time_delayed = time_snapshot_root-delayed_time_rom:t_sam:time_snapshot_root-t_sam;
                time_ext = cat(2,time_delayed,time_snapshot_root:t_sam:tsim);
                pred_mdot_c_leak_arr_ext = cat(2, outRoot.mdot_c_leak_signal.Data(end-length(time_delayed):end-1,1)',pred_mdot_c_leak_arr);

                pred_delta_tw_max_delayed = outRoot.tw_max_signal.Data(end-length(time_delayed):end,1) - ...
                                            outrom.tw_max_signal.Data(end-length(time_delayed):end,1);

                pred_delta_tw_max_arr_ext = cat(2,pred_delta_tw_max_delayed',pred_delta_tw_max_arr);

                %delayed taken from control_series matrix
                %control_series = time, u_mpc, active_mpc
                if snapshot_id == 1
                    [control_ind, ~] = find(control_series(:,1) == time_delayed); 
                    u_trial_ext = cat(2, control_series(control_ind,2)', u_trial);
    
                    d_mpc_ext = cat(2, control_series(control_ind,3)', d_mpc);

                else
                    historical_times = time_snapshot_root-delayed_time_rom:t_sam:time_snapshot_root-t_sam;
                    [control_ind, ~] = find(control_series(:,1) == historical_times);
                    u_trial_ext = cat(2, control_series(control_ind,2)', u_trial);

                    d_mpc_ext = cat(2, control_series(control_ind,3)', d_mpc);
                end

                inport_data = [time_ext', pred_mdot_c_leak_arr_ext', pred_delta_tw_max_arr_ext', u_trial_ext', d_mpc_ext'];
                simIn = simIn.setExternalInput(inport_data);
    
                % set sampled ACS velocity
                simIn = setVariable(simIn, 'vel_ACS', sampled_vel_ACS, 'Workspace', name_Wks);
    
                % RESAMPLING MTTF AND MTTR OF PS1 COMPONENTS 
                % Getting the variables from the model workspace
                %mdl_Wks = get_param(name_root,'ModelWorkspace');
    
                %Ncomps = getVariable(mdl_Wks,'Ncomps');
                %Lambda_s_resample_lb = getVariable(mdl_Wks,'Lambda_s_resample_lb');
                %Lambda_s_resample_ub = getVariable(mdl_Wks,'Lambda_s_resample_ub');
                %Mu_s_resample_lb = getVariable(mdl_Wks,'Mu_s_resample_lb');
                %Mu_s_resample_ub = getVariable(mdl_Wks,'Mu_s_resample_ub');
    
                % Get target distributions parameters
                %Lambda_t = getVariable(mdl_Wks,'Lambda_t');
                %Mu_t = getVariable(mdl_Wks,'Mu_t');
    
                %sampled_vel_ACS = getVariable(mdlWks,'vel_ACS');
                %disp([string(sampled_vel_ACS), string(tree_par_index)]);
    
                cell_trans = {};
                Lambda_s = [];
                Mu_s = [];
    
                % Resampling failure and repair rates
                for i = 1:Ncomps
                    %Lambda_s(i) = Lambda_s_resample_lb(i) + (Lambda_s_resample_ub(i) - Lambda_s_resample_lb(i))*rand();
                    %Mu_s(i) = Mu_s_resample_lb(i) + (Mu_s_resample_ub(i) - Mu_s_resample_lb(i))*rand();
    
                    %!!! v(10): Multiple IS ensemble
                    %candidates for the importance distribution
    
                    Lambda_s_ensemble = [Lambda_t(i) Lambda_t(i) Lambda_t(i)];
                    Mu_s_ensemble = [Mu_t(i) Mu_t(i) Mu_t(i)];
    
                    Lambda_s = Lambda_s_ensemble(randi(3));
                    Mu_s = Mu_s_ensemble(randi(3));
    
                    cell_trans{i} = [0 Lambda_s; Mu_s 0];
                end
    
                % Setting the resampled variables again into the model workspace
                mat_trans = cell2mat(cell_trans);
                mat_trans = reshape(mat_trans,2,2,5);
                mat_trans = reshape(mat_trans, 1, 20);
                
                simIn = setVariable(simIn, 'Lambda_s', Lambda_s, 'Workspace', name_root);
                simIn = setVariable(simIn, 'Mu_s', Mu_s, 'Workspace', name_root);
                simIn = setVariable(simIn, 'mat_trans', mat_trans, 'Workspace', name_root);
                simIn = setVariable(simIn, 'parent_last_time', time_snapshots(snapshot_id), 'Workspace', name_root);
                simIn = setVariable(simIn, 'initial_branch_state', finalOP2_Raiz.loggedStates{23}.Values.Data, 'Workspace', name_root);
    
                % Add these configuration options to improve performance:
                %simIn = setModelParameter(simIn, 'SimulationMode', 'accelerator');
                %simIn = simIn.setModelParameter('SaveState', 'on');
                
                % Add these configuration options to improve parsim performance:
                out = sim(simIn);
                % out = parsim(simIn,'ShowProgress', false, ...        % <-- Turn off console progress
                %                    'ShowSimulationManager', 'off' ...% <-- Hide the Simulation Manager GUI
                %                     );

                nombre = 'rom_MCDET_v2_2022b_control_CMEX';
    
                %close_system('SMDFR_rom_MCDET_v2_2022b_control.slx')
    
                % Saving information of current branch
                % oper_point = get_param(nombre,'CurrentOperatingPoint')
                init_state = branch{3};
                finalOP = out.finalOP;
                % dsm_thr / lambda_dsm ride along on the ROOT node so
                % fitness_cost can read the penalty settings even when it is
                % called inside luus_jaakola_det's parfor, where client globals
                % are not visible. Empty when unset -> historical behaviour.
                tree_data{end+1} = struct('id',current_node,'data',out,'final_OP',finalOP,'init_state',init_state,'Lambda_t',Lambda_t,'Mu_t', Mu_t, ...
                    'dsm_thr', dsm_thr_cfg, 'lambda_dsm', lambda_dsm_cfg);
                
                % Preparing new branch(es) to be stored in the queue
                time_snapshot = finalOP.snapshotTime;
                
                % End condition
                parent_last_state = finalOP2_Raiz.loggedStates{23}.Values.Data;
                new_queue_state = finalOP.loggedStates{12}.Values.Data;
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
                %parent_final_state = parent_finalOP.loggedStates{13}.Values.Data;
                parent_final_state = search_final_state(tree_data,parent_node,Ncomps);
    
                % Setting complementary branch event value, ON parent_finalOP
                compl_branch_event_value = search_node(tree,current_node).compl_branch_event;
                
                %parent_finalOP.loggedStates{21}.Values.Data = compl_branch_event_value;
    
                load_system('rom_MCDET_v2_2022b_control_CMEX.slx');
                
                % SET SIMULATION PARAMETERS
                nombre = 'rom_MCDET_v2_2022b_control_CMEX';
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
                    % FIX 2026-09-03: use isequal for the [demand_ps, ava_ps] pair.
                    % `if parent_last_state == [1 0]` compares element-wise and an
                    % `if` over a vector requires ALL elements true, so [1 0] vs
                    % [1 1] silently fell through to the next test. It happened to
                    % reach the intended branch here, but [0 0] vs [0 1] and any
                    % future state pair would mis-route. The ACS branch semantics
                    % (activate or not; if active, available or failed) depend on
                    % this pair being matched exactly.
                    if isequal(parent_last_state(:)', [1 0])
                        initial_branch_state = [1 1];
                        simIn2 = setVariable(simIn2,'initial_branch_state', initial_branch_state, 'Workspace', nombre);
                        simIn2 = setVariable(simIn2,'complementary_branch_flag', compl_branch_event_value, 'Workspace', nombre);
                    elseif isequal(parent_last_state(:)', [1 1])
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
                
                %!!! Updating states of curreng branch into parent_finalOP
                % modifying parent_finalOP
                % branch_state = branch{3};
                % parent_finalOP.loggedStates{13}.Values.Data = [0,1];
                
                % Setting intial state of SimulationInput object
                % Setting step_counters for branching
                parent_finalOP.loggedStates{13}.Values.Data = 0;
                %parent_finalOP.loggedStates{16}.Values.Data = 0;
    
                %
                simIn2 = setInitialState(simIn2,parent_finalOP);
    
                % Resampling MTTF and MTTR of PS1 components (can be done rdly)
                if current_node >= 1
                    % Getting the variables from the model workspace
                    %mdlWks = get_param(nombre,'ModelWorkspace');
    
                    %Ncomps = getVariable(mdlWks,'Ncomps');
                    %Lambda_s_resample_lb = getVariable(mdlWks,'Lambda_s_resample_lb');
                    %Lambda_s_resample_ub = getVariable(mdlWks,'Lambda_s_resample_ub');
                    %Mu_s_resample_lb = getVariable(mdlWks,'Mu_s_resample_lb');
                    %Mu_s_resample_ub = getVariable(mdlWks,'Mu_s_resample_ub');
    
                    % Get target distributions parameters
                    %Lambda_t = getVariable(mdlWks,'Lambda_t');
                    %Mu_t = getVariable(mdlWks,'Mu_t');
    
                    %sampled_vel_ACS = getVariable(mdlWks,'vel_ACS');
                    %disp([string(sampled_vel_ACS), string(tree_par_index)]);
    
                    cell_trans = {};
                    Lambda_s = [];
                    Mu_s = [];
    
                    % Resampling failure and repair rates
                    for i = 1:Ncomps
                        %Lambda_s(i) = Lambda_s_resample_lb(i) + (Lambda_s_resample_ub(i) - Lambda_s_resample_lb(i))*rand();
                        %Mu_s(i) = Mu_s_resample_lb(i) + (Mu_s_resample_ub(i) - Mu_s_resample_lb(i))*rand();
    
                        %!!! v(10): Multiple IS ensemble
                        %candidates for the importance distribution
    
                        Lambda_s_ensemble = [Lambda_t(i) Lambda_t(i) Lambda_t(i)];
                        Mu_s_ensemble = [Mu_t(i) Mu_t(i) Mu_t(i)];
    
                        Lambda_s = Lambda_s_ensemble(randi(3));
                        Mu_s = Mu_s_ensemble(randi(3));
    
                        cell_trans{i} = [0 Lambda_s; Mu_s 0];
                    end
    
                    % Setting the resampled variables again into the model workspace
                    mat_trans = cell2mat(cell_trans);
                    mat_trans = reshape(mat_trans,2,2,5);
                    mat_trans = reshape(mat_trans, 1, 20);
    
                    simIn2 = setVariable(simIn2,'Lambda_s', Lambda_s, 'Workspace', nombre);
                    simIn2 = setVariable(simIn2,'Mu_s', Mu_s, 'Workspace', nombre);
                    simIn2 = setVariable(simIn2,'mat_trans', mat_trans, 'Workspace', nombre);
                    
                end
                
                %! set predicted ts in matlab, INPORT
                %concatenating initial value + predicted ts of BB
                pred_mdot_c_leak_arr = cat(2, outRoot.mdot_c_leak_signal.Data(end,1), Y_denorm_50(:,1)');
                %pred_mdot_c_leak_arr = cat(2, outRoot.mdot_c_leak_signal.Data(end,1), Y_denorm(tree_par_index,:,1));
    
                %TESTING ANTICIPATED ROM INITALIZATION, AVOIDING TRANSIENT
                time_delayed = time_snapshot_root-delayed_time_rom:t_sam:time_snapshot_root-t_sam;
                time_ext = cat(2,time_delayed,time_snapshot_root:t_sam:tsim);
                pred_mdot_c_leak_arr_ext = cat(2, outRoot.mdot_c_leak_signal.Data(end-length(time_delayed):end-1,1)',pred_mdot_c_leak_arr);

                pred_delta_tw_max_delayed = outRoot.tw_max_signal.Data(end-length(time_delayed):end,1) - ...
                                            outrom.tw_max_signal.Data(end-length(time_delayed):end,1);

                pred_delta_tw_max_arr_ext = cat(2,pred_delta_tw_max_delayed',pred_delta_tw_max_arr);

                %delayed taken from control_series matrix
                %control_series = time, u_mpc, active_mpc
                if snapshot_id == 1
                    [control_ind, ~] = find(control_series(:,1) == time_delayed); 
                    u_trial_ext = cat(2, control_series(control_ind,2)', u_trial);
    
                    d_mpc_ext = cat(2, control_series(control_ind,3)', d_mpc);

                else
                    historical_times = time_snapshot_root-delayed_time_rom:t_sam:time_snapshot_root-t_sam;
                    [control_ind, ~] = find(control_series(:,1) == historical_times);
                    u_trial_ext = cat(2, control_series(control_ind,2)', u_trial);

                    d_mpc_ext = cat(2, control_series(control_ind,3)', d_mpc);
                end

                inport_data = [time_ext', pred_mdot_c_leak_arr_ext', pred_delta_tw_max_arr_ext', u_trial_ext', d_mpc_ext'];
                simIn2 = simIn2.setExternalInput(inport_data);
    
                % set sampled ACS velocity
                simIn2 = setVariable(simIn2,'vel_ACS', sampled_vel_ACS, 'Workspace', name_Wks);
                
                %sampled_mdot_c_drift = leakage_rate_data{tree_par_index};
                %simIn2 = setVariable(simIn2,'mdot_c_drift', sampled_mdot_c_drift, 'Workspace', name_root);
                
                %DEACTIVATED UNDER PERSISTENT SIMULATION
                %simIn2 = setModelParameter(simIn2,"SaveFinalState", "on", "SaveOperatingPoint", "on", "FinalStateName", "finalOP");
                
                % FIX 2026-09-03: this call was commented out while out2 is
                % consumed below (line ~462) and saved into tree_data, so any
                % child node that DID get queued would fail at the save step with
                % an undefined-variable error. Verified by
                % diagnose_child_branch_1170.m: the child branch simulates
                % cleanly (~10-26 s) and yields a valid finalOP.
                %
                % Scope note: this is a latent defect, not the reason the tree
                % often holds only the root. Expansion is gated upstream by
                % `time_snapshot < tsim` (line ~243/472), which is the "was this
                % run cut short by an ACS event?" test. When no transition time
                % falls inside the prediction window, no child is queued and this
                % branch never executes -- that is the intended behaviour.
                out2 = sim(simIn2);
                % Alternative for parallel execution:
                % out2 = parsim(simIn2,'ShowProgress', false, ...        % <-- Turn off console progress
                %                    'ShowSimulationManager', 'off' ...% <-- Hide the Simulation Manager GUI
                %                     );

                nombre = 'rom_MCDET_v2_2022b_control_CMEX';
    
                %close_system('SMDFR_rom_MCDET_control.slx');
                
                % Saving information of current branch
                % oper_point = get_param(nombre,'CurrentOperatingPoint')
                init_state = branch{3};
                finalOP = out2.finalOP;
                % same extra fields as the root node, so tree_data stays a
                % homogeneous struct array
                tree_data{end+1} = struct('id',current_node,'data',out2,'final_OP',finalOP,'init_state',init_state,'Lambda_t',Lambda_t,'Mu_t',Mu_t, ...
                    'dsm_thr', dsm_thr_cfg, 'lambda_dsm', lambda_dsm_cfg);
                
                % Preparing new branch(es) to be stored in the queue
                time_snapshot = finalOP.snapshotTime;
                
                % End condition
                new_queue_state = finalOP.loggedStates{12}.Values.Data;
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
                            new_queue_state = finalOP.loggedStates{12}.Values.Data;
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
        %close_system('rom_MCDET_v2_2022b_control_CMEX.slx');
        %time_snapshot_root = time_snapshot_root + dt_root;
        % END OF OUTER LOOP
    end
end