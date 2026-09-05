%% Analysis - initial nodes

Ncomps = 5; 

root_nodes = {};

for i = 1:length(forest_roots)
    node_data_extrc = forest_roots{i};
    root_nodes{i} = node_data_extrc;
end

%% PROBABILISTIC SAFETY MARGINS
% Calculating maximum values 

max_parents = zeros(1,1);
max_values = zeros(1,1);
max_indexes = zeros(1,1);

max_matrices = {};
for i = 1:length(forest_data)
    max_matrices{i} = zeros(1,2);
end

for k = 1:length(forest_data)
    forest_k = forest_data{k};
    tree_data_unzip = {};
    tree_unzip = {};
    % unzipping forest one by one
    for l = 1:length(forest_k)
        tree_data_unzip{end+1} = forest_data{k}(l).tree_data;
        tree_unzip{end+1} = forest_data{k}(l).tree;
    end
    %!
    root_data_unzip = forest_roots{k};
    terminal_nodes = search_terminal_nodes(tree_unzip);

    % for each forest, unzipping terminal nodes
    for i = 1:length(terminal_nodes)
        parents_list = search_parents(tree_unzip,terminal_nodes{i});
        branch_list = {};
        branch_list{end+1} = terminal_nodes{i};
        for j = 1:length(parents_list)
            %!!! REMOVING PARENT NODE FOR MAXIMUM
            %if parents_list{j} ~= 0
                branch_list{end+1} = parents_list{j};
            %end
        end
        max_parent_value = [];
        max_parent_index = [];
        for j = 1:length(branch_list)
            node_data = search_tree_data(tree_data_unzip,branch_list{j});
            %!!! RELEASING FIRST PEAKED VALUE
            %if branch_list{j} == 0
                %pred_twmax_signal = node_data.data.pred_tw_max_signal.Data(:,:);
            %else
                pred_twmax_signal = node_data.data.pred_tw_max_signal.Data(:,:);
            %end
            [max_r, max_c] = find(pred_twmax_signal == max(pred_twmax_signal));
            max_parent_value(j) = pred_twmax_signal(max_r(1));
            %display(max_values);
            max_parent_index(j) = max_r(1);
        end
        max_twmax_index = find(max_parent_value == max(max_parent_value));
        max_twmax_i = max_parent_value(max_twmax_index(1));
        max_parent_i = branch_list(max_twmax_index(1));
        max_parent_index_i = max_parent_index(max_twmax_index(1));
        %storing
        max_parents(end+1,1) = max_parent_i{1};
        max_values(end+1,1) = max_twmax_i;
        max_indexes(end+1,1) = max_parent_index_i;

        %for later
        max_matrices{k}(end+1,:) = [terminal_nodes{i}; max_twmax_i];
    end
end

max_parents(1,:) = []; max_values(1,:) = []; max_indexes(1,:) = []; 

for i = 1:length(forest_data)
    max_matrices{i}(1,:) = [];
end

%% Saving max_matrices 
% format_string_1 = 'data__max_matrices_%d.xls';
% string_1 = sprintf(format_string_1,ddbb_id);
% 
% for i = 1:length(max_matrices)
%     table = array2table(max_matrices{i},'VariableNames',{'Terminal Nodes', 'P_(N,max)'});
%     writetable(table, string_1,'Sheet',i)
% end

%% Saving parent nodes

% for k = 1:length(forest_data)
%     forest_k = forest_data{k};
%     tree_data_unzip = {};
%     tree_unzip = {};
%     cell_table = {};
%     % unzipping forest one by one
%     for l = 1:length(forest_k)
%         tree_data_unzip{end+1} = forest_data{k}(l).tree_data;
%         tree_unzip{end+1} = forest_data{k}(l).tree;
%     end
%     root_data_unzip = forest_roots{k};
%     terminal_nodes = search_terminal_nodes(tree_unzip);
% 
%     % for each forest, unzipping terminal nodes
%     for i = 1:length(terminal_nodes)
%         parents_list = search_parents(tree_unzip,terminal_nodes{i});
%         branch_list = {};
%         %branch_list{end+1} = terminal_nodes{i};
%         for j = 1:length(parents_list)
%             branch_list{end+1} = parents_list{j};
%         end
%         cell_table{i,1} = terminal_nodes{i};
%         cell_table{i,2} = string(branch_list);
%     end
% table = cell2table(cell_table);
% writetable(table,'data__parent_nodes.xls','Sheet',k);
% end

%% RESULTS - PSMs WITH BRUTE MC SAMPLING
% display(max_values);
% display(max_indexes);
% 
% % Using bracketing OS approach
% max_matrices_total = cat(1,max_matrices{:});
% sorted_samples = sort(max_matrices_total(:,2));
% yhat_max = sorted_samples(end);
% 
% % For upper safety margins
% N_b = length(max_matrices_total);
% % gamma_1: estimated 95-th percentile of maximum values
% gamma_1 = 0.95;
% 
% % U_y: upper threshold
% U_y = 1.015;
% 
% % y_ref: reference value
% y_ref = 1;
% 
% % calculating the upper safety margin
% [M_U_ya, beta_1] = SM_Up(U_y, y_ref, yhat_max, N_b, gamma_1);
% display(M_U_ya);
% display(beta_1);

%% Adjusted PSMs with Importance Sampling

%weight_nodes: structure
%weight_nodes{i} = [node_index IS_weight_per_node P_N,max_per_node P{IS}_{N,max}_per_node]
%what if we just count every time we weight the node? (to compare)
weight_nodes = {};
%weight_nodes{i} = [node_index, IS_weight, IS_weight_norm]
for i = 1:length(forest_data)
    weight_nodes{i} = zeros(1,8);
end

% storing max IS weight value per node
weight_matrices = {};
for i = 1:length(forest_data)
    weight_matrices{i} = zeros(1,4);
end

for l = 1:length(forest_data)
    forest_l = forest_data{l};
    tree_data_unzip = {};
    tree_unzip = {};
    % unzipping forest one by one
    for m = 1:length(forest_l)
        tree_data_unzip{end+1} = forest_data{l}(m).tree_data;
        tree_unzip{end+1} = forest_data{l}(m).tree;
    end
    %root_data_unzip = forest_roots{l};
    terminal_nodes = search_terminal_nodes(tree_unzip);
    romRoot_init_state = forest_l(1).tree_data.init_state;

    if all(romRoot_init_state == [0 1]) || all(romRoot_init_state == [0 0])
        %PS not requested. assume hot sb and include node 0 and
        %the previous simulation time up to requested PS instant?
        for i = 1:length(terminal_nodes)
            parents_list = search_parents(tree_unzip,terminal_nodes{i});
            %
            parents_list{end+1} = terminal_nodes{i};
            %parents_list = cell2mat(parents_list); parents_list = sort(parents_list);
            counter_probs_sample = [];
            counter_probs_target = [];
            for j = 1:length(parents_list)
                %display(parents_list{j});
                parent_node = parents_list{j};
                % Note: not considering the state of ps before activation 
                % i.e., skipping root node
                if parent_node ~= 0
                    parent_data = search_tree_data(tree_data_unzip,parent_node);
                    total_time = parent_data.data.tout;
                    t_counter = 1;
                    t = total_time(t_counter);
                    t_init = total_time(1);
                    t_finish = total_time(end-1);
                    %t_finish = total_time(end-2);
                    t_next_event = 0;
    
                    while t < t_finish 
                        trans_sys = parent_data.data.trans_sys.Data(t_counter+1,:);
                        mat_ttf = parent_data.data.mat_ttf.Data(:,:,t_counter+1);
                        trans_c = parent_data.data.trans_c.Data(t_counter+1);
                        % updating t_next_event
                        t_next_event = parent_data.data.new_ttf.Data(t_counter+1);
                        % Note: not considering the state of ps before activation 
                        % Classifying the transition
                        if any(trans_sys > 1) %working
                            mat_ttf_c = mat_ttf(trans_c,:);
                            [row_c, col_c] = find(mat_ttf_c == min(mat_ttf_c(mat_ttf_c > 0)));
                            new_ttf = mat_ttf_c(col_c);
                            t_next_event = new_ttf;
                            previous_state = ones(Ncomps,1); 
                            prob_vec_sample = zeros(1,Ncomps);
                            prob_vec_target = zeros(1,Ncomps);
                            mat_trans_sample = parent_data.data.mat_trans.Data;
                            for k = 1:Ncomps
                                mat_ttf_k = mat_ttf(k,:);
                                [row_k,trans_state_k] = find(mat_ttf_k == min(mat_ttf_k(mat_ttf_k > 0)));
                                transition_rate_k_sample = mat_trans_sample(previous_state(k),trans_state_k,k);
                                transition_rate_k_target = parent_data.Lambda_t(k);
                                delta_t = min(new_ttf,t_finish) - t;
                                if k ~= trans_c
                                    prob_vec_sample(1,k) = 1 - CDF(transition_rate_k_sample, delta_t);
                                    prob_vec_target(1,k) = 1 - CDF(transition_rate_k_target, delta_t);
                                    %prob_vec_sample(1,k) = PDF(transition_rate_k_sample, delta_t);
                                    %prob_vec_target(1,k) = PDF(transition_rate_k_target, delta_t);
                                else
                                    if new_ttf <= tsim
                                        prob_vec_sample(1,k) = PDF(transition_rate_k_sample, delta_t);
                                        prob_vec_target(1,k) = PDF(transition_rate_k_target, delta_t);
                                    else
                                        prob_vec_sample(1,k) = 1 - CDF(transition_rate_k_sample, delta_t);
                                        prob_vec_target(1,k) = 1 - CDF(transition_rate_k_target, delta_t);
                                    end
    
                                    % hot stand-by assumption.
                                    % node 0 does not fail.
                                    % not failed, not inspected,
                                    % starts at t_init. added to prob of
                                    % child node
                                    if search_parent_node(tree_unzip,parent_node) == 0 %
                                        %prob_vec_sample(1,k) = prob_vec_sample(1,k) * CDF(transition_rate_k_sample, t_init);
                                        %prob_vec_target(1,k) = prob_vec_target(1,k) * CDF(transition_rate_k_target, t_init);
                                        prob_node0_target = 1 - CDF(sum(parent_data.Lambda_t),t_init);
                                        prob_vec_target(1,k) = prob_vec_target(1,k) * prob_node0_target;
                                        prob_node0_sample = 1 - CDF(sum(mat_trans_sample(1,2,:)),t_init);
                                        prob_vec_sample(1,k) = prob_vec_sample(1,k) * prob_node0_sample;
                                    end
                                end
                            end
                            counter_probs_sample(end+1) = prod(prob_vec_sample);
                            counter_probs_target(end+1) = prod(prob_vec_target);
                        else %under repair
                            mat_ttf_c = mat_ttf(trans_c,:);
                            [row_c, col_c] = find(mat_ttf_c == min(mat_ttf_c(mat_ttf_c > 0)));
                            new_ttf = mat_ttf_c(col_c);
                            t_next_event = new_ttf;
                            previous_state = ones(Ncomps,1); 
                            previous_state(trans_c) = 2;
                            prob_vec_sample = zeros(1,Ncomps);
                            prob_vec_target = zeros(1,Ncomps);
                            mat_trans_sample = parent_data.data.mat_trans.Data;
                            for k = 1:Ncomps
                                transition_rate_k_sample = mat_trans_sample(previous_state(k),trans_sys(k),k);
                                transition_rate_k_target = parent_data.Mu_t(k);
                                delta_t = min(new_ttf,t_finish) - t;
                                if k ~= trans_c
                                    %hot standby?
                                    prob_vec_sample(1,k) = 1 - CDF(transition_rate_k_sample, delta_t);
                                    prob_vec_target(1,k) = 1 - CDF(transition_rate_k_target, delta_t);
                                    %prob_vec_sample(1,k) = 1;
                                    %prob_vec_target(1,k) = 1;
                                else
                                    if new_ttf <= tsim
                                        prob_vec_sample(1,k) = PDF(transition_rate_k_sample, delta_t);
                                        prob_vec_target(1,k) = PDF(transition_rate_k_target, delta_t);
                                    else
                                        prob_vec_sample(1,k) = 1 - CDF(transition_rate_k_sample, delta_t);
                                        prob_vec_target(1,k) = 1 - CDF(transition_rate_k_target, delta_t);
                                    end

                                    % hot stand-by assumption. only failed node right after
                                    % node 0 includes CDF(trans_c).
                                    % Failed, not inspected, intervention
                                    % starts at t_init, added to the child
                                    % node
                                    if search_parent_node(tree_unzip,parent_node) == 0 %
                                        %prob_vec_sample(1,k) = prob_vec_sample(1,k) * CDF(transition_rate_k_sample, t_init);
                                        %prob_vec_target(1,k) = prob_vec_target(1,k) * CDF(transition_rate_k_target, t_init);
                                        prob_node0_target = (parent_data.Lambda_t(k)/(sum(parent_data.Lambda_t)))*CDF(sum(parent_data.Lambda_t),t_init);
                                        prob_vec_target(1,k) = prob_vec_target(1,k) * prob_node0_target;
                                        prob_node0_sample = (mat_trans_sample(1,2,k)/(sum(mat_trans_sample(1,2,:))))*CDF(sum(mat_trans_sample(1,2,:)),t_init);
                                        prob_vec_sample(1,k) = prob_vec_sample(1,k) * prob_node0_sample;
                                    end
                                end
                            end
                            counter_probs_sample(end+1) = prod(prob_vec_sample);
                            counter_probs_target(end+1) = prod(prob_vec_target);
                        end
                        t = min(t_next_event,t_finish);
                        [time_r, time_c] = find(total_time > t);
                        t_counter = time_r(1);
                    end
                    parent_data = search_tree_data(tree_data_unzip,parent_node);
                    twmax_node = max(parent_data.data.pred_tw_max_signal.Data);
                    weight_node = prod(prob_vec_target)/prod(prob_vec_sample);
                    weight_nodes{l}(end+1,:) = [parent_node, terminal_nodes{i}, prod(prob_vec_target), prod(prob_vec_sample), weight_node, 0, twmax_node, 0];
                % else: probability of ocurrence = 1 (no prot system activated)
                end
            end
            %display(counter_probs_sample)
            %display(counter_probs_target)
    
            % storing UN-normalized weights for importance sampling
            %weight_nodes(end+1,:) = [prod(counter_probs_sample) prod(counter_probs_target)];
    
            % storing UN-normalized weights for importance sampling
            weight_matrices{l}(end+1,:) = [terminal_nodes{i}; prod(counter_probs_sample); prod(counter_probs_target); prod(counter_probs_target)/prod(counter_probs_sample)];
        end
    elseif all(romRoot_init_state == [1 1])
        %PS requested and available. only include node 0
        %(memoryless?)
        for i = 1:length(terminal_nodes)
            parents_list = search_parents(tree_unzip,terminal_nodes{i});
            %
            parents_list{end+1} = terminal_nodes{i};
            %parents_list = cell2mat(parents_list); parents_list = sort(parents_list);
            counter_probs_sample = [];
            counter_probs_target = [];
            for j = 1:length(parents_list)
                %display(parents_list{j});
                parent_node = parents_list{j};
                % Note: not considering the state of ps before activation 
                % i.e., including root node for IS weighting
                %if parent_node ~= 0
                    parent_data = search_tree_data(tree_data_unzip,parent_node);
                    total_time = parent_data.data.tout;
                    t_counter = 1;
                    t = total_time(t_counter);
                    t_init = total_time(1);
                    %t_finish = total_time(end-1);
                    t_finish = total_time(end-2);
                    t_next_event = 0;
    
                    while t < t_finish 
                        trans_sys = parent_data.data.trans_sys.Data(t_counter+2,:);
                        mat_ttf = parent_data.data.mat_ttf.Data(:,:,t_counter+2);
                        trans_c = parent_data.data.trans_c.Data(t_counter+2);
                        % updating t_next_event
                        t_next_event = parent_data.data.new_ttf.Data(t_counter+2);
                        % Note: not considering the state of ps before activation 
                        % Classifying the transition
                        if any(trans_sys > 1) %working
                            mat_ttf_c = mat_ttf(trans_c,:);
                            [row_c, col_c] = find(mat_ttf_c == min(mat_ttf_c(mat_ttf_c > 0)));
                            new_ttf = mat_ttf_c(col_c);
                            t_next_event = new_ttf;
                            previous_state = ones(Ncomps,1); 
                            prob_vec_sample = zeros(1,Ncomps);
                            prob_vec_target = zeros(1,Ncomps);
                            mat_trans_sample = parent_data.data.mat_trans.Data;
                            for k = 1:Ncomps
                                mat_ttf_k = mat_ttf(k,:);
                                [row_k,trans_state_k] = find(mat_ttf_k == min(mat_ttf_k(mat_ttf_k > 0)));
                                transition_rate_k_sample = mat_trans_sample(previous_state(k),trans_state_k,k);
                                transition_rate_k_target = parent_data.Lambda_t(k);
                                delta_t = min(new_ttf,t_finish) - t;
                                if k ~= trans_c
                                    prob_vec_sample(1,k) = 1 - CDF(transition_rate_k_sample, delta_t);
                                    prob_vec_target(1,k) = 1 - CDF(transition_rate_k_target, delta_t);
                                    %prob_vec_sample(1,k) = PDF(transition_rate_k_sample, delta_t);
                                    %prob_vec_target(1,k) = PDF(transition_rate_k_target, delta_t);
                                else
                                    if new_ttf <= tsim %all but j-th are right-censored
                                        prob_vec_sample(1,k) = PDF(transition_rate_k_sample, delta_t);
                                        prob_vec_target(1,k) = PDF(transition_rate_k_target, delta_t);
                                    else %all comps are right-censored
                                        prob_vec_sample(1,k) = 1 - CDF(transition_rate_k_sample, delta_t);
                                        prob_vec_target(1,k) = 1 - CDF(transition_rate_k_target, delta_t);
                                    end
                                end
                            end
                            counter_probs_sample(end+1) = prod(prob_vec_sample);
                            counter_probs_target(end+1) = prod(prob_vec_target);
                        else %under repair
                            mat_ttf_c = mat_ttf(trans_c,:);
                            [row_c, col_c] = find(mat_ttf_c == min(mat_ttf_c(mat_ttf_c > 0)));
                            new_ttf = mat_ttf_c(col_c);
                            t_next_event = new_ttf;
                            previous_state = ones(Ncomps,1); 
                            previous_state(trans_c) = 2;
                            prob_vec_sample = zeros(1,Ncomps);
                            prob_vec_target = zeros(1,Ncomps);
                            mat_trans_sample = parent_data.data.mat_trans.Data;
                            for k = 1:Ncomps
                                transition_rate_k_sample = mat_trans_sample(previous_state(k),trans_sys(k),k);
                                transition_rate_k_target = parent_data.Mu_t(k);
                                delta_t = min(new_ttf,t_finish) - t;
                                if k ~= trans_c 
                                    %hot standby?
                                    prob_vec_sample(1,k) = 1 - CDF(transition_rate_k_sample, delta_t);
                                    prob_vec_target(1,k) = 1 - CDF(transition_rate_k_target, delta_t);
                                    %prob_vec_sample(1,k) = 1;
                                    %prob_vec_target(1,k) = 1;
                                else % transitioned component
                                    if new_ttf <= tsim %all but j-th are right-censored
                                        prob_vec_sample(1,k) = PDF(transition_rate_k_sample, delta_t);
                                        prob_vec_target(1,k) = PDF(transition_rate_k_target, delta_t);
                                    else  %all comps are right-censored
                                        prob_vec_sample(1,k) = 1 - CDF(transition_rate_k_sample, delta_t);
                                        prob_vec_target(1,k) = 1 - CDF(transition_rate_k_target, delta_t);
                                    end
                                end
                            end
                            counter_probs_sample(end+1) = prod(prob_vec_sample);
                            counter_probs_target(end+1) = prod(prob_vec_target);
                        end
                        t = min(t_next_event,t_finish);
                        [time_r, time_c] = find(total_time > t);
                        t_counter = time_r(1);
                    end
                    parent_data = search_tree_data(tree_data_unzip,parent_node);
                    twmax_node = max(parent_data.data.pred_tw_max_signal.Data);
                    weight_node = prod(prob_vec_target)/prod(prob_vec_sample);
                    weight_nodes{l}(end+1,:) = [parent_node, terminal_nodes{i}, prod(prob_vec_target), prod(prob_vec_sample), weight_node, 0, twmax_node, 0];
                % else: probability of ocurrence = 1 (no prot system activated)
                %end
            end
            %display(counter_probs_sample)
            %display(counter_probs_target)
    
            % storing UN-normalized weights for importance sampling
            %weight_nodes(end+1,:) = [prod(counter_probs_sample) prod(counter_probs_target)];
    
            % storing UN-normalized weights for importance sampling
            weight_matrices{l}(end+1,:) = [terminal_nodes{i}; prod(counter_probs_sample); prod(counter_probs_target); prod(counter_probs_target)/prod(counter_probs_sample)];
        end
    %elseif all(romRoot_init_state == [1 0])
        %nothing so far. unlikely scenario for the HRM
    end

end

%weight_nodes(1,:) = [];
for l = 1:length(forest_data)
    weight_matrices{l}(1,:) = [];
    weight_nodes{l}(1,:) = [];
end

%% Normalizing weights (per set of parent nodes, for each terminal branch)
%weight_norm_matrices_v2 = {};
% max_matrices_IS_v2 = {};
% 
% for l = 1:length(forest_data)
%     %finding terminal nodes: ia
%     [C,ia] = unique(weight_nodes{l}(:,2));
%     terminal_nodess = weight_nodes{l}(ia,2);
%     weight_norm_matrix_v2 = zeros(length(ia),2);
%     max_matrix_IS_v2 = zeros(length(ia),2);
%     k = 1;
%     for t_node = terminal_nodess'
%         % finding terminal nodes
%         [rows_terminal, col_terminal] = find(weight_nodes{l}(:,2) == t_node);
%         % getting the rows of parent nodes
%         weight_node_temp = weight_nodes{l}(rows_terminal,:);
%         weight_nodes_norm_vec = weight_node_temp(:,5) ./ sum(weight_node_temp(:,5));
%         % passing the new normalized weights to the original matrix
%         weight_nodes{l}(rows_terminal,6) = weight_nodes_norm_vec;
%         twmax_node_IS = weight_node_temp(:,7) .* weight_nodes_norm_vec;
%         % passing the IS (with norm weights) weighted P_norm to the original matrix
%         weight_nodes{l}(rows_terminal,8) = twmax_node_IS;
% 
%         % getting maximum norm power adjusted by IS, for all parent nodes
%         twmax_IS_v2 = max(twmax_node_IS);
%         max_matrix_IS_v2(k,:) = [t_node, twmax_IS_v2];
%         k = k + 1;
%     end
%     max_matrices_IS_v2{l} = max_matrix_IS_v2;
% end

%% RESULTS - CORRECTED BY IS

%% Normalizing the weights, FOR EACH TREE!

% weight_norm_matrices = {};
% 
% for i = 1:length(forest_data)
%     weight_i = weight_matrices{i}(:,4);
%     % Normalizing the weights
%     weight_i_norm = weight_i ./ sum(weight_i);
%     % Not normalizing
%     %weight_i_norm = weight_i;
%     weight_norm_matrices{i} = horzcat(weight_matrices{i}(:,1), weight_i_norm);
% end

%% Saving weight_norm_matrices
% format_string_2 = 'data__weight_norm_matrices_%d.xls';
% string_2 = sprintf(format_string_2,ddbb_id);
% 
% for i = 1:length(weight_norm_matrices)
%     table = array2table(weight_norm_matrices{i},'VariableNames',{'Terminal Nodes', 'w_(i)'});
%     writetable(table, string_2, 'Sheet',i)
% end

%% Assigning max weight value per parent node

IS_prec_matrices = {};
IS_maxweight_vec = {};
for i = 1:length(forest_data)
    IS_prec_matrices{i} = -1*ones(length(forest_data{i}),length(forest_data{i}));
end

for l = 1:length(forest_data)
    forest_l = forest_data{l};
    tree_data_unzip = {};
    tree_unzip = {};
    % unzipping forest one by one
    for m = 1:length(forest_l)
        tree_data_unzip{end+1} = forest_data{l}(m).tree_data;
        tree_unzip{end+1} = forest_data{l}(m).tree;
    end
    root_data_unzip = forest_roots{l};
    terminal_nodes = search_terminal_nodes(tree_unzip);
    %terminal_nodes_vec = cell2mat(terminal_nodes)';
    for i = 1:length(terminal_nodes)
        % computing IS values for an squared matrix with the nodes
        parents_list = search_parents(tree_unzip,terminal_nodes{i});
        parents_list{end+1} = terminal_nodes{i};
        parents_list_vec = cell2mat(parents_list)';
        col_IS_matrix = zeros(terminal_nodes{i}+1, 1);
        col_IS_matrix(parents_list_vec + 1, 1) = weight_matrices{l}(i,4);
        IS_prec_matrices{l}(1:terminal_nodes{i}+1,terminal_nodes{i}+1) = col_IS_matrix;
    end
end

% vector of maximum weights per parent node, per tree
for i = 1:length(forest_data)
    IS_maxweight_vec{i} = max(IS_prec_matrices{i},[],2);
end

%%
max_matrices_IS = {};

for i = 1:length(forest_data)
    max_matrices_IS{i} = zeros(size(max_matrices{i}));
    max_matrices_IS{i}(:,1) = max_matrices{i}(:,1);
    %max_matrices_IS{i}(:,2) = P_0 + weight_norm_matrices{i}(:,2) .* (max_matrices{i}(:,2) - P_0);
    max_matrices_IS{i}(:,3) = weight_matrices{i}(:,4);
    max_matrices_IS{i}(:,2) = max_matrices{i}(:,2);
end

max_matrix_IS_cat = cat(1,max_matrices_IS{:});

%% Saving max_matrices_IS
% format_string_3 = 'data__max_matrices_IS_%d.xls';
% string_3 = sprintf(format_string_3,ddbb_id);
% 
% for i = 1:length(max_matrices_IS)
%     table = array2table(max_matrices_IS{i},'VariableNames',{'Terminal Nodes', 'P_(N_,max,IS)'});
%     writetable(table,string_3,'Sheet',i)
% end
% 
% format_string_4 = 'data__y_max_IS__l_%d.xls';
% string_4 = sprintf(format_string_4,ddbb_id);
% 
% table = array2table(y_max_IS__l,'VariableNames',{'Forest ID', 'y_max_IS__l'});
% writetable(table,string_4,'Sheet',i)

%% Computing IS and MC based safety margins

% %% WILKS FORMULA
% Using bracketing OS approach

%max_matrix_MC_index = find(max_matrix_IS_cat(:,1) == 0);
%max_matrix_MC_cat = max_matrix_IS_cat(max_matrix_MC_index,:);

max_matrix_MC_cat = zeros(1,2);
max_matrix_MC_cat(1,:) = max_matrix_IS_cat(1,1:2)';
row_count = 2;

max_node = max(max_matrix_IS_cat(:,1));
for row = 2:size(max_matrix_IS_cat,1)
    node_value = max_matrix_IS_cat(row,1);
    if node_value == 0 || node_value == 2
        %append that row into the matrix
        max_matrix_MC_cat(row_count,1:2) = max_matrix_IS_cat(row,1:2);
        %restart counter
        max_node = node_value;
        row_count = row_count+1;
    else 
        if mod(node_value,2) == 0
            if node_value <= max_node
                %previous was a maximum
                %append that row into the matrix
                max_matrix_MC_cat(row_count,1:2) = max_matrix_IS_cat(row,1:2);
                %restart counter
                max_node = node_value;
                row_count = row_count+1;
            else
                %keep increasing
                max_node = node_value;
            end
        else
            %keep increasing
            max_node = node_value;
        end
    end
end

yhat_max_MC = max(max_matrix_MC_cat(:,2));

%%

% For upper safety margins
N_l_MC = length(max_matrix_MC_cat(:,2));
% gamma_1: estimated 95-th percentile of maximum values
gamma = 0.95;

% U_y: upper threshold
U_y = 1244;

% y_ref: reference value
y_ref = 1144;

% calculating the upper safety margin
[M_U_ya_MC, beta_1_MC] = SM_Up(U_y, y_ref, yhat_max_MC, N_l_MC, gamma);
display(M_U_ya_MC);
display(beta_1_MC);

%% IS-BASED UPPER CONFIDENCE INTERVAL

% gamma percentile estimations based on IS
gamma=0.95;
N_l_IS = length(max_matrix_IS_cat(:,2));
[y_max_gamma_alpha_1n, y_max_gamma_alpha_2n, y_max_gamma_alpha_3n, y_max_gamma_alpha_4n] = perc_quantile_base(max_matrix_IS_cat, N_l_IS, gamma);

% upper confidence interval (alpha) of gamma percentile, based on IS
alpha = .95;
ci_IS_2n = os_ci_IS_2n(max_matrix_IS_cat, y_max_gamma_alpha_2n, N_l_IS, alpha, gamma);

% calculating the upper safety margin
[M_U_ya_IS, beta_1_IS] = SM_Up(U_y, y_ref, ci_IS_2n, N_l_IS, gamma);
display(M_U_ya_IS);
display(beta_1_IS);

%% Functions

%%%%%%%%%%%
% SAFETY MARGINS
%%%%%%%%%%%

function [M_U_ya, beta_1] = SM_Up(U_y, y_ref, yhat_max, N, gamma_1)
    beta_1 = 1-gamma_1^N;
    if  (yhat_max < U_y) && (yhat_max > y_ref)
        M_U_ya =  ( U_y - yhat_max ) / (U_y - y_ref);
    elseif yhat_max >= U_y
        M_U_ya = 0;
    elseif ymax_hat <= y_ref
        M_U_ya = 0;
    end
end


%%%%%%%%%%
% IS-BASED QUANTILE ESTIMATOR
% UPPER CONFIDENCE INTERVAL OF PERCENTILE BASED ON IS
%%%%%%%%%%

%
function [y_max_gamma_alpha_1n, y_max_gamma_alpha_2n, y_max_gamma_alpha_3n, y_max_gamma_alpha_4n] = perc_quantile_base(y_max_vec_IS_2, num_trials_IS, gamma)
    
    %y_IS_max = max(y_max_IS_estimations);
    
    %sorting
    y_max_vec_IS_2 = sortrows(y_max_vec_IS_2,2);
    
    %adding up the likelihood (Glynn, 1996, alpha_1n)
    y_max_vec_IS_2(:,4) = cumsum(y_max_vec_IS_2(:,3));

    alpha_1n_index = find(y_max_vec_IS_2(:,4) >= (gamma)*num_trials_IS);
    if ~ isempty(alpha_1n_index)
        y_max_gamma_alpha_1n = y_max_vec_IS_2(alpha_1n_index(1), 2);
    else
        y_max_gamma_alpha_1n = 0;
    end
    
    alpha_3n_index = find(y_max_vec_IS_2(:,4) >= gamma*sum(y_max_vec_IS_2(:,3)));
    if ~ isempty(alpha_3n_index)
        y_max_gamma_alpha_3n = y_max_vec_IS_2(alpha_3n_index(1), 2);
    else
        y_max_gamma_alpha_3n = 0;
    end
    
    %adding up the likelihood (Glynn, 1996, alpha_2n; Nakayama)
    y_max_vec_IS_2(:,5) = flip(cumsum(flip(y_max_vec_IS_2(:,3))));
    
    %
    alpha_2n_index = find(y_max_vec_IS_2(:,5) >= (1-gamma)*num_trials_IS);
    if ~ isempty(alpha_2n_index)
        y_max_gamma_alpha_2n = y_max_vec_IS_2(alpha_2n_index(end), 2);
    else
        y_max_gamma_alpha_2n = 0;
    end
    
    %alpha_nak_index = find(y_max_vec_IS_2(:,5) <= (1-gamma)*num_trials_IS);
    %y_max_gamma_nak = y_max_vec_IS_2(alpha_nak_index(end), 2);

    alpha_4n_index = find(y_max_vec_IS_2(:,5) >= num_trials_IS - gamma*sum(y_max_vec_IS_2(:,3)));
    if ~ isempty(alpha_4n_index)
        y_max_gamma_alpha_4n = y_max_vec_IS_2(alpha_4n_index(end), 2);
    else
        y_max_gamma_alpha_4n = 0;
    end
    
end
%

function ci_IS_2n = os_ci_IS_2n(y_max_vec_IS_2, y_max_gamma_alpha_2n, num_trials_IS, alpha, gamma)
    % estimate phi(function of c, w)
    %w_1 = 15/3; w_2 = -12/3;
    %testing w
    %w_1 = 1.3; w_2 = -0.3;
    %w_1 = -1/3; w_2 = 4/3;
    w_1 = 5/3; w_2 = -2/3;

    % constraints for upper bound (ub [c_ub_min, c_ub_max]) interval of c
    c_ub_max = (1-gamma)*((num_trials_IS)^(1/2));
    c_ub_min = (1)/((num_trials_IS)^(1/2));

    % options for c_1
    c1_1 = c_ub_max;
    c1_2 = c_ub_min;

    % qm: optimal values of c?
    %assuming c_1 as c_ub_max...
    %c_2 = c1_1; 
    %c_1 = c_2 / 2;

    %assuming c_1 as c_ub_min...
    %c_1 = c1_2;
    %c_2 = 2*c_1;

    %taking the extremes
    c_1 = 0.05;
    %c_1 = 0.2;
    %c_2 = 0.06;
    c_2 = c1_2;

    upper_perc = gamma + c_1*((num_trials_IS)^(-1/2));
    lower_perc = gamma - c_2*((num_trials_IS)^(-1/2));

    %
    [~, est_upper_perc, ~, ~] = perc_quantile_base(y_max_vec_IS_2, num_trials_IS, upper_perc);
    [~, est_lower_perc, ~, ~] = perc_quantile_base(y_max_vec_IS_2, num_trials_IS, lower_perc);

    % c1>0
    phi_bar_pn1 = (est_upper_perc - y_max_gamma_alpha_2n) / (c_1*((num_trials_IS)^(-1/2))) ;

    % c2>0
    phi_bar_pn2 = (est_lower_perc - y_max_gamma_alpha_2n) / ( c_2*((num_trials_IS)^(-1/2))) ;
    %phi_bar_pn2 = (est_upper_perc - est_lower_perc) / ( 2*c_2*((num_trials_IS)^(-1/2))) ;

    % estimate phi: weighted sum of finite difference estimators 
    phi_bar = w_1*phi_bar_pn1 + w_2*phi_bar_pn2;

    % estimate psi
    One_func_index = find(y_max_vec_IS_2(:,2) >= y_max_gamma_alpha_2n);
    weights = y_max_vec_IS_2(One_func_index, 3);
   
    psi_est = ((1 / num_trials_IS) * (sum((weights.^2)))) - (1 - gamma)^2;

    % estimate kappa = function(phi_bar and psi_est)
    kappa_est = phi_bar * psi_est;

    % calculate one-sided confidence interval, with conf. level = alpha
    z_gamma = norminv(alpha,0,1);

    ci_IS_2n = y_max_gamma_alpha_2n + z_gamma*(kappa_est * ((num_trials_IS)^(-1/2)));

end

%%%%%%%%%%
% TREE-BASED DATA
%%%%%%%%%%

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

function node_data = search_tree_data(tree_data,node_id)
for i = 1:length(tree_data)
    current_node = tree_data{i};
    if current_node.id == node_id
        node_data = current_node;
    end
end
end

function parent_node = search_parent_node(tree,value)
    if value == 0
        parent_node = 'R';
    else
        for i = 1:length(tree)
            current_node = tree{i};
            if current_node.id == value
                parent_node = current_node.left;
            end
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

% New functions
function terminal_nodes_list = search_terminal_nodes(tree)
    terminal_nodes_list = {};
    for i = 1:length(tree)
    current_node = tree{i};
        if isempty(current_node.right)
            terminal_nodes_list{end+1} = current_node.id;
        end
    end
end

function parents_list = search_parents(tree,value)
    %using "backprop" to find the parents
    parents_list = {};
    current_parent = search_parent_node(tree,value);
    if value ~= 0
        while current_parent ~= 0
            %add the current one
            parents_list{end+1} = current_parent;
            %find the next parent
            current_parent = search_parent_node(tree,current_parent);
        end
        %add the root node
        parents_list{end+1} = current_parent;
    else
        %do nothing
    end
end

function pdf = PDF(lambda,t)
    pdf = lambda*exp(-lambda*t);
end

function cdf = CDF(lambda,t)
    cdf = 1 - exp(-lambda*t);
end




