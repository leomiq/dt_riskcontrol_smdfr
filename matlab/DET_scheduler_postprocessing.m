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
            pred_twmax_signal = node_data.data.pred_tw_max_signal.Data(:,:);
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
                t_finish = total_time(end-2);
                %t_finish = total_time(end-1);
                t_next_event = 0;

                while t < t_finish 
                    trans_sys = parent_data.data.trans_sys.Data(t_counter+2,:);
                    mat_ttf = parent_data.data.mat_ttf.Data(:,:,t_counter+2);
                    trans_c = parent_data.data.trans_c.Data(t_counter+2);
                    % updating t_next_event
                    t_next_event = parent_data.data.new_ttf.Data(t_counter+2);
                    % Note: not considering the state of ps before activation 
                    % Classifying the transition
                    if any(trans_sys > ones(1,Ncomps)) %working
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
                                prob_vec_sample(1,k) = PDF(transition_rate_k_sample, delta_t);
                                prob_vec_target(1,k) = PDF(transition_rate_k_target, delta_t);

                                % hot stand-by assumption.
                                % node 0 does not fail.
                                % not failed, not inspected,
                                % starts at t_init
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
                                %prob_vec_sample(1,k) = 1 - CDF(transition_rate_k_sample, delta_t);
                                %prob_vec_target(1,k) = 1 - CDF(transition_rate_k_target, delta_t);
                                prob_vec_sample(1,k) = 1;
                                prob_vec_target(1,k) = 1;
                            else
                                prob_vec_sample(1,k) = PDF(transition_rate_k_sample, delta_t);
                                prob_vec_target(1,k) = PDF(transition_rate_k_target, delta_t);
                                
                                % hot stand-by assumption. only failed node right after
                                % node 0 includes CDF(trans_c).
                                % Failed, not inspected, intervention
                                % starts at t_init
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
                twmax_node = max(parent_data.data.tw_max_signal.Data);
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

max_matrix_MC_index = find(max_matrix_IS_cat(:,1) == 0);
max_matrix_MC_cat = max_matrix_IS_cat(max_matrix_MC_index,:);
yhat_max_MC = max(max_matrix_MC_cat(:,2))

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

%% Plotting histograms Wilks formula against IS based upper tolerance limit

figure

hold on
h1 = histogram(max_matrix_MC_cat(:,2), 'BinWidth', .05);
h1.FaceColor = 'blue';
h1.EdgeColor = 'none';
h1.FaceAlpha = 0.3;
xline(yhat_max_MC,'LineWidth', 2, 'Color', 'blue');

h4 = histogram(max_matrix_IS_cat(:,2), 'BinWidth', .05);
h4.FaceColor = 'red';
h4.EdgeColor = 'none';
h4.FaceAlpha = 0.3;
xline(ci_IS_2n,'LineWidth', 2, 'Color', 'red')
xlabel('${T}_{\mathrm{p,max}}$', 'Interpreter','latex','FontSize', 14);
ylabel('$\mathrm{frequency}$', 'Interpreter','latex');

xline(y_ref,'LineWidth', 1, 'Color', 'black')
str = '${y}_{\mathrm{ref}}$';
h5 = text(y_ref + 3, 20, str,'Interpreter','latex','FontSize', 14);

xline(U_y,'LineWidth', 1, 'Color', 'black')
str = '$U_{y}$';
h6 = text(U_y - 7, 20, str,'Interpreter','latex','FontSize', 14);
title(strcat('$\mathrm{SM}_{T_\mathrm{p}}(t=650s)$','=',num2str(M_U_ya_IS)), 'Interpreter', 'latex', 'FontSize', 14);


xlim([1160, 1166])

hlg = legend({'${T}_{\mathrm{p,max}}$ MC','${\hat{T}}_{\mathrm{p,max}}^{\gamma, \beta}$, Wilks',...
    '${T}_{\mathrm{p,max}}$ IS','${\hat{T}}_{\mathrm{p,max}}^{\gamma, \beta}$, IS'},...
    'Location','northeast','Interpreter','latex', 'FontSize', 10,'NumColumns', 2);

hold off

%% Plotting histograms Wilks formula against IS based upper tolerance limit

figure

subplot(1,2,1)
hold on
h1 = histogram(max_matrix_MC_cat(:,2), 'BinWidth', 0.5);
h1.FaceColor = 'blue';
h1.EdgeColor = 'none';
h1.FaceAlpha = 0.3;
xline(yhat_max_MC,'LineWidth', 2, 'Color', 'blue');
xlabel('${T}_{\mathrm{p,max}}$', 'Interpreter','latex','FontSize', 14);
ylabel('$\mathrm{frequency}$', 'Interpreter','latex');

hlg = legend({'${T}_{\mathrm{p,max}}$ distribution, MC-based sampling','${\hat{T}}_{\mathrm{p,max}}^{\gamma, \beta}$, Wilks formula'},...
    'Location','southoutside','Interpreter','latex', 'FontSize', 13);

hold off

subplot(1,2,2)
h4 = histogram(max_matrix_IS_cat(:,2), 'BinWidth', 1);
h4.FaceColor = 'red';
h4.EdgeColor = 'none';
h4.FaceAlpha = 0.3;
xline(ci_IS_2n,'LineWidth', 2, 'Color', 'red')
xlabel('${T}_{\mathrm{p,max}}$', 'Interpreter','latex','FontSize', 14);
ylabel('$\mathrm{frequency}$', 'Interpreter','latex');

hlg = legend({'${T}_{\mathrm{p,max}}$ distribution, IS method','${\hat{T}}_{\mathrm{p,max}}^{\gamma, \beta}$, IS method'},...
    'Location','southoutside','Interpreter','latex', 'FontSize', 13);

hold off

%%
%figure
%[h3_weights, h3_intervals] = histwc(max_matrix_IS_cat(:,2), max_matrix_IS_cat(:,3), 100);
%bar(h3_intervals, h3_weights);

figure
scatter(max_matrix_IS_cat(:,2), max_matrix_IS_cat(:,3));

%% Table with maximum values and weights

max_matrix_IS_2table = zeros(size(max_matrix_IS_cat,1),4);
l_index = 0;
for i = 1:size(max_matrix_IS_cat,1)
    if max_matrix_IS_cat(i,1) == 1 || max_matrix_IS_cat(i,1) == max_matrix_IS_cat(i-1,1)
        l_index = l_index + 1;
    end
    max_matrix_IS_2table(i,1) = l_index;
end

max_matrix_IS_2table(:,2:4) = max_matrix_IS_cat(:,:);

max_matrix_IS_table = array2table(max_matrix_IS_2table, 'VariableNames',{'l','b_l','Tb_{p,max}','w_b'});
% 
format_string = 'data__table_max_matrix_IS_2.xls';
writetable(max_matrix_IS_table,format_string);

%% Table for crude MC

%max_matrix_IS_2table_MC = max_matrix_IS_2table(find(max_matrix_IS_2table(:,2) == 1),:);

%max_matrix_IS_table_MC = array2table(max_matrix_IS_2table_MC(:,[1,3]), 'VariableNames', {'l','T_{p,max}'});

%format_string = 'data__table_max_matrix_MC.xls';
%writetable(max_matrix_IS_table_MC,format_string);

%% Plotting a colorscale plot, based on IS weights

% Assigning a color to each node, according to max IS for child
%weight_all = IS_prec_matrices{:};
max_norm_nodes = {};
min_norm_nodes = {};
for i = 1:length(forest_data)
    vec_w = IS_maxweight_vec{i};
    max_norm_weight = max(vec_w(vec_w > 0));
    max_norm_nodes{i} = max_norm_weight;
    min_norm_weight = min(vec_w(vec_w >= 0));
    min_norm_nodes{i} = min_norm_weight;
end

max_norm_nodes = cell2mat(max_norm_nodes);
min_norm_nodes = cell2mat(min_norm_nodes);

%%
cmap = flipud(colormap('winter')); % take your pick (doc colormap)
zmap = linspace(min(min_norm_nodes), max(max_norm_nodes), length(cmap));
alphamap = linspace(0.4, 1, length(cmap));

%
fig = figure();
hold on;
grid on;
box on;

xlim([410,1500])
ylim([1100,1240])
xline(500)

% for forest_data
for i = 8%:length(forest_data)
    forest_i = forest_data{i};
    tree_data_unzip = {};
    tree_unzip = {};
    % unzipping forest one by one
    for j = 1:length(forest_i)
        tree_data_unzip{end+1} = forest_data{i}(j).tree_data;
        tree_unzip{end+1} = forest_data{i}(j).tree;
    end
    root_data_unzip = forest_roots{i};
    terminal_nodes = search_terminal_nodes(tree_unzip);
    IS_prec_matrix = IS_prec_matrices{i};
    eps = 0.0006;
    for k = 1:size(IS_prec_matrix,1)
        node_data = search_tree_data(tree_data_unzip,k-1);
        max_weight_node = IS_maxweight_vec{i}(k,1);
        if max_weight_node == 0
            %display(node_data.id);
            %color = interp1(zmap, cmap, max_weight_node);
            %plot(node_data.data.npower.Time, node_data.data.npower.Data, 'Color', 'b', 'LineWidth', 1);
            %text(node_data.data.npower.Time(end-10),node_data.data.npower.Data(end-10),string(node_data.id),'Color','k');
        else
            if max_weight_node  == 1
                %color = interp1(zmap, cmap, max_weight_node);
                %alpha = interp1(zmap,alphamap,max_weight_node);
                color_alpha = 'r';
            else
                color = interp1(zmap, cmap, max_weight_node);
                alpha = interp1(zmap,alphamap,max_weight_node);
                color_alpha = cat(2,color,[alpha]);
            end
            plot(node_data.data.pred_tw_max_signal.Time, node_data.data.pred_tw_max_signal.Data(1,:), 'Color', color_alpha, 'LineWidth', 1)
            scatter(node_data.data.npower.Time(end),node_data.data.npower.Data(end),'k*','LineWidth',1.2);
            text(node_data.data.npower.Time(end)+5,node_data.data.npower.Data(end)+5,string(node_data.id),'Color','k');
            if node_data.id ~= 0
                %text(node_data.data.npower.Time(end)+2,node_data.data.npower.Data(end)+3,string(node_data.id),'Color','k','Interpreter','latex','FontSize',12);
            end
            if any(ismember(cell2mat(terminal_nodes),k-1))
                %if node_data.id == 1
                    %text(node_data.data.npower.Time(end)+5,node_data.data.npower.Data(end)+5,string(node_data.id),'Color','k');
                %eps = eps + 10;
                %end
            end
        end
    end
end

% for root_data
for i = 8%:length(forest_data)
    forest_i = forest_data{i};
    tree_data_unzip = {};
    tree_unzip = {};
    % unzipping forest one by one
    for j = 1:length(forest_i)
        tree_data_unzip{end+1} = forest_data{i}(j).tree_data;
        tree_unzip{end+1} = forest_data{i}(j).tree;
    root_data_unzip = forest_roots{i};
    terminal_nodes = search_terminal_nodes(tree_unzip);
    end
    %
    plot(root_data_unzip.data.tw_max_signal.Time, root_data_unzip.data.tw_max_signal.Data, 'Color', 'b', 'LineWidth', 1);
end

%hold off
%alpha(.5)

yline(1180)
text(350,1180 + 5,'$d(t)=1$','Color','black','Interpreter','Latex','FontSize',12)
yline(1220)
text(350,1220 + 5,'$d(t)=2$','Color','black','Interpreter','Latex','FontSize',12)
yline(1144)
text(350,1144 + 10,'$d(t)=2$','Color','black','Interpreter','Latex','FontSize',12)

xlabel('time [s]', 'Interpreter', 'Latex','FontSize',15)
ylabel('$T_{\mathrm{p}}$ [K]', 'Interpreter', 'Latex','FontSize',15)

colormap(cmap)
%colorbar
cb = colorbar('AxisLocation','out');
cb.Label.String = 'weight [-]';
cb.Label.FontSize = 15;
cb.Label.Interpreter = 'Latex';
%colormap flipud('winter')
%caxis([min(min_norm_nodes), max(max_norm_nodes)]);

%print(fig, '-dpdf', 'parloop_mcdet_forest_14_1.pdf');
%exportgraphics(fig, 'results.png', 'Resolution',800)

%% Scatter plot 
% 
% %cmap = colormap('winter'); % take your pick (doc colormap)
% %zmap = linspace(min(min_norm_nodes), max(max_norm_nodes), length(cmap));
% %alphamap = linspace(0.22, 1, length(cmap));
% 
% grouped_scatter = zeros(length(cat(1,max_matrices_IS{:})),4);
% row_counter = 1;
% for k = 1:length(forest_data)
%     max_matrix_IS = max_matrices_IS{k};
%     weight_norm_matrix = weight_norm_matrices{k};
%     for row = 1:length(max_matrix_IS)
%         grouped_scatter(row_counter,:) = [k, weight_norm_matrix(row,1), weight_norm_matrix(row,2), max_matrix_IS(row,2)];
%         row_counter = row_counter + 1;
%         %color = interp1(zmap, cmap, max_weight_node);
%         %scatter(max_matrix_IS(row,2), weight_norm_matrix(row,2), 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'LineWidth', 1)
%     end
% end
% 
% 
% %
% %figure();
% hold on;
% grid on;
% box on;
% 
% xlim([1,1.007])
% ylim([-0.1,1.1])
% %xline(510)
% 
% line=xline(1.006);
% %line.DisplayName= 'S_{P_n}';
% %line.legend = 'S_{P_n}';
% %line.Interpreter = 'Latex';
% text(1.006 + 0.00020, 0.5, '$S_{Tw_{max}}=1$','Color','black','Interpreter','Latex');
% 
% clr = lines(length(forest_data));
% hgs=gscatter(grouped_scatter(:,4), grouped_scatter(:,3), grouped_scatter(:,1),clr);
% %hgs(:).LineWidth = 6*ones(5);
% for i =1:5
%     hgs(i).LineWidth = 1.5;
%     hgs(i).Marker = 'o';
%     hgs(i).MarkerSize = 6;
%     hgs(i).MarkerFaceColor = 'auto';
% end
% hLG = legend([hgs],'l=1','l=2','l=3','l=4','l=5');
% hLG.Interpreter = 'Latex';
% hLG.Location = 'northeastoutside';
% 
% hold off
% %alpha(.5)
% 
% %xline(1.012)
% %text(1.012 + 0.00050, 0.5, '$S_{P_N}=2$','Color','black','Interpreter','Latex')
% 
% xlabel('$tw_{max}^{IS} [-]$', 'Interpreter', 'Latex');
% ylabel('$w_i [-]$', 'Interpreter', 'Latex');
% 
% %cb = colorbar('AxisLocation','out');
% %cb.Label.String = 'weight [-]';
% %cb.Label.FontSize = 12;
% %cb.Label.Interpreter = 'Latex';
% %colormap winter
% %caxis([min(min_norm_nodes), max(max_norm_nodes)]);
% 
% %% Scatter plot (v2)
% 
% cmap = flipud(colormap('winter')); % take your pick (doc colormap)
% zmap = linspace(min(min_norm_nodes), max(max_norm_nodes), length(cmap));
% %alphamap = linspace(0.22, 1, length(cmap));
% 
% %
% %figure();
% hold on;
% grid on;
% box off;
% 
% xlim([1,1.015])
% ylim([-0.1,1.1])
% %xline(510)
% 
% line=xline(1.006);
% %line.DisplayName= 'S_{P_n}';
% %line.legend = 'S_{P_n}';
% %line.Interpreter = 'Latex';
% text(1.006 + 0.00020, 0.5, '$S_{P_N}=1$','Color','black','Interpreter','Latex');
% 
% line=xline(1.012);
% text(1.012 - 0.00020, 0.5, '$S_{P_N}=2$','Color','black','Interpreter','Latex');
% 
% grouped_scatter = zeros(length(cat(1,max_matrices_IS{:})),4);
% row_counter = 1;
% for k = 1:length(forest_data)
%     max_matrix_IS = max_matrices_IS{k};
%     %max_matrix = max_matrices{k};
%     weight_norm_matrix = weight_norm_matrices{k};
%     for row = 1:length(max_matrix_IS)
%         %grouped_scatter(row_counter,:) = [k, weight_norm_matrix(row,1), weight_norm_matrix(row,2), max_matrix_IS(row,2)];
%         row_counter = row_counter + 1;
%         color = interp1(zmap, cmap, weight_norm_matrix(row,2));
%         scatter(max_matrix_IS(row,2), weight_norm_matrix(row,2), 'MarkerEdgeColor', color, 'LineWidth', 1.5)
%     end
% end
% 
% xlabel('$P_{N,max}^{IS} [-]$', 'Interpreter', 'Latex');
% ylabel('$w_i [-]$', 'Interpreter', 'Latex');
% 
% colormap(cmap)
% cb = colorbar('AxisLocation','out');
% cb.Label.String = 'weight [-]';
% %cb.Label.FontSize = 12;
% %cb.Label.Interpreter = 'Latex';
% %colormap winter
% caxis([min(min_norm_nodes), max(max_norm_nodes)]);
% 
% %hold off;
% 
% ax=axes;
% set(ax,'units','normalized','position',[0.5,0.2,0.2,0.2])
% %box(ax,'on')
% grid(ax,'on')
% scatter(grouped_scatter(:,4), grouped_scatter(:,3),'g','parent',ax)
% set(ax,'xlim',[1,1.0005],'ylim',[-0.0001,0.001])
% set(ax,'box','on')
% set(ax,'Xgrid','on')
% set(ax,'Ygrid','on')
% 
% %%
% % LAST CHECKPOINT (02-03)
% % load('WORKSPACE_TEST_3400-800.mat');
% %%___________
% 
% % Another useful plots
% 
% % Histograms, predicted values distribution
% subplot(2,1,1);
% %yyaxis left
% hold on
% histogram(max_matrices_total(:,2),50,'facealpha',.5,'facecolor','blue','edgecolor','black')
% %histfit(max_matrices_total(:,2),50)
% xline(yhat_max,'-','Color','b', 'LineWidth', 1.25)
% txt = {'${\hat{P}}_{\mathrm{N,max}}$'};
% text(yhat_max - 0.0015, 6, txt, 'Interpreter', 'Latex', 'Color', 'b', 'FontSize',13)
% xline(U_y_IS,'-','Color','k', 'LineWidth', 1.25)
% txt4 = {'${{U}}_{P_N}$'};
% text(U_y - 0.001, 6, txt4, 'Interpreter', 'Latex', 'Color', 'k', 'FontSize',13)
% xline(1,'-','Color','k', 'LineWidth', 1.5)
% txt5 = {'$P_{N,\mathrm{ref}}$'};
% text(1 + 0.0005, 6, txt5, 'Interpreter', 'Latex', 'Color', 'k', 'FontSize',13)
% ylabel('Frequency [-]', 'Interpreter', 'Latex')
% xlim([0.9995, 1.0155])
% title('Maximum values obtained with brute MC', 'Interpreter', 'Latex','FontSize',13)
% hold off
% 
% % Histograms, maximum values
% subplot(2,1,2)
% hold on
% %yyaxis right
% histogram(max_matrices_IS_total(:,2),50,'facealpha',.5,'facecolor','red','edgecolor','black')
% %histfit(max_matrices_IS_total(:,2),50,'lognormal')
% xline(yhat_max,'-','Color','b', 'LineWidth', 1.25)
% txt2 = {'${\hat{P}}_{\mathrm{N,max}}$'}
% text(yhat_max - 0.0015, 50, txt2, 'Interpreter', 'Latex', 'Color', 'b', 'FontSize',13)
% xline(yhat_max_IS,'-','Color','r', 'LineWidth', 1.25)
% txt3 = {'${\hat{P}}_{N,\mathrm{max}}^{\mathrm{IS}}$'}
% text(yhat_max_IS - 0.002, 50, txt3, 'Interpreter', 'Latex', 'Color', 'r', 'FontSize',13)
% xline(U_y_IS,'-','Color','k', 'LineWidth', 1.25)
% text(U_y - 0.001, 50, txt4, 'Interpreter', 'Latex', 'Color', 'k', 'FontSize',13)
% xlim([0.9995, 1.0155])
% xline(1,'-','Color','k', 'LineWidth', 1.5)
% txt5 = {'$P_{N,\mathrm{ref}}$'};
% text(1 + 0.0005, 50, txt5, 'Interpreter', 'Latex', 'Color', 'k', 'FontSize',13)
% ylabel('Frequency [-]','Interpreter','Latex')
% hold off
% %eps_histogram = 10E-5;
% % Adjusting estimator and label
% 
% %text(yhat_max - eps_histogram, 50, txt, 'Interpreter', 'Latex', 'Color', 'b', 'FontSize',11)
% 
% xlabel('${P}_{\mathrm{N,max}}$', 'Interpreter', 'Latex','FontSize',13)
% title('Maximum values corrected with IS', 'Interpreter', 'Latex')

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
    w_1 = 15/3; w_2 = -12/3;
    %testing w
    %w_1 = 1.3; w_2 = -0.3;
    %w_1 = -1/3; w_2 = 4/3;
    %w_1 = 5/3; w_2 = -2/3;

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
    c_1 = 0.71;
    c_2 = 0.06;

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




