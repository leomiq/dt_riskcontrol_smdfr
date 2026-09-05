%% GBM with IS-MCDET. Plot without zoom
    
figure;
grid on
hold on

delayed_indexes = delayed_time_rom / t_sam;

% Plot the main signal in black
plot(forest_roots_hrm{1}.data.tw_max_signal.Time(1:end), ...
     forest_roots_hrm{1}.data.tw_max_signal.Data(1:end), ...
     'Color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);

plot(forest_roots_hrm{1}.data.tw_max_signal.Time(end-60:end), ...
      forest_roots_hrm{1}.data.tw_max_signal.Data(end-60:end), ...
      'Color', 'b', 'LineStyle', '-', 'LineWidth', 1.5);

% Loop over the parent elements
for i = 1:N_par
    % Plot the corresponding roots in blue dashed line
    plot(forest_roots{i}.data.tw_max_signal.Time(:), ...
         forest_roots{i}.data.tw_max_signal.Data(:), ...
         'Color', 'b', 'LineStyle', '--', 'LineWidth', 1.5);
     
    % % Loop over forest_data elements
    % for j = 1:length(forest_data{i})
    %     if j == 1
    %        plot(forest_data{i}(j).tree_data.data.pred_tw_max_signal.Time(end-N_pred:end), ...
    %             forest_data{i}(j).tree_data.data.pred_tw_max_signal.Data(end-N_pred:end), ...
    %             'Color', [0 0 1 0.5], 'LineWidth', 1.5);
    %        plot(forest_data_opt{i}(j).tree_data.data.pred_tw_max_signal.Time(end-N_pred:end), ...
    %             forest_data_opt{i}(j).tree_data.data.pred_tw_max_signal.Data(end-N_pred:end), ...
    %             'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
    %        plot(forest_data_best{i}(j).tree_data.data.pred_tw_max_signal.Time(end-N_pred:end), ...
    %             forest_data_best{i}(j).tree_data.data.pred_tw_max_signal.Data(end-N_pred:end), ...
    %             'Color', 'm', 'LineStyle', '--', 'LineWidth', 1.5);
    %     elseif j == 2 
    %        plot(forest_data{i}(j).tree_data.data.pred_tw_max_signal.Time, ...
    %             forest_data{i}(j).tree_data.data.pred_tw_max_signal.Data, ...
    %             'Color', [0 0 1 0.5], 'LineWidth', 1.5);
    %        plot(forest_data_opt{i}(j).tree_data.data.pred_tw_max_signal.Time, ...
    %             forest_data_opt{i}(j).tree_data.data.pred_tw_max_signal.Data, ...
    %             'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
    %        plot(forest_data_best{i}(j).tree_data.data.pred_tw_max_signal.Time, ...
    %             forest_data_best{i}(j).tree_data.data.pred_tw_max_signal.Data, ...
    %             'Color', 'm', 'LineStyle', '--', 'LineWidth', 1.5);
    %     elseif j == 3 
    %        plot(forest_data{i}(j).tree_data.data.pred_tw_max_signal.Time, ...
    %             forest_data{i}(j).tree_data.data.pred_tw_max_signal.Data, ...
    %             'Color', [0 0 1 0.5], 'LineWidth', 1.5);
    %        plot(forest_data_opt{i}(j).tree_data.data.pred_tw_max_signal.Time, ...
    %             forest_data_opt{i}(j).tree_data.data.pred_tw_max_signal.Data, ...
    %             'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
    %        plot(forest_data_best{i}(j).tree_data.data.pred_tw_max_signal.Time, ...
    %             forest_data_best{i}(j).tree_data.data.pred_tw_max_signal.Data, ...
    %             'Color', 'm', 'LineStyle', '--', 'LineWidth', 1.5);
    %     end
    % end
end

% Add extra lines and texts
xline(500)
xline(1000)
xline(1160)
yline(1180)
text(430, 1180 + 2, '$d(t)=1$', 'Color', 'black', 'Interpreter', 'Latex')
yline(1220)
yline(target)
text(430, 1220 + 2, '$d(t)=2$', 'Color', 'black', 'Interpreter', 'Latex')

% Set labels
xlabel('Time [s]', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('$T_\mathrm{p}$ [K]', 'Interpreter', 'Latex', 'FontSize', 15);

% Set axes limits
xlim([1120, 1200]);
ylim([1175, 1190]);

% Add a legend with LaTeX notation
legend({'$T_{\rm p}, \mathrm{PO}$', ...
        %'$T_{\rm p}, \mathrm{GBM}$, no PCS action', ...
        %'$T_{\rm p}, \mathrm{GBM}$, best found solution (LJ)', ...
        %'$T_{\rm p}, \mathrm{GBM}$, previous PCS action', ...
        '$T_{\rm p}$, optimal control trajectory (LJ)', ...
        },...
       'Interpreter', 'Latex', 'FontSize', 12, 'Location', 'northeast');

hold off;

%xlim([400, 1400]);
%ylim([1050, 1225]);

%% GBM with IS-MCDET. Plot without zoom
    
figure;
grid on
hold on

delayed_indexes = delayed_time_rom / t_sam;
    

% Plot the main signal in black
plot(forest_roots_hrm{1}.data.mdot_c_net_signal.Time(:), ...
     forest_roots_hrm{1}.data.mdot_c_net_signal.Data(:), ...
     'Color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);

% Loop over the parent elements
for i = 1:N_par
    % Plot the corresponding roots in blue dashed line
    plot(forest_roots{i}.data.mdot_c_net_signal.Time(:), ...
         forest_roots{i}.data.mdot_c_net_signal.Data(:), ...
         'Color', 'b', 'LineStyle', '--', 'LineWidth', 1.5);
     
    % Loop over forest_data elements
    for j = 1:length(forest_data{i})
        if j == 1
           plot(forest_data{i}(j).tree_data.data.mdot_c_net_signal.Time(end-N_pred:end), ...
                forest_data{i}(j).tree_data.data.mdot_c_net_signal.Data(end-N_pred:end), ...
                'Color', [0 0 1 0.5], 'LineWidth', 1.5);
           plot(forest_data_opt{i}(j).tree_data.data.mdot_c_net_signal.Time(end-N_pred:end), ...
                forest_data_opt{i}(j).tree_data.data.mdot_c_net_signal.Data(end-N_pred:end), ...
                'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
           plot(forest_data_best{i}(j).tree_data.data.mdot_c_net_signal.Time(end-N_pred:end), ...
                forest_data_best{i}(j).tree_data.data.mdot_c_net_signal.Data(end-N_pred:end), ...
                'Color', 'm', 'LineStyle', '--', 'LineWidth', 1.5);
        elseif j == 2 
           plot(forest_data{i}(j).tree_data.data.mdot_c_net_signal.Time, ...
                forest_data{i}(j).tree_data.data.mdot_c_net_signal.Data, ...
                'Color', [0 0 1 0.5], 'LineWidth', 1.5);
           plot(forest_data_opt{i}(j).tree_data.data.mdot_c_net_signal.Time, ...
                forest_data_opt{i}(j).tree_data.data.mdot_c_net_signal.Data, ...
                'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
           plot(forest_data_best{i}(j).tree_data.data.mdot_c_net_signal.Time, ...
                forest_data_best{i}(j).tree_data.data.mdot_c_net_signal.Data, ...
                'Color', 'm', 'LineStyle', '--', 'LineWidth', 1.5);
        elseif j == 3 
           plot(forest_data{i}(j).tree_data.data.mdot_c_net_signal.Time, ...
                forest_data{i}(j).tree_data.data.mdot_c_net_signal.Data, ...
                'Color', [0 0 1 0.5], 'LineWidth', 1.5);
           plot(forest_data_opt{i}(j).tree_data.data.mdot_c_net_signal.Time, ...
                forest_data_opt{i}(j).tree_data.data.mdot_c_net_signal.Data, ...
                'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
           plot(forest_data_best{i}(j).tree_data.data.mdot_c_net_signal.Time, ...
                forest_data_best{i}(j).tree_data.data.mdot_c_net_signal.Data, ...
                'Color', 'm', 'LineStyle', '--', 'LineWidth', 1.5);
        end
    end
end

% Add extra lines and texts
xline(500)
xline(1000)
xline(1160)
yline(1180)
text(430, 1180 + 2, '$d(t)=1$', 'Color', 'black', 'Interpreter', 'Latex')
yline(1220)
yline(target)
text(430, 1220 + 2, '$d(t)=2$', 'Color', 'black', 'Interpreter', 'Latex')

% Set labels
xlabel('t [s]', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('${\dot{m}}_{\kappa,\mathrm{net}}$ [kg/s]', 'Interpreter', 'Latex', 'FontSize', 15);

% Set axes limits
xlim([1120, 1200]);
ylim([-5, 5]);

% Add a legend with LaTeX notation
legend({'${\dot{m}}_{\kappa,\mathrm{net}}, \mathrm{PO}$', ...
         '',...
        '${\dot{m}}_{\kappa,\mathrm{net}}$, GBM, no PCS action', ...
        '${\dot{m}}_{\kappa,\mathrm{net}}$, GBM, best found solution (LJ)', ...
        '${\dot{m}}_{\kappa,\mathrm{net}}$, GBM, previous PCS action', ...
        },...
       'Interpreter', 'Latex', 'FontSize', 12, 'Location', 'northeast');

hold off;

%xlim([400, 1400]);
%ylim([1050, 1225]);

%% GBM with IS-MCDET. Plot with zoom

figure;
grid off
hold on

% Main plot
plot(forest_roots_hrm{1}.data.tw_max_signal.Time(:), forest_roots_hrm{1}.data.tw_max_signal.Data(:), 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);
for i = 1:N_par
    plot(forest_roots{i}.data.tw_max_signal.Time(:), forest_roots{i}.data.tw_max_signal.Data(:), 'Color', 'b', 'LineStyle', '--', 'LineWidth', 1.5);
    for j = 1:length(forest_data{i})
        if j == 1
            plot(forest_data{i}(j).tree_data.data.pred_tw_max_signal.Time(1:end), forest_data{i}(j).tree_data.data.pred_tw_max_signal.Data(1:end), 'Color', [0 0 1 0.5], 'LineWidth', 1.5);
        elseif j == 2 
           plot(forest_data{i}(j).tree_data.data.pred_tw_max_signal.Time(1:end), forest_data{i}(j).tree_data.data.pred_tw_max_signal.Data(1:end),'Color',[1 0 0 0.5],'LineWidth',1.5);
        elseif j == 3 
           plot(forest_data{i}(j).tree_data.data.pred_tw_max_signal.Time(1:end), forest_data{i}(j).tree_data.data.pred_tw_max_signal.Data(1:end),'Color',[0 0 1 0.5],'LineWidth',1.5);
        end
    end
end

xline(500);
xline(1300);
yline(1180);
text(430, 1180 + 2, '$d(t)=1$', 'Color', 'black', 'Interpreter', 'Latex');
yline(1220);
text(430, 1220 + 2, '$d(t)=2$', 'Color', 'black', 'Interpreter', 'Latex');

% Add labels
xlabel('t [s]', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('$T_\mathrm{p}$ [K]', 'Interpreter', 'Latex', 'FontSize', 15);

% Add a legend with LaTeX notation
legend({'$T_p, \mathrm{PO}$', ...
        '$T_p, \mathrm{WBM}$'}, ...
       'Interpreter', 'Latex', 'FontSize', 12, 'Location', 'northeast');

xlim([400, 1450]);
ylim([1050, 1230]);

hold off 

% Add inset zoomed plot
axes('Position', [0.5 0.2 0.35 0.35]); % Adjust position and size
box on; % Add box around inset plot
hold on;
grid on;

% Plot zoomed region
plot(forest_roots_hrm{1}.data.tw_max_signal.Time(:), forest_roots_hrm{1}.data.tw_max_signal.Data(:), 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);
for i = 1:N_par
    plot(forest_roots{i}.data.tw_max_signal.Time(:), forest_roots{i}.data.tw_max_signal.Data(:), 'Color', 'b', 'LineStyle', '--', 'LineWidth', 1.5);
    for j = 1:length(forest_data{i})
        if j == 1
            plot(forest_data{i}(j).tree_data.data.pred_tw_max_signal.Time(1:end), forest_data{i}(j).tree_data.data.pred_tw_max_signal.Data(1:end), 'Color', [0 0 1 0.5], 'LineWidth', 1.5);
        elseif j == 2 
           plot(forest_data{i}(j).tree_data.data.pred_tw_max_signal.Time(1:end), forest_data{i}(j).tree_data.data.pred_tw_max_signal.Data(1:end),'Color',[1 0 0 0.5],'LineWidth',1.5);
        elseif j == 3 
           plot(forest_data{i}(j).tree_data.data.pred_tw_max_signal.Time(1:end), forest_data{i}(j).tree_data.data.pred_tw_max_signal.Data(1:end),'Color',[0 0 1 0.5],'LineWidth',1.5);
        end
    end
end

% Add inset labels
xlabel('t [s]', 'Interpreter', 'Latex', 'FontSize', 8);
ylabel('$T_\mathrm{p}$ [K]', 'Interpreter', 'Latex', 'FontSize', 8);

% Set x-axis limits for zoom
xlim([1295, 1405]);
ylim([1150, 1175]);

hold off;

%%

figure;
hold on
for i = 1:N_par
    %plot(forest_roots{i}.data.mdot_c_leak_signal.Time(:), forest_roots{i}.data.mdot_c_leak_signal.Data(:));
    %plot(forest_data{i}.tree_data.data.mdot_c_leak_signal.Time(:), forest_data{i}.tree_data.data.mdot_c_leak_signal.Data(:));
    %plot(forest_data{i}.tree_data.data.pred_delta_tw_max_signal.Data(:))
end

%% MC dropout generated sequences. snapshot: t=1200 s
% see MATLAB_RMM_pilot_4_1

figure;

% Subplot 1: Plotting Y_denorm for the second column (index 2)
subplot(1, 2, 1); % 1 row, 2 columns, first subplot
hold on;
for i = 1:N_par
    plot(Y_denorm(i,:,1), 'Color', [0 0 1 0.5]); % Plot the second column
end
%title('Subplot 1: array\_pred\_M(:, :, 2)');
xlabel('Sequence index [-]', 'Interpreter', 'Latex', 'FontSize', 12);
ylabel('$\dot{m}_{c}$', 'Interpreter', 'Latex', 'FontSize', 12);
grid on;
hold off;

% Subplot 2: Plotting Y_denorm for the first column (index 1)
subplot(1, 2, 2); % 1 row, 2 columns, second subplot
hold on;
for i = 1:N_par
    plot(Y_denorm(i,:,2), 'Color', [1 0 0 0.5]); % Plot the first column
end
%title('Subplot 2: array\_pred\_M(:, :, 1)');
xlabel('Sequence index [-]', 'Interpreter', 'Latex', 'FontSize', 12);
ylabel('$\Delta T_{\mathrm{p}}$', 'Interpreter', 'Latex', 'FontSize', 12);
grid on;
hold off;

%% Dynamic safety margins calculated in three snapshots

% Data from the table
t = [1000; 1100; 1300];           % Time in seconds
y1 = [ 0.6677; 0.6240; 0.6976];   % y1 values: MC-based
y2 = [ 0.6677; 0.6042; 0.6977];  % y2 values: IS-based

% Define the range for the fitted curve (only between t=950 and t=1200)
t_fit = linspace(1000, 1300, 100);

% Quadratic fit for y1
p1 = polyfit(t, y1, 2);         % Fit quadratic polynomial for y1
y1_fit = polyval(p1, t_fit);    % Evaluate polynomial between 950 and 1200

% Quadratic fit for y2
p2 = polyfit(t, y2, 2);         % Fit quadratic polynomial for y2
y2_fit = polyval(p2, t_fit);    % Evaluate polynomial between 950 and 1200

% Create the plot
figure;
hold on; grid off;

% Plot scatter points
scatter(t, y1, 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineWidth', 1.5);
scatter(t, y2, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'LineWidth', 1.5);

% Plot quadratic regression curves with alpha transparency
h1 = plot(t_fit, y1_fit, 'b-', 'LineWidth', 1.5); % Blue curve for y1
h2 = plot(t_fit, y2_fit, 'r-', 'LineWidth', 1.5); % Red curve for y2

% Apply transparency (alpha = 0.3)
h1.Color(4) = 0.3; % Set alpha for y1 curve
h2.Color(4) = 0.3; % Set alpha for y2 curve

% Add solid black lines with transparency connecting scatter points to the x-axis
for i = 1:length(t)
    % Solid black line with alpha transparency
    plot([t(i), t(i)], [0, y1(i)], '-', 'Color', [0 0 0 0.3], 'LineWidth', 1);
    plot([t(i), t(i)], [0, y2(i)], '-', 'Color', [0 0 0 0.3], 'LineWidth', 1);
end

% Set axes limits and labels
xlim([900 1400]);
ylim([0.5 0.8]); % Adjust as necessary
xlabel('snapshot time [s]', 'Interpreter', 'Latex', 'FontSize', 12);
ylabel('$\mathrm{PSM} \, [-]$', 'Interpreter', 'Latex', 'FontSize', 12);

% Add legend
legend({'$\mathrm{PSM}$, MC-based (Wilks formula)', '$\mathrm{PSM}$, IS-based (Glynn formula)', ...
        '$\mathrm{PSM}$ MC, quadratic fit', '$\mathrm{PSM}$ IS, quadratic fit'}, ...
       'Interpreter', 'Latex', 'FontSize', 12, 'Location', 'southoutside');

% Add title
title('Evolution of PSMs', 'Interpreter', 'Latex', 'FontSize', 14);

hold off;
