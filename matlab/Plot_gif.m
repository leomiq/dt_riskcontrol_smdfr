%% GBM with IS-MCDET. Plot without zoom

t_sam = 0.5;
time_snapshot_init = 1160;
time_snapshots = time_snapshot_init:t_sam:time_snapshot_init + 30;
delayed_time_rom = 10;

% Name of the GIF filea
gifFileName = 'LJ_animation_1180_tail.gif';

% Pre-calculate frame delay (in seconds) for the GIF
delayTime = 0.2;

for idx = 1:length(time_snapshots)
    snapshot = time_snapshots(idx);
    
    % Load global workspace if needed
    load(['LJ_test_sampled_vel_ACS_2_loop_' num2str(snapshot) '.mat'])
    
    delayed_indexes = delayed_time_rom / t_sam;
    
    % Create a new figure for each snapshot (or clear the current figure)
    figure('visible', 'off');  % set invisible if you don't want to display each frame
    grid on
    hold on
    
    % Plot the main signal in black
    plot(forest_roots_hrm{1}.data.tw_max_signal.Time(:), ...
         forest_roots_hrm{1}.data.tw_max_signal.Data(:), ...
         'Color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);
    
    % Loop over the parent elements
    for i = 1:N_par
        % Plot the corresponding roots in blue dashed line
        plot(forest_roots{i}.data.tw_max_signal.Time(:), ...
             forest_roots{i}.data.tw_max_signal.Data(:), ...
             'Color', 'b', 'LineStyle', '--', 'LineWidth', 1.5);
         
        % Loop over forest_data elements
        for j = 1:length(forest_data{i})
            if j == 1
               plot(forest_data{i}(j).tree_data.data.pred_tw_max_signal.Time(end-N_pred:end), ...
                    forest_data{i}(j).tree_data.data.pred_tw_max_signal.Data(end-N_pred:end), ...
                    'Color', [0 0 1 0.5], 'LineWidth', 1.5);
               plot(forest_data_opt{i}(j).tree_data.data.pred_tw_max_signal.Time(end-N_pred:end), ...
                    forest_data_opt{i}(j).tree_data.data.pred_tw_max_signal.Data(end-N_pred:end), ...
                    'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
               plot(forest_data_best{i}(j).tree_data.data.pred_tw_max_signal.Time(end-N_pred:end), ...
                    forest_data_best{i}(j).tree_data.data.pred_tw_max_signal.Data(end-N_pred:end), ...
                    'Color', 'm', 'LineStyle', '--', 'LineWidth', 1.5);
            elseif j == 2 
               plot(forest_data{i}(j).tree_data.data.pred_tw_max_signal.Time, ...
                    forest_data{i}(j).tree_data.data.pred_tw_max_signal.Data, ...
                    'Color', [0 0 1 0.5], 'LineWidth', 1.5);
               plot(forest_data_opt{i}(j).tree_data.data.pred_tw_max_signal.Time, ...
                    forest_data_opt{i}(j).tree_data.data.pred_tw_max_signal.Data, ...
                    'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
               plot(forest_data_best{i}(j).tree_data.data.pred_tw_max_signal.Time, ...
                    forest_data_best{i}(j).tree_data.data.pred_tw_max_signal.Data, ...
                    'Color', 'm', 'LineStyle', '--', 'LineWidth', 1.5);
            elseif j == 3 
               plot(forest_data{i}(j).tree_data.data.pred_tw_max_signal.Time, ...
                    forest_data{i}(j).tree_data.data.pred_tw_max_signal.Data, ...
                    'Color', [0 0 1 0.5], 'LineWidth', 1.5);
               plot(forest_data_opt{i}(j).tree_data.data.pred_tw_max_signal.Time, ...
                    forest_data_opt{i}(j).tree_data.data.pred_tw_max_signal.Data, ...
                    'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
               plot(forest_data_best{i}(j).tree_data.data.pred_tw_max_signal.Time, ...
                    forest_data_best{i}(j).tree_data.data.pred_tw_max_signal.Data, ...
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
    xlabel('Time [s]', 'Interpreter', 'Latex', 'FontSize', 15);
    ylabel('$T_\mathrm{p}$ [K]', 'Interpreter', 'Latex', 'FontSize', 15);

    % Add a legend with LaTeX notation
    % Add a legend with LaTeX notation
    legend({'$T_{\rm p}, \mathrm{PO}$', ...
        '',...
        '$T_{\rm p}, \mathrm{GBM}$, no PCS action', ...
        '$T_{\rm p}, \mathrm{GBM}$, best found solution (LJ)', ...
        '$T_{\rm p}, \mathrm{GBM}$, previous PCS action', ...
        },...
       'Interpreter', 'Latex', 'FontSize', 12, 'Location', 'northeast');
    
    % Set axes limits
    xlim([1120, 1200]);
    ylim([1175, 1190]);
    
    hold off;
    
    % Capture the plot as an image frame
    drawnow;
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % Write to the GIF File
    if idx == 1
        % For the first frame, create the file with loop count set to infinite
        imwrite(imind, cm, gifFileName, 'gif', 'DelayTime', delayTime, 'Loopcount', inf);
    else
        % For subsequent frames, append to the file
        imwrite(imind, cm, gifFileName, 'gif', 'DelayTime', delayTime, 'WriteMode', 'append');
    end
    
    % Close the figure if you don't want to display it during processing
    close(gcf);
end

disp(['GIF animation created: ', gifFileName])