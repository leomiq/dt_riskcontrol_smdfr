figure;

plot(control_series(end-60+2:end,1), control_series(end-60+2:end,2))

xlabel('Time [s]', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel(['Control sequence,', char(10), 'induced current [A]'], 'Interpreter', 'Latex', 'FontSize', 15);
%%
time_vector = zeros(length(time_snapshots),1);
%%
for idx = 1:length(time_snapshots)
    snapshot = time_snapshots(idx);
    
    % Load global workspace if needed
    load(['PID_test_sampled_vel_ACS_2_loop_' num2str(snapshot) '.mat'])

    time_vector(idx,1) = sum(t_vector);
end
%%
figure;
hist(time_vector);

ylabel('frequency', 'Interpreter', 'Latex', 'FontSize', 15);
xlabel('computing time per step $t_{\rm sim}$ [s] ', 'Interpreter', 'Latex', 'FontSize', 15);

%%
figure;

plot(time_snapshot_root*ones(1,26) + (0.5:1:25.5), u_opt);

xlabel('$t [s]$', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel(['Control sequence,', char(10), 'induced current [A]'], 'Interpreter', 'Latex', 'FontSize', 15);