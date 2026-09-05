function run_phase2c_full_loop_pilot()
% Validation pilot: one full LJ realization and one full GA realization
% (seed 1) to confirm the closed-loop pipeline and measure wall-clock time
% before deciding the number of realizations for the KPI statistics.
    t0 = tic;
    kpiLJ = run_phase2c_full_loop('LJ', 1, 'phase2c_fullloop_pilot');
    fprintf('=== PILOT LJ KPI ===\n');
    disp(kpiLJ);
    kpiGA = run_phase2c_full_loop('GA', 1, 'phase2c_fullloop_pilot');
    fprintf('=== PILOT GA KPI ===\n');
    disp(kpiGA);
    save('phase2c_fullloop_pilot_summary.mat', 'kpiLJ', 'kpiGA');
    fprintf('PILOT COMPLETO en %.1f min\n', toc(t0)/60);
end
