function pilot_single_run()
% Single LJ run at the fixed query state t_q = 1180 s to validate the
% experiment folder (models, mex binaries, accelerator cache) before the
% full 20+20 recorded-seed run.
    S = load('LJ_test_sampled_vel_ACS_2_loop_1180.mat');
    det_params = refresh_operating_points(S);
    if ~isfield(det_params, 'Y_denorm_50')
        det_params.Y_denorm_50 = det_params.Y_denorm;
    end
    rng(1, 'twister');
    tWall = tic;
    [Tp, u_opt, f_opt, hist, tvec] = luus_jaakola_det( ...
        S.u0, S.u_bounds, S.max_passes, S.M, S.r_init, S.gamma, S.R, ...
        det_params, S.target, S.d_mpc, ...
        S.time_snapshot_root, S.delayed_time_rom, S.t_sam);
    fprintf('PILOT OK f=%.6f u1=%.6f minDSM=%.3f runtime=%.2f passes=%d wall=%.2f s\n', ...
        f_opt, u_opt(1), 1244 - max(Tp(:,1)), sum(tvec), nnz(tvec), toc(tWall));
    save('pilot_result.mat', 'Tp', 'u_opt', 'f_opt', 'hist', 'tvec');
end
