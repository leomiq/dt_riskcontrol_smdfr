function det_params = refresh_operating_points(S)
% REFRESH_OPERATING_POINTS  Regenerate the Simulink operating points for the
% fixed query state stored in the archived workspace S.
%
% The archived det_params contains OperatingPoint objects saved in Oct 2025.
% Simulink refuses to load them into the current R2024a model ("model was
% changed after the operating point was saved"), so the operating points are
% rebuilt deterministically with the current model, using the archived
% scenario parameters, control history and leakage signal as inputs. This
% preserves the fixed-state experimental design (same accident seed, same
% archived scenario data) while making the OPs loadable.
    % Silence the interactive Simulink dialogs (model-workspace save,
    % fast-restart recompile, cross-release operating points) so this runs
    % unattended under matlab -batch.
    try, suppress_simulink_prompts(); catch, end

    name_hrm = 'SMDFR_hrm_MCDET_2022b_control_CMEX';
    name_rom = 'rom_MCDET_v2_2022b_control_CMEX';
    name_wks = 'simulation_core_dualfluid_FI3_sim_rom_extract_MCDET_v2';
    t0      = double(S.time_snapshot_root);
    ts      = double(S.t_sam);
    delayed = double(S.delayed_time_rom);

    %% --- HRM root at t0 (finalOP_Raiz) ---
    load_system(name_hrm);
    simRoot = Simulink.SimulationInput(name_hrm);
    simRoot = setVariable(simRoot, 'B', 10, 'Workspace', name_hrm);
    simRoot = setVariable(simRoot, 'A', double(S.sampled_A), 'Workspace', name_hrm);
    simRoot = setVariable(simRoot, 'mdot_c_drift', 10, 'Workspace', name_hrm);
    simRoot = setVariable(simRoot, 'vel_ACS', double(S.sampled_vel_ACS), 'Workspace', name_hrm);
    simRoot = setModelParameter(simRoot, 'StopTime', string(t0));
    simRoot = setModelParameter(simRoot, 'SaveFinalState', 'on', ...
        'SaveOperatingPoint', 'on', 'FinalStateName', 'finalOPRoot');
    simRoot = simRoot.setExternalInput(S.control_series(2:end, :));
    outRootNew = sim(simRoot);
    finalOP_Raiz = outRootNew.finalOPRoot;

    %% --- HRM root2 at t0 - delayed (initial state for ROM) ---
    simRoot2 = simRoot;
    simRoot2 = setModelParameter(simRoot2, 'StopTime', string(t0 - delayed));
    simRoot2 = setModelParameter(simRoot2, 'SaveState', 'off', ...
        'SaveFinalState', 'on', 'SaveOperatingPoint', 'on', ...
        'FinalStateName', 'finalOPRoot2');
    outRoot2 = sim(simRoot2);
    finalOP2_Raiz = outRoot2.finalOPRoot2;
    % Pass 0 (discard changes): without it, close_system opens the interactive
    % "Save contents of model workspace before closing?" dialog, because the
    % ROM/HRM ModelWorkspace has a FILE data source
    % (simulation_core_dualfluid_FI3_sim_rom_extract_MCDET_v2.m). That dialog
    % blocks any -batch run.
    close_system(name_hrm, 0);

    %% --- ROM at t0 - delayed (romOP / outrom) ---
    load_system([name_rom '.slx']);
    simromRoot = Simulink.SimulationInput(name_rom);
    simromRoot = setModelParameter(simromRoot, 'FastRestart', 'off');
    simromRoot = setModelParameter(simromRoot, 'SaveFinalState', 'on', ...
        'SaveOperatingPoint', 'on', 'FinalStateName', 'finalOP');
    simromRoot = setModelParameter(simromRoot, 'StopTime', string(t0 - delayed));
    tt = transpose(0:ts:t0);
    n = numel(tt);
    inport_data = [tt, S.outRoot.mdot_c_leak_signal.Data(:), ...
        zeros(n, 1), zeros(n, 1), zeros(n, 1)];
    simromRoot = simromRoot.setExternalInput(inport_data);
    simromRoot = setVariable(simromRoot, 'vel_ACS', double(S.sampled_vel_ACS), ...
        'Workspace', name_wks);
    outrom = sim(simromRoot);
    romOP = outrom.finalOP;
    close_system(name_rom, 0);   % 0 = discard, avoids the workspace-save dialog

    %% --- Update det_params ---
    det_params = S.det_params;
    det_params.finalOP_Raiz  = finalOP_Raiz;
    det_params.finalOP2_Raiz = finalOP2_Raiz;
    det_params.romOP         = romOP;
    det_params.outrom        = outrom;
    % Keep the archived outRoot as the signal source (same scenario data);
    % also store the regenerated HRM output for provenance.
    det_params.outRoot       = S.outRoot;
    det_params.outRootNew    = outRootNew;
    fprintf('Operating points regenerated: HRM t=%.1f, HRM t=%.1f, ROM t=%.1f\n', ...
        t0, t0 - delayed, t0 - delayed);
end
