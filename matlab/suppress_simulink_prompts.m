function old = suppress_simulink_prompts()
% SUPPRESS_SIMULINK_PROMPTS  Turn off the interactive Simulink dialogs that
% block an automated (-batch) run of these experiments.
%
% Three prompts show up in this workflow:
%
%  1. "Save contents of model workspace before closing?" -- raised by
%     close_system on a model whose ModelWorkspace has a FILE data source
%     (here: simulation_core_dualfluid_FI3_sim_rom_extract_MCDET_v2.m). Any
%     close_system called WITHOUT the second argument asks. The fix is to pass
%     0 (discard) explicitly; this function also clears the dirty flag so the
%     question does not arise.
%
%  2. Fast-restart / recompile confirmations -- raised when FastRestart is
%     toggled on a model that must be recompiled, or when an operating point
%     from another release is loaded.
%
%  3. NonCurrentReleaseOperatingPointMsg -- the archived operating points were
%     saved under a different release; left as 'warning' rather than 'error' so
%     a run does not abort.
%
% Call at the top of any batch script, and pass the returned struct to
% RESTORE_SIMULINK_PROMPTS at the end if the session is interactive.
%
%   old = suppress_simulink_prompts();
%   ... run experiments ...
%   restore_simulink_prompts(old);

    old = struct();

    % --- warnings that surface as dialogs / noise ---
    ids = { 'Simulink:Engine:ModelWorkspaceNotSaved', ...
            'Simulink:Commands:CloseSystemUnsavedChanges', ...
            'Simulink:Harness:WarnABoutSaveWithHarness', ...
            'Simulink:Engine:NonCurrentReleaseOperatingPoint', ...
            'Simulink:Engine:OperatingPointRelease', ...
            'MATLAB:MKDIR:DirectoryExists' };
    old.warnstates = struct('id', {}, 'state', {});
    for i = 1:numel(ids)
        try
            s = warning('query', ids{i});
            old.warnstates(end+1) = struct('id', ids{i}, 'state', s.state);
            warning('off', ids{i});
        catch
        end
    end
    % blanket: torch/py FutureWarnings and Simulink find_system variant notices
    old.all = warning('query', 'all');

    % --- per-model configuration ---
    models = { 'rom_MCDET_v2_2022b_control_CMEX', ...
               'SMDFR_hrm_MCDET_2022b_control_CMEX' };
    old.models = struct('name', {}, 'params', {});
    for i = 1:numel(models)
        nm = models{i};
        if ~bdIsLoaded(nm), continue; end
        p = struct();
        % Accept operating points saved by another release with a warning
        % instead of an error-level dialog.
        for prm = {'NonCurrentReleaseOperatingPointMsg'}
            try
                p.(prm{1}) = get_param(nm, prm{1});
                set_param(nm, prm{1}, 'warning');
            catch
            end
        end
        % Clear the dirty flag so close_system has nothing to ask about.
        try
            p.Dirty = get_param(nm, 'Dirty');
            set_param(nm, 'Dirty', 'off');
        catch
        end
        old.models(end+1) = struct('name', nm, 'params', p);
    end

    fprintf('[suppress_simulink_prompts] %d warning ids silenced, %d model(s) configured\n', ...
        numel(old.warnstates), numel(old.models));
end
