function check_torch_env()
% CHECK_TORCH_ENV  Verify the MATLAB Python environment needed to re-run the
% BiLSTM MC-dropout inference (required to regenerate a severe accident case
% from the HRM rather than varying ROM parameters, which has no effect --
% vel_ACS is an HRM scenario parameter and the ROM receives the ACS action
% through its u_ACS inport).

    fprintf('\n===== MATLAB Python / torch environment =====\n');

    try
        pe = pyenv;
        fprintf('Version    : %s\n', string(pe.Version));
        fprintf('Executable : %s\n', string(pe.Executable));
        fprintf('Library    : %s\n', string(pe.Library));
        fprintf('Status     : %s\n', string(pe.Status));
        fprintf('ExecutionMode: %s\n', string(pe.ExecutionMode));
    catch ME
        fprintf('pyenv unavailable: %s\n', ME.message);
        return;
    end

    fprintf('\n-- torch --\n');
    try
        t = py.importlib.import_module('torch');
        try
            fprintf('torch version   : %s\n', string(py.getattr(t, char([95 95 'version' 95 95]))));
        catch
            fprintf('torch version   : <unreadable>\n');
        end
        cu = t.cuda.is_available();
        fprintf('cuda available  : %s\n', string(cu));
        if cu
            fprintf('device count    : %d\n', double(t.cuda.device_count()));
            fprintf('device name     : %s\n', string(t.cuda.get_device_name(int32(0))));
        end
    catch ME
        fprintf('torch import FAILED: %s\n', ME.message);
        fprintf('\n>> Without torch the BiLSTM MC-dropout inference cannot run,\n');
        fprintf('   so a new accident scenario cannot be generated from the HRM.\n');
        return;
    end

    fprintf('\n-- mc_inference module --\n');
    try
        m = py.importlib.import_module('mc_inference');
        fprintf('mc_inference imported OK\n');
        try
            mdl = m.model;
            fprintf('m.model present : %s\n', string(class(mdl)));
        catch
            fprintf('m.model not exposed as an attribute\n');
        end
        d = py.dir(m);
        names = {};
        for i = 1:double(py.len(d))
            nm = char(string(d{i}));
            if isempty(nm) || nm(1) ~= char(95)
                names{end+1} = nm; %#ok<AGROW>
            end
        end
        fprintf('public attributes: %s\n', strjoin(names, ', '));
    catch ME
        fprintf('mc_inference import FAILED: %s\n', ME.message);
        return;
    end

    fprintf('\n-- smoke test: one MC-dropout call --\n');
    try
        X = rand(200, 2);
        m = py.importlib.import_module('mc_inference');
        Xpy = py.numpy.array(X);
        tic;
        r = m.mc_dropout_trajectories(m.model, Xpy, py.int(5));
        el = toc;
        Y = double(py.numpy.array(r));
        fprintf('OK in %.2f s -> output size %s\n', el, mat2str(size(Y)));
        fprintf('>> The inference path works. A severe HRM scenario can be generated.\n');
    catch ME
        fprintf('smoke test FAILED: %s\n', ME.message);
    end
    fprintf('=============================================\n\n');
end
