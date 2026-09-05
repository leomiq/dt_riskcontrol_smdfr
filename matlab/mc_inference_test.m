% (A) Import your mc_inference module
mod = py.importlib.import_module('mc_inference');
py.importlib.reload(mod);

% (B) Prepare some data for inference (shape: [batch_size, seq_len, 2])
% Let's do a single sample with seq_len=50
seq_len = 50;
x_input_mat = rand(1, seq_len, 2);  % [1,50,2] double in MATLAB

% Convert to Python
x_input_py = py.numpy.array(x_input_mat);
num_trajs = py.int(59); % 3 trajectories

% measure computational time 
tic;
result_py = mod.mc_dropout_trajectories(mod.model, x_input_py, num_trajs, py.str('cuda'));
toc;

% Convert python list -> MATLAB double
result_mat = double(result_py);
