%% DATA EXTRACTION FROM PREVIOUS SIMULATED ACCIDENTAL SCENARIOS

% Folder holding the failure-injection database. These sweeps are raw run
% output and are not redistributed with this repository; set the environment
% variable FAILURE_INJECTION_DIR to their location, or edit the fallback below.
folderPath = getenv('FAILURE_INJECTION_DIR');
if isempty(folderPath)
    folderPath = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                          'data', 'failure_injection');
end
if ~endsWith(folderPath, filesep)
    folderPath = [folderPath filesep];
end

% Specify the name of the .mat file you want to load
fileName = 'database';  % Change to your file's name

%database: index of simulation, time, mdot_c, T_p
db_Tp_mdotc = zeros(1,4);

index_nsim = 1;

for i = 1:5
    text_chain = [folderPath fileName num2str(i) '.mat'];

    load(text_chain);

    for j = 1:N_par
        fault_time_index = find(array_data_tpeak_hrm{j}(:,1)==round(list_fault_time(j,2)));
        time_indexes = fault_time_index:5:length(array_data_tpeak_hrm{j}(:,1));
        data_delta_tpeak_2 = array_data_tpeak_hrm{j}(time_indexes,3) - array_data_tpeak_rom_v2{j}(time_indexes,3);
        db_Tp_mdotc = cat(1,db_Tp_mdotc, [index_nsim*ones(1,length(time_indexes))' ...
                                        time_indexes' ...
                                        array_data_tpeak_rom_v2{j}(time_indexes,2) ...
                                        data_delta_tpeak_2]);
        index_nsim = index_nsim + 1;
    end

end

db_Tp_mdotc = db_Tp_mdotc(2:end,:);

x = db_Tp_mdotc(:,3);
y = db_Tp_mdotc(:,4);

scatter(x,y)


%% TESTING FIT (V2)

% For example, extracting from a matrix as:
x = db_Tp_mdotc(:,3);
y = db_Tp_mdotc(:,4);

% Fit a third-order polynomial to the data
fit_poly = fit(x, y, 'poly4');

% Plot the fitted curve along with the data points
h = plot(fit_poly, x, y);

% Increase the line width of the fitted curve
% Note: In the output of plot with a fit object, the fitted curve is usually the first handle,
% and the data points are the subsequent handles. Adjust if necessary.
set(h(1), 'Marker','+')
set(h(2), 'LineWidth', 3);  % Adjust the value (3) as needed for desired thickness

% Add labels with LaTeX interpreter
xlabel('$\dot{m}_{\kappa, \rm leak}$ [kg/s]', 'Interpreter', 'latex','FontSize',14);
ylabel('$\Delta T_{\rm p}$ [K]', 'Interpreter', 'latex','FontSize',14);

