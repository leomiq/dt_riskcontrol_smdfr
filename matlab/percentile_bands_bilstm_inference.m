% Assuming Y_denorm is a 3D array: (n_trials, sequence_length, n_features)
% Replace with your actual Y_denorm data
% Example dimensions: [59 trials, 100 timesteps, 2 features]
[n_trials, sequence_length, n_features] = size(Y_denorm);

% Create a figure for the plots
figure;
for feature_idx = 1:n_features
    % Extract data for the current feature
    feature_data = squeeze(Y_denorm(:, :, feature_idx)); % (n_trials x sequence_length)
    
    % Calculate statistics
    mean_values = mean(feature_data, 1);
    lower_percentile = prctile(feature_data, 5, 1); % 5th percentile
    upper_percentile = prctile(feature_data, 95, 1); % 95th percentile
    
    % Create a subplot for each feature
    subplot(n_features, 1, feature_idx);
    hold on;
    
    % Plot the mean
    plot(1:sequence_length, mean_values, 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
    
    % Plot the 5th and 95th percentiles as a shaded region
    fill([1:sequence_length, fliplr(1:sequence_length)], ...
         [upper_percentile, fliplr(lower_percentile)], ...
         'r', 'FaceAlpha', 0.9, 'EdgeColor', 'b', 'DisplayName', '5th-95th Percentile');
    
    % Add labels and title
    title(['Feature ', num2str(feature_idx)]);
    xlabel('Timestep');
    ylabel('Value');
    legend('show');
    grid on;
    
    hold off;
end

% Adjust layout for better visualization
sgtitle('Mean, 5th and 95th Percentiles for Each Feature');

%%
