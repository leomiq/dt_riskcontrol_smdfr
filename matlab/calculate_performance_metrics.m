function metrics = calculate_performance_metrics(time, PCT, control_signal, params)
    % Calculate comprehensive performance metrics
    % Inputs:
    %   time: time vector [s]
    %   PCT: Peak Cladding Temperature vector [K]
    %   control_signal: control input vector [A]
    %   params: structure with parameters
    
    metrics = struct();
    
    % Extract parameters
    PCT_ref = params.PCT_ref;
    PCT_max = params.PCT_max;
    tolerance = params.tolerance; % e.g., 0.02 for 2%
    
    % Error calculation
    error = PCT_ref - PCT;
    
    % 1. MINIMUM SAFETY MARGIN
    safety_margin = PCT_max - PCT;
    metrics.min_safety_margin = min(safety_margin);
    metrics.safety_violation = metrics.min_safety_margin < 0;
    
    % 2. MAXIMUM PCT OVERSHOOT
    metrics.min_overshoot = min(PCT - PCT_ref);
    metrics.percent_overshoot = (metrics.min_overshoot / PCT_ref) * 100;
    
    % 3. SETTLING TIME (with tolerance)
    settling_threshold = tolerance;
    settled_indices = find(abs(error) <= settling_threshold);
    
    if ~isempty(settled_indices)
        % Find the point after which it never leaves the tolerance band
        last_unsettled = find(abs(error) > settling_threshold, 1, 'last');
        if isempty(last_unsettled)
            metrics.settling_time = time(settled_indices(1));
        elseif last_unsettled + 1 <= numel(time)
            metrics.settling_time = time(last_unsettled + 1);
        else
            % FIX 2026-09-04: the trajectory is still outside the tolerance band
            % at the LAST sample, so there is no settling instant within the
            % recorded horizon and time(last_unsettled+1) is out of range.
            % Previously this threw MATLAB:badsubscript and destroyed the KPI
            % computation after a completed 61-query realization (the widened-
            % bound penalised run ends over-cooled and never re-enters the band).
            metrics.settling_time = inf;
        end
    else
        metrics.settling_time = inf;
    end
    
    % 4. INTEGRAL ABSOLUTE ERROR (IAE)
    metrics.IAE = trapz(time, abs(error));
    
    % 5. CUMULATIVE ERROR (Integral of Error)
    metrics.cumulative_error = trapz(time, error);
    
    % Additional useful metrics
    metrics.steady_state_error = mean(error(round(0.8*end):end));
    metrics.peak_PCT = max(PCT);
    
    % Control effort metrics
    control_changes = diff(control_signal);
    metrics.total_variation = sum(abs(control_changes));
    metrics.max_control_rate = max(abs(control_changes)) / mean(diff(time));
    
    % Display results
    fprintf('\n--- PERFORMANCE METRICS ---\n');
    fprintf('Minimum Safety Margin: %.2f K %s\n', metrics.min_safety_margin);
    fprintf('Minimum Overshoot: %.2f K (%.1f%%)\n', metrics.min_overshoot, metrics.percent_overshoot);
    fprintf('Settling Time (<%.2f K): %.2f s\n', tolerance, metrics.settling_time);
    fprintf('Integral Absolute Error (IAE): %.2f K·s\n', metrics.IAE);
    fprintf('Cumulative Error: %.2f K·s\n', metrics.cumulative_error);
    fprintf('Steady-State Error: %.2f K\n', metrics.steady_state_error);
    fprintf('Peak PCT: %.2f K\n', metrics.peak_PCT);
    fprintf('Control Total Variation: %.2f A\n', metrics.total_variation);
end