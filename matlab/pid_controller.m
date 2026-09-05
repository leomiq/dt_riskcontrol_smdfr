% pid_controller.m
function [optimal_current, pid_info] = pid_controller(current_state, current_control, dt, params)
    % PID controller for PCT regulation
    
    % Extract PCT from state
    
    PCT = current_state;
    
    % PID parameters
    Kp = params.Kp;
    Ki = params.Ki; 
    Kd = params.Kd;
    
    % Reference and error
    PCT_ref = params.PCT_ref;
    error = PCT_ref - PCT;
    
    % Persistent variables for integral and derivative terms
    persistent integral_term last_error;
    if isempty(integral_term)
        integral_term = 0;
        last_error = 0;
    end
    
    % PID terms
    proportional = Kp * error;
    integral_term = integral_term + Ki * error * dt;
    derivative = Kd * (error - last_error) / dt;
    
    % PID output
    pid_output = proportional + integral_term + derivative;
    
    % Convert to current (you may need scaling)
    optimal_current = current_control + pid_output;
    
    % Anti-windup and saturation
    optimal_current = max(min(optimal_current, params.I_max), params.I_min);
    
    % Rate limiting
    max_change = params.dI_max * dt;
    optimal_current = current_control + ...
        sign(optimal_current - current_control) * ...
        min(abs(optimal_current - current_control), max_change);
    
    % Update persistent variables
    last_error = error;
    
    % Store PID info
    pid_info.error = error;
    pid_info.proportional = proportional;
    pid_info.integral = integral_term;
    pid_info.derivative = derivative;
    pid_info.PCT_current = PCT;
end

