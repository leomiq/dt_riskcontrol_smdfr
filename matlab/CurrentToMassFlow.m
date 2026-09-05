%CurrentToMassFlow.m
function m_dot = CurrentToMassFlow(I_input)
    % currentToMassFlow calculates the coolant mass flow rate [kg/s] 
    % as a function of the pump input current [A] using the ECM model
    % from Nowak et al. (2022) with an adapted balanced square channel design.
    % DOI: https://doi.org/10.1016/j.anucene.2022.109142

    % Input:
    %   I_input - Pump input current [A]
    %
    % Output:
    %   m_dot   - Coolant mass flow rate [kg/s]
    %
    % Assumptions and Parameters:
    %   Coolant: Lead (Pb) at 800°C
    %   Density (rho): 10449 kg/m^3
    %   Channel: Square design with height (Hd) = width (Wd) = 0.515 m
    %   Channel length (L): 1.0 m
    %   Wall thickness (d): 0.001 m
    %   Magnetic induction (Bm): 1 T
    %   Loop pressure drop (DeltaP): 15000 Pa
    %   Friction factor (f): 0.02
    %   Resistivity of Pb (q0_LM): 0.118 (ECM units)
    %   Effective resistivity of pump wall (q0_w): 1.2
    %   Fringe factor (K): 1

    % Coolant properties
    rho = 10449;  % [kg/m^3]

    % Geometric parameters (balanced square channel)
    Hd = 0.515;   % Channel height [m]
    Wd = 0.515;   % Channel width [m]
    S  = Wd * Hd; % Cross-sectional area [m^2]
    L  = 1.0;     % Channel length [m]
    d  = 0.001;   % Wall thickness [m]

    % Operating parameters
    Bm     = 1;       % Magnetic induction [T]
    DeltaP = 15000;   % Pressure drop [Pa]
    f      = 0.02;    % Friction factor

    % Electrical properties (ECM parameters)
    q0_LM = 0.118;   % Resistivity of Pb (ECM units)
    q0_w  = 1.2;     % Resistivity of pump wall (assumed for molybdenum)
    K     = 1;       % Fringe factor

    % Compute resistances in the ECM model:
    R_LM = q0_LM * (Wd) / (Hd * L);
    R_w  = q0_w  * (Wd) / (2 * d * L);
    R_b  = q0_LM * (K^2 * Hd) / d;
    R_ver = (R_w * R_b) / (R_w + R_b);  % Parallel combination of R_w and R_b
    R_tot = R_LM + R_ver;
    
    % Compute ECM coefficients:
    % Pump gain coefficient (C) [Pa/A]
    C = (Bm * R_ver * Hd^2) / (R_tot * d);
    
    % Electromagnetic loss coefficient (B) [Pa/(m^3/s)]
    B = (Bm^2 * R_tot * Hd^2) / d;
    
    % Hydrodynamic loss coefficient (A) [Pa/(m^3/s)^2]
    A = (f * rho * L * (Wd + Hd)) / (4 * (S^3));

    % Given the ECM pressure-balance:
    %   DeltaP = C * I_input - B * Q - A * Q^2
    % Rearranged as a quadratic in Q:
    %   A*Q^2 + B*Q + (DeltaP - C * I_input) = 0
    % Solve for Q (volumetric flow rate in m^3/s)
    disc = B^2 - 4 * A * (DeltaP - C * I_input);
    if disc < 0
        error('Negative discriminant: No real solution for I_input = %.2f A', I_input);
    end
    % Use the positive root as the physically relevant solution
    Q = (-B + sqrt(disc)) / (2 * A);
    %physical constraints
    Q = max(Q,0);

    % Calculate coolant mass flow rate [kg/s]
    m_dot = rho * Q;
    
    %physically meaningful
    m_dot = max(m_dot,0);
    
    % Optional: Display the computed values for verification
    %fprintf('Input current: %.2f A\n', I_input);
    %fprintf('Computed volumetric flow rate: %.6f m^3/s\n', Q);
    %fprintf('Coolant mass flow rate: %.2f kg/s\n', m_dot);
end
