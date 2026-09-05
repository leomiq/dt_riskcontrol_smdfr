%CalculateTs.m
function Ts_required = CalculateTs(t_open, T_ave)
% calculateTs computes the required fuel boundary temperature Ts
% to achieve a given fuse plug opening time.
%
% Inputs:
%   t_open - Desired opening time (s), e.g., 100 seconds
%   T_ave  - Average initial fuse plug temperature (K), e.g., 545 K
%
% Output:
%   Ts_required - Required fuel boundary temperature (K)
%
% Constants
T_mp = 733;       % K, melting point of FLiBe (molten salt)
L    = 0.025;     % m, fuse plug thickness (25 mm)

% Geometry (SM-DFR)
D_inner = 0.05;   % m, inner diameter of the tube
D_outer = 0.1;    % m, outer diameter of the tube
h_plug  = L;      % m, fuse plug length (equal to thickness)
A_salt  = pi * (D_inner^2) / 4;             % m², cross-sectional area for salt
A_wall  = pi * ((D_outer^2 - D_inner^2)) / 4; % m², cross-sectional area for the pipe wall

% Material properties for FLiBe (fuse plug)
rho_plug = 2055;  % kg/m³, density of FLiBe
c_plug   = 2414;  % J/(kg·K), specific heat of FLiBe
Lh       = 424700; % J/kg, latent heat of fusion of FLiBe

% Material properties for the pipe wall (Hastelloy-N)
rho_pipe = 8987;  % kg/m³
c_pipe   = 577.8; % J/(kg·K)

% Thermal conductivities (assumed at the melting point for salt)
k_salt = 0.629 + 0.0005 * T_mp;  % W/m·K, thermal conductivity of salt
k_wall = 21.96;                % W/m·K, thermal conductivity of Hastelloy-N

% Compute volumes
V_plug = A_salt * h_plug;   % m³, volume of the fuse plug
V_pipe = A_wall * h_plug;   % m³, volume of the pipe wall section

% Total energy required to heat and melt the fuse plug and heat the pipe wall.
% For the fuse plug: sensible heating from T_ave to T_mp plus latent heat.
Q_plug = rho_plug * V_plug * ( c_plug * (T_mp - T_ave) + Lh );
% For the pipe wall: sensible heating from T_ave to T_mp.
Q_pipe = rho_pipe * V_pipe * ( c_pipe * (T_mp - T_ave) );
Q_total = Q_plug + Q_pipe;  % Total energy required (J)

% Conduction term (effective heat transfer coefficient)
conduction_term = k_wall * A_wall + k_salt * A_salt;  % W/K

% Rearranged equation:
% t_open = Q_total / ( (conduction_term*(T_s - T_mp)) / L )
% Solve for T_s:
Ts_required = T_mp + (Q_total * L) / (t_open * conduction_term);

end
