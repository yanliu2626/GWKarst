function [Qriv] = kinematic_wave(mat_alpha,mat_x,mat_q,Q_ini,Q_up)
%solve kinematic wave propogation using explicit scheme ode45, which uses
%adaptive time step to control the stability (Courant–Friedrichs–Lewy condition)

%   channel settings
%   mat_alpha: kinematic wave parameter, Q and A relationship
%   mat_x: spacing between two cross sections [m]

%   constant parameter
%   m: kinematic wave parameter, Q and A relationship

%   boundary conditions
%   mat_q: lateral flow [m3/m] negative: losing water, positive: gain water
%   Q0: discharge at the starting point of the river [m3/s]

m = 3/5;

% ode settings
dt_out = 86400;
options=odeset('reltol',1e-5,'abstol',1e-5,'maxstep',dt_out);

% initial condition 
% A_ini: the cross sectional area of the last time step

% boundary condition
% A_up

% call ode
% only calculate one day, since the interaction between river bed and
% aquifer (corresponding cell) will be updated every time step
[~,Q_all] = ode45(@ode_kinematic_wave,0:43200:86400,Q_ini,options);  
% A_update = Q_all(end,:)';

% calculate Q
Qriv = Q_all(end,:)';

% ode function for cross sectional area
    function dQdt = ode_kinematic_wave(~,Q)
        % reconstruct A to include the first cross section A0
        Q_new = [Q_up;Q];
        
        % calculate dQdt
        dQdt = (mat_q-(Q_new(2:end,1)-Q_new(1:end-1,1))./mat_x)./(mat_alpha.*m.*((Q_new(2:end,1)+Q_new(1:end-1,1))/2).^(m-1));
        
    end
end

