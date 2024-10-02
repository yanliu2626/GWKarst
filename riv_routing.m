function [Q,h_update] = riv_routing(channel_info,Qex,QnK,Qkarst,QrunoffK, ...
    riv_k,grid_nK_riv_k,grid_K_riv_k,bnd_gridK_ind_riv,grid_sprriv_riv_k,Q_ini,Q_tri_mat)
%call the kenimatic wave to solve flow in the river channel
%   Detailed explanation goes here
mat_alpha = channel_info.mat_alpha;
mat_x = channel_info.mat_x;
chn_w = channel_info.chn_w;
n_reach = channel_info.n_reach;

Q_rivCon_grid_nK  = zeros(n_reach,1);
Q_rivCon_gridK = zeros(n_reach,1);
Q_rivCon_sprriv = zeros(n_reach,1);
Q_rivCon_runoffK = zeros(n_reach,1);

for i=1:n_reach
    % river connected to non-karst grids
    Q_rivCon_grid_nK(i,1) = sum(QnK(grid_nK_riv_k==riv_k(i),1));

    % river connected to karst aquifer
    Q_rivCon_gridK(i,1) = sum(Qex(bnd_gridK_ind_riv==riv_k(i),1));

    % river connected to karst spring rivers
    Q_rivCon_sprriv(i,1) = sum(Qkarst(grid_sprriv_riv_k==riv_k(i),1));
    
    % surface runoff from karst grids
    Q_rivCon_runoffK(i,1) = sum(QrunoffK(grid_K_riv_k==riv_k(i),1));

end

Q_rivCon = Q_rivCon_grid_nK+Q_rivCon_gridK+Q_rivCon_sprriv+Q_rivCon_runoffK+Q_tri_mat;

mat_q = Q_rivCon./mat_x; % m2/s
% Q_up = Q_rivCon(1); % to be modified
Q_up = 0; % to be modified
[Q] = kinematic_wave(mat_alpha,mat_x,mat_q,Q_ini,Q_up);
A_update = mat_alpha.*Q.^0.6;
h_update = A_update./chn_w;
end

