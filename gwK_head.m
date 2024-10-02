function [Qex,Qkarst,Qmatrix,Qconduit,QexcessK,hK,hK_conduit] = gwK_head(hK_ini, ...
    hK_conduit_ini,rech_matrix,rech_conduit,h_riv_ini,gridK_info,sprriv_info, ...
    river_info_one,S,nK)
%update the storage of karst grids and calculate karst spring discharge and
%the exchange between karst aquifer and the river

% Input variables
% hK_ini: the initial hydraulic head of karst grids (m.a.s.l). [nK,1]
% rechK_t: recharge depth (m) of karst grids. [nK,1]
% h_riv_ini: the initial river stage (water depth m). [n_riv,1]
% gridK_info: all information about karst grids
%   gridK_info.ne_K: effective porosity of karst grids (-). [nK,1]
%   gridK_info.bnd_grid_ind: the index of connected karst cells (l,t,b,r). [nK,4]
%   gridK_info.gridK_ind: the index (1d) of karst grid in the study domain
%   (including 0,1,2 for out of domain,non-karst,karst)
%   gridK_info.A_gridK: the area of karst grid (km2). [nK,1]
%   gridK_info.hK_bot: the bottom of karst aquifer (m.a.s.l). [nK,1]
% river_info: all information about river reaches
%   river_info.bnd_ind_riv: the index of the cell that river is connected. [n_riv,1]
%   river_info.C_rb: conductance of river bed sediment (m2/s). [n_riv,1]
%   river_info.R_rb: bottom elevation of river bed (m.a.s.l). [n_riv,1]
%   river_info.A_rbK: bottom elevation of river bed (m.a.s.l). [n_riv,1]
% S: calculated as ne.*A (m2). [nK+n_riv,1]
% nK: the number of karst grids
% grid_domain_value: the index (0,1,2) of all grids

% output variables
% Qex: exchange dischange between river and grids. [n_riv,1]
% hK: updated hydraulic head. [nK,1]

% effective porosity of karst
ne_K = gridK_info.ne_K;
k_C = gridK_info.k_C;
K_k = gridK_info.K_k; %m/s
A_gridK = gridK_info.A_gridK;
hK_bot = gridK_info.hK_bot;
hK_top = gridK_info.hK_top;

% first update storage adding recharge
hK_ini = hK_ini + rech_matrix./ne_K;
hKover = max(0,hK_ini - hK_top);
hK_ini = min(hK_ini,hK_top);
hK_conduit_ini = hK_conduit_ini + rech_conduit;
hK_conduit = hK_conduit_ini;

% calculte matrix flow and flow direction for conduit flow
% calculate the between-grid flow and exchange with rivers
% configuration of karst grid properties

bnd_grid_ind = gridK_info.bnd_grid_ind;
bnd_ind_sprriv = sprriv_info.bnd_ind_sprriv;
bnd_ind_riv = river_info_one.bnd_ind_riv;
uni_bnd_ind_riv = river_info_one.uni_bnd_ind_riv;
uni_bnd_ind_riv_mat = river_info_one.uni_bnd_ind_riv_mat;
n_rivseg_uni = river_info_one.n_rivseg_uni;

R_sprriv = sprriv_info.R_sprriv; % elevation of the spring river bottom
C_sprriv = sprriv_info.C_sprriv; % conductance of spring river
n_sprriv = sprriv_info.n_sprriv; % number of cells with spring river

C_rb = river_info_one.C_rb;
R_rb = river_info_one.R_rb;
bnd_ind_l = bnd_grid_ind(:,1);
bnd_ind_t = bnd_grid_ind(:,2);
bnd_ind_b = bnd_grid_ind(:,3);
bnd_ind_r = bnd_grid_ind(:,4);

% ode settings
dt_out = 1;
options=odeset('reltol',1e-4,'abstol',1e-4,'NonNegative',1,'maxstep',dt_out);

% initial head of the karst storage and river
h_sprriv_ini = zeros(n_sprriv,1);
h_ini = [hK_ini-hK_bot;h_sprriv_ini;h_riv_ini];

l_k = ones(nK,1)*2500;% m 
w_k_l = ones(nK,1)*5000; %m
w_k_t = w_k_l;
w_k_b = w_k_l;
w_k_r = w_k_l;

% call ode
% only calculate one day, since the interaction between river bed and
% aquifer (corresponding cell) will be updated every time step
% h_all is head above the bottom of aquifer and river bottom
[~,h_all] = ode45(@ode_gw_riv,0:0.5:1,h_ini,options);  
haboveK = h_all(end,1:nK)';
hK = haboveK+hK_bot;

% head in karst grids cannot be larger than its surface elevation
habovesurK = max(0,hK - hK_top);
hK = min(hK,hK_top);

% river water cannot be smaller than zero
h_sprriv = h_all(end,nK+1:nK+n_sprriv)';
riv_stage = h_all(end,nK+n_sprriv+1:end)';
S_sprrivK = S(nK+1:nK+n_sprriv,1);
S_riv = S(nK+n_sprriv+1:end,1);
Qmatrix = max(0,(h_sprriv-h_sprriv_ini).*S_sprrivK/86400); % [m3/s];
Qex = (riv_stage-h_riv_ini).*S_riv/86400; % [m3/s];

% calculate conduit flow and update conduit storage, assume conduit flow
% distribution fraction is only based on matrix gradient
hK_diff = nan(nK,4);
hK_neighbor_l = nan(nK,1);
hK_neighbor_t = nan(nK,1);
hK_neighbor_b = nan(nK,1);
hK_neighbor_r = nan(nK,1);

hK_neighbor_bot_l = nan(nK,1);
hK_neighbor_bot_t = nan(nK,1);
hK_neighbor_bot_b = nan(nK,1);
hK_neighbor_bot_r = nan(nK,1);

hK_neighbor_l(bnd_ind_l>0,1) = hK(bnd_ind_l(bnd_ind_l>0),1);
hK_neighbor_t(bnd_ind_t>0,1) = hK(bnd_ind_t(bnd_ind_t>0),1);
hK_neighbor_b(bnd_ind_b>0,1) = hK(bnd_ind_b(bnd_ind_b>0),1);
hK_neighbor_r(bnd_ind_r>0,1) = hK(bnd_ind_r(bnd_ind_r>0),1);

hK_neighbor_bot_l(bnd_ind_l>0,1) = hK_bot(bnd_ind_l(bnd_ind_l>0),1);
hK_neighbor_bot_t(bnd_ind_t>0,1) = hK_bot(bnd_ind_t(bnd_ind_t>0),1);
hK_neighbor_bot_b(bnd_ind_b>0,1) = hK_bot(bnd_ind_b(bnd_ind_b>0),1);
hK_neighbor_bot_r(bnd_ind_r>0,1) = hK_bot(bnd_ind_r(bnd_ind_r>0),1);

% hK_diff >0 indicating upstream (flow to downstream), <0 downstream, nan no flow
hK_diff_l = hK-hK_neighbor_l;
hK_diff_l(hK_neighbor_bot_l>hK,1) = hK_neighbor_bot_l(hK_neighbor_bot_l>hK,1)-hK_neighbor_l(hK_neighbor_bot_l>hK,1);
hK_diff_l(hK_bot>hK_neighbor_l,1) = hK(hK_bot>hK_neighbor_l,1)-hK_bot(hK_bot>hK_neighbor_l,1);

hK_diff_t = hK-hK_neighbor_t;
hK_diff_t(hK_neighbor_bot_t>hK,1) = hK_neighbor_bot_t(hK_neighbor_bot_t>hK,1)-hK_neighbor_t(hK_neighbor_bot_t>hK,1);
hK_diff_t(hK_bot>hK_neighbor_t,1) = hK(hK_bot>hK_neighbor_t,1)-hK_bot(hK_bot>hK_neighbor_t,1);

hK_diff_b = hK-hK_neighbor_b;
hK_diff_b(hK_neighbor_bot_b>hK,1) = hK_neighbor_bot_b(hK_neighbor_bot_b>hK,1)-hK_neighbor_b(hK_neighbor_bot_b>hK,1);
hK_diff_b(hK_bot>hK_neighbor_b,1) = hK(hK_bot>hK_neighbor_b,1)-hK_bot(hK_bot>hK_neighbor_b,1);

hK_diff_r = hK-hK_neighbor_r;
hK_diff_r(hK_neighbor_bot_r>hK,1) = hK_neighbor_bot_r(hK_neighbor_bot_r>hK,1)-hK_neighbor_r(hK_neighbor_bot_r>hK,1);
hK_diff_r(hK_bot>hK_neighbor_r,1) = hK(hK_bot>hK_neighbor_r,1)-hK_bot(hK_bot>hK_neighbor_r,1);

hK_diff(:,1) = hK_diff_l;
hK_diff(:,2) = hK_diff_t;
hK_diff(:,3) = hK_diff_b;
hK_diff(:,4) = hK_diff_r;

grid_flowdir = hK_diff;
grid_flowdir(grid_flowdir>0) = 1;
grid_flowdir(grid_flowdir<0) = 2;
grid_flowdir(isnan(grid_flowdir)) = 0;

% separation factor calculation
Qdown_fsep_h = zeros(nK,4);
Qdown_fsep_h(hK_diff>0) = hK_diff(hK_diff>0);
Qdown_fsep = Qdown_fsep_h./repmat(sum(Qdown_fsep_h,2),1,4);

Qcon = k_C .* hK_conduit .* A_gridK/86400; % m3/s
QexcessK = (hKover+ habovesurK).*ne_K .* A_gridK/86400; % m3/s

hK_conduit = max(0,hK_conduit-k_C .* hK_conduit);

% calculate separation factor for downstream cells
[Qcond] = Qconduit_redis(Qcon,grid_flowdir,Qdown_fsep,bnd_grid_ind);
Qconduit = Qcond(bnd_ind_sprriv);

% distribute Q in the conduit that does not flow to sprriv, check within
% each karst aquifer
% set cells already connected to sprriv to 0
Qcond(bnd_ind_sprriv) = 0;

sprriv_aquifer = sprriv_info.sprriv_aquifer;
grid_karst_aquifer = sprriv_info.grid_karst_aquifer;
unique_aquiferK = sprriv_info.unique_aquiferK;
nK_aquiferK = sprriv_info.nK_aquiferK;

for ii=1:nK_aquiferK
    num_aquifer = unique_aquiferK(ii);
    ind_aquifer = (grid_karst_aquifer == num_aquifer);
    Qremain = sum(Qcond(ind_aquifer));
    
    ind_sprriv_aquifer = (sprriv_aquifer==num_aquifer);
    n_sprriv_aquifer = sum(ind_sprriv_aquifer);
    Qconduit(ind_sprriv_aquifer) = Qconduit(ind_sprriv_aquifer)+Qremain/n_sprriv_aquifer;
    
end

% discharge from karst
Qkarst = Qmatrix + Qconduit;

% matrix flow from karst to spring is represented by assuming flow for
% matrix to rivers (one direction), no flow from river to matrix which is
% the same as only possible of water from matrix to spring, the other
% direction is no occuring
% hini for spring river set to a fixed value and update to this value every
% time step

% ode function head change
    function dhdt = ode_gw_riv(~,h)
        habove_k = h(1:nK,1);
        h_k = habove_k+hK_bot;
        h_riv = h(nK+n_sprriv+1:end,1);
                
        % connection with left
        dh_left = zeros(nK,1);
        C_l = zeros(nK,1);
        
        thicknessn_l = habove_k(bnd_ind_l>0,1);% thickness of aquifer
        thicknessm_l = habove_k(bnd_ind_l(bnd_ind_l>0),1);
        hn_l = h_k(bnd_ind_l>0,1);
        hm_l = h_k(bnd_ind_l(bnd_ind_l>0),1);
        hn_bot_l = hK_bot(bnd_ind_l>0,1);
        hm_bot_l = hK_bot(bnd_ind_l(bnd_ind_l>0),1);
        ln_l = l_k(bnd_ind_l>0,1);
        lm_l = l_k(bnd_ind_l(bnd_ind_l>0),1);
        Kn_l  = K_k(bnd_ind_l>0,1);
        Km_l  = K_k(bnd_ind_l(bnd_ind_l>0),1);
        w_nm_l = w_k_l(bnd_ind_l>0,1);
        [C_nm_l] = conductance_karst(w_nm_l,Kn_l,Km_l,ln_l,lm_l,thicknessn_l,thicknessm_l);
        C_l(bnd_ind_l>0,1) = C_nm_l;
        dh_l = hm_l - hn_l;
        dh_l(hm_bot_l>hn_l,1) = hm_l(hm_bot_l>hn_l,1) - hm_bot_l(hm_bot_l>hn_l,1);
        dh_l(hn_bot_l>hm_l,1) = hn_bot_l(hn_bot_l>hm_l,1) - hn_l(hn_bot_l>hm_l,1);
        dh_left(bnd_ind_l>0,1) = dh_l;
        dQdt_left = C_l .* dh_left;
        
        % connection with top
        dh_top = zeros(nK,1);
        C_t = zeros(nK,1);
        
        thicknessn_t = habove_k(bnd_ind_t>0,1);% thickness of aquifer
        thicknessm_t = habove_k(bnd_ind_t(bnd_ind_t>0),1);
        hn_t = h_k(bnd_ind_t>0,1);
        hm_t = h_k(bnd_ind_t(bnd_ind_t>0),1);
        hn_bot_t = hK_bot(bnd_ind_t>0,1);
        hm_bot_t = hK_bot(bnd_ind_t(bnd_ind_t>0),1);
        ln_t = l_k(bnd_ind_t>0,1);
        lm_t = l_k(bnd_ind_t(bnd_ind_t>0),1);
        Kn_t  = K_k(bnd_ind_t>0,1);
        Km_t  = K_k(bnd_ind_t(bnd_ind_t>0),1);
        w_nm_t = w_k_t(bnd_ind_t>0,1);
        [C_nm_t] = conductance_karst(w_nm_t,Kn_t,Km_t,ln_t,lm_t,thicknessn_t,thicknessm_t);
        C_t(bnd_ind_t>0,1) = C_nm_t;
        dh_t = hm_t - hn_t;
        dh_t(hm_bot_t>hn_t,1) = hm_t(hm_bot_t>hn_t,1) - hm_bot_t(hm_bot_t>hn_t,1);
        dh_t(hn_bot_t>hm_t,1) = hn_bot_t(hn_bot_t>hm_t,1) - hn_t(hn_bot_t>hm_t,1);
        dh_top(bnd_ind_t>0,1) = dh_t;
        dQdt_top = C_t .* dh_top;

               
        % connection with bottom
        dh_bottom = zeros(nK,1);
        C_b = zeros(nK,1);
        
        thicknessn_b = habove_k(bnd_ind_b>0,1);% thickness of aquifer
        thicknessm_b = habove_k(bnd_ind_b(bnd_ind_b>0),1);
        hn_b = h_k(bnd_ind_b>0,1);
        hm_b = h_k(bnd_ind_b(bnd_ind_b>0),1);
        hn_bot_b = hK_bot(bnd_ind_b>0,1);
        hm_bot_b = hK_bot(bnd_ind_b(bnd_ind_b>0),1);
        ln_b = l_k(bnd_ind_b>0,1);
        lm_b = l_k(bnd_ind_b(bnd_ind_b>0),1);
        Kn_b  = K_k(bnd_ind_b>0,1);
        Km_b  = K_k(bnd_ind_b(bnd_ind_b>0),1);
        w_nm_b = w_k_b(bnd_ind_b>0,1);
        [C_nm_b] = conductance_karst(w_nm_b,Kn_b,Km_b,ln_b,lm_b,thicknessn_b,thicknessm_b);
        C_b(bnd_ind_b>0,1) = C_nm_b;
        dh_b = hm_b - hn_b;
        dh_b(hm_bot_b>hn_b,1) = hm_b(hm_bot_b>hn_b,1) - hm_bot_b(hm_bot_b>hn_b,1);
        dh_b(hn_bot_b>hm_b,1) = hn_bot_b(hn_bot_b>hm_b,1) - hn_b(hn_bot_b>hm_b,1);
        dh_bottom(bnd_ind_b>0,1) = dh_b;
        dQdt_bottom = C_b .* dh_bottom;
        
        % connection with right
        dh_right = zeros(nK,1);
        C_r = zeros(nK,1);
        
        thicknessn_r = habove_k(bnd_ind_r>0,1);% thickness of aquifer
        thicknessm_r = habove_k(bnd_ind_r(bnd_ind_r>0),1);
        hn_r = h_k(bnd_ind_r>0,1);
        hm_r = h_k(bnd_ind_r(bnd_ind_r>0),1);
        hn_bot_r = hK_bot(bnd_ind_r>0,1);
        hm_bot_r = hK_bot(bnd_ind_r(bnd_ind_r>0),1);
        ln_r = l_k(bnd_ind_r>0,1);
        lm_r = l_k(bnd_ind_r(bnd_ind_r>0),1);
        Kn_r  = K_k(bnd_ind_r>0,1);
        Km_r  = K_k(bnd_ind_r(bnd_ind_r>0),1);
        w_nm_r = w_k_r(bnd_ind_r>0,1);
        [C_nm_r] = conductance_karst(w_nm_r,Kn_r,Km_r,ln_r,lm_r,thicknessn_r,thicknessm_r);
        C_r(bnd_ind_r>0,1) = C_nm_r;
        dh_r = hm_r - hn_r;
        dh_r(hm_bot_r>hn_r,1) = hm_r(hm_bot_r>hn_r,1) - hm_bot_r(hm_bot_r>hn_r,1);
        dh_r(hn_bot_r>hm_r,1) = hn_bot_r(hn_bot_r>hm_r,1) - hn_r(hn_bot_r>hm_r,1);
        dh_right(bnd_ind_r>0,1) = dh_r;
        dQdt_right = C_r .* dh_right;
        
        % karst matrix flow to spring river (represented by fixed head rivers)
        dQdt_sprriv = zeros(nK,1);
        hn_sprriv = h_k(bnd_ind_sprriv(bnd_ind_sprriv>0),1);
        hnbot_sprriv = hK_bot(bnd_ind_sprriv(bnd_ind_sprriv>0),1);
        Rbot_sprriv = max(R_sprriv,hnbot_sprriv);
        Rn_sprriv = max(Rbot_sprriv,hn_sprriv); 
        dh_sprriv = Rbot_sprriv - Rn_sprriv;
        flux_sprriv = C_sprriv .* dh_sprriv; % C_sprriv: conductance of river bed sediment   
        dQdt_sprriv(bnd_ind_sprriv(bnd_ind_sprriv>0),1) = flux_sprriv;
        % spring river matrix flow flux change
        dQdt_sprriv_flux = -flux_sprriv; % m3/s   
        
        % connection with river
        % not two rivers may connect to the same grid
        dQdt_river = zeros(nK,1);
        hn_riv = h_k(bnd_ind_riv(bnd_ind_riv>0),1);
        Rn_riv = max(R_rb,hn_riv); % R_rb is the bottom elevation of the river
        dh_river = (h_riv + R_rb - Rn_riv);
        flux_exch = C_rb .* dh_river; % C_rb: conductance of river bed sediment
        flux_exch_mat = repmat(flux_exch,1,n_rivseg_uni);
        flux_exch_gridK = sum(flux_exch_mat.*uni_bnd_ind_riv_mat)';      
        dQdt_river(uni_bnd_ind_riv,1) = flux_exch_gridK;
        % river flux change
        dQdt_riv_flux = -flux_exch; % m3/s   
        
        % sum up all changes
        dQdt_gridK = (dQdt_left + dQdt_top + dQdt_right + dQdt_bottom + dQdt_sprriv + dQdt_river);
        dQdt_all = [dQdt_gridK;dQdt_sprriv_flux;dQdt_riv_flux];
        dhdt = dQdt_all./S*86400; % m/d       
    end

end
