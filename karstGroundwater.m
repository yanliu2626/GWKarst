% function calculating head, groundwater flow, and river routing

function [result] = karstGroundwater(ndays,n_gridRiv, ...
    n_sprriv,nK,n_rivK,gridK_info,sprriv_info,channel_info_one,river_info, ...
    river_info_one,riv_network_info,S,rivgrid_k_ind,grid_nK_riv_k,grid_K_riv_k, ...
    bnd_gridK_ind_riv,grid_sprriv_riv_k,Rech_condiut_karst,Rech_matrix_karst, ...
    Runoff_all_karst,Q_nK,C_sprriv)

sprriv_info.C_sprriv = C_sprriv; % conductance of spring river

% initial condition
h_riv_ini = 0.1*ones(n_rivK,1); %[m above the river bottom]
hK_ini = gridK_info.hK_top-1; 
hK_conduit_ini = 0.001 * ones(size(hK_ini)); % m of conduit storage

% time-series of spring and rivers
Qriv_Tseries = zeros(ndays,n_gridRiv);
Qex_Tseries = zeros(ndays,n_rivK);
Qkarst_Tseries = zeros(ndays,n_sprriv);
Qmatrix_Tseries = zeros(ndays,n_sprriv);
Qconduit_Tseries = zeros(ndays,n_sprriv);
QsurfK_Tseries = zeros(ndays,nK);
Qrunoff_Tseries = zeros(ndays,nK);
Qexcess_Tseries = zeros(ndays,nK);
hK_Tseries = zeros(ndays,nK);
hK_Tseries(1,:) = hK_ini';
Q_ini = 2*ones(n_gridRiv,1); %m3/s

% run river in sequence
riv_run_name = {'danube_tri1','neckar_tri1','neckar_tri2','neckar_tri3','upper_danube','upper_neckar'};
n_rivgrid_cum = zeros(1,length(riv_run_name)); % river in all grids
n_rivgridK_cum = zeros(1,length(riv_run_name)); % river that in karst grids
for jj=1:length(riv_run_name)
    n_rivgrid_i = channel_info_one.(riv_run_name{1,jj}).n_reach;
    n_rivgrid_cum(jj) = n_rivgrid_i;
    n_rivgridK_cum(jj) = river_info.(riv_run_name{1,jj}).n_rivK;
end
n_rivgrid_cum = cumsum(n_rivgrid_cum);
n_rivgrid_sta_end = [0,n_rivgrid_cum]; % using 0 since start ind will be sta+1

n_rivgridK_cum = cumsum(n_rivgridK_cum);
n_rivgridK_sta_end = [0,n_rivgridK_cum];

% which grid that tributaries join the main river 
ind_tri = [37,46,51,52,nan,nan]; % the index of tributary to the main river
Q_tri = zeros(1,4);

A_gridK = gridK_info.A_gridK;

% loop calculating time-series Q
for i=2:ndays
    % run the karst grids (considering inter-catchment groundwater flow)
    % discharge for all the karst springs (Qspr) and the exchange between river
    % and karst aquifer (Qex)

    rech_cond_t =  Rech_condiut_karst(i,:)';
    rech_matr_t = Rech_matrix_karst(i,:)';
    runoff_t = Runoff_all_karst(i,:)';
    QrunoffK = runoff_t.*A_gridK/86400;
    
    [Qex,Qkarst,Qmatrix,Qconduit,QexcessK,hK,hK_conduit] = gwK_head(hK_ini, ...
        hK_conduit_ini,rech_matr_t,rech_cond_t,h_riv_ini,gridK_info, ...
        sprriv_info,river_info_one,S,nK);
    
    hK_ini = hK;
    hK_conduit_ini = hK_conduit;
    hK_Tseries(i,:) = hK';
    Qex_Tseries(i,:) = Qex';
    Qkarst_Tseries(i,:) = Qkarst';
    Qmatrix_Tseries(i,:) = Qmatrix';
    Qconduit_Tseries(i,:) = Qconduit';
    
    QrunoffKtot = QrunoffK+QexcessK;
    QsurfK_Tseries(i,:) = QrunoffKtot';
    Qrunoff_Tseries(i,:) = QrunoffK';
    Qexcess_Tseries(i,:) = QexcessK';
    
    % channel routing for the major rivers
    % first simulate tributaries and then simulate main rivers
    QnK = Q_nK(i,:)';
    n_riv_network = length(riv_run_name);

    for iriv=1:n_riv_network
        % river index
        riv_run_i = riv_run_name{1,iriv};
        channel_info = channel_info_one.(riv_run_i);
        riv_k_i = riv_network_info.(riv_run_i).riv_k;
        rivgrid_k_ind_i = rivgrid_k_ind(n_rivgrid_sta_end(iriv)+1:n_rivgrid_sta_end(iriv+1),1);
        Q_ini_i = Q_ini(n_rivgrid_sta_end(iriv)+1:n_rivgrid_sta_end(iriv+1),1);
        Q_tri_mat = zeros(n_rivgrid_sta_end(iriv+1)-n_rivgrid_sta_end(iriv),1);
        
        if iriv==5
            % tributaries for upper danube
            Q_tri_mat(ind_tri(1)) = Q_tri(1);
        elseif iriv==6
            % tributaries for upper neckar
            Q_tri_mat(ind_tri(2:4)) = Q_tri(2:4);
        end
        
        % river routing
        [Qriv_i,h_update_i] = riv_routing(channel_info,Qex,QnK,Qkarst, ...
            QrunoffKtot,riv_k_i,grid_nK_riv_k,grid_K_riv_k, ...
            bnd_gridK_ind_riv,grid_sprriv_riv_k,Q_ini_i,Q_tri_mat);
        
        if n_rivgridK_sta_end(iriv+1)>0
            h_riv_ini(n_rivgridK_sta_end(iriv)+1:n_rivgridK_sta_end(iriv+1),1) = h_update_i(rivgrid_k_ind_i>0);
        end
        
        Q_ini(n_rivgrid_sta_end(iriv)+1:n_rivgrid_sta_end(iriv+1),1) = Qriv_i;
        Qriv_Tseries(i,n_rivgrid_sta_end(iriv)+1:n_rivgrid_sta_end(iriv+1)) = Qriv_i';
        
        if iriv<=4
            % 4 tributaries for uppder danube and upper neckar
            Q_tri(1,iriv) = Qriv_i(end);
        end

    end
end

result.Q = Qriv_Tseries; % simulated river discharge
result.head = hK_Tseries; % head timeseries of karst matrix
end





