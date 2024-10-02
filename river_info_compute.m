function [river_info] = river_info_compute(riv_network_info,grid_domain_value,gridK_ind)
%obtain the river information for karst grid
%   Detailed explanation goes here

fields = fieldnames(riv_network_info);
n_riv = length(fields);

for iriv=1:n_riv

    riv_i = riv_network_info.(fields{iriv,1}).riv_i;
    riv_j = riv_network_info.(fields{iriv,1}).riv_j;
    A_rb = riv_network_info.(fields{iriv,1}).A_rb;
    depth_rb = riv_network_info.(fields{iriv,1}).depth_rb;
    elev_rb = riv_network_info.(fields{iriv,1}).elev_rb;
    K_rb = riv_network_info.(fields{iriv,1}).K_rb;

    imax_domain = size(grid_domain_value,1);
    n_reach = length(riv_i);
    bnd_gridK_ind_riv = nan(n_reach,1);
    R_rb = elev_rb;
    C_rb = K_rb.*A_rb./depth_rb; % m2/s

    for ii=1:n_reach
        if grid_domain_value(riv_i(ii),riv_j(ii)) == 2
            bnd_gridK_ind_riv(ii) = (riv_j(ii)-1)*imax_domain+riv_i(ii);
        end
    end

    nan_ind = isnan(bnd_gridK_ind_riv);
    rivgrid_k_ind = ~nan_ind;
    bnd_gridK_ind_riv = bnd_gridK_ind_riv(~nan_ind,1);
    [~,ib] = ismember(bnd_gridK_ind_riv,gridK_ind);
    bnd_ind_riv = ib;
    R_rb = R_rb(~nan_ind,1);
    C_rb = C_rb(~nan_ind,1);
    A_rbK = A_rb(~nan_ind,1);
    n_rivK = length(bnd_gridK_ind_riv);

    river_info.(fields{iriv,1}).bnd_ind_riv = bnd_ind_riv;
    river_info.(fields{iriv,1}).R_rb = R_rb;
    river_info.(fields{iriv,1}).C_rb = C_rb;
    river_info.(fields{iriv,1}).n_rivK = n_rivK;
    river_info.(fields{iriv,1}).A_rbK = A_rbK;
    river_info.(fields{iriv,1}).bnd_gridK_ind_riv = bnd_gridK_ind_riv;
    river_info.(fields{iriv,1}).rivgrid_k_ind = rivgrid_k_ind;
end
end

