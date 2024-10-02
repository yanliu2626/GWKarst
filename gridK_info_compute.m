function [gridK_info] = gridK_info_compute(grid_karst_list,ne_gridK,A_gridK,elev_karst,m_karst,k_C,K_k_gridK)
%obtain the karst grid information

gridK_ind = grid_karst_list(:,1);
bnd_l = grid_karst_list(:,4);
bnd_t = grid_karst_list(:,5);
bnd_b = grid_karst_list(:,6);
bnd_r = grid_karst_list(:,7);

[~,ib_l] = ismember(bnd_l,gridK_ind);
[~,ib_t] = ismember(bnd_t,gridK_ind);
[~,ib_b] = ismember(bnd_b,gridK_ind);
[~,ib_r] = ismember(bnd_r,gridK_ind);

bnd_grid_ind = [ib_l,ib_t,ib_b,ib_r];
hK_bot = elev_karst - m_karst;
hK_top = elev_karst;

gridK_info.ne_K = ne_gridK;
gridK_info.bnd_grid_ind = bnd_grid_ind;
gridK_info.gridK_ind = gridK_ind;
gridK_info.A_gridK = A_gridK;
gridK_info.hK_bot = hK_bot;
gridK_info.hK_top = hK_top;
gridK_info.k_C = k_C;
gridK_info.K_k = K_k_gridK;
end

