clear all

%% 1 Study domain
% 1.1 computational grids
% shapefile of the study domain, karst
file_karst = '.\shapefiles\karst.shp';
file_domain = '.\shapefiles\study_domain.shp';
file_surface = '.\shapefiles\surface_catchment.shp';
file_elev = '.\dem\elev_bound.tif';

shape_karst = shaperead(file_karst);
shape_domain = shaperead(file_domain);
shape_surface = shaperead(file_surface);

% elevation tiff
[elev_tiff,R_elev] = geotiffread(file_elev);
elev_tiff = double(elev_tiff);
elev_tiff(elev_tiff==-32768)=nan; % set no data points to nan, resolution: 1/1200 degree
lon_elev_lim = R_elev.LongitudeLimits;
lat_elev_lim = R_elev.LatitudeLimits;
res_elev_tiff = R_elev.CellExtentInLatitude; % resolution of tiff grid [degree]
lon_elev_tiff = (lon_elev_lim(1)+res_elev_tiff/2:res_elev_tiff:lon_elev_lim(2))';
lat_elev_tiff = (lat_elev_lim(2)-res_elev_tiff/2:-res_elev_tiff:lat_elev_lim(1))';

% generate computational grids
% grid_karst_list: ind,i,j,l,t,b,r
res_spatial = 3; %[minute]
[grid_domain_value,grid_surface,grid_karst_list,grid_karst_aquifer,elev_domain,lon_domain,lat_domain,area_grid_lat] = domain2grid(shape_karst,shape_domain,shape_surface,elev_tiff,lon_elev_tiff,lat_elev_tiff,res_spatial);

elev_karst = elev_domain(grid_karst_list(:,1));
A_gridK = area_grid_lat(grid_karst_list(:,2)); % area of karst m2
nK = size(grid_karst_list,1); % number of karst grids

gridK_lat = lat_domain(grid_karst_list(:,2),1);
gridK_lon = lon_domain(grid_karst_list(:,3),1);

% index of grid domain
n_row_domain = size(grid_domain_value,1);
n_col_domain = size(grid_domain_value,2);
ind_domain = (1:n_row_domain*n_col_domain)';
grid_oneCol = reshape(grid_domain_value,[n_row_domain*n_col_domain,1]);

%save grid information
% save('grid_studydomain.mat','grid_surface','lon_domain','lat_domain')

% grid karst location 
gridK_ijk_latlon = [grid_karst_list(:,2),grid_karst_list(:,3),grid_karst_list(:,1),gridK_lat,gridK_lon];
% save('gridK_location.mat','gridK_ijk_latlon');

% in domain grid lat lon
indomain = grid_domain_value>0;
londomain_mat = repmat(lon_domain',n_row_domain,1);
latdomain_mat = repmat(lat_domain,1,n_col_domain);
londomain_grid = londomain_mat(indomain);
latdomain_grid = latdomain_mat(indomain);

% save('latlondomain_grid.mat','latdomain_grid','londomain_grid')
fprintf('### finish domain configuration %s ###\n', datestr(now))

% 1.2 river networks
file_river_list = {'danube_tri1','neckar_tri1','neckar_tri2','neckar_tri3','upper_danube','upper_neckar'};
n_riv_network = length(file_river_list);

rivnet_folder = '.\rivernetworks\';
gridriv = grid_surface;
gridriv(gridriv == 0) = nan;
gridriv(gridriv == 1) = 0;

for ii=1:n_riv_network
    
    riv_name_num = file_river_list{ii};
    riv_dir = strcat(rivnet_folder,'gridriv_',file_river_list{ii},'.txt');
    rivnet_data = readmatrix(riv_dir,'NumHeaderLines',1,'Delimiter','\t');
    
    [chn_length,chn_w,elev_rb,chn_s,riv_i,riv_j,riv_k,A_rb] = compute_riv_network_info(rivnet_data,n_row_domain);    

    % length [m]
    riv_network_info.(riv_name_num).chn_length = chn_length;
    % width [m]
    riv_network_info.(riv_name_num).chn_w = chn_w;
    % river bed elevation [m.a.s.l]
    riv_network_info.(riv_name_num).elev_rb = elev_rb;
    % slope [-]
    riv_network_info.(riv_name_num).chn_s = chn_s;
    % river index i [-]
    riv_network_info.(riv_name_num).riv_i = riv_i;
    % river index j [-]
    riv_network_info.(riv_name_num).riv_j = riv_j;
    % river index k [-]
    riv_network_info.(riv_name_num).riv_k = riv_k;
    % area contacted with grid [m2]
    riv_network_info.(riv_name_num).A_rb = A_rb;   
    % reach of the river
    riv_network_info.(riv_name_num).n_gridRiv = size(chn_length,1);
    
    % roughness 
    riv_network_info.(riv_name_num).chn_n = 0.01 * ones(size(chn_length));
    % bed sediment depth [m]
    riv_network_info.(riv_name_num).depth_rb = 0.5 * ones(size(chn_length));
    % hydraulic conductivity of bed sediment [m/s]
    riv_network_info.(riv_name_num).K_rb = 10^-10 * ones(size(chn_length));
    
    gridriv(riv_k) = 1;
end

% save('gridriv.mat','gridriv');
fprintf('### finish river network configuration %s ###\n', datestr(now))

% 1.2 river - grid connections
% non-karst grid
grid_nK_ind = find(grid_domain_value == 1 & grid_surface==1);
grid_nK_i = grid_nK_ind-floor(grid_nK_ind/n_row_domain)*n_row_domain;
grid_nK_j = ceil(grid_nK_ind/n_row_domain);
grid_nK_k = (grid_nK_j-1)*n_row_domain+grid_nK_i;
grid_nK_lat = lat_domain(grid_nK_i);
grid_nK_lon = lon_domain(grid_nK_j);
A_grid_nK = area_grid_lat(grid_nK_i); % m2
n_nonK = length(grid_nK_ind);
grid_nK_ijk = [grid_nK_i,grid_nK_j,grid_nK_ind];

% non-karst grid connets to rivers
load('grid2river_surface.mat')
grid_nK_riv_i = grid2river_surface_row(grid_nK_k);
grid_nK_riv_j = grid2river_surface_col(grid_nK_k);
grid_nK_riv_k = (grid_nK_riv_j-1)*n_row_domain+grid_nK_riv_i;

% karst grid connets to rivers
grid_K_riv_i = grid2river_surface_row(grid_karst_list(:,1));
grid_K_riv_j = grid2river_surface_col(grid_karst_list(:,1));
grid_K_riv_k = (grid_K_riv_j-1)*n_row_domain+grid_K_riv_i;

fprintf('### finish non-/karst grid configuration %s ###\n', datestr(now))

% spring river info: grids with karst discharge river segments
sprriv_file = load('grid_sprriv_sprrivinfo.mat','grid_sprriv_info');
grid_sprriv_info = sprriv_file.grid_sprriv_info;

grid_sprriv_ij = grid_sprriv_info(:,1:2);
grid_sprriv_width = grid_sprriv_info(:,3);
grid_sprriv_length = grid_sprriv_info(:,4)*1000; % km to meter
grid_sprriv_elev = grid_sprriv_info(:,5);

n_sprriv = size(grid_sprriv_ij,1);

% karst grid connected to which river
grid_sprriv_k = (grid_sprriv_ij(:,2)-1)*n_row_domain+grid_sprriv_ij(:,1);
grid_sprriv_riv_i = grid2river_surface_row(grid_sprriv_k);
grid_sprriv_riv_j = grid2river_surface_col(grid_sprriv_k);
grid_sprriv_riv_k = grid2river_surface(grid_sprriv_k);

A_sprrivK = grid_sprriv_length.*grid_sprriv_width;
R_sprriv = grid_sprriv_elev;

[~,bnd_ind_sprriv] = ismember(grid_sprriv_k,grid_karst_list(:,1),'rows');
sprriv_aquifer = grid_karst_aquifer(bnd_ind_sprriv);

sprriv_info.bnd_ind_sprriv = bnd_ind_sprriv;
sprriv_info.n_sprriv = n_sprriv; % number of cells with spring river
sprriv_info.A_sprrivK = A_sprrivK;
sprriv_info.R_sprriv = R_sprriv; % elevation of the spring river bottom
sprriv_info.sprriv_aquifer = sprriv_aquifer;
sprriv_info.grid_karst_aquifer = grid_karst_aquifer;

unique_aquiferK = unique(grid_karst_aquifer);
nK_aquiferK = length(unique_aquiferK);
sprriv_info.unique_aquiferK = unique_aquiferK;
sprriv_info.nK_aquiferK = nK_aquiferK;

fprintf('### finish karst discharge grid configuration %s ###\n', datestr(now))

% all rivers to one
gridK_ind = grid_karst_list(:,1);
[river_info] = river_info_compute(riv_network_info,grid_domain_value,gridK_ind);
[channel_info_one] = channel_info_compute(riv_network_info);
n_rivK = 0;
A_rbK =[];
bnd_gridK_ind_riv =[];
rivgrid_k_ind = [];
bnd_ind_riv = [];
C_rb = [];
R_rb = [];
riv_info_fields = fieldnames(river_info);

for iriv=1:n_riv_network
    n_rivK_i = river_info.(riv_info_fields{iriv,1}).n_rivK;
    A_rbK_i = river_info.(riv_info_fields{iriv,1}).A_rbK; %m2
    bnd_gridK_ind_riv_i = river_info.(riv_info_fields{iriv,1}).bnd_gridK_ind_riv;
    rivgrid_k_ind_i = river_info.(riv_info_fields{iriv,1}).rivgrid_k_ind;
    bnd_ind_riv_i = river_info.(riv_info_fields{iriv,1}).bnd_ind_riv;
    C_rb_i = river_info.(riv_info_fields{iriv,1}).C_rb;
    R_rb_i = river_info.(riv_info_fields{iriv,1}).R_rb;

    n_rivK = n_rivK+n_rivK_i;
    A_rbK = [A_rbK;A_rbK_i];
    bnd_gridK_ind_riv = [bnd_gridK_ind_riv;bnd_gridK_ind_riv_i];
    rivgrid_k_ind = [rivgrid_k_ind;rivgrid_k_ind_i];
    bnd_ind_riv = [bnd_ind_riv;bnd_ind_riv_i];
    C_rb = [C_rb;C_rb_i];
    R_rb = [R_rb;R_rb_i];
end

n_gridRiv = length(rivgrid_k_ind);

% unique bnd_ind_riv, some grids have main and tri in one grid, sum other
% up to this grid.
uni_bnd_ind_riv = unique(bnd_ind_riv);
n_rivseg_uni = length(uni_bnd_ind_riv);
n_rivseg = length(bnd_ind_riv);
uni_bnd_ind_riv_mat = zeros(n_rivseg,n_rivseg_uni);
for jj=1:length(uni_bnd_ind_riv)
   uni_bnd_ind_riv_mat(:,jj) = (bnd_ind_riv==uni_bnd_ind_riv(jj));
end

river_info_one.bnd_ind_riv = bnd_ind_riv;
river_info_one.uni_bnd_ind_riv = uni_bnd_ind_riv;
river_info_one.uni_bnd_ind_riv_mat = uni_bnd_ind_riv_mat;
river_info_one.n_rivseg = n_rivseg;
river_info_one.n_rivseg_uni = n_rivseg_uni;
river_info_one.C_rb = C_rb;
river_info_one.R_rb = R_rb;
river_info_one.A_rbK = A_rbK;

fprintf('### finish river connection configuration %s ###\n', datestr(now))

%% 2 Nonkarst discharge and karst recharge
% simulation time period
t0 = datetime(1951,1,1);
te = datetime(1955,12,31);
tdate = (t0:te)';
ndays = length(tdate);

% flow from non-karst grid cells
Q_nK_sim = load('Q_nK.mat');
Q_nK = Q_nK_sim.Q_nK;

% karst recharge and overland flow
rech_karst = load ('rech_karst.mat');
Rech_condiut_karst = rech_karst.Rech_condiut_karst;
Rech_matrix_karst = rech_karst.Rech_matrix_karst;
Runoff_all_karst = rech_karst.Runoff_all_karst;
fprintf('### finish reading nonkarst discharge and karst recharge %s ###\n', datestr(now))

%% 3 Model run
% parameter
load gridK_attr
m_karst = 200; % karst aquifer thickness meter
depth_sprriv = 0.5;

% conduit recession coefficient
load gridKc.mat

% spring river hydraulic conductivity
[ia_sprriv,ib_sprriv] = ismember(grid_sprriv_k,gridK_ind);

% a systematic multiplier for a single ensemble member
para_cor = [1.0	1.0	1.0	1.0]; 

k_C_gridK = 1./gridKc;
K_sprriv = K_k_gridK(ib_sprriv);

k_C_gridK = k_C_gridK*para_cor(1);
K_sprriv = K_sprriv*para_cor(2);  
K_k_gridK_cor = K_k_gridK*para_cor(3);
ne_gridK_cor = ne_gridK*para_cor(4);

k_C_gridK(k_C_gridK>1) = 1;

[gridK_info] = gridK_info_compute(grid_karst_list,ne_gridK_cor,A_gridK,elev_karst,m_karst,k_C_gridK,K_k_gridK_cor);
C_sprriv = K_sprriv.*A_sprrivK/depth_sprriv;

S_sprriv = 1*A_sprrivK;
S_riv = 1*A_rbK;
S_gridK = ne_gridK_cor .* A_gridK;
S = [S_gridK;S_sprriv;S_riv];

fprintf('### GWKarst model running ...... %s ###\n', datestr(now))

[result] = karstGroundwater(ndays,n_gridRiv,n_sprriv, ...
    nK,n_rivK,gridK_info,sprriv_info,channel_info_one,river_info, ...
    river_info_one,riv_network_info,S,rivgrid_k_ind,grid_nK_riv_k, ...
    grid_K_riv_k,bnd_gridK_ind_riv,grid_sprriv_riv_k,Rech_condiut_karst, ...
    Rech_matrix_karst,Runoff_all_karst,Q_nK,C_sprriv);

result.t = tdate;

fprintf('### finish GWKarst model simulation %s ###\n', datestr(now))



