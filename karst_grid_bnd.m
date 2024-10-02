function [karst_grid,karst_grid_value] = karst_grid_bnd(bnd_shp_karst,lon_shp_karst,lat_shp_karst,grid_domain,lat_bound,lon_bound,bnd_shp_domain,area_grid_lat)
%generate the grid of karst, which is within the study domain. Meanwhile,
%determine the boundary grids due to the disconnectivity of karst aquifer.

% Input variables
% bnd_shp_karst: the extent (lon; lat) of karst aquifer
% lon_shp_karst: the lon of karst aquifer boundary
% lat_shp_karst: the lat of karst aquifer boundary
% grid_domain: the grid that indicates the grid cells falling within the
% study domain
% lat_bound: the lat of each grid bound
% lon_bound: the lon of each grid bound
% bnd_shp_domain: the extent (lon; lat) of study domain grid

% Output variables
% karst_grid.ind: 1 there are study grid overlap with the karst aquifer;
% 0: no overlap
% karst_grid.grid: the grids that are in the study domain and are karst
% karst_grid_value: same size as grid_domain, but only show the karst grid
% with 2, the others are 0

%% check the overlap of karst aquifer with study domain
% find the overlap regions between grid_domain and karst aquifer
x_min = max(bnd_shp_domain(1,1),bnd_shp_karst(1,1));
x_max = min(bnd_shp_domain(2,1),bnd_shp_karst(2,1));
y_min = max(bnd_shp_domain(1,2),bnd_shp_karst(1,2));
y_max = min(bnd_shp_domain(2,2),bnd_shp_karst(2,2));
%%
if (x_min<x_max) && (y_min<y_max)
    % find index for grid that will be checked for karst
    ind_lon_min = find(lon_bound>x_min,1,'first')-1;
    ind_lon_max = find(lon_bound>=x_max,1,'first')-1;
    ind_lat_min = find(lat_bound<y_max,1,'first')-1;
    ind_lat_max = find(lat_bound<=y_min,1,'first')-1;
    
    poly_karst = polyshape(lon_shp_karst,lat_shp_karst);
    karst_grid_value = zeros(size(grid_domain));
    karst_ratio = zeros(size(grid_domain));
    
    for j=ind_lon_min:ind_lon_max
        for i=ind_lat_min:ind_lat_max        
            x_poly = [lon_bound(j),lon_bound(j),lon_bound(j+1),lon_bound(j+1)];
            y_poly = [lat_bound(i),lat_bound(i+1),lat_bound(i+1),lat_bound(i)];
            poly1 = polyshape(x_poly,y_poly);
            polyout = intersect(poly1,poly_karst);

            if isempty(polyout.Vertices)
                karst_grid_value(i,j) = 0;
                karst_ratio(i,j)= 0;
            elseif isequal(polyout.Vertices,poly1.Vertices)
                karst_grid_value(i,j) = 2;
                karst_ratio(i,j) = 1;
            else
                area_intersect = sum(areaint(polyout.Vertices(:,2),polyout.Vertices(:,1),referenceEllipsoid('wgs84')),'all','omitnan');
                ratio_intersect = area_intersect/area_grid_lat(i);
                karst_ratio(i,j) = ratio_intersect;
                if ratio_intersect>0.5
                    karst_grid_value(i,j) = 2;
                else
                    karst_grid_value(i,j) = 0;
                end
            end           
        end
    end 
    
    karst_grid_value(grid_domain==0) = 0;
    
    if sum(karst_grid_value==2,'all')==0
        karst_grid.ind = 0;
    else
        karst_grid.ind = 1;      
        [grid_karst] = grid_config(karst_grid_value);
        karst_grid.grid = grid_karst;   
    end
else
    karst_grid.ind = 0;
    karst_grid_value = nan;
end

end

