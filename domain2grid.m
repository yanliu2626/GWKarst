function [grid_domain_value,grid_surface,grid_karst_list,grid_karst_aquifer,elev_domain,lon_domain,lat_domain,area_grid_lat] = domain2grid(shape_karst,shape_domain,shape_surface,elev_tiff,lon_elev_tiff,lat_elev_tiff,res_spatial)
%generate the grid of the study domain and indicate the karst grid

%input variables
% shape_karst, the shapefile the karst map [polygon]
% bound_domain, the boundary of the study domain (outline) [polygon]
% elev_tiff, tiff file of the elevation [tiff]
% res_spatial, the spatial resolution of the grid, [minutes]

%output variables
% grid_domain, the study grid, 0 out of domain, 1 non-karst, 2 karst
% elev_domain, the mean elevation of the grid
% coord_ll, the coordinates of the lower left grid

%% bounds of the domain
bound_domain = shape_domain.BoundingBox;
minute_domain = (bound_domain-floor(bound_domain))*60;
ll_minute = floor(minute_domain(1,:)/res_spatial)*res_spatial;
rt_minute = ceil(minute_domain(2,:)/res_spatial)*res_spatial;
xy_ll = floor(bound_domain(1,:)) + ll_minute/60;
xy_rt = floor(bound_domain(2,:)) + rt_minute/60;
xy_ll_minute = floor(bound_domain(1,:))*60 + ll_minute; % left low
xy_rt_minute = floor(bound_domain(2,:))*60 + rt_minute; % right top
coord_ll = xy_ll;
bnd_shp_domain = [xy_ll;xy_rt];

% number of rows and columns of the study domain
n_row = (xy_rt_minute(2) - xy_ll_minute(2))/res_spatial;
n_col = (xy_rt_minute(1) - xy_ll_minute(1))/res_spatial;

% the bound lon, lat of every grid
lon_bound = (xy_ll_minute(1):res_spatial:xy_rt_minute(1))'/60;
lat_bound = (xy_rt_minute(2):-res_spatial:xy_ll_minute(2))'/60;

lon_domain = (lon_bound(1:end-1)+lon_bound(2:end))/2;
lat_domain = (lat_bound(1:end-1)+lat_bound(2:end))/2;

% area of grid in different lat
area_grid_lat = nan(n_row,1);
coord_lon = [floor(coord_ll(1)),floor(coord_ll(1)),floor(coord_ll(1))+res_spatial/60,floor(coord_ll(1))+res_spatial/60];
for i=1:n_row
    coord_lat = [lat_bound(i+1),lat_bound(i),lat_bound(i),lat_bound(i+1)];   
    area_grid_lat(i,1) = areaint(coord_lat,coord_lon,referenceEllipsoid('wgs84'));
end

%% determine if grid is in the study domain (subsurface)
grid_domain = nan(n_row,n_col);
poly_domain = polyshape(shape_domain.X,shape_domain.Y);
ratio_domain = nan(n_row,n_col);
for i=1:n_row
    for j=1:n_col
        x_poly = [lon_bound(j),lon_bound(j),lon_bound(j+1),lon_bound(j+1)];
        y_poly = [lat_bound(i),lat_bound(i+1),lat_bound(i+1),lat_bound(i)];
        poly1 = polyshape(x_poly,y_poly);
        polyout = intersect(poly1,poly_domain);
        
        if isempty(polyout.Vertices)
            grid_domain(i,j) = 0;
            ratio_domain(i,j)= 0;
        elseif isequal(polyout.Vertices,poly1.Vertices)
            grid_domain(i,j) = 1;
            ratio_domain(i,j) = 1;
        else
            area_intersect = sum(areaint(polyout.Vertices(:,2),polyout.Vertices(:,1),referenceEllipsoid('wgs84')),'all','omitnan');
            ratio_intersect = area_intersect/area_grid_lat(i);
            ratio_domain(i,j) = ratio_intersect;
            if ratio_intersect>0.5
                grid_domain(i,j) = 1;
            else
                grid_domain(i,j) = 0;
            end
        end                
    end
end

%% determine if grid is in the surface domain
grid_surface = nan(n_row,n_col);
poly_surface = polyshape(shape_surface.X,shape_surface.Y);
ratio_surface = nan(n_row,n_col);
for i=1:n_row
    for j=1:n_col
        x_poly = [lon_bound(j),lon_bound(j),lon_bound(j+1),lon_bound(j+1)];
        y_poly = [lat_bound(i),lat_bound(i+1),lat_bound(i+1),lat_bound(i)];
        poly1 = polyshape(x_poly,y_poly);
        polyout = intersect(poly1,poly_surface);
        
        if isempty(polyout.Vertices)
            grid_surface(i,j) = 0;
            ratio_surface(i,j)= 0;
        elseif isequal(polyout.Vertices,poly1.Vertices)
            grid_surface(i,j) = 1;
            ratio_surface(i,j) = 1;
        else
            area_intersect = sum(areaint(polyout.Vertices(:,2),polyout.Vertices(:,1),referenceEllipsoid('wgs84')),'all','omitnan');
            ratio_intersect = area_intersect/area_grid_lat(i);
            ratio_surface(i,j) = ratio_intersect;
            if ratio_intersect>0.5
                grid_surface(i,j) = 1;
            else
                grid_surface(i,j) = 0;
            end
        end                
    end
end

%% calculate the average elevation of the grid
elev_domain = nan(size(grid_domain));

for i=1:n_row
    grid_lat_up = lat_bound(i);
    grid_lat_low = lat_bound(i+1);
    
    ind_lat_elev_sta = find(lat_elev_tiff<=grid_lat_up,1,'first');
    ind_lat_elev_end = find(lat_elev_tiff<grid_lat_low,1,'first')-1;
    for j=1:n_col
        grid_lon_low = lon_bound(j);
        grid_lon_up = lon_bound(j+1);
        
        ind_lon_elev_sta = find(lon_elev_tiff>=grid_lon_low,1,'first');
        ind_lon_elev_end = find(lon_elev_tiff>grid_lon_up,1,'first')-1; 
        
        if grid_domain(i,j)>0
            elev_domain(i,j) = mean(elev_tiff(ind_lat_elev_sta:ind_lat_elev_end,ind_lon_elev_sta:ind_lon_elev_end),'all','omitnan');
        end       
    end
end

%% determine karst grid
% IMPORTANT: careful with the connection for grids overlap with two karst
% aquifers
n_shape_karst = length(shape_karst);
grid_domain_value = grid_domain;
grid_karst_list = [];
grid_karst_aquifer = [];
ind_k = 0;
for i=1:n_shape_karst
    bnd_shp_karst = shape_karst(i).BoundingBox;
    lon_shp_karst = shape_karst(i).X;
    lat_shp_karst = shape_karst(i).Y;
    [karst_grid,karst_grid_value] = karst_grid_bnd(bnd_shp_karst,lon_shp_karst,lat_shp_karst,grid_domain,lat_bound,lon_bound,bnd_shp_domain,area_grid_lat);
    grid_domain_value(karst_grid_value==2) = 2;
    if karst_grid.ind ==1
        ind_k = ind_k+1;
        grid_karst_list = [grid_karst_list;karst_grid.grid];
        karst_aquifer = ones(length(karst_grid.grid),1)*ind_k;
        grid_karst_aquifer = [grid_karst_aquifer;karst_aquifer];
    end
end

end

