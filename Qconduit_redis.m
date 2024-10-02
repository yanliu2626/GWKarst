function [Qcond] = Qconduit_redis(Qcon,grid_flowdir,Qdown_fsep,bnd_grid_ind)
%compute the conduit flow to rivers
%   grid_flowdir: 1=upstream, 2=downstream, 
%   id_grids: the number of karst grids in order
%note: if all the conduit flow is totally distributed to river cells
%(representing virtual karst spring features), conduit flow contribute to
%rivers. If there are sinks of conduit flow in the cells, add that to
%conduit storage

% grids in the most upstream
% if is 0, no flow at the karst boundary
upstr_faces = (grid_flowdir == 2); % grids with upstream faces
downstr_faces = (grid_flowdir == 1); % grids with downstream faces
upstr_check = sum(upstr_faces,2,'omitnan');
downstr_check = sum(downstr_faces,2,'omitnan');

grids_down = (upstr_check==0 & downstr_check>0); % grids only connects to downstream, no upstream connection
n_upstr = sum(grids_down); % number of grids which are upstream cells
n_grid = size(grid_flowdir,1); % number of grids left after removing upstream grids

num_grids = (1:n_grid)';

id_grids_down = bnd_grid_ind(grids_down==1,:);
num_grids_down = num_grids(grids_down==1,:);

%% start while loop, if there are upstream grids and the total number of
% grids is larger than 1
while n_upstr>0 && n_grid>1
    
    %% add Q of upstream to downstream grids
    for i=1:n_upstr
        num_grid_i = num_grids_down(i);
        Q_fsep_i = Qdown_fsep(num_grid_i,:);
        Qcon_i = Qcon(num_grid_i);
        Q_sep_i = (Q_fsep_i*Qcon_i)';
        
        % set the upstream conduit flow cells to zero after separating to
        % downstream cells
        Qcon(num_grid_i) = 0;
        
        % now separate Q to downstream and to Q in that grid
        id_down_i = id_grids_down(i,:);
        Qcon(id_down_i(id_down_i>0),1) = Qcon(id_down_i(id_down_i>0),1)+Q_sep_i(id_down_i>0); 
        
        % set upstream grids to zeros to exclude these grids in the calculations
        grid_flowdir_i = grid_flowdir(num_grid_i,:);
        ind_down = (grid_flowdir_i==1);
        grid_flowdir_i(ind_down) = 0;
        grid_flowdir(num_grid_i,:) = grid_flowdir_i;
        
        % set the corresponding downstream cells' upstream grids to zeros
        if ind_down(1) % left
            grid_flowdir(id_down_i(1),4) = 0;
        end
        
        if ind_down(2) % top
            grid_flowdir(id_down_i(2),3) = 0;
        end
        
        if ind_down(3) % bottom
            grid_flowdir(id_down_i(3),2) = 0;
        end
        
        if ind_down(4) % right
            grid_flowdir(id_down_i(4),1) = 0;
        end     
    end

    %% update all information
    % grids in the most upstream
    upstr_faces = (grid_flowdir == 2); % grids with upstream
    downstr_faces = (grid_flowdir == 1); % grids with downstream
    upstr_check = sum(upstr_faces,2,'omitnan');
    downstr_check = sum(downstr_faces,2,'omitnan');

    grids_down = (upstr_check==0 & downstr_check>0); % grids only connects to downstream, no upstream connection
    n_upstr = sum(grids_down); % number of grids which are upstream connections
    n_grid = size(grid_flowdir,1); % number of grids left after removing upstream grids

    id_grids_down = bnd_grid_ind(grids_down==1,:);
    num_grids_down = num_grids(grids_down==1,:);

    %% if only grids with river inside left, break the while loop and return
    % the conduit flow
    if n_upstr == 0
        % cells with river
%         Qconduit = Qcon(bnd_ind_sprriv); % return the conduit flow in grids with rivers inside
        
        % cells with conduit sink, which are not connected to rivers
        % ## update conduit storage of this cell
        break; % break the while loop, since all concuit flow has been redistributed to grids with rivers
    end   

end
Qcond = Qcon;
end


