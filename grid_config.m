function [grid_karst] = grid_config(grid_value)
%assign the attributes to every grid: 0 out of study domain, 1 non-karst, 
%2 karst
% input parameters
% res_spatial: spatial resolution [arc degree]
% coord_ll: the coordinates of the lower left corner [arc degree]
% grid_value: 0,1,2 to indicate out of domain, non-karst, and karst [-]

% output variables
% grid_karst: a matrix consists of all karst grid and indicates the connected faces for lateral flow
% grid_karst: % ind,i,j,l,t,b,r
%%
% grid size
n_row = size(grid_value,1);
n_col = size(grid_value,2);

% generate karst grid networks for lateral flow computation
n_karst = sum(grid_value==2,'all');
grid_karst = nan(n_karst,7); % ind,i,j,l,t,b,r

k = 1;
% karst grid with index k
for j=1:n_col
    for i=1:n_row
        if grid_value(i,j)==2
            grid_karst(k,1) = (j-1)*n_row+i;
            grid_karst(k,2) = i;
            grid_karst(k,3) = j;
            k = k+1;
        end            
    end
end

%% define the connected grid
for k=1:n_karst
    karst_i = grid_karst(k,2);
    karst_j = grid_karst(k,3);
    
    % left
    i_l = karst_i;
    j_l = karst_j-1;
    
    [Lia_l,Locb_l] = ismember([i_l,j_l],grid_karst(:,2:3),'row');     
    
    if Lia_l
        grid_karst(k,4) = grid_karst(Locb_l,1);
    else
        grid_karst(k,4) = 0;
    end
    
    % top
    i_t = karst_i-1;
    j_t = karst_j;
    
    [Lia_t,Locb_t] = ismember([i_t,j_t],grid_karst(:,2:3),'row');     
    
    if Lia_t
        grid_karst(k,5) = grid_karst(Locb_t,1);
    else
        grid_karst(k,5) = 0;
    end
    
    % bottom
    i_b = karst_i+1;
    j_b = karst_j;
    
    [Lia_b,Locb_b] = ismember([i_b,j_b],grid_karst(:,2:3),'row');     
    
    if Lia_b
        grid_karst(k,6) = grid_karst(Locb_b,1);
    else
        grid_karst(k,6) = 0;
    end
    
    % right
    i_r = karst_i;
    j_r = karst_j+1;
    
    [Lia_r,Locb_r] = ismember([i_r,j_r],grid_karst(:,2:3),'row');     
    
    if Lia_r
        grid_karst(k,7) = grid_karst(Locb_r,1);
    else
        grid_karst(k,7) = 0;
    end   
           
end

end

