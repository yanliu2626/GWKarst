function [row_rivconnect,col_rivconnect] = grid_flow_river(n_gridrow,n_gridcol,grid_row,grid_col,grid_flowdir,gridriv)
%calculate grids connected to which river cell
row_rivconnect = 0;
col_rivconnect = 0;
rowfrom = grid_row;
colfrom = grid_col;
indriv_gridfrom = gridriv(rowfrom,colfrom);

if indriv_gridfrom == 1
    row_rivconnect = rowfrom;
    col_rivconnect = colfrom;
    
else
    cell_fldr = grid_flowdir(rowfrom,colfrom);    
    while row_rivconnect==0 && col_rivconnect==0
        [rowto,colto] = togrid(rowfrom,colfrom,cell_fldr);
        
        if rowto<1 || colto<1 || rowto>n_gridrow || colto>n_gridcol
            %connect to out of boundary
            row_rivconnect = -999;
            col_rivconnect = -999;
            break;
        elseif gridriv(rowto,colto) == 1
            %connect to river grids
            row_rivconnect = rowto;
            col_rivconnect = colto;
            break;
        elseif isnan(gridriv(rowto,colto))
            %connect to out of surface domain, possible outlets
            row_rivconnect = -99;
            col_rivconnect = -99;
            break;
        else
            %connect to other grids (not river grids) within domain
            rowfrom = rowto;
            colfrom = colto;
            cell_fldr = grid_flowdir(rowfrom,colfrom);
            
            if isnan(cell_fldr)
                %connect to out of surface domain, possible outlets
                row_rivconnect = -99;
                col_rivconnect = -99;
            break;
            end
        end      
    end
    
end

end

