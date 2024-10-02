function [chn_length,chn_w,elev_rb,chn_s,riv_i,riv_j,riv_k,A_rb] = compute_riv_network_info(rivnet_data,n_row_domain)
%compute river network info, especially careful for slope
% length [m]
chn_length = rivnet_data(:,4)*1000;
% width [m]
chn_w = rivnet_data(:,3);
% area contacted with grid [m2]
A_rb = chn_length .* chn_w;
% river index i [-]
riv_i = rivnet_data(:,1);
% river index j [-]
riv_j = rivnet_data(:,2);
% river index k [-]
riv_k = (riv_j-1)*n_row_domain+riv_i;

% river elevation before pre-processing
elev0 = rivnet_data(:,5);
n_reach = length(elev0);

% river bed elevation [m.a.s.l]
elev_rb = elev0;
% correct elev to form descending elev

for i = 2:n_reach
       
    elev_diff = elev_rb(i-1)-elev_rb(i);

    if elev_diff<=0
        indlow = find(elev_rb<elev_rb(i-1),1,'first');
        
        if ~isempty(indlow)
            elev_dd = elev_rb(i-1) - elev_rb(indlow);
            elev_dd_len = sum(chn_length(i-1:indlow))-0.5*(chn_length(i-1)+chn_length(indlow));
            
        else
            elev_dd = elev_rb(i-2) - elev_rb(i-1);
            elev_dd_len = 0.5*(chn_length(i-1)+chn_length(i-2));            
        end
        
        elev_grad = elev_dd/elev_dd_len;
        elev_rb(i) = elev_rb(i-1) - elev_grad*0.5*sum(chn_length(i-1:i)); 
    end
end
    


% slope [-]
slo = -diff(elev_rb) ./ ((chn_length(1:end-1)+chn_length(2:end))/2);
chn_s = [slo(1);slo];
end

