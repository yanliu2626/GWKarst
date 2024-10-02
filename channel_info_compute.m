function [channel_info] = channel_info_compute(riv_network_info)
%generate channel setting information for the kinematic wave calculations

fields = fieldnames(riv_network_info);
n_riv = length(fields);

for iriv=1:n_riv
    chn_w = riv_network_info.(fields{iriv,1}).chn_w;
    chn_s = riv_network_info.(fields{iriv,1}).chn_s;
    chn_n = riv_network_info.(fields{iriv,1}).chn_n;
    chn_length = riv_network_info.(fields{iriv,1}).chn_length;
    mat_alpha = (chn_n.*chn_w.^(2/3)./chn_s.^0.5).^0.6;
    mat_x = chn_length; % m
    channel_info.(fields{iriv,1}).mat_alpha = mat_alpha;
    channel_info.(fields{iriv,1}).mat_x = mat_x;
    channel_info.(fields{iriv,1}).chn_w = chn_w;
    channel_info.(fields{iriv,1}).n_reach = length(mat_x);
end
end

