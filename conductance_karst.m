function [C_nm] = conductance_karst(w_nm,Kn,Km,ln,lm,hn,hm)
%calculate the conductance for grid n with the neighbor m

% ## needed to check later
Tn = Kn.*hn;
Tm = Km.*hm;
C_nm = w_nm.*Tn.*Tm./(Tn.*lm+Tm.*ln);
end

