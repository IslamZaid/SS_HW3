%% 3- Earth Model
function [Jp] = earth_model(sigma,epslon_p,T_p,r_ear,r_sat_ear)
    Jp = (r_ear/(r_sat_ear))^2 * (epslon_p * sigma * T_p^4) ;
end