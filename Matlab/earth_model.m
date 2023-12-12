%% 3- Earth Model
function [Jp] = earth_model(sigma,T_p,r_ear,r_sat_ear)
    epslon_p = 1;       % Earth emmisivity 
    Jp = (r_ear/(r_sat_ear))^2 * (epslon_p * sigma * T_p^4) ;
end