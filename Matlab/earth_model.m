%% 3- Earth Model
function [Jp] = earth_model()
    epslon_p = 1;
    sigma = 5.67e-8;
    T_p = 254;    %(K)
    r_ear = 6.371e+6;
    r_sat_ear = 6.371e+6 + 2e+6;
    Jp = (r_ear/(r_sat_ear))^2 * (epslon_p * sigma * T_p^4) ;
end