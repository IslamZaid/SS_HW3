%% 1- Sun Model
function Q_s = sun_model(r_sun_earth, r_earth,R_ear_sun,R_sat_sun, A_sat )
%     R_sat_sun  : sun to satallite vector 
%     A_sat : Satallite area vector
    Ps = 3.856*10^26;    % (w) 
    J_s = Ps / (4*pi*norm(R_sat_sun)^2);    % (w/m^2)  
    A_sat_proj = abs(dot(A_sat, R_sat_sun)/norm(R_sat_sun));
    Q_s = A_sat_proj *J_s;
    
    ang = acos((dot(R_ear_sun,R_sat_sun))/(norm(R_ear_sun)*norm(R_ear_sun)));
    ang_min = atan(r_earth/r_sun_earth);

    if ang < ang_min 
        Q_s = 0;
    end
end