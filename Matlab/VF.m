%% 2- Planetary Albedo Model
 function [F   ] = VF(Normal_Vectors,centerNodes3,A_sat,  R_sat_ear,area_sat)
      % find the elements that see the sun and the radiant power emitted by sun per unit projected area
 
    n_nodes = size(centerNodes3,1);
 
    
    % calcualte the albedo
    F = zeros(n_nodes,1);  
    for i = 1 : n_nodes-1
        for j = 1 : n_nodes-1  % in y
            r_sat_elm = R_sat_ear -  squeeze(centerNodes3(i,j,1:3))';
            A = squeeze(Normal_Vectors(i,j,:))';
            B = r_sat_elm;
            cos_1 = dot(A,B)/(norm(A)*norm(B));

            A = A_sat;
            B = r_sat_elm ;
            cos_2 = dot(A,B)/(norm(A)*norm(B));
            area_elm = centerNodes3(i,j,4);
            d = norm(r_sat_elm);
            F(i,j) =    cos_1 * cos_2 *area_sat*area_elm/ d^2;
        end
    end
 
end