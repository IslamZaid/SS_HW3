%% 2- Planetary Albedo Model
function [Q_A,I,Q ] = albedo_model(Normal_Vectors,centerNodes3,A_sat, R_ear_sun,R_sat_ear,area_sat,Ps, a)
    % find the elements that see the sun and the radiant power emitted by sun per unit projected area
     n_nodes = size(centerNodes3,1);
    I = zeros(n_nodes);
    for i = 1 : n_nodes-1
        for j = 1 : n_nodes-1   % in y
            A = squeeze(Normal_Vectors(i,j,:))';
            B = R_ear_sun;
            dot_prod = dot(A,B);
            if dot_prod < 0
                r_elm_sun = R_ear_sun - squeeze(centerNodes3(i,j,1:3))';
                d = norm(r_elm_sun);
                I(i,j) = Ps / (4*pi*d^2);
            else
                I(i,j) = 0;
            end
        end
    end
    
    % calcualte the albedo
    Q = zeros(n_nodes,1);  
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
            Q(i,j) = a * I(i,j) * cos_1 * cos_2 *area_sat*area_elm/ d^2;
        end
    end
    Q(Q>0) = 0;
    Q = -1 * Q;
    Q_A = sum(sum(Q)); 
end