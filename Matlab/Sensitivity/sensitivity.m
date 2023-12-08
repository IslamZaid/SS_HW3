clc;
close all;
clear all;
 

b_ear_sun = 0;
b_sat_ear = pi;



%% Sphere Model Development

% Number of nodes
numNodes_list = [10, 20, 40, 80 ];
for nn = 1 : length(numNodes_list)
numNodes = numNodes_list(nn);
r_ear = 6.371e+6;  %(m);
[X,Y,Z] = sphere(numNodes);
X = r_ear * X;
Y = r_ear * Y;
Z = r_ear * Z;

% XYZ = 
centerNodes = zeros([size(X), 4]); % First column for labels, fifth column for lateral surface areas

Area = zeros(numNodes);

for i = 1 : numNodes
    for j = 1 : numNodes

        node1 = [X(i, j), Y(i, j), Z(i, j)];
        node2 = [X(i, j+1), Y(i, j+1), Z(i, j+1)];
        node3 = [X(i+1, j), Y(i+1, j), Z(i+1, j)];
        node4 = [X(i+1, j+1), Y(i+1, j+1), Z(i+1, j+1)];

        V1 = node2 - node1;
        V2 = node3 - node1;
        V3 = node4 - node1;

        Area(i,j) = 0.5 * ( norm(cross(V1,V2)) + norm(cross(V2,V3)) );

        % Save center node information
        centerNodes(i, j, :) = [  (X(i, j) + X(i, j+1) + X(i+1, j) + X(i+1, j+1)) / 4;
                                  (Y(i, j) + Y(i, j+1) + Y(i+1, j) + Y(i+1, j+1)) / 4;
                                  (Z(i, j) + Z(i, j+1) + Z(i+1, j) + Z(i+1, j+1)) / 4;
                                  Area(i,j)];
        % Add label to center node
        text(centerNodes(i, j, 1), ...
             centerNodes(i, j, 2), ...
             centerNodes(i, j, 3), ...
             sprintf('%d, %d', i, j), ...
             'FontSize', 8, 'Color', 'r');
    end
end

h = surf(X,Y,Z);
Normal_Vectors = Area.*h.FaceNormals;



 %% Radiation models
 % position vectors definition

% r_ear_sun = 1.496e+11;
r_sat_ear = 6.371e+6 + 2e+6 ; 
r_ear = 6.371e+6;
area_sat = 100;   %(m^2)
R_ear_sun = @(r,beta_e_sun) r * [ cos(beta_e_sun) , sin(beta_e_sun) ,0 ];          % position vector of the earth respect to the sun
R_sat_ear = @(beta_sat_ear) r_sat_ear * [ cos(beta_sat_ear) , sin(beta_sat_ear) ,0 ] ;        % position vector of the sat respect to the earth
R_sat_sun = @(r,beta_e_sun, beta_sat_e) R_ear_sun(r,beta_e_sun) + R_sat_ear(beta_sat_e);    % position vector of the sat respect to the sun

unit_vector_normalto_sat = @(beta_sat_ear) -1 * R_sat_ear(beta_sat_ear) / norm(R_sat_ear(beta_sat_ear));
A_sat = @(beta_sat_ear)  area_sat * unit_vector_normalto_sat(beta_sat_ear);







 % Total radiation
  
% Define time parameters
totaltime = 3.154e7; % Total simulation time (1 year in s)
% numsteps = 50;
% timestep = totaltime / numsteps;

timestep = 7620 /20;
numsteps=round(totaltime/timestep);
numsteps = 1;
%Time Implementation
% loop through time steps
n_earth = 2*pi/totaltime;  % rad/s
e = 0.0167;
G = 6.6371e-11; % Gravitational Constant in m^3/kg*s^2
mass_earth = 5.97e24; % in kg
mu = G*mass_earth; % Standard gravitational parameter in m^3/s^2
a = 1.49e11;   % m
h = sqrt( mu * a * (1-e^2) ) ;   % kg.m^2/s

alpha = 0.12; 
eps = 0.5; 
ae = alpha/eps;
sigma = 5.67*10^-8; % m^2/K.w

I_mat = zeros(numsteps,numNodes+1,numNodes+1);
Q_mat = zeros(numsteps,numNodes+1,numNodes);
Mat = zeros(numsteps,3);

% parfor step = 1:   10000%0.1*numsteps+1
step = 1;
t = (step - 1)*timestep;
M = n_earth*(t); 
theta  = M+(2*e-0.25*e^3)*sin(M)+(5/4)*e^2*sin(2*M)+(13/12)*e^3*sin(3*M);
r = (h^2/mu)/(1+e*cos(theta));   % earth sun dist

% 1- Sun Model
Q_sun = sun_model(r,r_ear,R_ear_sun(r,b_ear_sun),R_sat_sun(r,b_ear_sun, b_sat_ear), A_sat(b_ear_sun));
%2- Planetary Albedo Model
[Q_A,I,Q ]  = albedo_model(Normal_Vectors,centerNodes,A_sat(b_sat_ear), R_ear_sun(r,b_ear_sun),R_sat_ear(b_sat_ear),area_sat);
% 3- Earth Model
Q_p = earth_model();
Q_gen = 10000; 
Temp(nn) =( ((Q_gen/eps)+( Q_p + ae*(Q_sun + Q_A)))/(2*area_sat*sigma) )^0.25;

%     Mat(step,1) = Q_p ;
%     Mat(step,2) = Q_sun ;
%     Mat(step,3) = Q_A ;

theta_earth_mat(step) = theta;
theta_sat_mat(step) = b_sat_ear;

I_mat(step,:,:) = I; 
Q_mat(step,:,:) = Q; 
 
% end
Q_tot(nn) = Q_sun  + Q_A + Q_p ;

nn
Q_tot(nn)
 

end

%%
for i = 1 : length(Q_tot)-1

    Rel_Err(i) = 100* abs((Q_tot(i+1)-Q_tot(i))/Q_tot(i));

end


f = figure 

plot(numNodes_list , Q_tot, '-o','LineWidth',2)
hold on

xlabel("Number of elements in long and lat")
ylabel("Absorbed radiation (W)")
yyaxis right

plot(numNodes_list(2:end) , Rel_Err, '-s','LineWidth',2)

grid on 
grid minor 

% Increase the thickness of the ticks
ax = gca; % Get the current axis handle
ax.FontSize = 16; % Set the font size of the ticks
% ax.LineWidth = 2; % Set the thickness of the ticks
xticks(numNodes_list)
% You can also set other properties if needed, e.g., ax
ylabel("Relative Error (%)")
saveas(f,['fig_Q_sensitivity_',num2str(rad2deg(b_ear_sun)),'_',num2str(rad2deg(b_sat_ear))],"png")

%%

for i = 1 : length(Q_tot)-1

    Rel_Err_temp(i) = 100* abs((Temp(i+1)-Temp(i))/Temp(i));

end
f = figure 

plot(numNodes_list , Temp, '-o','LineWidth',2)
hold on

xlabel("Number of elements in long and lat")
ylabel("Satallite temperature (K)")
yyaxis right

plot(numNodes_list(2:end) , Rel_Err_temp, '-s','LineWidth',2)

grid on 
grid minor 

% Increase the thickness of the ticks
ax = gca; % Get the current axis handle
ax.FontSize = 16; % Set the font size of the ticks
% ax.LineWidth = 2; % Set the thickness of the ticks
xticks(numNodes_list)
% You can also set other properties if needed, e.g., ax
ylabel("Relative Error (%)")
saveas(f,['fig_temp_sensitivity_',num2str(rad2deg(b_ear_sun)),'_',num2str(rad2deg(b_sat_ear))],"png")
 
%% 
f= figure('Position',[100 100 1000 600]);
 
 
r_eartosun = 1.496e11;
f_sun = 1e+3;

pe = R_ear_sun(r_eartosun, b_ear_sun)/f_sun/2.5;
psat = R_sat_ear(b_sat_ear);

subplot(1,2,[1,2])
quiver3(0,0,0,pe(1), pe(2), pe(3)) 
hold on
quiver3(pe(1),pe(2),pe(3),psat(1), psat(2), psat(3))
% sun 
r_sun = 7e+9/f_sun;
[x,y,z] = sphere;
x = r_sun*x; y = r_sun*y; z = r_sun*z; 
surf(x,y,z, 10000*ones(size(z)),"EdgeColor","red","LineStyle","-.")
% earth
x = pe(1)+X; y = pe(2)+Y; z = pe(3)+Z; 
% surf(x,y,z, squeeze(I_mat(k,:,:)), 'EdgeAlpha',0.2)
Q_mat(:,:,end+1) = Q_mat(:,:,end);
surf(x,y,z, squeeze(Q_mat(1,:,:)), 'EdgeAlpha',0.2)
%satallite
pss = pe + psat ;
f_sat = 1e+6;
n_sat = 10;
[x,y,z] = ellipsoid(pss(1), pss(2), pss(3), 1*f_sat,1*f_sat,1*f_sat, n_sat);
% x = f_sat * x; y = f_sat* y; z = f_sat*z;
sat = surf(x,y,z, zeros(n_sat) );
% rotate(sat, [0,0,1],rad2deg(b_sat_ear) )
axis equal
view(-20, 65)
title("Sun radiation on the earth")
hold off
% color bar; 
% caxis([0 round(max(max(Q_mat)))]);

xlim([-1e7 8e7])
ylim([-1e7 1e7])
saveas(f,['fig_',num2str(rad2deg(b_ear_sun)),'_',num2str(rad2deg(b_sat_ear))],"png")












