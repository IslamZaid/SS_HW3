clc;
close all;
clear all;

%% Constants
e = 0.0167;
G = 6.6371e-11; % Gravitational Constant in m^3/kg*s^2
mass_earth = 5.97e24; % in kg
mu = G*mass_earth; % Standard gravitational parameter in m^3/s^2
a_SMA = 1.49e11;   % semi-major axis (m)
h = sqrt( mu * a_SMA * (1-e^2) ) ;   % kg.m^2/s


% Number of nodes
numNodes = 40;    % number of eleemtns in longitue and lattitde


% time incremintation 
totaltime = 3.154e7; % Total simulation time (1 year in s)
time_steps = 0 :60 : totaltime ;
numsteps=length(time_steps);


% Sun data 
Ps = 3.856*10^26;   % Total emmited power from the sun (w) 


% Earth data
sigma = 5.67e-8;    % boltz man constant
T_p = 254;          % Earth temperature (K)    
r_ear = 6.371e+6;   % Earth radius (m)
a = 0.39;           % fraction of solar radiation reflected from earth. the average varies from 0.31 to 0.39 


% Satellite data
r_sat_ear = 6.371e+6 + 2e+6;   % Satellite distance to the center of the earth (m)
area_sat = 100;     % Satelite Area (m^2)
alpha = 0.26;       % absorptance
eps = 0.83;          % Emmitance 
ae = alpha/eps;     % (alpha/eps)


% Please you the path in accordance with alpha and eps values (((IMPORTANT)))
path_img = 'images/WhitePaint_silicone' ;         % alpha = 0.26 , esp = 0.83
% path_img = 'images/WhitePaint_silicone_1000h' ;  % alpha = 0.29 , esp = 0.83
% path_img = 'images/WhitePaint_silicate' ;        % alpha = 0.12 , esp = 0.90       
% path_img = 'images/WhitePaint_silicate_1000h' ;  % alpha = 0.14 , esp = 0.90         
% path_img = 'images/Aluminized_kapton' ;          % alpha = 0.40 , esp = 0.63         
if ~exist(path_img, 'dir')           
   mkdir(path_img)
end
%% Sphere Model Development
% Number of nodes
st = 0;
[X,Y,Z] = sphere(numNodes);     % discretize the domain into sherical domain 
X = r_ear * X;
Y = r_ear * Y;
Z = r_ear * Z;
    
centerNodes = zeros([size(X), 4]); % First column for labels, forth column for lateral surface areas

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

        Area(i,j) = 0.5 * ( norm(cross(V1,V2)) + norm(cross(V2,V3)) );   % calculate the area of each element by dividing into two triangles

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

srfc = surf(X,Y,Z);
Normal_Vectors = Area.*srfc.FaceNormals;


%% Radiation models
%% position vectors definition

R_ear_sun = @(r,beta_ear_sun) r * [ cos(beta_ear_sun) , sin(beta_ear_sun) ,0 ];          % position vector of the earth wrt the sun
R_sat_ear = @(beta_sat_ear) r_sat_ear * [ cos(beta_sat_ear) , sin(beta_sat_ear) ,0 ] ;        % position vector of the sat wrt the earth
R_sat_sun = @(r,beta_ear_sun, beta_sat_ear) R_ear_sun(r,beta_ear_sun) + R_sat_ear(beta_sat_ear);    % position vector of the sat wrt the sun

unit_vector_normalto_sat = @(beta_sat_ear) -1 * R_sat_ear(beta_sat_ear) / norm(R_sat_ear(beta_sat_ear));
A_sat = @(beta_sat_ear)  area_sat * unit_vector_normalto_sat(beta_sat_ear);    % vector of the satellite area it's always pointing to the earth center

 

%% Total radiation

% calculations
n_earth = 2*pi/totaltime;                           %earth rotations around the sun (rad/s)
period_sat = 2*pi*sqrt(r_sat_ear^3/mu);              % satellite period around the earth (s)

% declare matrices where the data will be recorded
I_mat = zeros(numsteps,numNodes+1,numNodes+1);      % incident intensity of the soalr radiation on the earth (w/m^2)
Q_mat = zeros(numsteps,numNodes+1,numNodes);        % earth surface elements contribution to the satellite albedo radiation
Q_p_mat = zeros(numsteps,1);                        % planetary radiation to satellite(w)
Q_A_mat = zeros(numsteps,1);                        % Total Albedo radiation to satellite(w)
Q_sun_mat = zeros(numsteps,1);                      % Sun direct radiation to satellite (w)
r_EarToSun = zeros(numsteps,1);                     % earth to sun distance (m)
Temp = zeros(numsteps,1);                           % Satellite temperature (K)
theta_earth_mat = zeros(numsteps,1);                % angular position of the earth around the sum  (rad)
theta_sat_mat = zeros(numsteps,1);                  % angular position of the satellite around the earth (rad)

% loop through time steps
parfor step = 1: numsteps
    t = time_steps(step);
    M = n_earth*(t); % mean anomaly
    b_ear_sun  = M+(2*e-0.25*e^3)*sin(M)+(5/4)*e^2*sin(2*M)+(13/12)*e^3*sin(3*M); % the angular posion of the earth wrt sun (rad)
    r = (h^2/mu)/(1+e*cos(b_ear_sun));              % earth to sun distance (m)
    r_EarToSun(step) = r;                           % list of earth to sun distance (m)
    b_sat_ear = 2*pi*t / period_sat;                 % the angular posion of the satellte wrt earth (rad)

    % 1- Sun Model
    Q_sun = sun_model(r,r_ear,R_ear_sun(r,b_ear_sun),R_sat_sun(r,b_ear_sun, b_sat_ear), A_sat(b_ear_sun),Ps);
    
    %2- Planetary Albedo Model
    [Q_A,I,Q ]  = albedo_model(Normal_Vectors,centerNodes,A_sat(b_sat_ear), R_ear_sun(r,b_ear_sun),R_sat_ear(b_sat_ear),area_sat,Ps, a);
    
    % 3- Earth Model
    Q_p = earth_model(sigma,epslon_p,T_p,r_ear,r_sat_ear);
    Q_gen = 10000;                                  % generated heat due to the electronic devices inside the satellite 
    Temp(step) =( ((Q_gen/eps)+( Q_p + ae*(Q_sun + Q_A)))/(2*area_sat*sigma) )^0.25;  % temperature of the satellite (K)
    
    % record the data
    Q_p_mat(step) = Q_p ;
    Q_sun_mat(step) = Q_sun ;
    Q_A_mat(step) = Q_A ;

    theta_earth_mat(step) = b_ear_sun;
    theta_sat_mat(step) = b_sat_ear;

    I_mat(step,:,:) = I; 
    Q_mat(step,:,:) = Q; 

    if rem(step,100) == 0
        step
    end
end
 


%% plotting
cond = 1;
if cond == 1 

if size(Q_mat,2) ~= size(Q_mat,3)
    Q_mat(:,:,end+1) = Q_mat(:,:,end);
end
k = 80; %
f=figure();  %'Position', [0, 0, 1500, 1200],
r_eartosun = a_SMA;
f_sun = 1e+3;
b_ear_sun = theta_earth_mat(k);
b_sat_ear = theta_sat_mat(k);
pe = R_ear_sun(r_eartosun, b_ear_sun)/f_sun/2.5;
psat = R_sat_ear(b_sat_ear); 

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
surf(x,y,z, squeeze(Q_mat(k,:,:)), 'EdgeAlpha',0.2)
%satellite
pss = pe + psat ;
f_sat = 1e+6;
n_sat = 10;
[x,y,z] = ellipsoid(pss(1), pss(2), pss(3), 1*f_sat,1*f_sat,1*f_sat, n_sat);
% x = f_sat * x; y = f_sat* y; z = f_sat*z;
sat = surf(x,y,z, zeros(n_sat) );
% rotate(sat, [0,0,1],rad2deg(b_sat_ear) )

axis equal
view(-20, 45)
title("Earth surface elements contribution to Albedo")
hold off
c=colorbar;
c.LineWidth =2 ;
c.FontSize = 16 ;
clim([0 max(max(max(Q_mat)))]);

xlim([-1e7 8e7])
ylim([-1e7 2e7])
zlim([-2e7 1e7])
axes('position',[.3 .05 .5 .4])
box on
% earth
x = pe(1)+X; y = pe(2)+Y; z = pe(3)+Z;
% surf(x,y,z, squeeze(I_mat(k,:,:)), 'EdgeAlpha',0.2)
surf(x,y,z, squeeze(Q_mat(k,:,:)), 'EdgeAlpha',0.2)
view(-20, 45)
axis equal

%f.FontSize = 16;
exportgraphics(f,[path_img,'/albedo_elements.png'], Resolution=1200);

%%


if size(Q_mat,2) ~= size(Q_mat,3)
    Q_mat(:,:,end+1) = Q_mat(:,:,end);
end
k = 80 ; %
f=figure();  %'Position', [0, 0, 1500, 1200],
r_eartosun = a_SMA;
f_sun = 1e+3;
b_ear_sun = theta_earth_mat(k);
b_sat_ear = theta_sat_mat(k);
pe = R_ear_sun(r_eartosun, b_ear_sun)/f_sun/2.5;
psat = R_sat_ear(b_sat_ear); 

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
surf(x,y,z, squeeze(I_mat(k,:,:)), 'EdgeAlpha',0.2)
% surf(x,y,z, squeeze(Q_mat(k,:,:)), 'EdgeAlpha',0.2)
%satellite
pss = pe + psat ;
f_sat = 1e+6;
n_sat = 10;
[x,y,z] = ellipsoid(pss(1), pss(2), pss(3), 1*f_sat,1*f_sat,1*f_sat, n_sat);
% x = f_sat * x; y = f_sat* y; z = f_sat*z;
sat = surf(x,y,z, zeros(n_sat) );
% rotate(sat, [0,0,1],rad2deg(b_sat_ear) )

axis equal
view(-20, 45)
title("Incident radiation on the earth surface")
hold off
c=colorbar;
c.LineWidth =2 ;
c.FontSize = 16 ;
clim([0 max(max(max(I_mat)))]);

xlim([-1e7 8e7])
ylim([-1e7 2e7])
zlim([-2e7 1e7])
axes('position',[.3 .05 .5 .4])
box on
% earth
x = pe(1)+X; y = pe(2)+Y; z = pe(3)+Z;
surf(x,y,z, squeeze(I_mat(k,:,:)), 'EdgeAlpha',0.2)
% surf(x,y,z, squeeze(Q_mat(k,:,:)), 'EdgeAlpha',0.2)
view(-20, 45)
axis equal

%f.FontSize = 16;
exportgraphics(f,[path_img,'/inc_sun_rad.png'], Resolution=1200);

%%
 
f=figure('Position', [0, 0, 1800, 600]);  %'Position', [0, 0, 1500, 1200],
subplot(1,3,[1,2])
plot(Temp(1:60*24), 'LineWidth',2)
grid on 
grid minor 
title("Temperature Variation")
xlabel("Time (hours)")
ylabel("Temperature (K)")
xticks(0:6*60:24*60)
xticklabels(0:6:24) 
% ylim([170 300])
 subplot(1,3,3) 
plot(Temp(1 : period_sat/60), 'LineWidth',2)
title("Temperature - one cycle")
xlabel("Time (hours)")
xticks(0:0.5*60:2*60)
xticklabels(0:0.5:2)
grid on 
grid minor 
% ylim([170 300])

fontsize(gcf,16,"points")


exportgraphics(f,[path_img,'/temp_time.png'], Resolution=1200);
%%
 
f=figure('Position', [0, 0, 1800, 600]);  %'Position', [0, 0, 1500, 1200],
subplot(1,3,[1,2])
plot(Q_A_mat(1:60*24), 'LineWidth',2)
hold on
plot(Q_sun_mat(1:60*24), 'LineWidth',2)
plot(Q_p_mat(1:60*24), 'LineWidth',2)
grid on 
grid minor 
title("Radiation modes")
legend(["Albedo","Sun","Planetary"])
xlabel("Time (hours)")
ylabel("Radiation power (w)")
xticks(0:6*60:24*60)
xticklabels(0:6:24)  

subplot(1,3,3) 
plot(Q_A_mat( 1:period_sat/60), 'LineWidth',2)
hold on
plot(Q_sun_mat(1:period_sat/60), 'LineWidth',2)
plot(Q_p_mat(1:period_sat/60), 'LineWidth',2)
title("Radiation - one cycle")
xlabel("Time (hours)")
xticks(0:0.5*60:60*2 )
xticklabels(0:0.5:2)
grid on 
grid minor 
fontsize(gcf,16,"points")

exportgraphics(f,[path_img,'/QQQ_time.png'], Resolution=1200);


%%
 Q_tot = Q_A_mat+Q_sun_mat+Q_p_mat;
f=figure('Position', [0, 0, 1800, 600]);  %'Position', [0, 0, 1500, 1200],
subplot(1,3,[1,2])
plot(Q_tot(1:60*24), 'LineWidth',2)
grid on 
grid minor 
title("Total incident radiation")
xlabel("Time (hours)")
ylabel("Radiation power (w)")
xticks(0:6*60:24*60)
xticklabels(0:6:24) 
 subplot(1,3,3) 
plot(Q_tot(1:period_sat/60  ), 'LineWidth',2)
title("Total radiation - one cycle")
xlabel("Time (hours)")
xticks(0:0.5*60:2*60)
xticklabels(0:0.5:2)
grid on 
grid minor 
fontsize(gcf,16,"points")


exportgraphics(f,[path_img,'/Q_time.png'], Resolution=1200);









%% Initialize gif
ccond = 0;
if ccond == 1
if size(Q_mat,2) ~= size(Q_mat,3)
    Q_mat(:,:,end+1) = Q_mat(:,:,end);
end
timeloop = 1:128 ;

for k = timeloop
    f=figure('visible','off');  %'Position', [0, 0, 1500, 1200],
    r_eartosun = 1.496e11;
    f_sun = 1e+3;
    b_ear_sun = theta_earth_mat(k);
    b_sat_ear = theta_sat_mat(k);
    pe = R_ear_sun(r_eartosun, b_ear_sun)/f_sun/2.5;
    psat = R_sat_ear(b_sat_ear); 

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
    surf(x,y,z, squeeze(Q_mat(k,:,:)), 'EdgeAlpha',0.2)
    %satellite
    pss = pe + psat ;
    f_sat = 1e+6;
    n_sat = 10;
    [x,y,z] = ellipsoid(pss(1), pss(2), pss(3), 1*f_sat,1*f_sat,1*f_sat, n_sat);
    % x = f_sat * x; y = f_sat* y; z = f_sat*z;
    sat = surf(x,y,z, zeros(n_sat) );
    % rotate(sat, [0,0,1],rad2deg(b_sat_ear) )
 
    axis equal
    view(-20, 45)
    title("Albedo radiation reflected to satallite")
    hold off
    c=colorbar;
    c.LineWidth =2 ;
    c.FontSize = 16 ;
    clim([0 max(max(max(Q_mat)))]);
 
    xlim([-1e7 8e7])
    ylim([-1e7 2e7])
    zlim([-2e7 1e7])

%     f.FontSize = 16;
    exportgraphics(f,[path_img,'/Animate_sat.gif'], 'Resolution', 1200,  'BackgroundColor', 'none','Append',true);

    k
end
 
%%

for i = timeloop
    f = figure('Visible','off');
    plot(time_steps(1:i)/3600,Temp(1:i),"LineWidth",2)
    hold on 
    scatter(time_steps(i)/3600,Temp(i) )
    hold off
    i
    xlim([0 time_steps(timeloop(end))/3600])
    ylim([min(min(Temp))   max(max(Temp))])
    xlabel("Time (hours)")
    ylabel("Distance between sun and earth (m)")
    
    exportgraphics(f,[path_img,'/temp.gif'],'Append',true);
end





%% Initialize gif  (Year)

if size(Q_mat,2) ~= size(Q_mat,3)
    Q_mat(:,:,end+1) = Q_mat(:,:,end);
end
timeloop = 1:1000 : numsteps;

for k = timeloop
    f=figure('visible','off');  %'Position', [0, 0, 1500, 1200],
    r_eartosun = r_EarToSun(k);
    f_sun = 1e+3;
    b_ear_sun = theta_earth_mat(k);
    b_sat_ear = theta_sat_mat(k);
    pe = R_ear_sun(r_eartosun, b_ear_sun)/f_sun/2.5;
    psat = R_sat_ear(b_sat_ear); 

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
    surf(x,y,z, squeeze(Q_mat(k,:,:)), 'EdgeAlpha',0.2)
    %satellite
    pss = pe + psat ;
    f_sat = 1e+6;
    n_sat = 10;
    [x,y,z] = ellipsoid(pss(1), pss(2), pss(3), 1*f_sat,1*f_sat,1*f_sat, n_sat);
    % x = f_sat * x; y = f_sat* y; z = f_sat*z;
    sat = surf(x,y,z, zeros(n_sat) );
    % rotate(sat, [0,0,1],rad2deg(b_sat_ear) )
 
    axis equal
    view(-20, 45)
    title("Albedo radiation reflected to satallite")
    hold off
    c=colorbar;
    c.LineWidth =2 ;
    c.FontSize = 16 ;
    clim([0 max(max(max(Q_mat)))]);
 
    xlim([-8e7 8e7])
    ylim([-8e7 8e7])
    zlim([-2e7 1e7])

%     f.FontSize = 16;
    exportgraphics(f,[path_img,'/Animate_year.gif'],'Append',true);

    k
end
 
%% 



for i = 1 :1000: numsteps
    f = figure('Visible','off');
    plot(time_steps(1:i)/max(time_steps),r_EarToSun(1:i),"LineWidth",2)
    hold on 
    scatter(time_steps(i)/max(time_steps),r_EarToSun(i) )
    hold off
    i
    xlim([0 1])
    ylim([min(min(r_EarToSun))   max(max(r_EarToSun))])
    xlabel("Time (years)")
    ylabel("Distance between sun and earth (m)")
    
    exportgraphics(f,[path_img,'/r_earth.gif'],'Append',true);
end

end
end

