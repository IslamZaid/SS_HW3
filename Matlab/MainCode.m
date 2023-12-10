clc;
close all;
clear all;
 %% 

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



b_ear_sun = pi/2;
% r_ear_sun = 1.496e+11;
r_sat_ear = 6.371e+6 + 2e+6 ; 
r_ear = 6.371e+6;
area_sat = 100;   %(m^2)


%% Sphere Model Development
% Number of nodes
numNodes = 40;
st = 0;
r_ear = 6.371e+6;  %(m);
[X,Y,Z] = sphere(numNodes);
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

srfc = surf(X,Y,Z);
Normal_Vectors = Area.*srfc.FaceNormals;



%% Radiation models
%% position vectors definition

R_ear_sun = @(r,beta_e_sun) r * [ cos(beta_e_sun) , sin(beta_e_sun) ,0 ];          % position vector of the earth respect to the sun
R_sat_ear = @(beta_sat_ear) r_sat_ear * [ cos(beta_sat_ear) , sin(beta_sat_ear) ,0 ] ;        % position vector of the sat respect to the earth
R_sat_sun = @(r,beta_e_sun, beta_sat_e) R_ear_sun(r,beta_e_sun) + R_sat_ear(beta_sat_e);    % position vector of the sat respect to the sun

unit_vector_normalto_sat = @(beta_sat_ear) -1 * R_sat_ear(beta_sat_ear) / norm(R_sat_ear(beta_sat_ear));
A_sat = @(beta_sat_ear)  area_sat * unit_vector_normalto_sat(beta_sat_ear);

 

%% Total radiation
%% %%%%%%%%%%%%%%%%%%%
tic
% Define time parameters
totaltime = 3.154e7; % Total simulation time (1 year in s)
time_steps = 0 :60 : totaltime ;
numsteps=length(time_steps);

%Time Implementation
% loop through time steps
n_earth = 2*pi/totaltime;  % rad/s
 
I_mat = zeros(numsteps,numNodes+1,numNodes+1);
Q_mat = zeros(numsteps,numNodes+1,numNodes);
Q_p_mat = zeros(numsteps,1);
Q_A_mat = zeros(numsteps,1);
Q_sun_mat = zeros(numsteps,1);
parfor step = 1:   numsteps
    t = time_steps(step);
    M = n_earth*(t); 
    theta  = M+(2*e-0.25*e^3)*sin(M)+(5/4)*e^2*sin(2*M)+(13/12)*e^3*sin(3*M);
    r = (h^2/mu)/(1+e*cos(theta));   % earth sun dist
    r_EarToSun(step) = r;
    b_ear_sun = theta;
    b_sat_ear = 2*pi*t / (127*60); %%%%%!!!!!!!
    % 1- Sun Model
    Q_sun = sun_model(r,r_ear,R_ear_sun(r,b_ear_sun),R_sat_sun(r,b_ear_sun, b_sat_ear), A_sat(b_ear_sun));
    %2- Planetary Albedo Model
    [Q_A,I,Q ]  = albedo_model(Normal_Vectors,centerNodes,A_sat(b_sat_ear), R_ear_sun(r,b_ear_sun),R_sat_ear(b_sat_ear),area_sat);
    % 3- Earth Model
    Q_p = earth_model();
    Q_gen = 10000; 
    Temp(step) =( ((Q_gen/eps)+( Q_p + ae*(Q_sun + Q_A)))/(2*area_sat*sigma) )^0.25;
    
    Q_p_mat(step) = Q_p ;
    Q_sun_mat(step) = Q_sun ;
    Q_A_mat(step) = Q_A ;

    theta_earth_mat(step) = theta;
    theta_sat_mat(step) = b_sat_ear;

    I_mat(step,:,:) = I; 
    Q_mat(step,:,:) = Q; 

    if rem(step,100) == 0
        step
    end
end
 


toc
  
 %%
% figure
% f_sun = 5e+3;
% pe = R_ear_sun(r_ear, b_ear_sun)/f_sun;
% psat = R_sat_ear(b_sat_ear);
% 
% subplot(2,2,1)
% quiver3(0,0,0,pe(1), pe(2), pe(3)) 
% hold on
% quiver3(pe(1),pe(2),pe(3),psat(1), psat(2), psat(3))
% % sun 
% r_sun = 7e+9/f_sun;
% [x,y,z] = sphere;
% x = r_sun*x; y = r_sun*y; z = r_sun*z; 
% surf(x,y,z, ones(size(z)))
% % earth
% x = pe(1)+X; y = pe(2)+Y; z = pe(3)+Z; 
% surf(x,y,z, I, 'EdgeAlpha',0.2)
% %satallite
% pss = pe + psat ;
% f_sat = 1e+6;
% n_sat = 10;
% [x,y,z] = ellipsoid(pss(1), pss(2), pss(3), 1*f_sat,1*f_sat,1*f_sat, n_sat);
% % x = f_sat * x; y = f_sat* y; z = f_sat*z;
% sat = surf(x,y,z, zeros(n_sat) )
% % rotate(sat, [0,0,1],rad2deg(b_sat_ear) )
% axis equal
% view(0, 90)
% title("Sun radiation on the earth")
% 
% 
% 
% 
% subplot(2,2,2)
% quiver3(0,0,0,pe(1), pe(2), pe(3)) 
% hold on
% quiver3(pe(1),pe(2),pe(3),psat(1), psat(2), psat(3))
% % sun 
% r_sun = 7e+9/f_sun;
% [x,y,z] = sphere;
% x = r_sun*x; y = r_sun*y; z = r_sun*z; 
% surf(x,y,z, ones(size(z)))
% % earth
% x = pe(1)+X; y = pe(2)+Y; z = pe(3)+Z; 
% surf(x,y,z, Q, 'EdgeAlpha',0.2)
% %satallite
% pss = pe + psat ;
% f_sat = 1e+6;
% n_sat = 10;
% [x,y,z] = ellipsoid(pss(1), pss(2), pss(3), 1*f_sat,1*f_sat,1*f_sat, n_sat);
% % x = f_sat * x; y = f_sat* y; z = f_sat*z;
% sat = surf(x,y,z, zeros(n_sat) )
% % rotate(sat, [0,0,1],rad2deg(b_sat_ear) )
% axis equal
% view(90+rad2deg(b_sat_ear), 50)
% title("Sun radiation on the earth")
%  
% 
% 
% 
% 
% % subplot(2,2,3)
% % quiver3(0,0,0,pe(1), pe(2), pe(3)) 
% % hold on
% % quiver3(pe(1),pe(2),pe(3),psat(1), psat(2), psat(3))
% % % sun 
% % r_sun = 7e+9/f_sun;
% % [x,y,z] = sphere;
% % x = r_sun*x; y = r_sun*y; z = r_sun*z; 
% % surf(x,y,z, ones(size(z)))
% % % earth
% % x = pe(1)+X; y = pe(2)+Y; z = pe(3)+Z; 
% % surf(x,y,z, I)
% % axis equal
% % 
% % subplot(2,2,4)
% % quiver3(0,0,0,pe(1), pe(2), pe(3)) 
% % hold on
% % quiver3(pe(1),pe(2),pe(3),psat(1), psat(2), psat(3))
% % % sun 
% % r_sun = 7e+9/f_sun;
% % [x,y,z] = sphere;
% % x = r_sun*x; y = r_sun*y; z = r_sun*z; 
% % surf(x,y,z, ones(size(z)))
% % % earth
% % x = pe(1)+X; y = pe(2)+Y; z = pe(3)+Z; 
% % surf(x,y,z, I)
% % axis equal
% % 
% % 
% % 
% % 
% % 

%% 
%    figure
% h = animatedline;
% axis([0,4*pi,-1,1])
% 
% 
% for k = 1:1000
%   r_eartosun = 1.496e11;
% f_sun = 1e+3;
% b_ear_sun = theta_earth_mat(k);
% b_sat_ear=theta_sat_mat(k);
% pe = R_ear_sun(r_eartosun, b_ear_sun)/f_sun/2.5;
% psat = R_sat_ear(b_sat_ear);
% 
% % subplot(2,2,1)
% quiver3(0,0,0,pe(1), pe(2), pe(3)) 
% hold on
% quiver3(pe(1),pe(2),pe(3),psat(1), psat(2), psat(3))
% % sun 
% r_sun = 7e+9/f_sun;
% [x,y,z] = sphere;
% x = r_sun*x; y = r_sun*y; z = r_sun*z; 
% surf(x,y,z, 10000*ones(size(z)),"EdgeColor","red","LineStyle","-.")
% % earth
% x = pe(1)+X; y = pe(2)+Y; z = pe(3)+Z; 
% % surf(x,y,z, squeeze(I_mat(k,:,:)), 'EdgeAlpha',0.2)
% Q_mat(:,:,33) = Q_mat(:,:,end);
% surf(x,y,z, squeeze(Q_mat(k,:,:)), 'EdgeAlpha',0.2)
% %satallite
% pss = pe + psat ;
% f_sat = 1e+6;
% n_sat = 10;
% [x,y,z] = ellipsoid(pss(1), pss(2), pss(3), 1*f_sat,1*f_sat,1*f_sat, n_sat);
% % x = f_sat * x; y = f_sat* y; z = f_sat*z;
% sat = surf(x,y,z, zeros(n_sat) );
% % rotate(sat, [0,0,1],rad2deg(b_sat_ear) )
% axis equal
% view(0, 45)
% title("Sun radiation on the earth")
% hold off
% colorbar; 
% caxis([0 10000]);
% 
% xlim([-1e7 8e7])
% ylim([-1e7 1e7])
%     drawnow
% 
% end

cond = 0;
if cond == 1 

if size(Q_mat,2) ~= size(Q_mat,3)
    Q_mat(:,:,end+1) = Q_mat(:,:,end);
end
k = 80 %
f=figure();  %'Position', [0, 0, 1500, 1200],
r_eartosun = a;
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
caxis([0 max(max(max(Q_mat)))]);

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
exportgraphics(f,'images/albedo_elements.png', Resolution=1200);

%%


if size(Q_mat,2) ~= size(Q_mat,3)
    Q_mat(:,:,end+1) = Q_mat(:,:,end);
end
k = 80 %
f=figure();  %'Position', [0, 0, 1500, 1200],
r_eartosun = a;
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
caxis([0 max(max(max(I_mat)))]);

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
exportgraphics(f,'images/inc_sun_rad.png', Resolution=1200);

%%
 
f=figure();  %'Position', [0, 0, 1500, 1200],
plot(Temp(1:60*24), 'LineWidth',2)
grid on 
grid minor 
title("Temperature variation over 24 hours")
xlabel("Time (hours)")
ylabel("Temperature (K)")
xticks(0:6*60:24*60)
xticklabels(0:6:24)
ax=axes(Position=[0.5 0.25 0.4 0.4 ]);
plot(Temp(0*60+1:60*2+1), 'LineWidth',2)
xticks(0:0.5*60:2*60)
xticklabels(0:0.5:2)
grid on 
grid minor 
fontsize(gcf,16,"points")

exportgraphics(f,'images/temp_time.png', Resolution=1200);


%%
 
f=figure('Position', [0, 0, 1200, 600]);  %'Position', [0, 0, 1500, 1200],
period = 0*24*60+1 : 1*24*60;
tickss = 0:6*60:24*60;
ticksslabels = 0:6:24;

plot(Q_A_mat(period), 'LineWidth',2)
hold on 
plot(Q_p_mat(period), 'LineWidth',2)
plot(Q_sun_mat(period), 'LineWidth',2)
grid on 
grid minor 
title("Temperature variation over 24 hours")
xlabel("Time (hours)")
ylabel("Temperature (K)")
xticks(tickss)
xticklabels(ticksslabels)

% ax=axes(Position=[0.5 0.25 0.4 0.4 ]);
% plot(Temp(0*60+1:60*2+1), 'LineWidth',2)
% xticks(0:0.5*60:2*60)
% xticklabels(0:0.5:2)
% grid on 
% grid minor 

fontsize(gcf,16,"points")

exportgraphics(f,'images/QQQ_time.png', Resolution=1200);

%% Initialize gif
if size(Q_mat,2) ~= size(Q_mat,3)
    Q_mat(:,:,end+1) = Q_mat(:,:,end);
end
timeloop = 1:120 ;

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
    caxis([0 max(max(max(Q_mat)))]);
 
    xlim([-1e7 8e7])
    ylim([-1e7 2e7])
    zlim([-2e7 1e7])

%     f.FontSize = 16;
    exportgraphics(f,'images/Animate_sat.gif', 'Resolution', 1200,  'BackgroundColor', 'none','Append',true);

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
    
    exportgraphics(f,'images/temp.gif','Append',true);
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
    caxis([0 max(max(max(Q_mat)))]);
 
    xlim([-8e7 8e7])
    ylim([-8e7 8e7])
    zlim([-2e7 1e7])

%     f.FontSize = 16;
    exportgraphics(f,'images/Animate_year.gif','Append',true);

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
    
    exportgraphics(f,'images/r_earth.gif','Append',true);
end


end

