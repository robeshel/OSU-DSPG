clc 
clear 
close all

%-------------------------------------------------------------------------%
%                         INITIALIZATION    
%-------------------------------------------------------------------------%
% E = Total electric field matrix using Poisson's equation
% V = Potential matrix
% Nx = Number of grid points in X- direction
% Ny = Number of grid points in Y-Direction
%-------------------------------------------------------------------------%

V = csvread('VM_withPlasmaPot.csv');

size = 5*10e-3/(length(V)-1); % It is the grid size (e.g. 0.1 m)

Nx = length(V);
Ny = length(V);

mpx = ceil(Nx/2); % Mid-point of x
mpy = ceil(Ny/2); % Mid point of y

A = gradient(V);
[Ex,Ey] = gradient(V);
 Ex = -Ex;
 Ey = -Ey;


% Electric field Magnitude
E = sqrt(Ex.^2+Ey.^2);
%[Ex, Ey] = gradient(E);
x = (1:Nx);
y = (1:Ny);

% Contour Display for electric potential
% figure(1)
% contour_range_V = -201:0.5:1505;
% contour(x,y,V,contour_range_V,'linewidth',0.05);
% axis([min(x) max(x) min(y) max(y)]);
% colorbar('location','eastoutside','fontsize',10);
% xlabel('x-axis in meters','fontsize',10);
% ylabel('y-axis in meters','fontsize',10);
% title('Electric Potential distribution, V(x,y) in volts','fontsize',10);
% h1=gca;
% set(h1,'fontsize',10);
% fh1 = figure(1);
% set(fh1, 'color', 'white')

% Quiver Display for electric field Lines
figure(2)
contour(x,y,E,'linewidth',0.5);
hold on, quiver(x,y,Ex,Ey,2)
title('Electric field Lines, E (x,y) in V/m','fontsize',14);
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',14);
xlabel('x-axis in meters','fontsize',14);
ylabel('y-axis in meters','fontsize',14);
h3=gca;
set(h3,'fontsize',14);
fh3 = figure(2);
set(fh3, 'color', 'white')

% figure(3)
% contour(Ex,'linewidth',0.5)
% axis([min(x) max(x) min(y) max(y)]);
% colorbar('location','eastoutside','fontsize',14);
% xlabel('x-axis in meters','fontsize',14);
% ylabel('y-axis in meters','fontsize',14);
% h3=gca;
% set(h3,'fontsize',14);
% fh3 = figure(3);
% set(fh3, 'color', 'white')
% 
% figure(4)
% contour(Ey,'linewidth',0.5)

Nj = 100000;
%-------------------------------------------------------------------------%
%                       Positon Simulation
%-------------------------------------------------------------------------%
NPos_x = zeros(25,Nj);
NPos_y = zeros(25,Nj);
for t = 1:26
    r = 126:10:376;
    [NPos_x(t,:), NPos_y(t,:)] = Simulation(110, r(1,t), Ex, Ey);
end
% x1 = 200;
% y1 = 250;
% [NPos_x, NPos_y] = Simulation(x1, y1, Ex, Ey);



    
    
% qi = 1;
% e_charge = 1.602e-19;
% m_Xe = 2.1801714e-25; %kg
% del_t = 2.5e-7; % seconds
% Nj = 10000;  % Number of iterations for the position sim
% 
% acc_x = zeros(1,Nj);
% acc_y = zeros(1,Nj);
% Ext= zeros(1,Nj);
% Eyt = zeros(1,Nj);
% vel_x =  zeros(1,Nj);
% vel_y =  zeros(1,Nj);
% t = zeros(1,Nj);
% 
% Pos_x = zeros(1,Nj);
% Pos_y = zeros(1,Nj);
% Pos_x(1,1) = 130;
% Pos_y(1,1) = 213;
% min_p = 0.5*10e-3; %m
% max_p = 1*10e-3; %m
% Vl_x = 0 ;  % m/s
% Vl_y = 0;
% disp_x = zeros(1,Nj);
% disp_y = zeros(1,Nj);
% xt_ceil = zeros(1,Nj);
% yt_ceil = zeros(1,Nj);
% xt_floor = zeros(1,Nj);
% yt_floor = zeros(1,Nj);
% Fx = zeros(1,Nj);
% Fy = zeros(1,Nj);
% 
% 
% for z = 1:Nj-1    % Number of iterations
%    
%    xt_ceil(1,z) = ceil(Pos_x(1,z));
%    yt_ceil(1,z) = ceil(Pos_y(1,z));
%   
%    xt_floor(1,z) = floor(Pos_x(1,z));
%    yt_floor(1,z) = floor(Pos_y(1,z));
%    
%    if xt_ceil(1,z) == 1 || yt_ceil(1,z)== 1 || xt_ceil(1,z)== 501 || yt_ceil(1,z)== 501 || xt_floor(1,z)== 1 || yt_floor(1,z)== 1 || xt_floor(1,z)== 501 || yt_floor(1,z)== 501
%        break
%    else
%    Ext(1,z) = Ex(yt_floor(1,z),xt_floor(1,z))+((Ex(yt_ceil(1,z),xt_ceil(1,z))-Ex(yt_floor(1,z),xt_floor(1,z))))*(Pos_x(1,z)-xt_floor(1,z));
%    Fx(1,z) = Ext(1,z) * qi * e_charge;
%    
%    Eyt(1,z) = Ey(yt_floor(1,z),xt_floor(1,z))+((Ey(yt_ceil(1,z),xt_ceil(1,z))-Ey(yt_floor(1,z),xt_floor(1,z))))*(Pos_y(1,z)-yt_floor(1,z));
%    Fy(1,z) = Eyt(1,z) * qi * e_charge;
% 
%     
%     acc_x(1,z) = ((Fx(1,z))./m_Xe); %m/s^2
%     acc_y(1,z) = ((Fy(1,z))./m_Xe); %m/s^2
%     vel_x(1,z) = Vl_x + (acc_x(1,z) * del_t);
%     vel_y(1,z) = Vl_y + (acc_y(1,z) * del_t);
%     disp_x(1,z) = (vel_x(1,z) + (acc_x(1,z)*(del_t)))*del_t;
%     disp_y(1,z) = (vel_y(1,z) + (acc_y(1,z)*(del_t)))*del_t; 
%   
%         if abs(sqrt((disp_x(1,z))^2 + (disp_y(1,z))^2)) < min_p
%             del_t = 1.05 * del_t;
%              t(z) = t(z)+del_t;
%              Pos_x(1,z+1) = Pos_x(1,z);
%              Pos_y(1,z+1) = Pos_y(1,z);
%         
%         elseif abs(sqrt((disp_x(1,z))^2 + (disp_y(1,z))^2)) > max_p
%             del_t = 0.95 * del_t;
%              t(z) = t(z)+del_t;
%              Pos_x(1,z+1) = Pos_x(1,z);
%              Pos_y(1,z+1) = Pos_y(1,z);
% 
%         else
%              Vl_x =  vel_x(1,z);
%              Vl_y =  vel_y(1,z);
%              t(z) = t(z)+del_t;
%              Pos_x(1,z+1) = Pos_x(1,z) +  disp_x(1,z);
%              Pos_y(1,z+1) = Pos_y(1,z) +  disp_y(1,z);
% 
%          end
%    end
% end