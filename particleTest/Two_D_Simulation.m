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

qi= -1;
e_charge = 1.602e-19;
eps = 8.854*(10^-12);
size = 0.1; % this variable is needed for interpolation, it is the grid size (e.g. 0.1 mm)
m_Xe = 2.1801714e-25; %kg
Nx = 501;     % Number of X-grids
Ny = 501;     % Number of Y-grids
mpx = ceil(Nx/2); % Mid-point of x
mpy = ceil(Ny/2); % Mid point of y
% Take transpose for proper x-y orientation


V = csvread('Voltage_mat.csv');
A = gradient(V);
[Ex,Ey]=gradient(V);


% Electric field Magnitude
E = (qi/(4*pi*eps))* (1./Ex.^2+Ey.^2);  
x = (1:Nx);
y = (1:Ny);

% Contour Display for electric potential
figure(1)
contour_range_V = -201:0.5:1505;
contour(x,y,V,contour_range_V,'linewidth',0.05);
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',10);
xlabel('x-axis in meters','fontsize',10);
ylabel('y-axis in meters','fontsize',10);
title('Electric Potential distribution, V(x,y) in volts','fontsize',10);
h1=gca;
set(h1,'fontsize',10);
fh1 = figure(1);
set(fh1, 'color', 'white')

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


%-------------------------------------------------------------------------%
% Positon Simulation
%-------------------------------------------------------------------------%


del_t = 0.000001;
Nj = 10000;  % Number of iterations for the position sim

acc_x = zeros(1,Nj+1);
acc_y = zeros(1,Nj+1);
vel_x =  zeros(1,Nj+1);
vel_y =  zeros(1,Nj+1);
time = zeros(1,Nj+1);
t = zeros(1,Nj+1);


Pos_y(1,1) = 153; %ALERT THESE ARE TRANSPOSED PLEASE CHECK, x is y y is x?
Pos_x(1,1) = 250;
min_p = 0.005; %mm
max_p = 0.01; %mm
Vl_x = 0 * 10e-9; % 
Vl_y = 0 * 10e-9;
disp_x = zeros(1,Nj);
disp_y = zeros(1,Nj);
for z = 1:Nj    % Number of iterations
%     while i<51
%         while j<51    % Number of iterations
           
%     Fx = Ex(Pos_x(1,z),Pos_y(1,z)) * qi * e_charge;% need to interpolate the cahrge between points
%     Fy = Ey(Pos_x(1,z),Pos_y(1,z)) * qi * e_charge;
   

   xt_ceil = ceil(Pos_x(1,z));
   yt_ceil = ceil(Pos_y(1,z));
   xt_floor = floor(Pos_x(1,z));
   yt_floor = floor(Pos_y(1,z));
   Pos_x(1,z);
   
   Ext = Ex(xt_floor,yt_floor)+((Ex(xt_ceil,yt_ceil)-Ex(xt_floor,yt_floor))/size)*(Pos_x(1,z)-xt_floor);
   Fx = Ext * qi * e_charge;% need to interpolate the cahrge between points
   
   Eyt = Ey(xt_floor,yt_floor)+((Ey(xt_ceil,yt_ceil)-Ey(xt_floor,yt_floor))/size)*(Pos_y(1,z)-yt_floor);
   Fy = Eyt * qi * e_charge;% need to interpolate the cahrge between points
    
    acc_x(1,z) = (Fx./m_Xe) * 10e-9; %mm/ms^2
    acc_y(1,z) = (Fy./m_Xe) * 10e-9; %mm/ms^2
    vel_x(1,z) = Vl_x + (acc_x(1,z) * del_t);
    vel_y(1,z) = Vl_y + (acc_x(1,z) * del_t);
    disp_x(1,z) = (vel_x(1,z+1) + (acc_x(1,z)*(del_t)))*del_t;
    disp_y(1,z) = (vel_y(1,z+1) + (acc_y(1,z)*(del_t)))*del_t; 
    
    Pos_x(1,z+1) = Pos_x(1,z) +  disp_x(1,z);
    Pos_y(1,z+1) = Pos_y(1,z) +   disp_y(1,z);
   
        if abs(disp_x(1,z)) < min_p || abs(disp_y(1,z)) < min_p
            del_t = 2 * del_t;
             Vl_x = vel_x(1,z);
             Vl_y = vel_y(1,z);
             t(z) = t(z)+del_t;
             i = Pos_x(1, z+1);
             j = Pos_y(1, z+1);
            continue
        end
        
        if abs(disp_x(1,z)) > max_p || abs(disp_y(1,z)) > max_p
            del_t = del_t/2;
             Vl_x = vel_x(1,z);
             Vl_y = vel_y(1,z);
             t(z) = t(z)+del_t;
             i = Pos_x(1, z+1);
             j = Pos_y(1, z+1);
            continue
        end
    
end


%title('Acceleration','fontsize',14);
% hold on
% figure(5)
% plot(vel)
% %title('Velocity','fontsize',14, 'r');
% figure(6)
% plot(state_mat)
% %title('Position','fontsize', 14, 'k--');
% hold off
%  % this is a new model