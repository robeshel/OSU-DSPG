clc
close all;
clear all;
tic
%-------------------------------------------------------------------------%
%                         INITIALIZATION    
%-------------------------------------------------------------------------%
% E = Total electric field matrix using Poisson's equation
% V = Potential matrix
% Nx = Number of grid points in X- direction
% Ny = Number of grid points in Y-Direction
%-------------------------------------------------------------------------%
% Enter the dimensions
qi= -1;
e_charge = 1.602e-19;
eps = 8.854*(10^-12);
res = 20;
size = 0.1; % this variable is needed for interpolation, it is the grid size (e.g. 0.1 mm)
m_Xe = 2.1801714e-25; %kg
Nx = 401;     % Number of X-grids
Ny = 401;     % Number of Y-grids
mpx = ceil(Nx/2); % Mid-point of x
mpy = ceil(Ny/2); % Mid point of y
Ni = 50;  % Number of iterations for the position sim
V = zeros(Nx,Ny);   % Potential (Voltage) matrix
T = 0;            % Top-wall potential
B = 0;            % Bottom-wall potential
L = 30;            % Left-wall potential
R = 0;            % Right-wall potential
%-------------------------------------------------------------------------%
% Initializing edges potentials
%-------------------------------------------------------------------------%
V(1,:) = L;
V(Nx,:) = R;
V(:,1) = B;
V(:,Ny) = T;
%-------------------------------------------------------------------------%
% Initializing Corner potentials
%-------------------------------------------------------------------------%
V(1,1) = 0.5*(V(1,2)+V(2,1));
V(Nx,1) = 0.5*(V(Nx-1,1)+V(Nx,2));
V(1,Ny) = 0.5*(V(1,Ny-1)+V(2,Ny));
V(Nx,Ny) = 0.5*(V(Nx,Ny-1)+V(Nx-1,Ny));
%-------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
% Initializing Particle Properties
% -------------------------------------------------------------------------%
mat_size = 51; % particle matrix size
Position_mat = zeros(mat_size,mat_size);
Time_mat = zeros(mat_size,mat_size);%initializing Particle Matrix
% -------------------------------------------------------------------------%
res = 4;
t1 = 10*res; %Thickness of screen Grid
t2 = 16*res; %Thickness of Acc Grid
g = 18*res; % Gap Between Screen And Acc Grid
r_s = 32*res; % Radius of Screen Grid
r_a = 14*res; % Radius of Acc Grid
lp_s = 28*res;   % Length of plate in terms of number of grids
lp_a = 36*res;   % Length of plate in terms of number of grids  
pp_s = 26*res; %Position of plate_1 on x axis
pp_a = pp_s + t1 + g; %Position of plate_2 on x axis


for z = 1:Ni    % Number of iterations
    for i=2:Nx-1
        for j=2:Ny-1
% -------------------------------------------------------------------------%
                V(1:mpy - r_s, pp_s:pp_s+t1) = 1500;
                V(mpy + r_s:Ny, pp_s:pp_s+t1) = 1500;
                V(1:mpy - r_a, pp_a:pp_a+t2) = -200;
                V(mpy + r_a:Ny,  pp_a:pp_a+t2) = -200;
% -------------------------------------------------------------------------%
                V(i,j)=0.25*(V(i+1,j)+V(i-1,j)+V(i,j+1)+V(i,j-1));
        end
    end      
end
% Take transpose for proper x-y orientation
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
Nj = 500;  % Number of iterations for the position sim

acc_x = zeros(1,Nj+1);
acc_y = zeros(1,Nj+1);
vel_x =  zeros(1,Nj+1);
vel_y =  zeros(1,Nj+1);
time = zeros(1,Nj+1);
t = zeros(1,Nj+1);


Pos_x(1,1) = 103; %ALERT THESE ARE TRANSPOSED PLEASE CHECK, x is y y is x?
Pos_y(1,1) = 229;
min_p = 0.01; %mm
max_p = 0.02; %mm
Vl_x = 0 * 10e-9; % 
Vl_y = 0 * 10e-9;
for z = 1:Nj    % Number of iterations
%     while i<51
%         while j<51    % Number of iterations
           
%     Fx = Ex(Pos_x(1,z),Pos_y(1,z)) * qi * e_charge;% need to interpolate the cahrge between points
%     Fy = Ey(Pos_x(1,z),Pos_y(1,z)) * qi * e_charge;
   

   xt_ceil = ceil(Pos_x(1,z))
   yt_ceil = ceil(Pos_y(1,z))
   xt_floor = floor(Pos_x(1,z))
   yt_floor = floor(Pos_y(1,z))
   Pos_x(1,z)
   
   Ext = Ex(xt_floor,yt_floor)+((Ex(xt_ceil,yt_ceil)-Ex(xt_floor,yt_floor))/size)*(Pos_x(1,z)-xt_floor);
   Fx = Ext * qi * e_charge;% need to interpolate the cahrge between points
   
   Eyt = Ey(xt_floor,yt_floor)+((Ey(xt_ceil,yt_ceil)-Ey(xt_floor,yt_floor))/size)*(Pos_y(1,z)-yt_floor);
   Fy = Eyt * qi * e_charge;% need to interpolate the cahrge between points
   
   z
    acc_x(1,z) = (Fx./m_Xe) * 10e-9; %mm/ms^2
    acc_y(1,z) = (Fy./m_Xe) * 10e-9; %mm/ms^2
    vel_x(1,z) = Vl_x + (acc_x(1,z) * del_t);
    vel_y(1,z) = Vl_y + (acc_x(1,z) * del_t);
    Pos_x(1,z+1) = Pos_x(1,z) + (vel_x(1,z+1) + (acc_x(1,z)*(del_t)))*del_t;
    Pos_y(1,z+1) = Pos_y(1,z) + (vel_y(1,z+1) + (acc_y(1,z)*(del_t)))*del_t;  
   
        if abs(Pos_x(1,z) - Pos_x(1,z+1)) < min_p || abs(Pos_y(1,z) - Pos_y(1,z+1)) < min_p
            del_t = 2 * del_t;
             Vl_x = vel_x(1,z);
             Vl_y = vel_y(1,z);
             t(z) = t(z)+del_t;
             i = Pos_x(1, z+1);
             j = Pos_y(1, z+1);
            continue
        end
        if abs(Pos_x(1,z) - Pos_x(1,z+1)) > max_p || abs(Pos_y(1,z) - Pos_y(1,z+1)) > max_p
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