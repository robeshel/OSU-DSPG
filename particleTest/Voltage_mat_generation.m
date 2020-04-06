clc
close all;
clear;
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
qi= 1;
e_charge = 1.602e-19;
eps = 8.854*(10^-12);
size = 0.1; % this variable is needed for interpolation, it is the grid size (e.g. 0.1 mm)
m_Xe = 2.1801714e-25; %kg
Nx = 1001;     % Number of X-grids
Ny = 1001;     % Number of Y-grids
mpx = ceil(Nx/2); % Mid-point of x
mpy = ceil(Ny/2); % Mid point of y
Ni = 500;  % Number of iterations
V = zeros(Nx,Ny);   % Potential (Voltage) matrix

% -------------------------------------------------------------------------%
% Initializing Particle Properties
% -------------------------------------------------------------------------%
res = 10;
t1 = 10*res; %Thickness of screen Grid
t2 = 16*res; %Thickness of Acc Grid
g = 18*res; % Gap Between Screen And Acc Grid
r_s = 25*res; % Radius of Screen Grid
r_a = 14*res; % Radius of Acc Grid
lp_s = 28*res;   % Length of plate in terms of number of grids
lp_a = 36*res;   % Length of plate in terms of number of grids  
pp_s = 26*res; %Position of plate_1 on x axis
pp_a = pp_s + t1 + g; %Position of plate_2 on x axis
eff = zeros(1,Ni);


for z = 1:Ni    % Number of iterations
    V_old = V(mpx,mpy);
    for i=2:Nx-1
        for j=2:Ny-1
% -------------------------------------------------------------------------%
                V(1:mpy - r_s, pp_s:pp_s+t1) = 1500;
                V(mpy + r_s:Ny, pp_s:pp_s+t1) = 1500;
                V(1:mpy - r_a, pp_a:pp_a+t2) = -200;
                V(mpy + r_a:Ny,  pp_a:pp_a+t2) = -200;
% -------------------------------------------------------------------------%
                
                V(i,j)=0.25*(V(i+1,j)+V(i-1,j)+V(i,j+1)+V(i,j-1));
                
%-------------------------------------------------------------------------%
% Edges potentials
%-------------------------------------------------------------------------%
                
                V(1,j) = (V(1,j-1) + V(2,j) + V(1,j+1))./3 ;
                V(Nx,j) = (V(Nx,j-1) + V(Nx-1,j) + V(Nx,j+1))./3;
                V(i,1) = (V(i-1,1) + V(i,2) + V(i+1,1))./3;
                V(i,Ny) = (V(i-1,Ny) + V(i, Ny-1) + V(i+1, Ny))./3;
                
%-------------------------------------------------------------------------%
%Corner potentials
%-------------------------------------------------------------------------%
                
                V(1,1) = 0.5*(V(1,2)+V(2,1));
                V(Nx,1) = 0.5*(V(Nx-1,1)+V(Nx,2));
                V(1,Ny) = 0.5*(V(1,Ny-1)+V(2,Ny));
                V(Nx,Ny) = 0.5*(V(Nx,Ny-1)+V(Nx-1,Ny));
                
%-------------------------------------------------------------------------%
        end
    end
    V_new = V(mpx,mpy);
    eff(1,z) = (V_new - V_old);
end
toc

writematrix(V, 'Voltage_mat10_5_2p8.csv') 