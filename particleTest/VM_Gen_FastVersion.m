%function Voltage_mat = Voltage_mat_generation(t, u, v, w, x, y, z)
clc
clear
%-------------------------------------------------------------------------%
%                         INITIALIZATION    
%-------------------------------------------------------------------------%
% E = Total electric field matrix using Poisson's equation
% V = Potential matrix
% Nx = Number of grid points in X- direction
% Ny = Number of grid points in Y-Direction
%-------------------------------------------------------------------------%
% Enter the dimensions
Input_params = table2array(readtable('Grid_parameter_new.xlsx', 'Range','B1:B15'));
D1 = round(table2array(readtable('Input_Data_new.xlsx','Sheet',1, 'Range','C1:G3')));
D2 = round(table2array(readtable('Input_Data_new.xlsx', 'Sheet',2,'Range','C1:G3')));
%D = [D1' D2'];

Nx = round(Input_params(10,1)*Input_params(12,1))+1;     % Number of X-grids
Ny = round(Input_params(10,1)*Input_params(11,1))+1;     % Number of Y-grids
% mpx = ceil(Nx/2); % Mid-point of x
% mpy = ceil(Ny/2); % Mid point of y
Ni = 300;  % Number of iterations
V_m = zeros(Nx,Ny);   % Potential (Voltage) matrix


%-------------------------------------------------------------------------%
% Initializing Particle Properties
%-------------------------------------------------------------------------%
[V_m] = voltageCalcFunction(Ni,Nx,Ny,D1,D2,V_m);

writematrix(V_m, 'VoltageDistribution_17Aug.csv')  % writes the generated voltage matrix to a given name


A = gradient(V_m); % Gradient Matrix
[Ex,Ey] = gradient(V_m); % Gradient in X and Y axis
Ex = -Ex;
Ey = -Ey;


% Electric field Magnitude
E = sqrt(Ex.^2+Ey.^2);

x = (1:Nx);
y = (1:Ny);

% Contour Display for electric potential
figure(1)
contour_range_V = -1501:0.5:1501;
contour(x,y,V_m,contour_range_V,'linewidth',0.5);
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',14);
xlabel('x-axis in mesh nodes','fontsize',14);
ylabel('y-axis in mesh nodes','fontsize',14);
title('Electric Potential distribution, V(x,y) in volts','fontsize',14);
h1=gca;
set(h1,'fontsize',14);
fh1 = figure(1); 
set(fh1, 'color', 'white')
