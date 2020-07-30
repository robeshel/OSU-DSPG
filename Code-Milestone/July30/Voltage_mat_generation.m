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
Input_params = table2array(readtable('Grid_parameter.xlsx', 'Range','B1:B15'));
D1 = round(table2array(readtable('Input_Data_multi.xlsx','Sheet',1, 'Range','C1:L21')));
D2 = round(table2array(readtable('Input_Data_multi.xlsx', 'Sheet',2,'Range','C1:L21')));
%D = [D1' D2'];

Nx = round(Input_params(10,1)*Input_params(12,1));     % Number of X-grids
Ny = round(Input_params(10,1)*Input_params(11,1));     % Number of Y-grids
% mpx = ceil(Nx/2); % Mid-point of x
% mpy = ceil(Ny/2); % Mid point of y
Ni = 300;  % Number of iterations
V_m = zeros(Nx,Ny);   % Potential (Voltage) matrix


%-------------------------------------------------------------------------%
% Initializing Particle Properties
%-------------------------------------------------------------------------%

tic
for f = 1:Ni    % Number of iterations
    for i=2:Nx-1
        for j=2:Ny-1
%-------------------------------------------------------------------------%
                V_m(i,j)=0.25*(V_m(i+1,j)+V_m(i-1,j)+V_m(i,j+1)+V_m(i,j-1)); % approximates the potential over the space
%-------------------------------------------------------------------------%
            for r = 1:20
                V_m(D1(r,4):D1(r,2), D1(r,1):D1(r,3)) = D1(r,5); % defines the potential on screen grid
                V_m(D1(r,9):D1(r,7), D1(r,6):D1(r,8)) = D1(r,10); % defines the potential on screen grid
                V_m(D2(r,4):D2(r,2), D2(r,1):D2(r,3)) = D2(r,5); % defines the potential on accelerator grid
                 V_m(D2(r,9):D2(r,7), D2(r,6):D2(r,8)) = D2(r,10); % defines the potential on accelerator grid
            end
               
%-------------------------------------------------------------------------%
% Edges potentials calculations
%-------------------------------------------------------------------------%
                
                V_m(1,j) = (V_m(1,j-1) + V_m(2,j) + V_m(1,j+1))./3 ;
                V_m(Nx,j) = (V_m(Nx,j-1) + V_m(Nx-1,j) + V_m(Nx,j+1))./3;
                V_m(i,1) = (V_m(i-1,1) + V_m(i,2) + V_m(i+1,1))./3;
                V_m(i,Ny) = (V_m(i-1,Ny) + V_m(i, Ny-1) + V_m(i+1, Ny))./3;
                
%-------------------------------------------------------------------------%
%Corner potentials calculations
%-------------------------------------------------------------------------%
                
                V_m(1,1) = 0.5*(V_m(1,2)+V_m(2,1));
                V_m(Nx,1) = 0.5*(V_m(Nx-1,1)+V_m(Nx,2));
                V_m(1,Ny) = 0.5*(V_m(1,Ny-1)+V_m(2,Ny));
                V_m(Nx,Ny) = 0.5*(V_m(Nx,Ny-1)+V_m(Nx-1,Ny));
                
%-------------------------------------------------------------------------%
        end
    end
end
toc
writematrix(V_m, 'VoltageDistribution_multiG.csv')  % writes the generated voltage matrix to a given name


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
