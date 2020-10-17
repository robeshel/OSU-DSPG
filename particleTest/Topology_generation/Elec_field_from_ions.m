%-------------------------------------------------------------------------%
%       Electric Fields due to ion in a 2-D plane using the Coulomb's Law
%-------------------------------------------------------------------------%
clc
close all; 
clear all;




%-------------------------------------------------------------------------%
%                   SYMBOLS USED IN THIS CODE                             
%-------------------------------------------------------------------------%
% E = Total electric field
% Ex = X-Component of Electric-Field
% Ey = Y-Component of Electric-Field
% Nx = Number of grid points in X- direction
% Ny = Number of grid points in Y-Direction
% eps_r = Relative permittivity
% r = distance between a selected point and the location of charge
% ex = unit vector for x-component electric field
% ey = unit vector for y-component electric field
e_charge = 1.602e-19; % amount of charge on a electron Coloumb
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%                           LOADING THE INPUT FILES
%-------------------------------------------------------------------------%
load('Test') %Variable number for the test iteration
date_yes = char(datetime(2020,9,10));
folder_name = ['Test_Data\' date_yes];
Test = Test -1;
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%                         INITIALIZATION                                  
%          Here, all the grid, size, charges, etc. are defined
%-------------------------------------------------------------------------%
% Reading the charge distribution matrix
cd_mat = csvread([folder_name '\Test' num2str(Test) '_Chrg_distrb_Mat.csv']);
% Constant 1/(4*pi*epsilon_0) = 9*10^9
k = 9*10^9;
% Enter the Relative permittivity
eps_r = 1;
charge_order = 1; % milli, micro, nano etc..
const = k*charge_order/eps_r;
% Enter the dimensions
Nx = size(cd_mat,1);
Ny = size(cd_mat,2);

% Electric fields Initialization
E_f = zeros(Nx,Ny);
Ex = E_f;
Ey = E_f;
% Vectors initialization
ex = E_f;
ey = E_f;
r = E_f;
r_square = E_f;


%-------------------------------------------------------------------------%
%                   COMPUTATION OF ELECTRIC FIELDS
%-------------------------------------------------------------------------%
%  Repeat for all the 'n' charges
for m = 1:Nx
    for n = 1:Ny
        q = cd_mat(m,n)* 10^-10;
        

        % Compute the unit vectors
        for i=m-15:m+15
            for j=n-15:n+15
                if i == m ||  j == n || i < 1 || j < 1 || i > 500 || j > 500
                  continue
                else
                r_square(i,j) = ((i-m)^2+(j-n)^2) * 10^-10;
%                 nhh = r_square(i,j)
%                 r(i,j) = abs(sqrt(r_square(i,j)));
                Ex(m,n) = Ex(m,n)+(q.*const./r_square(i,j));
%                 exx = ex(i,j)
                Ey(m,n) = Ey(m,n)+(q.*const./r_square(i,j));
%                 ayy = ey(i,j)
                end
            end
        end

    end
end

%-------------------------------------------------------------------------%
%                        SAVING THE TRAJECTORIES
%-------------------------------------------------------------------------%
if isfolder(folder_name)
    writematrix(Ex, [folder_name '\Test' num2str(Test) '_Ex_ion.csv'])  % writes the generated Electric Field matrix to a given name
    writematrix(Ey, [folder_name '\Test' num2str(Test) '_Ey_ion.csv'])  % writes the generated Electric Field matrix to a given name
else
    mkdir(fullfile('Test_Data\', date))
    writematrix(Ex, [folder_name '\Test' num2str(Test) '_Ex_ion.csv'])  % writes the generated Electric Field matrix to a given name
    writematrix(Ey, [folder_name '\Test' num2str(Test) '_Ey_ion.csv'])  % writes the generated Electric Field matrix to a given name
end


%-------------------------------------------------------------------------%
%                   PLOT THE RESULTS
%-------------------------------------------------------------------------%
% x_range = (1:Nx)-ceil(Nx/2);
% y_range = (1:Ny)-ceil(Ny/2);
% contour_range = -8:0.02:8;
contour(E_f', 'linewidth',0.7);
colorbar('location','eastoutside','fontsize',12);
xlabel('x ','fontsize',14);
ylabel('y ','fontsize',14);
title('Electric field distribution, E (x,y) in V/m','fontsize',14);

