clc 
clear


load('Test') %Variable number for the test iteration'
%Test = Test - 2;
SG_data = round(table2array(readtable('Input_Data\Input_Data_new.xlsx', 'Sheet', 1, 'Range','C1:G3')));
pt1 = [SG_data(1,2),SG_data(1,3)+1];
pt2 = [SG_data(1,2),SG_data(2,1)-1];

[h,xn,yn] = CircleSegment(pt1,pt2);
size_mat = size(yn, 2);

date_yes = char(datetime(2020,10,09));
folder_name = ['Test_Data\' date_yes];
NPos_x = csvread([folder_name '\Test' num2str(Test) '_newT_NPos_x.csv']);  % writes the generated Trajectory matrix to a given name
NPos_x = reshape(NPos_x,[size(NPos_x,1)/size_mat, size_mat])'; 
NPos_y = csvread([folder_name '\Test' num2str(Test) '_newT_NPos_y.csv']);
NPos_y = reshape(NPos_y,[size(NPos_y,1)/size_mat, size_mat])';

V = csvread([folder_name '\Test' num2str(Test) '_VtgDistMat_newT.csv']);

% size1 = 6*10e-3/size(V, 1); % this variable is needed for interpolation, it is the grid size (e.g. 0.1 m)

Nx = size(V,1);
Ny  = size(V,2);
Ex = zeros(Nx,Ny);
Ey = zeros(Nx,Ny);

for i = 3:Nx-2
    for j = 3:Ny-2
    Ex(i,j) = -(-V(i, j+2) + 8*V(i,j+1) - 8*V(i,j-1) + V(i,j-2))/12;
    Ey(i,j) = -(-V(i+2, j) + 8*V(i+1,j) - 8*V(i-1,j) + V(i-2,j))/12;
    end
end

p_dim = size(NPos_x, 1);
pos_cellX = {};
pos_cellY = {};

for s = 1:p_dim
    [pos_cellX{s}, pos_cellY{s}] = Post_process(NPos_x(s,:), NPos_y(s,:));
end

% Electric field Magnitude
E = sqrt(Ex.^2+Ey.^2);
 
x = (1:Ny);
y = (1:Nx);

% Quiver Display for electric field Lines
figure(1)
contour(x,y,E,'linewidth',0.5);
hold on, quiver(x,y,Ex,Ey,2)
title('Electric field Lines, E (x,y) in V/m','fontsize',14);
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',14);
xlabel('x-axis in meters','fontsize',14);
ylabel('y-axis in meters','fontsize',14);
h3=gca;
set(h3,'fontsize',14);


