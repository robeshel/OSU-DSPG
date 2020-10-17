%function Primary_Trajectory_Generation
clc
clear all
%-------------------------------------------------------------------------%
%                         INITIALIZATION    
%-------------------------------------------------------------------------%
% E = Total electric field matrix using Poisson's equation
% V = Potential matrix
% Nx = Number of grid points in X- direction
% Ny = Number of grid points in Y-Direction
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%                           LOADING THE INPUT FILES
%-------------------------------------------------------------------------%
load('Test') %Variable number for the test iteration'
Test = Test;
date_yes = char(datetime(2020,10,09));
folder_name = ['Test_Data\' date_yes];
%-------------------------------------------------------------------------%
qi = 100;
e_charge = 1.602e-19;
eps = 8.854e-12;
m_Xe = 2.1801714e-25; %kg
f_bohm = 1;
v_bohm = f_bohm * sqrt(100* 1.60217662e-19 *3/(100 * 2.18017e-25));


V = csvread([folder_name '\Test' num2str(Test) '_VtgDistMat_newT.csv']);
GR_Data =  round(table2array(readtable('Input_Data\Grid_parameter_new.xlsx', 'Range','B11:B12' )));
SG_data = round(table2array(readtable('Input_Data\Input_Data_new.xlsx', 'Sheet', 1, 'Range','C1:G3')));
Node_factor = 1e-3 / GR_Data(1,1);

Nx = size(V,1);
Ny  = size(V,2);
Ex = zeros(Nx,Ny);
Ey = zeros(Nx,Ny);

for i = 3:Nx-2
    for j = 3:Ny-2
    Ex(i,j) = -(-V(i, j+2) + 8*V(i,j+1) - 8*V(i,j-1) + V(i,j-2))/(12*Node_factor);
    Ey(i,j) = -(-V(i+2, j) + 8*V(i+1,j) - 8*V(i-1,j) + V(i-2,j))/(12*Node_factor);
    end
end

% Electric field Magnitude
E = sqrt(Ex.^2+Ey.^2);
 
x = (1:Ny);
y = (1:Nx);

Nj = 100000;
%-------------------------------------------------------------------------%
%                       Positon Simulation
%-------------------------------------------------------------------------%
% xn = SG_data(1,2); % Initial X position
% yn = SG_data(1,3)+1 : SG_data(2,1)-1; % Initial Y position


pt1 = [SG_data(1,2),SG_data(1,3)+1];
pt2 = [SG_data(1,2),SG_data(2,1)-1];

[h,xn,yn] = CircleSegment(pt1,pt2);
size_mat = size(yn, 2);

NPos_x = zeros(size_mat,Nj); % X position matrix to multiple trajectories
NPos_y = zeros(size_mat,Nj); % Y position matrix to multiple trajectories
time_step = zeros(size_mat,Nj);
Vx_in= v_bohm;
Vy_in = 0;
Vx_new = zeros(size_mat,Nj);
Vy_new = zeros(size_mat,Nj);



for itr = 1:size(yn, 2)
    [NPos_x(itr,:), NPos_y(itr,:), Vx_new(itr, :), Vy_new(itr, :),  time_step(itr,:)] = Path_Calculation(xn(1,itr), yn(1,itr), Ex, Ey, Vx_in, Vy_in, Nj);
end
p_dim = size(NPos_x, 1);
pos_cellX = {};
pos_cellY = {};

for s = 1:p_dim
    [pos_cellX{s}, pos_cellY{s}] = Post_process(NPos_x(s,:), NPos_y(s,:));
end

%-------------------------------------------------------------------------%
%                        SAVING THE TRAJECTORIES
%-------------------------------------------------------------------------%


if isfolder(folder_name)
    writematrix(NPos_x, [folder_name '\Test' num2str(Test) '_newT_NPos_x.csv'])  % writes the generated Trajectory matrix to a given name
    writematrix(NPos_y, [folder_name '\Test' num2str(Test) '_newT_NPos_y.csv'])  % writes the generated Trajectory matrix to a given name
    writematrix(time_step, [folder_name '\Test' num2str(Test) '_newT_TimeStep.csv'])  % writes the generated Time Step matrix to a given name
    writematrix(Vx_new, [folder_name '\Test' num2str(Test) '_newT_Vx.csv'])  % writes the generated Velocity matrix to a given name
else
    mkdir(fullfile('Test_Data\', date))
    writematrix(NPos_x, [folder_name '\Test' num2str(Test) '_newT_NPos_x.csv']) % writes the generated Trajectory matrix to a given name
    writematrix(NPos_y, [folder_name '\Test' num2str(Test) '_newT_NPos_y.csv'])  % writes the generated Trajectory matrix to a given name
    writematrix(time_step, [folder_name '\Test' num2str(Test) '_newT_TimeStep.csv'])  % writes the generated Time Step matrix to a given name
    writematrix(Vx_new, [folder_name '\Test' num2str(Test) '_newT_Vx.csv'])  % writes the generated Velocity matrix to a given name
end

