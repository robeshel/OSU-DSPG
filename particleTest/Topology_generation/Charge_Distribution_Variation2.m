clc 
clear

%-------------------------------------------------------------------------%
%                           LOADING THE INPUT FILES
%-------------------------------------------------------------------------%
load('Test') %Variable number for the test iteration
date_yes = char(datetime(2020,09,10));
folder_name = char(['Test_Data\' date_yes]);
% Test = Test -2;
%-------------------------------------------------------------------------%

V = csvread([folder_name '\Test' num2str(Test) '_VtgDistMat.csv']);
SG_data = round(table2array(readtable('Input_Data\Input_Data_new.xlsx', 'Sheet', 1, 'Range','C1:G3')));


Nx = size(V, 1);
Ny = size(V, 2);
cd_mat = zeros(Nx, Ny);

%-------------------------------------------------------------------------%
% xn = SG_data(1,2); % Initial X position
% yn = SG_data(1,3)+1 : SG_data(2,1)-1; % Initial Y position
% Nj = 1;
pt1 = [SG_data(1,2),SG_data(1,3)+1];
pt2 = [SG_data(1,2),SG_data(2,1)-1];

[h,xn,yn] = CircleSegment(pt1,pt2);
size_mat = size(yn, 2);
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
f_bohm = 1;
v_bohm = f_bohm * sqrt(100* 1.60217662 * (10^-19) *3/(100 * 2.18017 * 10^-25));
tar_bm_crt = 0.001; %Ampere
beam_current = sum(cd_mat(:, end));
Clmb_chrg = 1.60217662 * 10^-19;
[Ex,Ey] = gradient(V);
del_x = 0.005 / Nx; %m
del_y = del_x; %m
A_cell = del_x * del_y; %m*m
Vx_in = v_bohm;
Vy_in = 0;


%-------------------------------------------------------------------------%
a = 0;
b = 0;
c = 0;
d = 0;
while beam_current < tar_bm_crt
    for p = 1: size(yn, 2)
        y1 = yn(1,p)-0.5;
        x1 = xn-0.5;
        Vx_in = v_bohm;
        Vy_in = 0;
        while x1 < Nx-1 || y1 < Ny-1 || x1 > 0 || y1 > 0
                [x, y, Vx_new, Vy_new, del_t] = Path_Calculation2(x1, y1, Ex, Ey, Vx_in,Vy_in, Nj);
                x2 = ceil(x);
                y2 = floor(y);
                x3 = ceil(x);
                y3 = floor(y);
                i = ceil(y1);
                j = ceil(x1);
                cd_mat(i,j) = cd_mat(i, j) + (100 * Clmb_chrg/A_cell) - (a + b + c + d);
                if del_t == 0
                    continue
                else
                    J_part =  cd_mat(i, j) / del_t;
                    a = cd_mat(i, j) + abs((100 * J_part * (del_t) * 10^-16 * (x2- x)*(y2-y)/(A_cell * del_x * del_y)));
                    b = cd_mat(i, j) + abs((100 * J_part * (del_t) * 10^-16 * (x3- x)*(y3-y)/(A_cell * del_x * del_y)));
                    c = cd_mat(i, j) + abs((100 * J_part * (del_t) * 10^-16 * (x2- x)*(y3-y)/(A_cell * del_x * del_y)));
                    d = cd_mat(i, j) + abs((100 * J_part * (del_t) * 10^-16 * (x3- x)*(y2-y)/(A_cell * del_x * del_y)));
                    cd_mat(i, j) = a + b + c + d;
                end
                
                if abs(sqrt((x1-x)^2 + (y1-y)^2)) < 0.0001
                    break
                else                    
                    x1 = x;
                    y1 = y;
                    Vx_in = Vx_new;
                    Vy_in = Vy_new;
                    beam_current = beam_current + (cd_mat(i, end-1) * Vx_new);
                end
        end
    end
end

%-------------------------------------------------------------------------%


x = (1:Nx);
y = (1:Ny);

% Contour Display Charge Distribution
figure(2)
contour_range_V = -0.2:0.05:0.2;
contour(x,y,cd_mat,contour_range_V,'linewidth',0.5);
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',14);
h1=gca;
set(h1,'fontsize',14);
fh1 = figure(1); 
set(fh1, 'color', 'white')

%-------------------------------------------------------------------------%
%                       Writing Output File
%-------------------------------------------------------------------------%
if isfolder(folder_name)
    writematrix(cd_mat, [folder_name '\Test' num2str(Test) '_Chrg_distrb_Mat_V2.csv'])  % writes the generated Charge distribution matrix to a given name
else
    mkdir(fullfile('Test_Data\', date))
    writematrix(cd_mat, [folder_name '\Test' num2str(Test) '_Chrg_distrb_Mat_V2.csv'])  % writes the generated Charge distribution matrix to a given name
end

%-------------------------------------------------------------------------%

