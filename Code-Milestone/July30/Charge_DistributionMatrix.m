clc 
clear
V = csvread('Voltage_mat.csv');
Nx = size(V, 1);
Ny = size(V,2);
del_t = 2.5e-7;
cd_mat = zeros(Nx-1, Ny-1);

V_mesh = zeros(Nx, Ny);
Vx = zeros(Nx-1, Ny-1);
f_bohm = 1;
v_bohm = f_bohm * sqrt(100*3/(100 * 2.18017 * 10^-25));
tar_bm_crt = 0.001; %Ampere
beam_current = sum(cd_mat(:, end) * v_bohm);
[Ex,Ey] = gradient(V);
xn = [150:10:350];
yn = 135;
del_x = 0.005 / (Nx-1);
del_y = del_x;
A_cell = del_x * del_y;
Nj = 1;
tic
while beam_current < tar_bm_crt
    cd_mat(150:350, 135) = cd_mat(135, 150:350) +  100;
    for p = 1: size(xn, 2)
        x1 = xn(1,p);
        y1= yn;
        for i = 1:Nx-1
            for j = 1:Ny-1
                [x, y, Vx] = Simulation(x1, y1, Ex, Ey, v_bohm, Nj)
                x2 = ceil(x);
                y2 = ceil(y);
                x3 = ceil(x);
                y3 = floor(y);
                J_part =  cd_mat(i, j) / del_t;
                cd_mat(i, j) = cd_mat(i, j) + (100 * J_part * (del_t) * (x2- x)*(y2-y)/(A_cell * del_x * del_y));
                cd_mat(i, j) = cd_mat(i, j) + (100 * J_part * (del_t) * (x3- x)*(y3-y)/(A_cell * del_x * del_y));
                x1 = x;
                y1 = y;
            end 
            beam_current = beam_current + (cd_mat(i, end) * Vx)
        end
         
    end
end
toc

    