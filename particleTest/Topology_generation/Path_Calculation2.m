function [Pos_x, Pos_y, vel_x, vel_y, del_t] = Path_Calculation2(x1, y1, Ex, Ey, Vx_in, Vy_in, N_itr)
qi = 100; % Intial charge on the particle e.g. single charge or a double charge
e_charge = 1.602e-19; % amount of charge on a electron
m_Xe = 2.1801714e-25; %Mass of a single Xenon atom
del_t = 2.5e-7; % Time delta in seconds
Nj = N_itr;  % Number of iterations for the position sim
mat_size_x = size(Ex, 1)-1;
mat_size_y = size(Ex, 2)-1;

acc_x = zeros(1,Nj); % Acceleration in X direction
acc_y = zeros(1,Nj); % Acceleration in Y direction
Ext= zeros(1,Nj); % Value of gradient for the position in X direction
Eyt = zeros(1,Nj); % Value of gradient for the position in Y direction
vel_x =  zeros(1,Nj); % Velocity in X direction
vel_y =  zeros(1,Nj); % Velocity in Y direction
t = zeros(1,Nj); % Time matrix for all iterartion, stores time in each iteration

Pos_x = zeros(1,Nj);  % Position in X direction
Pos_y = zeros(1,Nj); % Position in Y direction
Pos_x(1,1) = x1; % Initial Position in X direction
Pos_y(1,1) = y1; % Initial Position in Y direction
min_p = 1*10e-7; %m  % Minimum position delta 
max_p = 2*10e-7; %m % Maximum position delta 
Vl_x = Vx_in; % m/s % Initial Velocity in X direction
Vl_y = Vy_in; % Initial Velocity in X direction
disp_x = zeros(1,Nj); % Storing all displacement in X direction
disp_y = zeros(1,Nj); % Storing all displacement in Y direction
xt_ceil = zeros(1,Nj); % Storing upper bound in X direction for interpolation
yt_ceil = zeros(1,Nj); % Storing upper bound in Y direction for interpolation
xt_floor = zeros(1,Nj); % Storing lower bound in X direction for interpolation
yt_floor = zeros(1,Nj); % Storing Lower bound in Y direction for interpolation
Fx = zeros(1,Nj); %Force Matrix in X direction
Fy = zeros(1,Nj); % Force Matrix in Y direction
tot_disp_1 = 0;

for z = 1: Nj
while tot_disp_1> max_p || tot_disp_1 < min_p     % Number of iterations
   
   xt_ceil(1,z) = ceil(Pos_x(1,z));
   yt_ceil(1,z) = ceil(Pos_y(1,z));
  
   xt_floor(1,z) = floor(Pos_x(1,z));
   yt_floor(1,z) = floor(Pos_y(1,z));
   
   if xt_ceil(1,z) == 1 || yt_ceil(1,z)== 1 || xt_ceil(1,z) == mat_size_y || yt_ceil(1,z)== mat_size_x || xt_floor(1,z)== 1 || yt_floor(1,z)== 1 || xt_floor(1,z)== mat_size_y || yt_floor(1,z)== mat_size_x
       break % This will stop the loop if the position goes beyond the size of the consideration
   else
   Ext(1,z) = Ex(yt_floor(1,z),xt_floor(1,z))+((Ex(yt_ceil(1,z),xt_ceil(1,z))-Ex(yt_floor(1,z),xt_floor(1,z))))*(Pos_x(1,z)-xt_floor(1,z)); %Interpolation for electric field in X direction
   Fx(1,z) = Ext(1,z) * qi * e_charge;  %Force Calculation in X direction
   
   Eyt(1,z) = Ey(yt_floor(1,z),xt_floor(1,z))+((Ey(yt_ceil(1,z),xt_ceil(1,z))-Ey(yt_floor(1,z),xt_floor(1,z))))*(Pos_y(1,z)-yt_floor(1,z)); %Interpolation for electric field in Y direction
   Fy(1,z) = Eyt(1,z) * qi * e_charge; %Force Calculation in Y direction

    
    acc_x(1,z) = ((Fx(1,z))./m_Xe); %m/s^2
    acc_y(1,z) = ((Fy(1,z))./m_Xe); %m/s^2
    vel_x(1,z) = Vl_x + (acc_x(1,z) * del_t); %m/s
    vel_y(1,z) = Vl_y + (acc_y(1,z) * del_t); %m/s
    disp_x(1,z) = (vel_x(1,z) + (acc_x(1,z)*(del_t)))*del_t; %m
    disp_y(1,z) = (vel_y(1,z) + (acc_y(1,z)*(del_t)))*del_t; %m
    tot_disp_1 = abs(sqrt((disp_x(1,z))^2 + (disp_y(1,z))^2));
  
        if tot_disp_1 < min_p % This modifies the time delta so that the displacement is greater than the min allowed displacement
            del_t = 1.02 * del_t;
             t(z) = t(z)+del_t;
             Pos_x(1,z) = Pos_x(1,z);
             Pos_y(1,z) = Pos_y(1,z);
        
        elseif tot_disp_1 > max_p % This modifies the time delta so that the displacement is smaller than the max allowed displacement
            del_t = 0.98 * del_t;
             t(z) = t(z)+del_t;
             Pos_x(1,z) = Pos_x(1,z);
             Pos_y(1,z) = Pos_y(1,z);

        else 
             Vl_x =  vel_x(1,z);
             Vl_y =  vel_y(1,z);
             t(z) = t(z)+del_t;
             Pos_x(1,z) = Pos_x(1,z) + (disp_x(1,z)*10^5);
             Pos_y(1,z) = Pos_y(1,z) + (disp_y(1,z)*10^5);

        end
   end
end
end