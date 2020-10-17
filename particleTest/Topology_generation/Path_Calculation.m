function [Pos_x, Pos_y, vel_x, vel_y, time_step] = Path_Calculation(x1, y1, Ex, Ey, Vx_in, Vy_in, Nj)
qi = 100; % Intial charge on the particle e.g. single charge or a double charge
e_charge = 1.602e-19; % amount of charge on a electron Coloumb
m_Xe = 2.1801714e-25; %Mass of a single Xenon atom in Kg
del_t = 2.5e-9; % Time delta in seconds
mat_size_x = size(Ex, 1);

acc_x = zeros(1,Nj); % Acceleration in X direction
acc_y = zeros(1,Nj); % Acceleration in Y direction
Ext = zeros(1,Nj); % Value of gradient for the position in X direction
Eyt = zeros(1,Nj); % Value of gradient for the position in Y direction
vel_x =  zeros(1,Nj); % Velocity in X direction
vel_y =  zeros(1,Nj); % Velocity in Y direction
t = zeros(1,Nj); % Time matrix for all iterartion, stores time in each iteration


Pos_x = zeros(1,Nj);  % Position in X direction
Pos_y = zeros(1,Nj); % Position in Y direction
time_step = zeros(1,Nj); %Matrix to save the time steps 
Pos_x(1,1) = x1; % Initial Position in X direction
Pos_y(1,1) = y1; % Initial Position in Y direction
min_p = 5*10e-7; %m  % Minimum position delta 
max_p = 1*10e-6; %m  % Maximum position delta 
Vl_x = Vx_in;  % m/s % Initial Velocity in X direction
Vl_y = Vy_in;  % m/s % Initial Velocity in Y direction
disp_x = zeros(1,Nj); % Storing all displacement in X direction
disp_y = zeros(1,Nj); % Storing all displacement in Y direction
xt_ceil = zeros(1,Nj); % Storing upper bound in X direction for interpolation
yt_ceil = zeros(1,Nj); % Storing upper bound in Y direction for interpolation
xt_floor = zeros(1,Nj); % Storing lower bound in X direction for interpolation
yt_floor = zeros(1,Nj); % Storing Lower bound in Y direction for interpolation
Fx = zeros(1,Nj); %Force Matrix in X direction
Fy = zeros(1,Nj); % Force Matrix in Y direction



for z = 2:Nj   % Number of iterations
   
   xt_ceil(1,z) = ceil(Pos_x(1,z-1));
   yt_ceil(1,z) = ceil(Pos_y(1,z-1));
   xt_floor(1,z) = floor(Pos_x(1,z-1));
   yt_floor(1,z) = floor(Pos_y(1,z-1));
   Co_ordMat = [xt_ceil(1,z) yt_ceil(1,z) xt_floor(1,z) yt_floor(1,z)];
   
   if any(Co_ordMat(:) < 1) || any(Co_ordMat(:) > mat_size_x)
%        disp('1')
       continue % This will stop the loop if the position goes beyond the size of the consideration

   else
       Ext(1,z) = Ex(yt_floor(1,z),xt_floor(1,z)) + ((Ex(yt_ceil(1,z),xt_ceil(1,z))-Ex(yt_floor(1,z),xt_floor(1,z))))*(Pos_x(1,z)-xt_floor(1,z)); %Interpolation for electric field in X direction
       Fx(1,z) = Ext(1,z) * qi * e_charge;  %Force Calculation in X direction

       Eyt(1,z) = Ey(yt_floor(1,z),xt_floor(1,z))+((Ey(yt_ceil(1,z),xt_ceil(1,z))-Ey(yt_floor(1,z),xt_floor(1,z))))*(Pos_y(1,z)-yt_floor(1,z)); %Interpolation for electric field in Y direction
       Fy(1,z) = Eyt(1,z) * qi * e_charge; %Force Calculation in Y direction


        acc_x(1,z) = (Fx(1,z)./m_Xe); %m/s^2
        acc_y(1,z) = (Fy(1,z)./m_Xe); %m/s^2
        vel_x(1,z) = Vl_x + (acc_x(1,z) * del_t);
        vel_y(1,z) = Vl_y + (acc_y(1,z) * del_t);
        disp_x(1,z) = (vel_x(1,z) + (acc_x(1,z)*(del_t)))* del_t;
        disp_y(1,z) = (vel_y(1,z) + (acc_y(1,z)*(del_t))) * del_t;
        tot_disp = abs(sqrt(disp_x(1,z)^2 + disp_y(1,z)^2));
    
  
        if tot_disp < min_p % This modifies the time delta so that the displacement is greater than the min allowed displacement
             del_t = 1.10 * del_t;
             t(1, z) = t(1,z)+del_t;
             Pos_x(1,z) = Pos_x(1,z-1);
             Pos_y(1,z) = Pos_y(1,z-1);
        
        elseif tot_disp > max_p % This modifies the time delta so that the displacement is smaller than the max allowed displacement
            del_t = 0.90 * del_t;
             t(z) = t(z)+del_t;
             Pos_x(1,z) = Pos_x(1,z-1);
             Pos_y(1,z) = Pos_y(1,z-1);

        else 
             Vl_x =  vel_x(1,z);
             Vl_y =  vel_y(1,z);
             t(z) = t(z)+del_t;
             Pos_x(1,z) = Pos_x(1,z-1) +  disp_x(1,z) * 10^5;
%              Ds_x = disp_x(1,z)
%              Ps_x = Pos_x(1,z+1)
             Pos_y(1,z) = Pos_y(1,z-1) +  disp_y(1,z) * 10^5;
%              Ps_y = Pos_y(1,z+1)
%              Ds_y = disp_y(1,z)
             time_step(1,z) = del_t;

        end
   end
end