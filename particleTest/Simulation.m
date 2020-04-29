function [Pos_x, Pos_y] = Simulation(x1, y1, Ex, Ey)
qi = 1;
e_charge = 1.602e-19;
m_Xe = 2.1801714e-25; %kg
del_t = 2.5e-7; % seconds
Nj = 100000;  % Number of iterations for the position sim

acc_x = zeros(1,Nj);
acc_y = zeros(1,Nj);
Ext= zeros(1,Nj);
Eyt = zeros(1,Nj);
vel_x =  zeros(1,Nj);
vel_y =  zeros(1,Nj);
t = zeros(1,Nj);

Pos_x = zeros(1,Nj);
Pos_y = zeros(1,Nj);
Pos_x(1,1) = x1;
Pos_y(1,1) = y1;
min_p = 0.5*10e-3; %m
max_p = 1*10e-3; %m
Vl_x = 0;  % m/s
Vl_y = 0;
disp_x = zeros(1,Nj);
disp_y = zeros(1,Nj);
xt_ceil = zeros(1,Nj);
yt_ceil = zeros(1,Nj);
xt_floor = zeros(1,Nj);
yt_floor = zeros(1,Nj);
Fx = zeros(1,Nj);
Fy = zeros(1,Nj);


for z = 1:Nj-1    % Number of iterations
   
   xt_ceil(1,z) = ceil(Pos_x(1,z));
   yt_ceil(1,z) = ceil(Pos_y(1,z));
  
   xt_floor(1,z) = floor(Pos_x(1,z));
   yt_floor(1,z) = floor(Pos_y(1,z));
   
   if xt_ceil(1,z) == 1 || yt_ceil(1,z)== 1 || xt_ceil(1,z)== 501 || yt_ceil(1,z)== 501 || xt_floor(1,z)== 1 || yt_floor(1,z)== 1 || xt_floor(1,z)== 501 || yt_floor(1,z)== 501
       break
   else
   Ext(1,z) = Ex(yt_floor(1,z),xt_floor(1,z))+((Ex(yt_ceil(1,z),xt_ceil(1,z))-Ex(yt_floor(1,z),xt_floor(1,z))))*(Pos_x(1,z)-xt_floor(1,z));
   Fx(1,z) = Ext(1,z) * qi * e_charge;
   
   Eyt(1,z) = Ey(yt_floor(1,z),xt_floor(1,z))+((Ey(yt_ceil(1,z),xt_ceil(1,z))-Ey(yt_floor(1,z),xt_floor(1,z))))*(Pos_y(1,z)-yt_floor(1,z));
   Fy(1,z) = Eyt(1,z) * qi * e_charge;

    
    acc_x(1,z) = ((Fx(1,z))./m_Xe); %m/s^2
    acc_y(1,z) = ((Fy(1,z))./m_Xe); %m/s^2
    vel_x(1,z) = Vl_x + (acc_x(1,z) * del_t);
    vel_y(1,z) = Vl_y + (acc_y(1,z) * del_t);
    disp_x(1,z) = (vel_x(1,z) + (acc_x(1,z)*(del_t)))*del_t;
    disp_y(1,z) = (vel_y(1,z) + (acc_y(1,z)*(del_t)))*del_t; 
  
        if abs(sqrt((disp_x(1,z))^2 + (disp_y(1,z))^2)) < min_p
            del_t = 1.05 * del_t;
             t(z) = t(z)+del_t;
             Pos_x(1,z+1) = Pos_x(1,z);
             Pos_y(1,z+1) = Pos_y(1,z);
        
        elseif abs(sqrt((disp_x(1,z))^2 + (disp_y(1,z))^2)) > max_p
            del_t = 0.95 * del_t;
             t(z) = t(z)+del_t;
             Pos_x(1,z+1) = Pos_x(1,z);
             Pos_y(1,z+1) = Pos_y(1,z);

        else
             Vl_x =  vel_x(1,z);
             Vl_y =  vel_y(1,z);
             t(z) = t(z)+del_t;
             Pos_x(1,z+1) = Pos_x(1,z) +  disp_x(1,z);
             Pos_y(1,z+1) = Pos_y(1,z) +  disp_y(1,z);

         end
   end
end