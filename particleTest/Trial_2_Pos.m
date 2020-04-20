clc 
clear 
close all

%-------------------------------------------------------------------------%
%                         INITIALIZATION    
%-------------------------------------------------------------------------%
% E = Total electric field matrix using Poisson's equation
% V = Potential matrix
% Nx = Number of grid points in X- direction
% Ny = Number of grid points in Y-Direction
%-------------------------------------------------------------------------%

qi = 1;
e_charge = 1.602e-19;
eps = 8.854*(10^-12);
m_Xe = 2.1801714e-25; %kg

V = csvread('Voltage_mat510_5_2p8.csv');

size = 5*10e-3/(length(V)-1); % this variable is needed for interpolation, it is the grid size (e.g. 0.1 m)

Nx = length(V);
Ny = length(V);

mpx = ceil(Nx/2); % Mid-point of x
mpy = ceil(Ny/2); % Mid point of y

A = gradient(V);
[Ex,Ey] = gradient(V);
Ex = -Ex/size;
Ey = -Ey/size;


% Electric field Magnitude
E = sqrt(Ex.^2+Ey.^2)
 
x = (1:Nx);
y = (1:Ny);

% Contour Display for electric potential
% figure(1)
% contour_range_V = -201:0.5:1505;
% contour(x,y,V,contour_range_V,'linewidth',0.05);
% axis([min(x) max(x) min(y) max(y)]);
% colorbar('location','eastoutside','fontsize',10);
% xlabel('x-axis in meters','fontsize',10);
% ylabel('y-axis in meters','fontsize',10);
% title('Electric Potential distribution, V(x,y) in volts','fontsize',10);
% h1=gca;
% set(h1,'fontsize',10);
% fh1 = figure(1);
% set(fh1, 'color', 'white')

% Quiver Display for electric field Lines
figure(2)
contour(x,y,E,'linewidth',0.5);
hold on, quiver(x,y,-Ex,-Ey,2)
title('Electric field Lines, E (x,y) in V/m','fontsize',14);
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',14);
xlabel('x-axis in meters','fontsize',14);
ylabel('y-axis in meters','fontsize',14);
h3=gca;
set(h3,'fontsize',14);
fh3 = figure(2);
set(fh3, 'color', 'white')


%-------------------------------------------------------------------------%
% Positon Simulation
%-------------------------------------------------------------------------%


del_t = 0.151; % seconds
Nj = 5;  % Number of iterations for the position sim

acc_x = zeros(1,Nj+1);
acc_y = zeros(1,Nj+1);
Ext= zeros(1,Nj+1);
Eyt = zeros(1,Nj+1);
vel_x =  zeros(1,Nj+1);
vel_y =  zeros(1,Nj+1);
t = zeros(1,Nj+1);

Pos_x = zeros(1,Nj);
Pos_y = zeros(1,Nj);
Pos_x(1,1) = 122;
Pos_y(1,1) = 344;
min_p = 0.5*10e-3; %m
max_p = 1*10e-3; %m
Vl_x = 0 * 10e-9;  
Vl_y = 0 * 10e-9;
disp_x = zeros(1,Nj);
disp_y = zeros(1,Nj);
xt_ceil = zeros(1,Nj);
yt_ceil = zeros(1,Nj);
xt_floor = zeros(1,Nj);
yt_floor = zeros(1,Nj);

for z = 1:Nj    % Number of iterations
   Pos_x(1,z)
   Pos_y(1,z)
   
   Ex_arrayy = [E(Pos_x(1,z)+1,Pos_y(1,z)+1) E(Pos_x(1,z)+1,Pos_y(1,z)) E(Pos_x(1,z)+1,Pos_y(1,z)-1) E(Pos_x(1,z),Pos_y(1,z)+1)...
       E(Pos_x(1,z),Pos_y(1,z)-1) E(Pos_x(1,z)-1,Pos_y(1,z)+1) E(Pos_x(1,z)-1,Pos_y(1,z)) E(Pos_x(1,z)-1,Pos_y(1,z)-1)] 
   [Ext(1,z), ind_ex] = min(Ex_arrayy);
   if ind_ex == 1
       Pos_x(1, z+1) = Pos_x(1,z)+1;
       Pos_y(1, z+1) = Pos_y(1,z)+1;
   elseif ind_ex == 2
       Pos_x(1, z+1) = Pos_x(1,z)+1;
       Pos_y(1, z+1) = Pos_y(1,z);
    elseif ind_ex == 3
       Pos_x(1, z+1) = Pos_x(1,z)+1;
       Pos_y(1, z+1) = Pos_y(1,z)-1;
   elseif ind_ex ==4
       Pos_x(1, z+1) = Pos_x(1,z);
       Pos_y(1, z+1) = Pos_y(1,z)+1;
   elseif ind_ex ==5
       Pos_x(1, z+1) = Pos_x(1,z);
       Pos_y(1, z+1) = Pos_y(1,z)-1;
   elseif ind_ex ==6
       Pos_x(1, z+1) = Pos_x(1,z)-1;
       Pos_y(1, z+1) = Pos_y(1,z)+1;
   elseif ind_ex ==7
       Pos_x(1, z+1) = Pos_x(1,z)-1;
       Pos_y(1, z+1) = Pos_y(1,z);
   else
       Pos_x(1, z+1) = Pos_x(1,z)-1;
       Pos_y(1, z+1) = Pos_y(1,z)-1;
   end 
   
   
%       Ey_arrayy = [Ey(Pos_x(1,z)+1,Pos_y(1,z)+1) Ey(Pos_x(1,z)+1,Pos_y(1,z)) Ey(Pos_x(1,z)+1,Pos_y(1,z)-1) Ey(Pos_x(1,z),Pos_y(1,z)+1)...
%        Ey(Pos_x(1,z),Pos_y(1,z)-1) Ey(Pos_x(1,z)-1,Pos_y(1,z)+1) Ey(Pos_x(1,z)-1,Pos_y(1,z)) Ey(Pos_x(1,z)-1,Pos_y(1,z)-1)]; 
%    [Eyt(1,z), ind_ex] = min(Ey_arrayy);
%    if ind_ex == 1
% %        Pos_x(1, z+1) = Pos_x(1,z)+1;
%        Pos_y(1, z+1) = Pos_y(1,z)+1;
%    elseif ind_ex == 2
% %        Pos_x(1, z+1) = Pos_x(1,z)+1;
%        Pos_y(1, z+1) = Pos_y(1,z);
%     elseif ind_ex == 3
% %        Pos_x(1, z+1) = Pos_x(1,z)+1;
%        Pos_y(1, z+1) = Pos_y(1,z)-1;
%    elseif ind_ex ==4
% %        Pos_x(1, z+1) = Pos_x(1,z);
%        Pos_y(1, z+1) = Pos_y(1,z)+1;
%    elseif ind_ex ==5
% %        Pos_x(1, z+1) = Pos_x(1,z);
%        Pos_y(1, z+1) = Pos_y(1,z)-1;
%    elseif ind_ex ==6
% %        Pos_x(1, z+1) = Pos_x(1,z)-1;
%        Pos_y(1, z+1) = Pos_y(1,z)+1;
%    elseif ind_ex ==7
% %        Pos_x(1, z+1) = Pos_x(1,z)-1;
%        Pos_y(1, z+1) = Pos_y(1,z);
%    else
% %        Pos_x(1, z+1) = Pos_x(1,z)-1;
%        Pos_y(1, z+1) = Pos_y(1,z)-1;
%    end 
%        
%    Fx = Ext(1,z) * qi * e_charge;
%    
%    Ey_array = [Ey(xt_floor(1,z),yt_floor(1,z)) Ey(xt_ceil(1,z),yt_floor(1,z)) Ey(xt_floor(1,z),yt_ceil(1,z)) Ey(xt_ceil(1,z),yt_ceil(1,z))];
%    Eyt(1,z) = min(Ey_array);
%    
% %    Eyt(1,z) = Ey(xt_floor(1,z),yt_floor(1,z))+((Ey(xt_ceil(1,z),yt_ceil(1,z))-Ey(xt_floor(1,z),yt_floor(1,z)))/size)*(Pos_y(1,z)-yt_floor(1,z));
% %    eeyy = Eyt(1,z);
% %       if V(xt_floor(1,z),yt_floor(1,z)) > 0
% %        Eyt(1,z) = -Eyt(1,z);
% %       end
%    Fy = Eyt(1,z) * qi * e_charge;
%     
%     acc_x(1,z) = (Fx./m_Xe); %m/s^2
%     acc_y(1,z) = (Fy./m_Xe); %m/s^2
%     vel_x(1,z) = Vl_x + (acc_x(1,z) * del_t);
%     vel_y(1,z) = Vl_y + (acc_x(1,z) * del_t);
%     disp_x(1,z) = (vel_x(1,z+1) + (acc_x(1,z)*(del_t)))*del_t;
%     disp_y(1,z) = (vel_y(1,z+1) + (acc_y(1,z)*(del_t)))*del_t; 
%   
%         if abs(sqrt((disp_x(1,z))^2 + (disp_y(1,z))^2)) < min_p
%             del_t = 1.02 * del_t;
%              t(z) = t(z)+del_t;
%              Pos_x(1,z+1) = Pos_x(1,z);
%              Pos_y(1,z+1) = Pos_y(1,z);
%         
%         elseif abs(sqrt((disp_x(1,z))^2 + (disp_y(1,z))^2)) > max_p
%             del_t = 0.98 * del_t;
%              t(z) = t(z)+del_t;
%              Pos_x(1,z+1) = Pos_x(1,z);
%              Pos_y(1,z+1) = Pos_y(1,z);
% 
%         else
%              Vl_x = vel_x(1,z);
%              Vl_y = vel_y(1,z);
%              t(z) = t(z)+del_t;
%              Pos_x(1,z+1) = Pos_x(1,z) +  disp_x(1,z);
%              psx = Pos_x(1,z+1);
%              Pos_y(1,z+1) = Pos_y(1,z) +  disp_y(1,z);
%              psy = Pos_y(1,z+1);
%         end
end
