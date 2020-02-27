clc
close all;
clear all;
%-------------------------------------------------------------------------%
%                         INITIALIZATION    
%-------------------------------------------------------------------------%
% E = Total electric field matrix using Poisson's equation
% V = Potential matrix
% Nx = Number of grid points in X- direction
% Ny = Number of grid points in Y-Direction
%-------------------------------------------------------------------------%
% Enter the dimensions
qi= 1;
e_charge = 1.602e-19;
eps = 8.854*(10^-12);
res = 20;
m_Xe = 2.1801714e-25; %kg
Nx = 51;     % Number of X-grids
Ny = 51;     % Number of Y-grids
mpx = ceil(Nx/2); % Mid-point of x
mpy = ceil(Ny/2); % Mid point of y
Ni = 75;  % Number of iterations for the Poisson solver
V = zeros(Nx,Ny);   % Potential (Voltage) matrix
% T = 0;            % Top-wall potential
% B = 0;            % Bottom-wall potential
% L = 0;            % Left-wall potential
% R = 0;            % Right-wall potential
% %-------------------------------------------------------------------------%
% % Initializing edges potentials
% %-------------------------------------------------------------------------%
% V(1,:) = L;
% V(Nx,:) = R;
% V(:,1) = B;
% V(:,Ny) = T;
% %-------------------------------------------------------------------------%
% % Initializing Corner potentials
% %-------------------------------------------------------------------------%
% V(1,1) = 0.5*(V(1,2)+V(2,1));
% V(Nx,1) = 0.5*(V(Nx-1,1)+V(Nx,2));
% V(1,Ny) = 0.5*(V(1,Ny-1)+V(2,Ny));
% V(Nx,Ny) = 0.5*(V(Nx,Ny-1)+V(Nx-1,Ny));
%-------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
% Initializing Particle Properties
% -------------------------------------------------------------------------%
mat_size = 51; % particle matriz size
Position_mat = zeros(mat_size,mat_size);
Time_mat = zeros(mat_size,mat_size);%initializing Particle Matrix
% -------------------------------------------------------------------------%

length_plate_s = 15;   % Length of plate in terms of number of grids 
length_plate_a = 17;   % Length of plate in terms of number of grids  
lp_s = floor(length_plate_s/2);
lp_a = floor(length_plate_a/2);
position_plate_s = 13; %Position of plate on x axis
position_plate_a = 24;
pp_s = position_plate_s;
pp_a= position_plate_a;
t1 = 5;
t2 = 16;
g = 5;
st = 16;

for z = 1:Ni    % Number of iterations
        
        for i=2:Nx-1
        for j=2:Ny-1      
            
            % The next two lines are meant to force the matrix to hold the 
            % potential values for all iterations
            
                V(st:st+t1, 36:51) = 1200;
                V(st:st+t1, 1:15) = 1200;
                V(st+t1+g:st+t1+g+t2, 1:17) = -250;
                V(st+t1+g:st+t1+g+t2, 34:51) = -250;
                
                V(i,j)=0.25*(V(i+1,j)+V(i-1,j)+V(i,j+1)+V(i,j-1));
        end
        end
        
end
% Take transpose for proper x-y orientation
V = V';
[Ex,Ey]=gradient(V);
Ex = -Ex;
Ey = -Ey;

% Electric field Magnitude
E = (qi/(4*pi*eps))* (1./Ex.^2+Ey.^2);  
x = (1:Nx)-mpx;
y = (1:Ny)-mpy;

% Contour Display for electric potential
figure(1)
contour_range_V = -1201:0.5:1201;
contour(x,y,V,contour_range_V,'linewidth',0.05);
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',10);
xlabel('x-axis in meters','fontsize',10);
ylabel('y-axis in meters','fontsize',10);
title('Electric Potential distribution, V(x,y) in volts','fontsize',10);
h1=gca;
set(h1,'fontsize',10);
fh1 = figure(1); 
set(fh1, 'color', 'white')


% Quiver Display for electric field Lines
figure(2)
contour(x,y,E,'linewidth',0.5);
hold on, quiver(x,y,Ex,Ey,2)
title('Electric field Lines, E (x,y) in V/m','fontsize',14);
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',14);
xlabel('x-axis in meters','fontsize',14);
ylabel('y-axis in meters','fontsize',14);
h3=gca;
set(h3,'fontsize',14);
fh3 = figure(2); 
set(fh3, 'color', 'white')
acc = zeros(1,Ny); 
vel =  zeros(1,Ny);
pos = zeros(1,Ny);

for i= 2:Nx    % Number of iterations
    F =  qi * e_charge * E(25,i);
    acc(1,i) = (F./m_Xe) * 10e-18; %Um/ns^2
    vel(1, i) = vel(1,i-1) + (acc(1,i)); % Um/ns;
    pos(1,i) = pos(1,i-1) + vel(1,i);%Um;
end

figure(4)
plot(acc)
title('Acceleration','fontsize',14);
figure(5)
plot(vel)
title('Velocity','fontsize',14);
 % this is a new model

