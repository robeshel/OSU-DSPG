function [V_m] = voltageCalcFunction(Ni,Nx,Ny,D1,D2,V_m)
%This function calculates the voltages of the grid for the thruster
tic
for f = 1:Ni    % Number of iterations
    for i=2:Nx-1
        for j=2:Ny-1
%-------------------------------------------------------------------------%
                V_m(i,j)=0.25*(V_m(i+1,j)+V_m(i-1,j)+V_m(i,j+1)+V_m(i,j-1)); % approximates the potential over the space
%-------------------------------------------------------------------------%
            for r = 1:2
                V_m(D1(r,4):D1(r,2), D1(r,1):D1(r,3)) = D1(r,5); % defines the potential on screen grid
                V_m(D2(r,4):D2(r,2), D2(r,1):D2(r,3)) = D2(r,5); % defines the potential on accelerator grid
            end
               V_m(1:10, 1: D1(1,1)) = D1(1,5);
               V_m(end-10:end, D1(1,1)) =  D1(1,5);
               V_m(:, 1:10) =  D1(1,5);
%-------------------------------------------------------------------------%
% Edges potentials calculations
%-------------------------------------------------------------------------%
                
                V_m(1,j) = (V_m(1,j-1) + V_m(2,j) + V_m(1,j+1))./3 ;
                V_m(Nx,j) = (V_m(Nx,j-1) + V_m(Nx-1,j) + V_m(Nx,j+1))./3;
                V_m(i,1) = (V_m(i-1,1) + V_m(i,2) + V_m(i+1,1))./3;
                V_m(i,Ny) = (V_m(i-1,Ny) + V_m(i, Ny-1) + V_m(i+1, Ny))./3;
                
%-------------------------------------------------------------------------%
%Corner potentials calculations
%-------------------------------------------------------------------------%
                
                V_m(1,1) = 0.5*(V_m(1,2)+V_m(2,1));
                V_m(Nx,1) = 0.5*(V_m(Nx-1,1)+V_m(Nx,2));
                V_m(1,Ny) = 0.5*(V_m(1,Ny-1)+V_m(2,Ny));
                V_m(Nx,Ny) = 0.5*(V_m(Nx,Ny-1)+V_m(Nx-1,Ny));
                
%-------------------------------------------------------------------------%
        end
    end
end
toc
end

