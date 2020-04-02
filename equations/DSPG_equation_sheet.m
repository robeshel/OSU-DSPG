%% Define Variables

% Ag = Active Grid Area
% alpha = Thrust correction factor for multiply charged ions
% da = diameter of acclelerator grid holes
% ds = diameter of screen holes
% E = Effectie electric field (Vt/Le)
% e = the charge of a sinlge electron
% Ea = Electric field at accel grid
% epsilon = Permittivity of free space
% Es = Electric field at screen grid
% Fi = Force on an ion
% Ft = Correction factor integral(2*pi*r*J(r)*cos(theta(r))dr / Ib
% gamma = Total thrust correction factor
% Ib = Current of ion bean
% Ii = Current of total ions
% Ip = singly charged ion current
% Ipp = doubly charged ion current
% J = Perveance
% J(r) = ion current density as a function of radius
% ld = Length from accelerator grid to Beam plasma potential surface
% le = length of langmuire sheath
% lg = length between screen and accelerator grid
% M = Mass of an ion
% mdot = mass flow rate
% mp = mass of propellant
% ni = ion number density in gap
% q = charge of an ion
% Q = Propellant particle flow rate
% R = ratio of net beam voltage to total voltage (Vb/Vt)
% rho = ion charge density in the gap
% sigma = screen grid surface charge density
% T = Thrust
% ta = Thickness of accelerator grid
% theta(r) = effective thrust vector angle as a function of radius
% Ts = Screen grid transparency
% ts = thickness of screen grid
% Va = Voltage accelerator grid
% Vb = net voltage through which ion is accelerated
% Vbp = Voltage Beam plasma potantial surface
% vi = velocity of ions
% Vp = discharge plasma sheath voltage
% Vs = Voltage screen grid
% Vt = Voltage total across grids (Vs+abs(Va))


% For Xenon 
    % sqrt(2*M/e) = 1.65e-3


%% Equations to define variables using other variables
le = sqrt((lg+ts)^2+ds^2/4)
Vt = Vs+abs(Va)
R = Vb/Vt
E = Vt/Le
alpha = (Ip+1/sqrt(2)*Ipp)/(Ip+Ipp)
gamma = alpha*Ft

%% Mass flow rate
mdot = Ib*mp/q

%% Ion beam current and voltage

Ib = 
%% Grid transparency
Ts = Ib/Ii

%% Sheath thickness
Le = sqrt((Lg+Ts)^2+(Ds^2/4))

%% Max Perveance
Jmax = 4*epsilon/9*sqrt(2*e/M)*Vt^(3/2)/Le^2

%% Ion velocity
vi = sqrt(2*e*Vb/M)

%% Thrust and force equations
Fs = 0.5*epsilon*Es^2
Fa = 0.5*epsilon*Ea^2
Fi = 0.5*epsilon*(Ea^2-Es^2)
T = -Fi
T = Fs+Fa
T = gamma*mdot*vi
T = sqrt(2*M/e)*Ib*sqrt(Vb) 
Tmax = 4/9*(epsilon*gamma*Ts/e)*sqrt(2*e/M)*(Vt^(3/2)/le^2)*M*sqrt(2*e*Vb/M)
T = 1.65*gamma*Ib*sqrt(Vb) %(mN)


%%Magnetron
%Heat of the anode of the magnetron. Raw temperature data.
T1= 1679.31+(-22252.17-1659.31)./(1+(x./1.881098e-39).^0.03126938);
