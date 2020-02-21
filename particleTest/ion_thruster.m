%2/21/20 - 9:01AM
%source: https://github.com/mccarts3/IonThruster/blob/master/IonThrusterVer3.m

clc
clear all
close all

%Constants
mol = 6.02*10^23; % A's Number atoms/mol
xenonAtomicMass = 131.29;   %g/mol
xenonIonMass = xenonAtomicMass / (mol*1000);    %kg/ions

%display xenon ion mass
fprintf('Xenon Ion Mass = %d\n', xenonIonMass);

%Note Ion if negatively charge have more mass than atoms

q = 1.602*10^-19;   %charge of electron || Assumed
speedOfLight = 3*10^8;

%values to change
rocketMass = 1999.8;    %kilograms
fuelMass = .2;  %kilograms
ionFlowRate = 1*10^10;  %(atom/ion)/sec  || Assumed
Voltage = 10000;

totalMass = rocketMass + fuelMass;

%End Constants
ionVelocity = sqrt(2*q*Voltage/xenonIonMass);

%displaying ion velocity
fprintf('Xenon Velocity = %d', ionVelocity);


%Array
timeArray = 1:100;
voltageArray = 1:100;
voltageArray = 1*10^4*voltageArray;     %starting from 10kV to 1MV

ThrustTimeXenon = 1:100;
velTimeXenon = 1:100;
ThrustVoltageXenon = 1:100;
accelTimeXenon = 1:100;
ThrustTimeXenondt = 1:100;


%Speed of Light
ThrustTimeXenonC = 1:100;
ThrustVoltageXenonC = 1:000;
accelTimeXenonC = 1:100;

k = 1;

for j = 1:15;   
   
    if j == 1 
        ionFlowRate = 1*10^9;
    elseif j ==2
        ionFlowRate = 1*10^10;
    elseif j == 3
        ionFlowRate = 1*10^11;
    elseif j == 4
        ionFlowRate = 1*10^12;
    elseif j == 5
        ionFlowRate = 1*10^13;
    elseif j == 6
        ionFlowRate = 1*10^14;
    elseif j == 7
        ionFlowRate = 1*10^15;
    elseif j == 8
        ionFlowRate = 1*10^16;
    elseif j == 9
        ionFlowRate = 1*10^17;
    elseif j == 10
        ionFlowRate = 1*10^18;
    elseif j == 11
        ionFlowRate = 1*10^19;
    elseif j == 12
        ionFlowRate = 1*10^20;
    elseif j == 13
        ionFlowRate = 1*10^21;
    elseif j == 14
        ionFlowRate = 1*10^22;
    else
        ionFlowRate = 1*10^23;
    end
    
for n = 1:100;
    if n == 1
        velTimeXenon(n) = sqrt(2*q*Voltage/xenonIonMass)*log((totalMass)/(totalMass - ionFlowRate*xenonIonMass*timeArray(n)));
        
    else
        velTimeXenon(n) = sqrt(2*q*Voltage/xenonIonMass)*log((totalMass)/(totalMass - ionFlowRate*xenonIonMass*timeArray(n))) + velTimeXenon(n-1);
        
    end
    
        ThrustVoltageXenon(n) = (totalMass - ionFlowRate*xenonIonMass*timeArray(2))*sqrt(2*q*voltageArray(n)/xenonIonMass)*log((totalMass - ionFlowRate*xenonIonMass*timeArray(1))/(totalMass - ionFlowRate*xenonIonMass*timeArray(2))); 
        
        ThrustTimeXenon(n) = (totalMass - ionFlowRate*xenonIonMass*timeArray(n))*sqrt(2*q*Voltage/xenonIonMass)*log((totalMass)/(totalMass - ionFlowRate*xenonIonMass*timeArray(n))); 
                
        ThrustVoltageXenonC(n) = (totalMass - ionFlowRate*xenonIonMass*timeArray(50))*speedOfLight*log((totalMass - ionFlowRate*xenonIonMass*timeArray(49))/(totalMass - ionFlowRate*xenonIonMass*timeArray(50))); 
        
    if n == 1
        ThrustTimeXenondt(n) = (totalMass - ionFlowRate*xenonIonMass*timeArray(n))*sqrt(2*q*Voltage/xenonIonMass)*log((totalMass)/(totalMass - ionFlowRate*xenonIonMass*timeArray(n))); 
        
        accelTimeXenon(n) = sqrt(2*q*Voltage/xenonIonMass)*log((totalMass)/(totalMass - ionFlowRate*xenonIonMass*timeArray(n))); 
        
        ThrustTimeXenonC(n) = (totalMass - ionFlowRate*xenonIonMass*timeArray(n))*speedOfLight*log((totalMass)/(totalMass - ionFlowRate*xenonIonMass*timeArray(n))); 
        
        accelTimeXenonC(n) = speedOfLight*log((totalMass)/(totalMass - ionFlowRate*xenonIonMass*timeArray(n))); 
    else
        ThrustTimeXenondt(n) = (totalMass - ionFlowRate*xenonIonMass*timeArray(n))*sqrt(2*q*Voltage/xenonIonMass)*log((totalMass - ionFlowRate*xenonIonMass*timeArray(n-1))/(totalMass - ionFlowRate*xenonIonMass*timeArray(n))); 

        accelTimeXenon(n) = sqrt(2*q*Voltage/xenonIonMass)*log((totalMass - ionFlowRate*xenonIonMass*timeArray(n-1))/(totalMass - ionFlowRate*xenonIonMass*timeArray(n))); 
       
        ThrustTimeXenonC(n) = (totalMass - ionFlowRate*xenonIonMass*timeArray(n))*speedOfLight*log((totalMass - ionFlowRate*xenonIonMass*timeArray(n-1))/(totalMass - ionFlowRate*xenonIonMass*timeArray(n))); 

        accelTimeXenonC(n) = speedOfLight*log((totalMass - ionFlowRate*xenonIonMass*timeArray(n-1))/(totalMass - ionFlowRate*xenonIonMass*timeArray(n))); 
    end
    
    if ionFlowRate*xenonIonMass*timeArray(n) > fuelMass
        ThrustTimeXenon(n) = 0;
        ThrustTimeXenondt(n) = 0;
        ThrustVoltageXenon(n) = 0;
        accelTimeXenon(n) = 0;
        ThrustTimeXenonC(n) = 0;
        ThrustVoltageXenonC(n) = 0;
        accelTimeXenonC(n) = 0;
    end
end

if (j-1)/2 == 1||(j-1)/2 == 2||(j-1)/2 == 3 || j == 1||(j-1)/2 == 4||(j-1)/2 == 5||(j-1)/2 == 6||(j-1)/2 == 7
    figure('units','normalized','outerposition',[0 0 1 1]) % Maximize plot window
    k = 1;
end

subplot(2,5, 5*k - 4);
plot(timeArray, ThrustTimeXenon);
title(sprintf('Thrust Vs Time|IonFlowRate = %.1e', ionFlowRate));
xlabel('Time(s)');
ylabel('Thrust(N)');
legend('Xenon');

subplot(2,5,5*k - 3);
plot(timeArray, ThrustTimeXenondt);
title(sprintf('Thrust Vs. Time|IonFlowRate = %.1e', ionFlowRate'));
xlabel('Time(s)');
ylabel('Thrust(N)/dt');
legend('Xenon');

subplot(2,5,5*k - 2);
plot(voltageArray, ThrustVoltageXenon);
title(sprintf('Thrust Vs. Voltage|Ion Flow Rate = %.1e', ionFlowRate'));
xlabel('Voltage(v)');
ylabel('Thrust(N)');
legend('Xenon');

subplot(2,5,5*k - 1);
plot(timeArray, ThrustTimeXenonC);
title(sprintf('Thrust Vs. Time|IonFlowRate = %.1e|10kV @c', ionFlowRate'));
xlabel('Time(s)');
ylabel('Thrust(N)/dt');
legend('Xenon');

subplot(2,5,5*k);
plot(timeArray, velTimeXenon);
title(sprintf('Thrust Vs. Time|IonFlowRate = %.1e', ionFlowRate'));
xlabel('Time(s)');
ylabel('Velocity (m/s^{2})');
legend('Xenon');

k = k + 1;
end