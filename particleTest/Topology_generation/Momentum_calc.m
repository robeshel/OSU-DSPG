clc
clear all
load('Test') %Variable number for the test iteration'
Test = Test - 2;
date_yes = char(datetime(2020,9,22));
folder_name = ['Test_Data\' date_yes];
NPos_x = csvread([folder_name '\Test' num2str(Test) '_RITEvo_NPos_x.csv']);  % writes the generated Trajectory matrix to a given name
NPos_x = reshape(NPos_x,[100000, 253])'; 
NPos_y = csvread([folder_name '\Test' num2str(Test) '_RITEvo_NPos_y.csv']);
NPos_y = reshape(NPos_y,[100000, 253])';
Vx = csvread([folder_name '\Test' num2str(Test) '_RITEvo_Vx.csv']);
Vx = reshape(Vx,[100000, 253])';


Vcell_x = {};
for s = 1:253
    row_mat = Vx(s,:);
    Vcell_x{s} = row_mat(row_mat~=0);
end
 
momentum = 0;


for s = 1:253
    row_mat = Vcell_x{s};
    momentum = momentum + row_mat(1,end) * 2.1801714e-25 * 100;
end

% V = csvread([folder_name '\Test' num2str(Test) '_VtgDistMat.csv']);

save('Momentum_RITEvo', "momentum")


