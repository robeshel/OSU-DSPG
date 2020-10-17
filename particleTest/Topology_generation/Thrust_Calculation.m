% Thrust Calculation
%function Thrust = Thrust_Calculation(Ib, Vb)
% 
clc
clear all
load('Test') %Variable number for the test iteration'
% Test = 17;
m_dot = 0.000000000199358;
date_yes = char(datetime(2020,10,09));
SG_data = round(table2array(readtable('Input_Data\Input_Data_new.xlsx', 'Sheet', 1, 'Range','C1:G3')));
AG_data = round(table2array(readtable('Input_Data\Input_Data_new.xlsx', 'Sheet', 2, 'Range','C1:G3')));
pt1 = [SG_data(1,2),SG_data(1,3)+1];
pt2 = [SG_data(1,2),SG_data(2,1)-1];

[h,xn,yn] = CircleSegment(pt1,pt2);
size_mat = size(yn, 2);
thrust_index = AG_data(1,4);


folder_name = ['Test_Data\' date_yes];
Vb = csvread([folder_name '\Test' num2str(Test) '_newT_Vx.csv']);
T_s = csvread([folder_name '\Test' num2str(Test) '_newT_TimeStep.csv']);
NPos_x = csvread([folder_name '\Test' num2str(Test) '_newT_NPos_x.csv']);  % writes the generated Trajectory matrix to a given name
NPos_x = reshape(NPos_x,[size(NPos_x,1)/size_mat, size_mat])';
Vb = (reshape(Vb,[size(Vb,1)/size_mat, size_mat]))';
T_s = (reshape(T_s,[size(T_s,1)/size_mat, size_mat]))';
sum_t = sum(T_s,2);
thr = 0;
for i = 1:size(Vb)
    row_vb = Vb(i,:);
    Vb_cell = row_vb(row_vb~=0);
    result = find(NPos_x(i,:)>thrust_index);% & NPos_x(i,:)<(thrust_index+1));
    if isempty(result)
        continue
        % thr = thr + (100 * Vb_cell(1,result) * 2.1801714e-25 * max(sum_t));
    else
        result = result(1,1);
        thr = thr + Vb_cell(1,result);
    end
end

vel_avg = thr ./  size_mat;

thrust_total = vel_avg * m_dot;
