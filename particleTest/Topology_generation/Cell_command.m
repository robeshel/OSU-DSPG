clc
clear all
fprintf('plot(');
SG_data = round(table2array(readtable('Input_Data\Input_Data_new.xlsx', 'Sheet', 1, 'Range','C1:G3')));
pt1 = [SG_data(1,2),SG_data(1,3)+1];
pt2 = [SG_data(1,2),SG_data(2,1)-1];

[h,xn,yn] = CircleSegment(pt1,pt2);
size_mat = size(yn, 2);

for i = 1:size_mat
    fprintf('pos_cellX{%d}, pos_cellY{%d}, ', i,i)
end

fprintf(')');