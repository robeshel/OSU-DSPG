function [N_Pos_x, N_Pos_y] = Post_process(Pos_x1, Pos_y1)
filter_mat = logical(Pos_x1) & logical(Pos_y1);
N_Pos_x = Pos_x1(filter_mat);
N_Pos_y = Pos_y1(filter_mat);
end

%%
%n c xc
 % xbc
 
 %%
 %d vjdklfjdjkf
 %cdxjknd