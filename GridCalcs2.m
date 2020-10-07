function [x,y] = GridCalcs2(apts,grids,gridThick,AptDiam,Gap,lengthX,lengthY,res, voltage)
%This function takes in specs for ion accelerator grids and apts
%And spits out the points where the apts and grids are located


%% X direction

m = 1; %Counter for the grids
for i = 1:grids*2-1
    x(1) = Gap.*res; %The first value is equal to the gap thickness
    if mod(i,2)==1 %If we just added a gap, then add a grid thickness
        x(i+1) = x(i)+(gridThick(m)).*res;
        m = m+1;
    else %If we just added a grid thickness, then add a gap thickness
        x(i+1) = x(i)+(Gap).*res;
    end
end

%% Y direction
%Calculate how much space is used by the apertures
%Use this to find the grid height
YaptsTotal = sum(AptDiam); %number of gaps times gap length
gridHeight = (lengthY-YaptsTotal)./(apts+1); %length of the X direction - the gap length / the number of grids

n = 1;
%Find the y coordinates
for j = 1:apts.*2-1
    
    %IS THERE A SPACE ABOVE THE FIRST GRID???? - No apertures at the top
    %and bottom
    
    y(1) = gridHeight.*res; %The first value is equal to the grid height
    if mod(j,2)==1 %If we just added an apt diameter, add a grid height
        y(j+1) = y(j)+ (gridHeight).*res;
    else %If we just added a grid height, add an apt diameter
        y(j+1) = y(j)+ (AptDiam(n)).*res;
        n= n + 1;
    end
end
%% output: Each grid has its own sheet
%This function also outputs a .xls file
%Each grid has its own sheet of coordinates
%Each row is the coordinates for one square
%The number of rows is the number of squares per grid
%The number of sheets is the number of grids
% x = flip(x);
% y = flip(y);

f = 1;
for l = 1:grids %amount of sheets
    sheet = l;
    
    
    %Titles for the columns
    names = {'X1','Y1','X2','Y2','voltage'};
    
    Y1 = y(1:2:end);                    %Column 2
    Y2 = y(2:2:end);                    %Column 4
    X1 = ones(1,length(Y1)).*x(f);      %Column 1
    X2 = ones(1,length(Y1)).*x(f+1);    %Column 3
    V = ones(1,length(Y1)).*voltage(l); %Column 5
    
    f = f+2;
    
    %Create a table of the values
    Tab = table(X1',Y1',X2',Y2',V','VariableNames',names);
    %Write the table to the .xls file
    writetable(Tab,'gridvals.xls','Sheet',sheet)
    
    k = length(Y1);
    %% Plot the points
    set(gca,'YDir','normal')
    for m = 1:length(Y1)
        
        plot(X1(k),Y1(k),'xb')
        hold on
        plot(X2(k),Y2(k),'xk')
        hold on
        
        k = k-1;
    end
    
    
end

end


