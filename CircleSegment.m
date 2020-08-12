function [h,x3,y3] = CircleSegment(pt1,pt2)
%% Calculated values
L = pt2(2)-pt1(2); %segment line length
r = sqrt((L.^2)./(2.*(1-cosd(90)))); %Full radius of the circle

%% Normal circle
%ranges from angle 135 to 225. This can be changed as needed
theta = 135:1:225; %range in degrees
x = r.*cosd(theta);
y = r.*sind(theta);
x2 = round(x + pt1(1)-x(1)); %Integer x vals
y2 = round(y + pt2(2)-y(1)); %Integer y vals

%% Eliminate duplicate points
[y3,index] = unique(y2); %Eliminate duplicate y values
x3 = x2(index); %Match the x values to the y values

%% Plot: this can be removed if you want
%If so, just put x3 and y3 as the outputs of the function, instead of h
h = plot(x3,y3);

%% Axis equal: this can be deleted 
axis equal
end

