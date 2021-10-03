% tradeOff.m - Andrew J. Cross , 2014
%
% A function to find the 'knee' of a trade-off curve. A trade-off curve is used. Add this code as a .m file in your file path.

% to visualise two populations with different magnitudes of a given
% variable. Trade-off curves are ideally shaped like an exponential function 
% - the 'knee' of the curve (the point furthest from a line connecting the 
% first and last points of the curve) provides the critical threshold value
% of the given variable, to separate the two populations. 
%
%
%         ^  
%         |                        o
%         |                        o
%         |                        o
%     y   |                       o 
%         |                      x  <== knee in the curve
%         |                   o
%         |    o    o   o
%         |     
%         +------------------------------>
%              number of measurements
%
% 
% Here, the knee in the curve is found by drawing a line between the first
% and last points in the curve; the point on the curve which is furthest
% from the line is the knee. 


function knee = tradeOff(y) 

y = sort(y)'; % Variable used to separate two populations
x = 1:length(y); % Number of measurements

% Find 'the line' passing through the first and last points 
% of the trade-off curve
line = polyfit([x(1) x(end)],[y(1) y(end)],1);

% Find distance between each point and the curve
dist = abs((line(1).*x)+(-1.*y)+line(2))/...
        sqrt(line(1)^2+(-1)^2);

% The 'knee' in the trade-off curve is at the point furthest from 'the line'
[c id] = max(dist);

% Find y value at this point
knee = y(id);
    
% Plot trade-off curve and the knee in the curve
figure, plot(x,y,'linewidth',3)
hold on
scatter(x(id),knee,100,'or','linewidth',2)
plot([0 x(id)],[knee knee],'--k','linewidth',2)

end
