function [bound] = GetRoadBoundary(road)
lines = zeros(2,2,3);
bound = lines;
y1 = 100;
y2 = 600;
tmp = road(y1,:);           % Get one row of image
tmp = find(tmp == 1);
% line 1, point 1
lines(1,1,1) = min(tmp);
lines(1,2,1) = y1;
% line 3, point 1
lines(1,1,3) = max(tmp);
lines(1,2,3) = y1;
% line 2, point 1
lines(1,1,2) = lines(1,1,1) + round((lines(1,1,3) - lines(1,1,1))/2);
lines(1,2,2) = y1;

tmp = road(y2,:);           % Get one row of image
tmp = find(tmp == 1);
%line 1, point 2
lines(2,1,1) = min(tmp);
lines(2,2,1) = y2;
% line 3, point 2
lines(2,1,3) = max(tmp);
lines(2,2,3) = y2;
% line 2, point 2
lines(2,1,2) = lines(2,1,1) + round((lines(2,1,3) - lines(2,1,1))/2);
lines(2,2,2) = y2;

for i = 1:3                 % calculate boundary of all road
    k = (y2 - y1)/(lines(2,1,i) - lines(1,1,i));    
    q = y2 - k*(lines(2,1,i));

    bound(1,1,i) = round(-q/k);
    bound(1,2,i) = 1;

    if ((1080 - q)/k) < 0   % if line have y = 0 and x < 0
        bound(2,1,i) = 1; 
        bound(2,2,i) = round(q);
    else
        bound(2,1,i) =round((1080 - q)/k);
        bound(2,2,i) = 1080;
    end
end


