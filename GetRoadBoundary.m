function [bound] = GetRoadBoundary(road)
% Function find and calculate boundary of two-lane road
% road - color background of traffic situation
% bound - return 3 lines - boundary of lane

imSize = size(road);
if size(imSize) ~= [1 3]
   error('imBcg is not color image.'); 
end

lines = zeros(2,2,3);
bound = lines;
y = zeros(2,1);
y(1) = 0.10 * imSize(1);     % 10 % of high image
y(2) = 0.55 * imSize(1);     % 55 % of high image 

for i=1:2                    % Get two poit for boundary of all road lane
    tmp = road(y(i),:);     % Get one row of image
    tmp = find(tmp == 1);
    % line 1, point 1
    lines(i,1,1) = min(tmp);
    lines(i,2,1) = y(i);
    % line 3, point 1
    lines(i,1,3) = max(tmp);
    lines(i,2,3) = y(i);
    % line 2, point 1
    lines(i,1,2) = lines(i,1,1) + round((lines(i,1,3) - lines(i,1,1))/2);
    lines(i,2,2) = y(i);
end

for i = 1:3                 % calculate boundary of all road lane
    k = (y(2) - y(1)) / (lines(2,1,i) - lines(1,1,i));    
    q = y(2) - k * (lines(2,1,i));

    % upper edge of the picture
    bound(1,1,i) = round(-q/k);
    bound(1,2,i) = 1;

    % lower edge of the picture
    if ((imSize(1) - q) / k) < 0   % if line have x = 0 and y > 0
        bound(2,1,i) = 1; 
        bound(2,2,i) = round(q);
    else                    % if line have y = high of image and x > 0 
        bound(2,1,i) =round((imSize(1) - q) / k);
        bound(2,2,i) = imSize(1);
    end
end


