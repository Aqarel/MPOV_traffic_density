function [roadLane] = GetTrafficLane(bcg,show)
% Function find traffic lanes in image
% bcg - color picture of background
% show - show histogram, original and segmented road 
% roadLane - structure:
%           .left(x,y,line) - boundary of left two-lane road, two points for every line 
%           .right(x,y,line) - boundary of right two-lane road, two points for every line
%           .surfLeft(:,:,1) - left two-lane road - mask of left lane
%           .surfLeft(:,:,2) - left two-lane road - mask of right lane
%           .surfRight(:,:,1) - right two-lane road - mask of left lane
%           .surfRight(:,:,2) - right two-lane road - mask of right lane

BTH = 1.15;                                     % Bulgarian constant for bigger threshold
                                                % , from analyze of an histogram and experience   
imSize = size(bcg);
if size(imSize) ~= [1 3]
   error('imBcg is not color image.'); 
end

roadLane = struct('left',[],'right',[],'surfLeft',[],'surfRight',[]);
bcgB = bcg(:,:,3);                              % Get blue composit
hist = ImageHistogramSum(bcgB);                 % Cumulative histogram
hist = round(255.*hist./max(hist));             % Normalize histogram

thRoad1 = 0;                                    % smaller threshold 
thRoad2 = 0;                                    % bigger threshold
thRoadFind = false;                             % Do it find bigger threshold?
cntr = 0;                                       
% Finding thresholds: 
% Algorithm used is based on an analysis of the
% histogram for our composition.
% For 00026.avi is thresholds cca: thRoad1 = 80; thRoad2 = 130;
for i=255:-1:2                                  
   if thRoadFind == false                       % finding bigger threshold
       if (hist(i) - hist(i - 1)) >= 2          % if first derivative >= 2 
           cntr = cntr + 1;
           if cntr == 3                         % 3x - first derivative >= 2 
               thRoad2 = round(i * BTH);        % Get bigger threshold for lines on road
               thRoadFind = true;
               cntr = 0;
           end
       else
           cntr = 0;
       end 
   else                                         % finding smaller threshold
       if (hist(i) - hist(i - 1)) <= 1          % if first derivative <= 1 
           cntr = cntr + 1;
           if cntr == 3                         % 3x - first derivative <= 1 
                thRoad1 = i;                    % Get smaller threshold for lines on road
                break;
           end
       else
           cntr = 0;
       end  
   end
end

if (thRoad1 == 0) || (thRoad2 == 0)                 % Find any threshold?
    error('Threshold not found.');
end

roadSeg = (bcgB > thRoad1) == (bcgB < thRoad2);     % Segmentation of road
roadFilter = strel('disk',10);                      
roadSeg = imopen(roadSeg, roadFilter);              % filter small object

regions = regionprops(roadSeg, 'Area');             % properties of picture
surfaces = cat(1, regions.Area);                    % size of surfaces
surfacesSort = sort(surfaces,'descend');            % sorting

idRoad1 = find(surfacesSort(1) == surfaces);        % find biggest surface
idRoad2 = find(surfacesSort(2) == surfaces);        % find second biggest surface

regions = regionprops(roadSeg, 'Centroid');         % properties of picture
centroid = cat(1, regions.Centroid);                % centroid of surface

if centroid(idRoad1,1) > centroid(idRoad2,1)        % Is finded left road really left road?
    tmp = idRoad1;                                  % Correct order of roads
    idRoad1 = idRoad2;
    idRoad2 = tmp;
end

idIm = bwlabel(roadSeg);                            % labeling of surface in image

surfLeft = idIm == idRoad1;                         % get left road
surfRight = idIm == idRoad2;                        % get right road

roadLane.left = GetRoadBoundary(surfLeft);          % get boundary of left road
roadLane.right = GetRoadBoundary(surfRight);        % get boundary of right road

% Only generated mask for every road lanes
imSize = size(bcgB);
roadLane.surfLeft(:,:,1) = poly2mask([roadLane.left(1,1,1) roadLane.left(1,1,2)...
                               roadLane.left(2,1,2) 0 roadLane.left(2,1,1)],...
                               [roadLane.left(1,2,1) roadLane.left(1,2,2)... 
                               roadLane.left(2,2,2) imSize(1) roadLane.left(2,2,1)],...
                               imSize(1),imSize(2));
roadLane.surfLeft(:,:,2) = poly2mask([roadLane.left(1,1,2) roadLane.left(1,1,3)...
                               roadLane.left(2,1,3) roadLane.left(2,1,2)],...
                               [roadLane.left(1,2,2) roadLane.left(1,2,3)... 
                               roadLane.left(2,2,3) roadLane.left(2,2,2)],...
                               imSize(1),imSize(2));
roadLane.surfRight(:,:,1) = poly2mask([roadLane.right(1,1,1) roadLane.right(1,1,2)...
                               roadLane.right(2,1,2) roadLane.right(2,1,1)],...
                               [roadLane.right(1,2,1) roadLane.right(1,2,2)... 
                               roadLane.right(2,2,2) roadLane.right(2,2,1)],...
                               imSize(1),imSize(2));
roadLane.surfRight(:,:,2) = poly2mask([roadLane.right(1,1,2) roadLane.right(1,1,3)...
                               roadLane.right(2,1,3) roadLane.right(2,1,2)],...
                               [roadLane.right(1,2,2) roadLane.right(1,2,3)... 
                               roadLane.right(2,2,3) roadLane.right(2,2,2)],...
                               imSize(1),imSize(2));

% Show my masterpieces %
if show == true
    figure();
    subplot(2,2,1:2);
    PlotHist(hist);
    hold on;
    plot([thRoad1 thRoad1],[0 255],'r');
    plot([thRoad2 thRoad2],[0 255],'r');
    subplot(2,2,3);
    imshow(bcgB,[]);
    subplot(2,2,4);
    imshow(roadSeg,[]);
end;
end