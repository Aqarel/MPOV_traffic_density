clc;
clear all; 
close all;

trafficObj = mmreader('../00026.avi'); %nactu video
bcg = get_background(trafficObj,50);
%%
close all;

bcgB = bcg(:,:,3);  % Get blue composit
%G = fspecial('gaussian',[60 60],50);
%bcgB = imfilter(bcgB,G,'same');
hist = ImageHistogramSum(bcgB);
hist = round(255.*hist./max(hist)); % Normalize histogram

thRoad1 = 0;
thRoad2 = 0;
thRoadFind = false;
for i=255:-1:2
   if thRoadFind == false
       if (hist(i) - hist(i - 1)) >= 3
           thRoad2 = round(i*1.05);                  % Get threshold for lines on road
           thRoadFind = true;
       end
   else
       if (hist(i) - hist(i - 1)) <= 1
           thRoad1 = i;                         % Get threshold for lines on road
           break;
       end  
   end
end

if (thRoad1 == 0) || (thRoad2 == 0)
    error('Threshold not found.');
end



roadSeg = (bcgB > thRoad1) == (bcgB < thRoad2);
roadFilter = strel('disk',10);
roadSeg = imopen(roadSeg, roadFilter);
%roadSeg = roadSeg ~= 1;
%roadFilter = strel('disk',20);
%roadSeg = imopen(roadSeg, roadFilter);
roadLightLines = bcgB > thRoad2;
%roadLightLines = imopen(roadLightLines, strel('disk',3));
roadLines = roadSeg == roadLightLines;

% Show my masterpieces
figure(1);
subplot(2,3,1);
imshow(bcgB,[]);
subplot(2,3,2:3);
PlotHist(hist);
subplot(2,3,4);
imshow(roadSeg,[]);
subplot(2,3,5);
imshow(roadLightLines,[]);
subplot(2,3,6);
imshow(roadLines,[]);






%% Old version using blur
clc;
clear all; 
close all;

trafficObj = mmreader('../00026.avi'); %nactu video
nframes = get(trafficObj, 'NumberOfFrames'); %pocet snimku ve videu

bcg = get_background(trafficObj,50);

bcg = rgb2gray(bcg);
summ = zeros(trafficObj.Height,trafficObj.Width);
se = strel('disk',5);
h = waitbar(0, 'Traffic lanes');
maxx = 200;
figure(1);
title('Traffic lanes');
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1]);
for i=1:round(nframes/500):nframes
    waitbar(i/nframes, h);
    I = rgb2gray(read(trafficObj, i));
    temp = imopen(double((I - bcg)+(bcg - I)),se) > 20;
    summ = summ + temp;
    
    subplot(1,2,1);
    imshow(I,[]);
    subplot(1,2,2);
    imshow(summ,[]);
end
close(h);
figure(2);
mesh(summ);


