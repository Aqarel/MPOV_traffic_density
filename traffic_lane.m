clc;
clear all; 
close all;

trafficObj = mmreader('00026.avi'); %nactu video
nframes = get(trafficObj, 'NumberOfFrames'); %pocet snimku ve videu

bcg = get_background(trafficObj,50);

%%
close all;


summ = zeros(trafficObj.Height,trafficObj.Width);
se = strel('disk',5);
h = waitbar(0, 'Traffic lanes');
maxx = 200;
figure(1);
title('Traffic lanes');
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1]);
for i=1:maxx
    waitbar(i/maxx, h);
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


