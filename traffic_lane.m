clc;
clear all; 
close all;

trafficObj = mmreader('../00025.avi'); %nactu video
bcg = get_background(trafficObj,50);
imwrite(bcg,'00025.png','png');

%%
% Load saved roadLane:
clear roadLane;
whos('-file','data/00012.mat');     % Information about file
roadLane = load('data/00012.mat');  % Load data from file


%%
% Testing of finding traffic Lane %
close all;
clear all;
clc;

IM_PATH = ['data/00012';'data/00013';'data/00023';...
            'data/00024';'data/00025';'data/00026'];

imgCount = size(IM_PATH);
for j = 1:imgCount(1)
    bcg = imread([IM_PATH(j,:) '.png']);
    roadLane = GetTrafficLane(bcg,false);
    save([IM_PATH(j,:) '.mat'],'roadLane');
    
    figure(2);
    subplot(3,4,j * 2 - 1);
    imagesc(roadLane.surfLeft(:,:,1)+roadLane.surfLeft(:,:,2)*2);
    colormap(bone);
    hold on;
    title(IM_PATH(j,:));
    for i=1:3
       plot(roadLane.left(:,1,i),roadLane.left(:,2,i),'LineWidth',3); 
    end

    subplot(3,4,j * 2);
    imagesc(roadLane.surfRight(:,:,1)+roadLane.surfRight(:,:,2)*2);
    colormap(bone);
    hold on;
    title(IM_PATH(j,:));
    for i=1:3
       plot(roadLane.right(:,1,i),roadLane.right(:,2,i),'LineWidth',3); 
    end

    figure(3);
    subplot(2,3,j);
    imshow(bcg);
    hold on;
    title(IM_PATH(j,:));
    for i=1:3
       plot(roadLane.left(:,1,i),roadLane.left(:,2,i),'LineWidth',3); 
       plot(roadLane.right(:,1,i),roadLane.right(:,2,i),'LineWidth',3); 
    end
end
