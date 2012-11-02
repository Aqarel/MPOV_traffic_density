clc;clear all; close all;

trafficObj = mmreader('00012.avi'); %nactu video
nframes = get(trafficObj, 'NumberOfFrames'); %pocet snimku ve videu
se = strel('disk',2);

width = 2; %velikost znacici kostky
bcg = get_background(trafficObj,50);
I = read(trafficObj, 1);
fig=figure;
rec = avifile('motion2.avi','Compression', 'FFDS', 'fps', get(trafficObj, 'FrameRate'  ));
%taggedCars = zeros([size(I,1) size(I,2) 3 nframes], class(I));
M = fspecial('laplacian', 0.2);

for i=2:nframes
    %imshow(I);
    title(sprintf('%3.2f%% done', i/nframes*100));
    I2 = I;
    I = rgb2gray(read(trafficObj, i));
    edg = imopen(abs(I-bcg),se)>20;
    imshow(edg);
    s = regionprops(edg, {'Centroid'});
    if ~isempty(s)
        centroids = cat(1,s.Centroid);
        hold on
        plot(centroids(:,1), centroids(:,2), 'b*');
        hold off
        rec = addframe(rec,  getframe(fig));
    end
end
close(rec)
clear all %vyprazdnit buffer


%frametare =  get(trafficObj, 'FrameRate'  );
%movie2avi(mov, 'motion.avi', 'compression', 'FFDS', 'fps', get(trafficObj, 'FrameRate'  ));
%frameRate = get(trafficObj,'FrameRate');
%implay(taggedCars,frameRate);

% hold on
% 
% figure(2) 
% imshow(regions)
