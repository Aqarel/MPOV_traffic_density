clc;clear all; close all;

trafficObj = mmreader('../00012.avi'); %nactu video
nframes = get(trafficObj, 'NumberOfFrames'); %pocet snimku ve videu
se = strel('disk',2);
width = 2; %velikost znacici kostky
H = fspecial('gaussian', 5, 0.7);
% rec = avifile('motion3.avi','Compression', 'FFDS', 'fps', get(trafficObj, 'FrameRate'  ));

try
    bcg= double(imread('bcg.bmp'));
catch Me
    bcg = get_background(trafficObj,50);
    imwrite(bcg, 'bcg.bmp');
end
bcg = double(bcg);
fig = figure;
%h = waitbar(0, 'processing');

for i=1:nframes
    %waitbar(i/nframes, h);
    I = double(read(trafficObj, i));
    diff = abs(bcg-I) + abs(I-bcg);
    ss = sqrt(sum(diff.^2,3));
    ed = imopen(ss>40, strel('disk',3)); %prahovani s filtraci
    bw = bwareaopen(ed, 1000); %vybereme pouze plochy, ktere maji vic, jak 1000px
    cc = bwconncomp(bw);
    s = regionprops(cc, {'Centroid', 'Area'});
    
    subplot(1,2,1);
    title(sprintf('snimek c.%d', i))
    imshow(uint8(I))
    subplot(1,2,2);
    imshow(bw,[]);
    
    if ~isempty(s)
        title(sprintf('%d detekovanych objektu', cc.NumObjects));
        centroids = cat(1,s.Centroid);
        hold on
        plot(centroids(:,1), centroids(:,2), 'b*');
        hold off
        %rec = addframe(rec,  getframe(fig));
    end
    
    %rec = addframe(rec, uint8(ed));
%     imshow(ed);
end
close(h)
% close(rec)
% clear(rec)