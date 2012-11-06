clc;clear all; close all;

trafficObj = mmreader('../00012.avi'); %nactu video
nframes = get(trafficObj, 'NumberOfFrames'); %pocet snimku ve videu
se = strel('disk',3); % fill a gap
width = 2; %velikost znacici kostky
%filtracni masky
H = fspecial('gaussian', 3, 0.7);
Sx = fspecial('prewit');
Sy = Sx';
P = fspecial('average',12);
% rec = avifile('motion3.avi','Compression', 'FFDS', 'fps', get(trafficObj, 'FrameRate'  ));

try
    bcg= double(imread('bcg.bmp'));
    edg_bcg= double(imread('edg_bcg.bmp'));
catch Me
    [bcg, edg_bcg] = get_background(trafficObj,50, Sx);
    imwrite(bcg, 'bcg.bmp');
    imwrite(uint8(edg_bcg), 'edg_bcg.bmp');
end
bcg = double(bcg);
%fig = figure;
%h = waitbar(0, 'processing');

for i=1:nframes
    %waitbar(i/nframes, h);
    I = double(read(trafficObj, i));
    I = imfilter(I,H); % odstraneni sumu
    
    diff = abs(bcg-I) + abs(I-bcg); % rozdiln smiku od pozadi
    ss = sqrt(sum(diff.^2,3)); % slouceni barev podle kvadratickeho prumeru 
    
    Ax = imfilter(ss, Sx);
    Ay = imfilter(ss, Sy);
    edg = sqrt(Ax.*Ax + Ay.*Ay)-edg_bcg; % hranovy filtr
    
    sa = imadd(ss,edg)>120; %slouceni prahu s rozdilem 
    ed = imopen(sa, strel('disk',1)); % prahovani / filtr
    edp = imfilter(ed,P); %prumerovani okoli
    edp = imfilter(edp,P); %druhe kolo prumerovani
    bw = bwareaopen(edp, 1000); % odstraneni malych ploch
    %bw2 = imfill(bw, 'holes'); %vyplneni der
    cc = bwconncomp(bw); % Find connected components in binary image
    s = regionprops(cc, {'Centroid', 'Area'});
    
    subplot(1,2,1);
    title(sprintf('snimek c.%d', i))
    imshow(uint8(I))
    subplot(1,2,2);
    imshow(bw,[]);
    
    if ~isempty(s)
        title(sprintf(' detekovanych objektu: %d', cc.NumObjects));
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