clc;clear all; close all;

trafficObj = mmreader('../00012.avi'); %nactu video
nframes = get(trafficObj, 'NumberOfFrames'); %pocet snimku ve videu
se = strel('disk',4); % fill a gap
width = 2; %velikost znacici kostky
%filtracni masky
G = fspecial('gaussian', 3, 0.7);
H = fspecial('log', 3, 0.7);
P = fspecial('average',5);
% rec = avifile('motion3.avi','Compression', 'FFDS', 'fps', get(trafficObj, 'FrameRate'  ));

try
    bcg= double(imread('bcg.bmp'));
    %edg_bcg= double(imread('edg_bcg.bmp'));
catch Me
    bcg = get_background(trafficObj,50);
    imwrite(uint8(bcg), 'bcg.bmp');
    %imwrite(uint8(edg_bcg), 'edg_bcg.bmp');
end

%fig = figure;
%h = waitbar(0, 'processing');

for i=1:nframes
    %waitbar(i/nframes, h);
    I = double(read(trafficObj, i));
    %I = imfilter(I,G); % odstraneni sumu
    
    D = I./bcg; % rozdiln smiku od pozadi
    M = D<1; % binarni maska rozdilnych pixelu
    bcg = bcg + (0.1*(1-M)+0.01*M).*D; %aktualizovat pozadi
    D = bgremove(I,bcg, 15);
   % bw = bwmorph(D,'erode',2); %erode to remove small moise
    bw = bwmorph(D,'close');
    bw2 = imfilter(bw, P); %spojeni prumerovanim - pomale
    bo = bwareaopen(bw2, 1000); % odstraneni malych ploch
    cc = bwconncomp(bo); % Find connected components in binary image
    s = regionprops(cc, {'Centroid'});
    
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