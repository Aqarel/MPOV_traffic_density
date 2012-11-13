clc;clear all; close all;
%check for presence of Intel Integrated Performance Primitives (Intel IPP) library
if ippl
    iptsetpref('UseIPPL', true)
    disp('IPPL library loaded')
else
    disp('IPP not found')
end

trafficObj = mmreader('../00012.avi'); %nactu video
nframes = get(trafficObj, 'NumberOfFrames'); %pocet snimku ve videu
duration = get(trafficObj, 'Duration'); % delka videa
H = fspecial('log',3,0.6);
se = strel('square',5); % fill a gap
width = 2; %velikost znacici kostky

% rec = avifile('motion3.avi', 'Compression', 'i420', 'fps', get(trafficObj, 'FrameRate'  ));
disp('getting bacground image...');
try
    bcg= double(imread('bcg.bmp'));
    %edg_bcg= double(imread('edg_bcg.bmp'));
catch Me
    bcg = get_background(trafficObj,50);
    imwrite(uint8(bcg), 'bcg.bmp');
    %imwrite(uint8(edg_bcg), 'edg_bcg.bmp');
end
disp('separating traffic lines...')
trafficLane = GetTrafficLane(bcg,0);
L = trafficLane.surfLeft(:,:,1) + trafficLane.surfLeft(:,:,2);
R = trafficLane.surfRight(:,:,1) + trafficLane.surfRight(:,:,2);
LR = L+R;
fig = figure(1);
subplot(1,2,1);
h = waitbar(0, 'processing');
disp('counting cars...')
for i=1:nframes
    tic
    waitbar(i/nframes, h);
    I = double(read(trafficObj, i));
    %I = imfilter(I,G); % odstraneni sumu
    I = imadd(imfilter(I,H),I); % posilit hrany
    D = uint8(I./bcg); % rozdiln smiku od pozadi
    M = uint8(D<0.9); % binarni maska rozdilnych pixelu
    bcg = bcg + double((0.2*(1-M)+0.02*M).*D); %aktualizovat pozadi
    D = bgremove(I,bcg, 30);
%     bw = bwmorph(D,'erode',2); %erode to remove small moise
    bw = bwmorph(D,'close');
    %bw2 = imfilter(bw, P); %spojeni prumerovanim - pomale
    
    %bw = imopen(D, se);
    
   % bo = imopen(bw);
    cc = bwconncomp(bw.*LR); % separovat prvni jizdni pruh
    L1 = regionprops(cc , {'Centroid', 'Area','BoundingBox'});
    idx = [L1.Area] > 1500; %vyprat pouze plochy s velkou plochou

    subplot(1,2,1);
    imshow(uint8(I));
    title(sprintf('snimek c.%d', i));
    subplot(1,2,2);
    imshow(bw,[]);
    
    if ~isempty(L1)
        title(sprintf(' detekovanych objektu: %d', sum(idx)));
        centroids = cat(1,L1.Centroid);
        boxes = cat(1,L1.BoundingBox);
        hold on
        subplot(1,2,1);
        for r=find(idx==1) % draw bounding boxes
            rectangle('Position',boxes(r,:));
        end
        subplot(1,2,2);
        plot(centroids(idx,1), centroids(idx,2), 'b*');
        hold off
        toc
    end
%     rec = addframe(rec,  getframe(fig));
%     rec = addframe(rec, im2frame(uint8(bo), cmap));
    %rec = addframe(rec, uint8(ed));
%     imshow(ed);
end
close(h)
% rec = close(rec);
clear(rec)