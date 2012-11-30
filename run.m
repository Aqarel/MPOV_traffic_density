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

[MR,MC] = size(bcg);

DOLNI_PRAH = 740;
HORNI_PRAH = 110;
COUNTED = 0;

% Kalman filter initialization
Rk=[[0.2845,0.0045]',[0.0045,0.0455]'];
Hk=[[1,0]',[0,1]',[0,0]',[0,0]'];
Q=0.01*eye(4);
P = 100*eye(4);
dt=1;
A=[[1,0,0,0]',[0,1,0,0]',[dt,0,1,0]',[0,dt,0,1]'];
g = 6; % pixels^2/time step
Bu = [0,0,0,g]';
kfinit=0;
x=zeros(1,4);
xp_init = [MC/2,MR/2,0,0]';
% strukt init

s_init = struct('R', Rk, 'Q', Q, 'P',P, 'x', x); % init struct
cars = s_init; % structure array
counted_cars = s_init;

disp('separating traffic lines...')
trafficLane = GetTrafficLane(bcg,0);
L = trafficLane.surfLeft(:,:,1) + trafficLane.surfLeft(:,:,2);
R = trafficLane.surfRight(:,:,1) + trafficLane.surfRight(:,:,2);
LR = L+R;

% generuje hranice, kde se vuz prestane sledovat
Lp = L;
Lp(1:DOLNI_PRAH,:)=0;
Rp = R;
Rp(HORNI_PRAH:end, :)=0;
LRp = Lp+Rp;

fig = figure(1);
subplot(1,2,1);
h = waitbar(0, 'processing');
disp('counting cars...')


for i=1:nframes
    tic
    waitbar(i/nframes, h);
    I = double(read(trafficObj, i));
    I = imadd(imfilter(I,H),I); % posilit hrany
    D = uint8(I./bcg); % rozdiln smiku od pozadi
    M = uint8(D<0.9); % binarni maska rozdilnych pixelu
    bcg = bcg + double((0.2*(1-M)+0.02*M).*D); %aktualizovat pozadi
    D = bgremove(I,bcg, 30);
    bw = bwmorph(D,'close'); %erode to remove small moise

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
        centroids = round(centroids(idx, 1:2)); % delete noise
        b = [];
        for j=1:size(centroids,1)  % delete centroids in off area
            if LRp(centroids(j,2),centroids(j,1)) == 1
                b = [b j];
            end
        end
        hold on
        co = centroids(b,:);
        centroids(b,:) = [];
        if size(co,1)>0
            plot(co(1), co(2), 'ro'); % mark centroids - not counted
        end 
        
        hold on
        boxes = cat(1,L1.BoundingBox);
        subplot(1,2,1);
        for r=find(idx==1) % draw bounding boxes
            rectangle('Position',boxes(r,:));
        end
        cc = centroids(:,1);
        cr = centroids(:,2);
        subplot(1,2,2);
        plot(cc, cr, 'b*');
        
        %kalman - predikce polohy vozu

        while size(cars,2)-sum(idx)-size(counted_cars,2)+1< 0% add new detected car
                cars(size(cars,2)+1) = s_init;
        end
                
        for j = 1:size(centroids,2)
            if cars(j).x(1) == 0
                xp = xp_init;
            else
                xp=A*cars(j).x(end,:)' + Bu;
                dx = cc - cars(j).x(end,1); % prirazeni spravneho bodu
                dy = cr - cars(j).x(end,2);
                d = dx.^2 + dy.^2;
                ck = find(d == min(d));
            end
            PP = A*cars(j).P*A' + cars(j).Q;
            Kk = PP*Hk'*1/(Hk*PP*Hk'+Rk);
            if cars(j).x(1) == 0 % aby vektor x nezacinal nulama, #XXX: asi by slo napsat lip
                cars(j).x(1,:) = (xp + Kk*([cc(1),cr(1)]' - Hk*xp))';
                ck = 1; % jenom pro pripad, ze by to naslo najednou vic aut
            else
                cars(j).x(end+1,:) = (xp + Kk*([cc(ck),cr(ck)]' - Hk*xp))';
            end
            cars(j).P = (eye(4)-Kk*Hk)*PP;
            
            cc(ck) = []; % remove used elements
            cr(ck) = [];
            
            a = round(cars(j).x(end,1:2));
            
            %zapocist, jestli bude v nasledujicim kroce za prahem
            if LRp(a(2),a(1)) == 1 
                 if counted_cars(size(counted_cars,1)).x(1) == 0
                    counted_cars(1) = cars(j);
                 else
                   counted_cars(size(counted_cars,1)) = cars(j);
                 end
                 cars(j) = [];
            end
        end
        
        
        toc
        
        %plot kalman prediction
        for j = 1:size(cars,2)
            plot(cars(j).x(:,1),cars(j).x(:,2), 'r-');
        end
        hold off
        
    end
%     rec = addframe(rec,  getframe(fig));
%     rec = addframe(rec, im2frame(uint8(bo), cmap));
    %rec = addframe(rec, uint8(ed));
%     imshow(ed);
end
close(h)
% rec = close(rec);
clear(rec)