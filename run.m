clc;clear all; close all;
%check for presence of Intel Integrated Performance Primitives (Intel IPP) library
if ippl
    iptsetpref('UseIPPL', true)
    disp('IPPL library loaded')
else
    disp('IPP not found')
end

trafficObj = mmreader('00012.avi'); %nactu video
nframes = get(trafficObj, 'NumberOfFrames'); %pocet snimku ve videu
duration = get(trafficObj, 'Duration'); % delka videa
H = fspecial('log',3,0.6);
se = strel('square',5); % fill a gap
width = 2; %velikost znacici kostky

% rec = avifile('motion3.avi', 'Compression', 'i420', 'fps', get(trafficObj, 'FrameRate'  ));
disp('getting bacground image...');
try
    bcg= double(imread('bcg.bmp'));
catch Me
    bcg = get_background(trafficObj,50);
    imwrite(uint8(bcg), 'bcg.bmp');
end

[MR,MC] = size(bcg);
tt=1;
DOLNI_PRAH = 750;
HORNI_PRAH = 80;
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

s_init = struct('R', Rk, 'Q', Q, 'P', P, 'x', x); % init struct
cars = struct('R', {}, 'Q', {}, 'P', {}, 'x', {}); % aktualne sledovane
counted_cars = cars; % ukonceno sledovani

disp('separating traffic lines...')
trafficLane = GetTrafficLane(bcg,0);
L = trafficLane.surfLeft(:,:,1) + trafficLane.surfLeft(:,:,2);
R = trafficLane.surfRight(:,:,1) + trafficLane.surfRight(:,:,2);
LR = L+R;

% generuje hranice, kde se vuz prestane sledovat
LRp = LR;
LRp(HORNI_PRAH:DOLNI_PRAH,:) = 0;

fig = figure(1);
subplot(1,2,1);
h = waitbar(0, 'processing');
disp('counting cars...')


for i=1:nframes
    tic
    waitbar(i/nframes, h, sprintf('EAT: %.2f minutes',(nframes-i)*tt/60));
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
        centroids = cat(1,L1.Centroid);
        centroids = round(centroids(idx, 1:2)); % delete noise
        b = [];
        for j=1:size(centroids,1)  % delete centroids in off area
            if LRp(centroids(j,2),centroids(j,1)) == 1
                b = [b j];
            end
        end
        title(sprintf(' sledovanych objektu: %d\n celkem vozidel: %d'...
                        , size(cars,2), COUNTED));
        hold on
        co = centroids(b,:);
%         centroids(b,:) = [];
        if size(co,1)>0 % mark centroids - not counted
            plot(co(1), co(2), 'rd'); 
        end 
        
        hold on
        boxes = cat(1,L1.BoundingBox);
        subplot(1,2,1);
        for r=find(idx==1) % draw bounding boxes
            rectangle('Position',boxes(r,:));
        end
        cc = centroids(:,1);
        cc(b) = [];
        cr = centroids(:,2);
        cr(b) = [];
        subplot(1,2,2);
        plot(cc, cr, 'b*');
        plot(centroids(b,1), centroids(b,2), 'bo');
        
        %kalman - predikce polohy vozu

        while size(cars,2)-sum(idx)+size(b,2)< 0% pridej vozidlo
            cars(size(cars,2)+1) = s_init;
            COUNTED = COUNTED+1;
        end
                
        to_remove = [];
        for j = 1:size(cars,2)
            if cars(j).x(1) == 0
                xp = xp_init;
                ck=1;
            else
                xp=A*cars(j).x(end,:)' + Bu;
                
                % prirazeni spravneho bodu
                % TODO: mame vic aut, nez centroidu -> hledat auto, ktere
                %       je v zakrytu
                dx = centroids(:,1) - cars(j).x(end,1); 
                dy = centroids(:,2) - cars(j).x(end,2);
                d = dx.^2 + dy.^2;
                ck = find(d == min(d)); % index hledaneho centroidu
                for bb = 1:size(co,1) % mame vozilo(a) mame mimo slodovane pruhy
                    if sum(centroids(ck,:) == co(bb,:)) == 2 % vozidlo se dostalo za prah
                        counted_cars(size(counted_cars,2)+1) = cars(j);
                        to_remove = [to_remove j];
                    end
                end
                
            end
            PP = A*cars(j).P*A' + cars(j).Q;
            Kk = PP*Hk'*1/(Hk*PP*Hk'+Rk);
            
            % predikce polohy
            if cars(j).x(1) == 0 % aby vektor x nezacinal nulama, #XXX: asi by slo napsat lip
                cars(j).x(1,:) = (xp + Kk*([cc(1),cr(1)]' - Hk*xp))'; % nove je vzdy na konci, proto zbyva posledni souradnice
                ck = 1; % jenom pro pripad, ze by to naslo najednou vic aut
                cc(1) = [];
                cr(1) = [];
            else
                cars(j).x(end+1,:) = (xp + Kk*([centroids(ck,1),centroids(ck,2)]' - Hk*xp))';
                
                ck2=find(centroids(ck,1) == cc);
                centroids(ck,:) = [];
                cc(ck2) = []; % odeber prirazene centroidy
                cr(ck2) = [];
            end %
            cars(j).P = (eye(4)-Kk*Hk)*PP;
            
        end % prediction loop
        
        if size(to_remove,2)>0
            % nelze odebirat primo ve smycce
            cars(to_remove) = [];
        end
        
        tt = toc;
        
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
