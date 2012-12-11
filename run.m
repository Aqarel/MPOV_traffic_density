clc;clear variables; close all;
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

% rec = avifile('motion3.avi', 'Compression', 'i420', 'fps', get(trafficObj, 'FrameRate'  ));
disp('getting bacground image...');
try
    bcg= double(imread('bcg.bmp'));
catch Me
    bcg = get_background(trafficObj,50);
    imwrite(uint8(bcg), 'bcg.bmp');
end

[MR,MC,z] = size(bcg);
tt=1;
DOLNI_PRAH = 760;
HORNI_PRAH = 80;
COUNTED = 0;

% Kalman filter initialization
Rk=[[0.0645,0.0045]',[0.0045,0.00445]'];
Hk=[[1,0]',[0,1]',[0,0]',[0,0]'];
Q=0.15*eye(4);
P = 14*eye(4);
dt=1;
A=[[1,0,0,0]',[0,1,0,0]',[dt,0,1,0]',[0,dt,0,1]'];
g = 6; % pixels^2/time step
Bu = [0,0,0,g]';
% kfinit=0;
x=zeros(1,4);
xp_init = [MC/2,MR/2,0,0]';
% strukt init

s_init = struct('R', Rk, 'Q', Q,  'P', P,  'x',  x, 'centroid', [0 0], 'bb', [0 0 0 0]); % init struct
cars   = struct('R', {}, 'Q', {}, 'P', {}, 'x', {}, 'centroid', {}, 'bb', {}); % aktualne sledovane
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
noCentroids = 0;

for i=440:nframes %250, 265, 450
    tic
    waitbar(i/nframes, h, sprintf('EAT: %.2f minutes',(nframes-i)*tt/60));
    I = double(read(trafficObj, i));
    I = imadd(imfilter(I,H),I); % posilit hrany
    D = uint8(I./bcg); % rozdiln smiku od pozadi
    M = uint8(D<0.99); % binarni maska rozdilnych pixelu
    bcg = bcg + double((0.15*(1-M)+0.03*M).*D); %aktualizovat pozadi
    D = bgremove(I,bcg, 30);
    bw = bwmorph(D,'close'); % erode to remove small moise
    bw = imfill(bw,'holes');

    ccbw = bwconncomp(bw.*LR); % separovat prvni jizdni pruh
    L1 = regionprops(ccbw , {'Centroid', 'Area','BoundingBox', 'FilledImage'});
    idx = [L1.Area] > 1900; % vyprat pouze plochy s velkou plochou
    
    subplot(1,2,1);
    imshow(uint8(I));
    title(sprintf('snimek c.%d', i));
    subplot(1,2,2);
    imshow(bw,[]);
    
    if ~isempty(L1)                 % object found
        centroids = cat(1,L1.Centroid);
        centroids = round(centroids(idx, 1:2)); 

        cOffArea = [];
        for j=1:size(centroids,1)  % delete centroids in off area
            if LRp(centroids(j,2),centroids(j,1)) == 1
                cOffArea = [cOffArea j]; %#ok<AGROW>
            end
        end
        
        hold on
        co = centroids(cOffArea,:);
        if size(co,1)>0 % mark centroids - not counted
            plot(co(1), co(2), 'rd'); % nekdy oznaci divnou chybu O_o
        end 
        line([0 1920],[HORNI_PRAH HORNI_PRAH],'color','g');
        line([0 1920],[DOLNI_PRAH DOLNI_PRAH],'color','g');
        
        hold on
        boxes = cat(1,L1.BoundingBox);
        subplot(1,2,1);
        for r=find(idx==1) % draw bounding boxes
            rectangle('Position',boxes(r,:));
        end
        cx = centroids(:,1);
        cx(cOffArea) = [];
        cy = centroids(:,2);
        cy(cOffArea) = [];
        subplot(1,2,2);
        plot(cx, cy, 'b*'); % centroidy ve sledovane oblasti
        plot(centroids(cOffArea,1), centroids(cOffArea,2), 'bo'); % centroidy mimo sedovanou oblast
        centroids(cOffArea,:) = [];
        %kalman - predikce polohy vozu

        while (sum(idx) - size(cOffArea,2) - size(cars,2) + noCentroids) > 0% pridej vozidlo, pocetAut - autMimoOblast - pocetSledAut + sledAutaBezCen
            cars(size(cars,2)+1) = s_init;
            COUNTED = COUNTED+1;
        end
        
        
        for j = 1:size(cars,2) % zresetuj posledni centroidy pro vsechna sledovana vozidla
            cars(j).centroid = 0;
        end
        
        for j = 1:size(centroids,1) %prirazeni nalezenych centroidu ke sledovanym vozidlum
            idx = 0;  % index hledaneho vozu
            min = 2300; % dosud nejmensi vzdalenost, vzdalenost z rohu do rohu
            for k = 1:size(cars,2)       
                if cars(k).centroid ~= 0 % toto vozidlo uz ma prideleny centroid
                    continue              % pokracuj na dalsi
                end
                if cars(k).x(1) == 0 && LRp(centroids(j,2), centroids(j,1)) == 1 
                    continue % nove sledovany vuz nedostane centroid mimo sledovanou oblast
                end
                
                % pokud pribyl novy objekt, pak by doslo k chybnemu vypoctu
                % vzdalenosti
                if (cars(k).x(end,1) == 0) && (cars(k).x(end,2) == 0)  
                    if min > 120     % Centroid je pøílš vzdálen od všech sledovaných aut  
                        idx = k;
                    end
                    cars(idx).centroid = centroids(j,:);
                    idx = 0;
                    break;
                end
                
                dx = centroids(j,1) - cars(k).x(end,1); 
                dy = centroids(j,2) - cars(k).x(end,2);
                d = sqrt(dx.^2 + dy.^2);
                if d<min
                    min = d;
                    idx = k;
                end
            end
            if idx % pro tento nectroid jsme nasli prislusne sledovane vozidlo
                if min <= 120
                    cars(idx).centroid = centroids(j,:);
                else % pokud se prida a zaroven odebere jedno sledovane auto, jinak by se nove auto priradilo ke staremu
                    if (centroids(j,2) < (HORNI_PRAH + 50) && centroids(j,2) > (DOLNI_PRAH - 50))
                        cars(end + 1) = s_init;
                        cars(end).centroid = centroids(j,:);
                        COUNTED = COUNTED+1;
                    else
                        cars(idx).centroid = centroids(j,:);
                        for q = 1:size(noCentroids,2)
                           if cars(end).x(1) == 0
                              cars(end) = []; 
                              COUNTED = COUNTED-1;
                           end
                        end
                    end
                end
            end
        end
        noCentroids = 0;
        to_remove = [];        
        for j = 1:size(cars,2)
            if cars(j).x(1) == 0
                xp = xp_init;
            else
                xp=A*cars(j).x(end,:)' + Bu; % kalman - 1. krok

                % pokud nebyl pridelen centroid
                if sum(cars(j).centroid == 0) % vozidlo se dostalo mimo sledovanou oblast, nebo do zakrytu jineho vozu 
                    match = 0;
                    for k = 1:size(co,1)
                        d = co(k,:) - cars(j).x(end,1:2); 
                        if sum(abs(d)<abs(3*cars(j).x(end,3:4))) == 2
                            match = 1;
                            break
                        end
                    end
                    if match % prida vuz na seznam pro odstraneni ze sledovani
                        to_remove = [to_remove j]; %#ok<AGROW>
                    else  % bude treba odhadovat pozici
                        if i == 312
                           xxxxxxx = 3; 
                        end
                        if size(cars(j).x(:,1),1) > 1 % interp1 potrebuje alespon 2 body
                            noCentroids = noCentroids + 1;
                            y = interp1(cars(j).x(:,1), cars(j).x(:,2), xp(1),'linear','extrap'); % zarucime, ze auto nebude menit smer
                             % Pouzijeme bod z predikce a z interpolace
                             % presuneme na primku 
                            cars(j).centroid = [xp(1) y];
                            if (xp(1) < 0) || (LRp(round(y), round(abs(xp(1)))) == 1) % v nasledujicim kroce bude vuz mimo oblast sledovani
                                to_remove = [to_remove j]; %#ok<AGROW>
                                noCentroids = noCentroids - 1;
                                if noCentroids < 0 
                                    noCentroids = 0;
                                end
                            end
                        end
                    end
                end  % nalezeni/odhad centroidu  
                
            end % endif prvni krok - odhad pozice
            PP = A*cars(j).P*A' + cars(j).Q; % aktualizace 
            Kk = PP*Hk'*1/(Hk*PP*Hk'+Rk); % aktualizace kalmanova zesileni
            
            % predikce polohy
            if cars(j).x(1) == 0 % aby vektor x nezacinal nulama, #XXX: asi by slo napsat lip
                cars(j).x(1,:) = (xp + Kk*(cars(j).centroid' - Hk*xp))'; % nove je vzdy na konci, proto zbyva posledni souradnice
            else
                cars(j).x(end+1,:) = (xp + Kk*(cars(j).centroid' - Hk*xp))'; % aktualizace polohy
            end %
            cars(j).P = (eye(4)-Kk*Hk)*PP;
            
        end % prediction loop
        
        for j = to_remove% smazat ze seznamu sledovanych
            counted_cars(size(counted_cars,2)+1) = cars(j); % ulozit pro naslednou analyzu
            cars(j) = []; 
        end
        
        tt = toc;
        
        %plot kalman prediction
        for j = 1:size(cars,2)
            plot(cars(j).x(:,1),cars(j).x(:,2), 'r-');
        end
        
        title(sprintf(' sledovanych objektu: %d\n celkem vozidel: %d'...
                        , size(cars,2), COUNTED));
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

