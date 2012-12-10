function [ bcg ] = get_background( video, n )
% Function is trying to get background of the video scene by averaging
% video - video 
% n - count of averaged frames

nframes = get(video, 'NumberOfFrames');
bcg = zeros(video.Height, video.Width, 3);
h = waitbar(0, 'Prumerovani pozadi');
for i=1:round(nframes / n):nframes
    waitbar(i / nframes, h);
    bcg = bcg + double(read(video, i));
end
    
bcg = bcg/n;
close(h)

end
