trafficObj = mmreader('../00012.avi'); %nactu video
nframes = get(trafficObj, 'NumberOfFrames'); %pocet snimku ve videu
duration = get(trafficObj, 'Duration'); % delka videa

trafficObj2 = trafficObj;
%figure(1)
k = zeros(240,320,3,nframes);
for i=1:nframes
    I = double(read(trafficObj, i));
    k(:,:,:,i) = imresize(I, [240 320], 'bilinear');
    
end

%%
% Example – to scale a given video to QVGA (320x240) resolution
fin = 'D:\Programing\Matlab\00012.avi';
fout = 'd:\out.avi';
fileinfo = aviinfo(fin);
nframes = fileinfo.NumFrames;
aviobj = avifile(fout, 'compression', 'none', 'fps',fileinfo.FramesPerSecond);
for i = 1:nframes
    %Read frames from input video
    mov_in = aviread(fin,i);
    im_in = frame2im(mov_in);
    %Do processing on each frame of the video
    %In this example - scaling
    [ht wd ch] = size(im_in);
    im_out = imresize(im_in, [240 320], 'bilinear');
    %Write frames to output video
    frm = im2frame(im_out);
    aviobj = addframe(aviobj,frm);
    i %Just to display frame number
end;
%Don't forget to close output file
aviobj = close(aviobj);
return;