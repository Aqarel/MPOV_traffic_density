function [ R ] = get_background( Obj, n )
%GET_BACKGROUND se snazi prumerovanim obrazu dostat pozadi sceny

nframes = get(Obj, 'NumberOfFrames');
R = zeros(Obj.Height,Obj.Width);
%S(:,:,1) = double(rgb2gray(read(Obj, 1)));
h = waitbar(0, 'Prumerovani pozadi');
for i=1:round(nframes/n):nframes
    waitbar(i/nframes, h);
    R = R + double(rgb2gray(read(Obj,i)));
end
    
R = uint8(R/n);
close(h)

end

