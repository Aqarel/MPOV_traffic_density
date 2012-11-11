function [ R ] = bgremove( I, bcg,d)
%   Detailed explanation goes here

[MR,MC,Dim] = size(I);

  % subtract background & select pixels with a big difference
  R = zeros(MR,MC);
  R = (abs(I(:,:,1 )- bcg(:,:,1)) > d) ...
    | (abs(I(:,:,2) - bcg(:,:,2)) > d) ...
    | (abs(I(:,:,3) - bcg(:,:,3)) > d);

    %eval(['imwrite(uint8(fore),''BGONE/nobg',int2str(index),'.jpg'',''jpg'')']);  
   
end

