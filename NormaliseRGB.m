function [ R ] = NormaliseRGB( I )
%NORMALISERGB tak touto metodou cesta nevede, ale mozna to bude uzitecne
s = sum(I,3);

R(:,:,1) = I(:,:,1)./s;
R(:,:,2) = I(:,:,2)./s;
R(:,:,3) = I(:,:,3)./s;

end

