function [H] = ImageHistogramSum(I)
%Kumulovaný histogram obrazku
H = cumsum(ImageHistogram(I));
%H = round((H.*255)./(max(H)));
end