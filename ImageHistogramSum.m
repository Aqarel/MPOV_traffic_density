function [H] = ImageHistogramSum(I)
% Cumulative histogram
% I - image
% H - vector - normalize cumulative histogram

H = cumsum(ImageHistogram(I));          % cumulative histogram
H = round(255.*H./max(H));              % Normalize histogram

end