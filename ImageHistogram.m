function [H] = ImageHistogram(I)
%Histogram obrazku
hist = zeros(1,256);
for j = 0:255
    hSize = size(find(I == j));
	hist(1,j+1) = hSize(1);
end
H = hist;