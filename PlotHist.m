function PlotHist(H)
%Vykresli histogram
for j = 0:255
   line([j j], [0 H(j+1)]);
end