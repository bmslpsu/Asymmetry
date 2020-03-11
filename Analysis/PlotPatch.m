function [h1,h2] = PlotPatch(dataMean,dataSTD,dataX,nSTD,nFly,lineColor,faceColor,alphaVal,lineWidth)
%% PlotPatch: Plots mean/med & STD as patch
%   INPUT:
%       dataMean    : mean/med data
%       dataSTD     : std of data
%       dataX       : x-axis points
%       nSTD        : # of STD's to plot for patch
%       nFly        : # of animals for standard error (nFly=1 for STD)
%       lineColor  	: color of mean/med data
%       edgeColor   : patch edge color
%       faceColor   : pathc area color
%       alphaVal    : patch transparancy (0-1)
%       lineWidth  	: linewidth of main curve
%   OUTPUT:
%       h1          : patch handle
%       h2          : mean curve handle
%

dataMean    = dataMean(:);
dataSTD     = dataSTD(:);
dataX       = dataX(:);

uE = dataMean + nSTD*dataSTD./sqrt(nFly);
lE = dataMean - nSTD*dataSTD./sqrt(nFly);
yP = [lE;flipud(uE)];
xP = [dataX;flipud(dataX)];

hold on
h1 = patch(xP,yP,1,'facecolor',faceColor,'edgecolor',lineColor);
h2 = plot(dataX,dataMean,'-','Color',lineColor,'LineWidth',lineWidth);
h1.EdgeColor = 'none';
alpha(h1,alphaVal)
end