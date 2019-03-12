function [h1,h2] = PlotPatch(dataMean,dataSTD,dataX,nSTD,nFly,lineColor,faceColor,alphaVal,lineWidth)
%% PlotPatch: Plots mean/med & STD as patch
%   INPUTS:
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
%   OUTPUTS
%       h1          : patch handle
%       h2          : mean curve handle
%---------------------------------------------------------------------------------------------------------------------------------      
uE = dataMean + nSTD*dataSTD./sqrt(nFly);
lE = dataMean - nSTD*dataSTD./sqrt(nFly);
yP = [lE;flipud(uE)];
xP = [dataX;flipud(dataX)];
h1 = patch(xP,yP,1,'facecolor',faceColor,'edgecolor',lineColor);
h2 = plot(dataX,dataMean,'Color',lineColor,'LineWidth',lineWidth);
alpha(h1,alphaVal)
end