%% FUNCTION:    PatchSTD
function [h1,h2] = PlotPatch(dataMean,dataSTD,dataX,nSTD,nFly,lineColor,faceColor,alphaVal,lineWidth)
%---------------------------------------------------------------------------------------------------------------------------------
% PatchSTD: Plots patch for data mean
    % INPUTS:
        % dataMean: mean of data
        % dataSTD: std of data
      	% dataX: x-axis points
        % edgeColor: patch edge color
        % faceColor: pathc area color
	% OUTPUTS
        % h: curve & patch handles
%---------------------------------------------------------------------------------------------------------------------------------      
    uE = dataMean + nSTD*dataSTD./sqrt(nFly);
    lE = dataMean - nSTD*dataSTD./sqrt(nFly);
    yP = [lE;flipud(uE)];
    xP = [dataX;flipud(dataX)];
    h1 = patch(xP,yP,1,'facecolor',faceColor,'edgecolor',lineColor);
    h2 = plot(dataX,dataMean,'Color',lineColor,'LineWidth',lineWidth);
    alpha(h1,alphaVal)
%---------------------------------------------------------------------------------------------------------------------------------
end