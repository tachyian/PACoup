function comodulogram( X, Y, Z, isContourPlot, rangeZ , labelsOn, tickOn,  colorbarOn )
%COMODULOGRAM Plots a two-dimensional color map
%   Inputs 
%       X,Y: coordinates along horizontal and vertical axes, respectively.
%       Use <meshgrid> function to generate these coordinates
%
%       Z: Coordinate to be represented with the color map along X and Y axes. 
%
%       isContourPlot: if set to 1 plots smooth contout plot
%
%       rangeZ: [minZ, maxZ] sets the lower (minY) and upper (maxZ) limits
%       of the Z range to be plotted.
%
%       labelsOn: Set to 1 to label the X and Y axes. Use a vector of the
%       form [xLabelOn, yLabelOn] where xLabelOn and yLabelOn are booleans
%       (0 or 1) that indicate whether to label X and Y axes, respectively.
%       
%       tickOn: Set to 0 to remove X and Y axes ticks. Use a vector of the
%       form [xTicksOn, yTicksOn] where xTicksOn and yTicksOn are booleans
%       (0 or 1) that indicate whether to include ticks in the X and Y axes
%       ,respectively.
%
%       colorbarOn: Set to 1 to display the colorbar on the right of the plot. 
%
%       Authors: David Escobar & Luke Johnson
%       

%%
fsize= 12;
if isContourPlot
    n_interpx_pac = 20; 
    n_interpy_pac = 20;
    stepx1 = diff(X(1,1:2))/n_interpx_pac;
    stepy1 = diff(Y(1:2,1))/n_interpy_pac;
    [plotx1, ploty1] = meshgrid(X(1,1) : stepx1 : X(1,end), Y(1,1): stepy1 : Y(end,1));
    plotz1 = interp2(X,Y,Z,plotx1,ploty1,'cubic');
    imagesc(plotx1(1,:),(ploty1(:,1)),(plotz1)); set(gca,'Ydir','normal');
else
    surf(X,Y,Z)
    view(2);
    axis tight;
end
xlabel(' ','fontsize',fsize);
%set(gca,'fontsize',14)
if length(labelsOn) == 1  && labelsOn ==1
    ylabel('Frequency for amplitude (Hz)' , 'fontsize',fsize)
    xlabel('Frequency for phase (Hz)' , 'fontsize',fsize)
    %colorlabel('Modulation index' , 'fontsize',fsize)
elseif length(labelsOn) == 2 && labelsOn(1)==1  && labelsOn(2)==0 
    xlabel('Frequency for phase (Hz)' , 'fontsize',fsize)
elseif length(labelsOn) == 2 && labelsOn(1)==0  && labelsOn(2)==1
    ylabel('Frequency for amplitude (Hz)' , 'fontsize',fsize) 
end


if length(tickOn) == 1  && tickOn ==0
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);    
elseif length(tickOn) == 2 && tickOn(1)==1  && tickOn(2)==0   
    set(gca,'YTick',[]);
elseif length(tickOn) == 2 && tickOn(1)==0  && tickOn(2)==1
    set(gca,'XTick',[]);    
elseif length(tickOn) == 2 && tickOn(1)==0  && tickOn(2)==0
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);    
end

if colorbarOn
    cb= colorbar;
    ylabel(cb, 'Index', 'fontsize',fsize)
end

contourcmap('jet')
if mean(rangeZ) == 0
    caxis('auto');
else
    caxis(rangeZ);
end

set(gca,  'fontsize',fsize) ;
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') - [0 4 0]);
xlabh = get(gca,'YLabel');
set(xlabh,'Position',get(xlabh,'Position') - [0 0 0]);

end