function [h] = DispImage(img, targetFig)

%
% DESCRIPTION
% – Display the image with fixed aspect ratio with double axes, one direct
% axes in px units and one axes with evaluated units (e.g. nm).
% – Part of Trace_y
%
% USAGE
% Standard method usage
% >> img.DispImage;
%
% INPUTS
% img  –  AFMimage object. Implicit with dot notation
% targetFig  –  Target figure window handle to use
%
% OUTPUTS
% h  –  Handle to the tiledlayout graphics object
% The image is ploted in existing or new figure. The figure index is
%   implicit and not fixed
%
% DEPENDENCIES
% – Method for Trace_y's @AFMimage/ object.
% – Used also by other routines
%
% AUTHORS
% Wei-Feng Xue
%
% HISTORY
% 2015.08  –  Imitial draft of simple DispImage.m
% 2024.08  –  WFX updated to plot image on double axes based on user 
% feedback that length scales are complicated to understand.
%



%
% Set up tile overlap to prep for double axes
%

% Init the image
if ~exist('targetFig', 'var') || isempty(targetFig)
    h_fig = figure(10);
    h_fig.OuterPosition(3:4) = [960 960];
else
    h_fig = figure(targetFig);
end

set(h_fig, 'NumberTitle', 'off', 'Name', 'Trace_y – DispImage');

%drawnow;
screenSize = get(0, 'ScreenSize');
if h_fig.OuterPosition(2)+h_fig.OuterPosition(4) > screenSize(4)
    h_fig.OuterPosition(2) = screenSize(4)-h_fig.OuterPosition(4)-10;
end

clf;
%cla;
h_tl = tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'tight');
title(h_tl, '   ', ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontSize', 12);
tx = xlabel(h_tl, ...
    sprintf('%s  %s', img.dataFile, img.imgType), ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontSize', 16);
set(tx, 'FontWeight', 'normal', 'FontAngle', 'italic', 'FontSize', 16);



%
% Set up the axes with secondary (indirect / evaluated) units
%
ax2 = axes(h_tl);
set(ax2, 'yDir', 'normal', ...
    'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1], ...
    'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1], ...
    'box', 'on', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 12);
xlabel(sprintf('x / %s', img.scanSizeUnit), ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 16);
ylabel(sprintf('y / %s', img.scanSizeUnit), ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 16);
% Turn off everything for these axes in the back layer
ax2.Toolbar.Visible = 'off';
ax2.Interactions = [];
ax2.TickDir = 'none';
% Initial axes lims
xlim([0 img.scanSizex]);
ylim([0 img.scanSizey]);



%
% Set up the axes with direct (px) units and put image into it
%
ax1 = axes(h_tl);
imagesc(ax1, 1:img.pixelPerLine, 1:img.nLines, img.z);
set(gca, 'yDir', 'normal', ...
    'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1], ...
    'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1], ...
    'box', 'on', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 12);

colormap(gca, img.cMap);
set(gca, 'CLim', img.cScale);
xlabel('x / pixel', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 16);
ylabel('y / pixel', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 16);

% Colourbar
cb = colorbar(ax1, 'eastoutside');
set(cb, 'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontSize', 12)
set(cb.Label, 'String', sprintf('%s / %s', img.imgType, img.zUnit), ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontSize', 16)
ax1.XAxisLocation = 'top';
ax1.YAxisLocation = 'right';
ax1.TickDir = 'none';

% Call back listeners
xRes = img.xResolution;
yRes = img.yResolution;
%xLimListener = addlistener(ax1, 'XLim', 'PostSet', @Update2ndAxes);
%yLimListener = addlistener(ax1, 'YLim', 'PostSet', @Update2ndAxes);
addlistener(ax1, 'XLim', 'PostSet', @Update2ndAxes);
addlistener(ax1, 'YLim', 'PostSet', @Update2ndAxes);
hz = zoom;
hz.ActionPostCallback = @Update2ndAxes;
hp = pan;
hp.ActionPostCallback = @Update2ndAxes;


%
% Function to update the evaluated axes
%
    %function [] = Update2ndAxes(src, evt)
    function [] = Update2ndAxes(~, ~)
        xlims1 = get(ax1, 'XLim');
        ylims1 = get(ax1, 'YLim');
        xlims2 = (xlims1-1).*xRes;
        ylims2 = (ylims1-1).*yRes;
        set(ax2, 'XLim', xlims2);
        set(ax2, 'YLim', ylims2);
    end



%
% Output handle if asked
%
if nargout == 1
    h = h_fig;
end



end