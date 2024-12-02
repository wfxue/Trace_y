function [] = DispLine(img, scan_line, slowAxis_line)

%
% DESCRIPTION
% – Display the image scan lines and the slow scan axis lines (horizontal
% and vertical lines).
% – Click on the image to select coordinates for the lines to show.
% – Part of Trace_y by WFX
%
% USAGE
% Standard method usage with filament number as input:
% >> img.DispFilament;
% Specify the initial coordinates for the lines to display:
% >> img.DispLine(scan_line, slowAxis_line)
%
% INPUTS
% img  –  AFMimage object. Implicit with dot notation
% scan_line  –  The index number (y coordinate in px) of the scan line to
% display.
% slowAxis_line  –  The index number (x coordinate in px) of the slow scan
% axis line to display.
%
% OUTPUTS
% Interactive plot of the image and the line profiles.
%
% DEPENDENCIES
% – Method for Trace_y's @AFMimage/ object
%
% AUTHORS
% Wei-Feng Xue
%
% HISTORY
% 2014.03  –  First draft
% 2024.09  –  WFX Trace_y update including updates to the visualisation.
% Here, also added better interactions with improved callback code.
%


%
% Defaults
%
if ~exist('scan_line', 'var')
    scan_line = round(img.nLines/2);
end
if ~exist('slowAxis_line', 'var')
    slowAxis_line = round(img.pixelPerLine/2);
end

%scan_line = round(scan_line);
%ll = scan_line;

%button = 0;

% Data
z = img.z;
pixelPerLine = img.pixelPerLine;
nLines = img.nLines;



%
% Init the figure window
%
h_fig = figure(101);
set(h_fig, 'NumberTitle', 'off', 'Name', 'Trace_y – DispLine');
h_fig.OuterPosition(3:4) = [805 880];

%drawnow;
screenSize = get(0, 'ScreenSize');
if h_fig.OuterPosition(2)+h_fig.OuterPosition(4) > screenSize(4)
    h_fig.OuterPosition(2) = screenSize(4)-h_fig.OuterPosition(4)-10;
end
clf;

%figPos = [1 300 1400 700];
%set(gcf, 'OuterPosition', figPos);
%ax1 = subplot(1, 2, 1);
%cla;
%img.DispImage;
%hold on

tl = tiledlayout(4, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tl, '   ', ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontSize', 16);
tx = xlabel(tl, ...
    sprintf('%s  %s', img.dataFile, img.imgType), ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontSize', 16);
set(tx, 'FontWeight', 'normal', 'FontAngle', 'italic', 'FontSize', 16);


h_ax_img = nexttile(tl, 5, [3 3]);
h_img = imagesc(img.z);
hold('on');
set(gca, 'yDir', 'normal', ...
    'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1], ...
    'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1], ...
    'box', 'on', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 12);
h_img.ButtonDownFcn = @SelectLines;
colormap(gca, img.cMap);
set(gca, 'CLim', img.cScale);
%colorbar;
xlabel('x / pixel', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 16);
ylabel('y / pixel', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 16);
%title(' ', ...
%    'interpreter', 'none', ...
%    'FontName', 'helvetica', 'FontWeight', 'bold', ...
%    'FontUnits', 'points', 'FontSize', 16);

% Listeners for updating line plot axes
addlistener(h_ax_img, 'XLim', 'PostSet', @UpdateLineAxes);
addlistener(h_ax_img, 'YLim', 'PostSet', @UpdateLineAxes);
hz = zoom;
hz.ActionPostCallback = @UpdateLineAxes;
hp = pan;
hp.ActionPostCallback = @UpdateLineAxes;


% Initial state
h_vl = [];
h_hl = [];
h_ax1 = [];
h_ax2 = [];
PlotLines(tl, scan_line, slowAxis_line);



%
% Callback function 
%
    function [] = SelectLines(~, evt)
        yy = round(evt.IntersectionPoint(2));
        xx = round(evt.IntersectionPoint(1));
        PlotLines(tl, yy, xx);
    end

    function [] = PlotLines(target, yy, xx)

        %img_xlims = h_ax_img.XLim;
        %img_ylims = h_ax_img.YLim;

        % Vertical slow scan axis line
        h_ax1 = nexttile(target, 8, [3, 1]);
        cla;
        plot(z(:, xx), 1:pixelPerLine, 'b-');
        set(h_ax1, 'YTickMode', 'manual');
        %ylim([1 pixelPerLine]);
        set(h_ax1, ...
            'box', 'on', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', ...
            'FontUnits', 'points', 'FontSize', 12);
        xlabel('z / nm', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', ...
            'FontUnits', 'pixels', 'FontSize', 16);
        legend(h_ax1, 'Slow scan axis', 'Location', 'southeast', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', 'FontSize', 10);

        delete(h_vl);
        h_vl = plot(h_ax_img, [xx, xx], ylim, 'c-', 'LineWidth', 0.5);


        % Horizontal scan line
        h_ax2 = nexttile(target, 1, [1, 3]);
        cla;
        plot(1:nLines, z(yy, :), 'b-');
        set(h_ax2, 'XTickMode', 'manual');
        %xlim([1 nLines]);
        set(h_ax2, ...
            'box', 'on', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', ...
            'FontUnits', 'points', 'FontSize', 12);
        ylabel('z / nm', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', ...
            'FontUnits', 'pixels', 'FontSize', 16);
        legend(h_ax2, 'Scan line', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', 'FontSize', 10);

        delete(h_hl);
        h_hl = plot(h_ax_img, xlim, [yy, yy], 'c-', 'LineWidth', 0.5);

        % Link axes to the image
        %linkaxes([h_ax_img h_ax1], 'y');
        %linkaxes([h_ax_img h_ax2], 'x');
        %h_ax_img.XLim = img_xlims;
        %h_ax_img.YLim = img_ylims;
        UpdateLineAxes;

        % Colorbar
        % Looks silly so not using
        %{
        h_cb = colorbar;
        h_cb.Layout.Tile = 'east';
        set(h_cb, 'FontName', 'helvetica', 'FontWeight', 'normal', ...
            'FontSize', 12)
        set(h_cb.Label, 'String', sprintf('%s / %s', img.imgType, img.zUnit), ...
            'FontName', 'helvetica', 'FontWeight', 'normal', ...
            'FontSize', 16)
        colormap(h_ax2, img.cMap);
        set(h_ax2, 'CLim', img.cScale);
        %}
    end

    function [] = UpdateLineAxes(~, ~)
        xlims = get(h_ax_img, 'XLim');
        ylims = get(h_ax_img, 'YLim');
        set(h_ax2, 'XLim', xlims);
        set(h_ax1, 'YLim', ylims);
    end



% This is old code using ginput and is obsolete
%{

hl = plot(xlim, [scan_line, scan_line], 'w-');

%hf2 = figure(2);
%figPos = get(gcf, 'OuterPosition');
%figPos(1) = figPos(3)+6;
%figPos = [701 300 700 700];
%set(gcf, 'OuterPosition', figPos);
ax2 = subplot(1, 2, 2);

while button <= 3
    % While mouse button is pressed in for ginput
    %figure(hf1);
    %subplot(ax1);
    
    if ll >= 1 && ll <= img.nLinesImg
        scan_line = round(ll);
        
        delete(hl);
        hl = plot(xlim, [scan_line, scan_line], 'w-', 'LineWidth', 2);
        
        %figure(hf2);
        subplot(ax2);
        plot(1:img.pixelPerLineImg, img.z(scan_line, :), 'b-');
        set(gca, 'XLim', [1 img.pixelPerLineImg]);
        xlabel('x / pixel', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', ...
            'FontUnits', 'pixels', 'FontSize', 14);
        ylabel('z / nm', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', ...
            'FontUnits', 'pixels', 'FontSize', 14);
        title(sprintf('Scan line y = %g, press ''esc'' to quit', scan_line), ...
            'interpreter', 'none', ...
            'FontName', 'helvetica', 'FontWeight', 'bold', ...
            'FontUnits', 'pixels', 'FontSize', 14);
        %figure(hf1);
        subplot(ax1);
    end

    [~, ll, button] = ginput(1);

end


% Finishing
subplot(ax2);
title(sprintf('Scan line y = %g', scan_line), ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontUnits', 'pixels', 'FontSize', 14);

% Close figure windows
%close(hf1);
%close(hf2);
%}



end