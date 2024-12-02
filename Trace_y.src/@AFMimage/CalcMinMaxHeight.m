function [z_max, z_min, periodicity, peak_loc, peak_heights, min_loc, min_heights] = ...
    CalcMinMaxHeight(img, filament)

%
% DESCRIPTION
% – CalcMinMaxHeight finds the "significant" peaks and troughs of a 
% periodic signal using no statistics but based on peak prominance and 
% loose guidelines determined by the periodicity and the 
% prominence/amplitude of the signal.
% – Alows the user to manually pick/edit the peaks. The fig shows selected
% peaks for visual inspection of data and allows the guideline freqency 
% parameters to be adjusted if necessary.
% – Used by MakeHelicalFilamentModel.m as first step in 3D reconstruction
% – Part of Trace_y
%
% USAGE
% Standard method usage
% >> [z_max, z_min, periodicity, peak_loc, peak_heights, min_loc, min_heights] = ...
% >>   img.DeleteFilament(filament);
%
% INPUTS
% img  –  AFMimage object. Implicit with dot notation
% filament  –  The index number of the traced filament to be analysed
%
% OUTPUTS
% z_max  –  Average height value of the peaks
% z_min  –  Average height value of the troughs
% periodicity  –  FFT based period length in pixels
% peak_loc  –  List of sub-pixel peak locations in pixel units
% peak_heights  –  The height value of the individual peaks in the list
% min_loc  –  List of sub-pixel trough locations in pixel units
% min_heights  –  The height value of the individual troughs in the list
%
% DEPENDENCIES
% – Uses Signal Processing Matlab toolbox
% – Method for Trace_y's @AFMimage/ object.
% – Used by @AFMImage/MakeHelicalFilamentModel.m
%
% AUTHORS
% Liisa Lutter, Wei-Feng Xue
%
% HISTORY
% 2019.02  –  Initial draft by LL edited by WFX to incorporated draft into 
% the AFMimage class methods.
% 2021.03  –  WFX edited to include smoothing using spap2 and ginput for 
% graphical user peak selection/editing
% 2024.09  –  WFX Trace_y update including updates to the visualisation
% code and replaced old ginput method with callback functions for the UI.
%



%
% Get data
%
%{
if ~exist('filament', 'var') || filament == 0
    z = img.features.lastFilament.z;
    l = img.features.lastFilament.l;
    filament = 0;
else
    z = img.features.filaments(filament).z;
    l = img.features.filaments(filament).l;
end
%}

if ~exist('filament', 'var') || isempty(filament)
    filament = 1;
end
z = img.features.filaments(filament).z;
l = img.features.filaments(filament).l;



%
% Some initial evaluations based on frequency analysis
%

% Calculating FFT to find periodicity
[ff, amp] = img.CalcFilamentFFT1D(filament);

% Finding the peaks on the FFT
[pks, locs] = findpeaks(amp(2:end), ff(2:end));

%[ps, is] = findpeaks(flip(amp(9:end)), flip(1./ff(9:end)),...
%    'Annotate','extents','WidthReference','halfheight');

% Sort peaks in the frequncy domain from high to low amplitude
[pks, idx] = sort(pks, 'descend');
locs = locs(idx);

%[~, i] = max(pks);
%periodicity = 1/locs(i);

% Find all possible periodicity
periodicity_all = 1./locs;



%
% Plotting these same peaks for visual inspection of accuracy and selection
%

% Old settings no longer needed
%img.DispFilament(filament);
%set(1, 'OuterPosition', [1 601 500 500]);
%set(2, 'OuterPosition', [1 51 500 500]);

h_fig = figure(13);
clf;
set(h_fig, 'NumberTitle', 'off', 'Name', 'Trace_y – CalcMinMaxHeight');
h_fig.OuterPosition(3:4) = [1100 880];

%drawnow;
screenSize = get(0, 'ScreenSize');
if h_fig.OuterPosition(2)+h_fig.OuterPosition(4) > screenSize(4)
    h_fig.OuterPosition(2) = screenSize(4)-h_fig.OuterPosition(4)-10;
end
clf;


% Set up button check
%quitEval = 0;
% Max peak by default
current_peak = 1;
% Fine grid at 0.1 pixels
l_fine = (min(l(:)):0.1:max(l(:)))';
z_smooth = [];
min_peaks = [];
peak_dist = [];

% Init variables for Callback functions
h_peaks = [];
h_peaks2 = [];
h_ax1 = [];
h_ax3 = [];
h_ax5 = [];
h_ax6 = [];
h_tt = [];
h_freq = [];
h_period = [];
doneBtn = [];

x_lims = [];
freq_lims = [];
period_lims = [];

UpdatePlots;

    function [] = UpdatePlots(~, ~)
        % While not pressed 'esc' (27) or 'q' (113)
        %while ~quitEval

        % LL introduced smoothness for smoothing filament model using spap2 to
        % define nknots based on pixel size.  Use here to refine peak picking.
        % WFX simplifid to just 1 based smoothness. 1 is not additonal
        % smoothing
        smoothness = max([periodicity_all(current_peak)/4 1.01]);
        % Apparently a period can be described by 4 spline pieses
        % smoothness cannot be less than 1 so add lower bound

        % Order: k = 4 for cubic
        k = 4;
        % Define knots
        l_step = smoothness;
        l_knt = [l(1):l_step:l(end) l(end)];
        knots = optknt(l_knt, k);

        % Make spline (B-form)
        sp = spap2(knots, k, l, z);
        z_smooth = fnval(sp, l_fine);


        % Determining loose guidelines to identify all relevant peaks on the
        % height graph. These conditions are less good for troughs as the data
        % can be noisy in those regions

        periodicity = 1/locs(current_peak);

        [~, peak_loc] = findpeaks(z_smooth, l_fine, ...
            'MinPeakDistance', periodicity*0.5);

        [min_peaks, min_loc] = findpeaks(-z_smooth, l_fine,...
            'MinPeakDistance', periodicity*0.5);


        %
        % Plot info
        %

        tl = tiledlayout(3, 2, 'Padding', 'loose', 'TileSpacing', 'compact');
        title(tl, '   ', ...
            'interpreter', 'none', ...
            'FontName', 'helvetica', 'FontWeight', 'bold', ...
            'FontSize', 16);
        tx = xlabel(tl, sprintf('%s, filament %d', img.dataFile, filament), ...
            'interpreter', 'none', ...
            'FontName', 'helvetica', 'FontWeight', 'bold', ...
            'FontSize', 16);
        set(tx, 'FontWeight', 'normal', 'FontAngle', 'italic', 'FontSize', 16);

        doneBtn = uicontrol('Style', 'pushbutton', ...
            'String', 'Done', ...
            'Units', 'pixels', 'Position', [25 20 150 25], ...
            'Callback', @Done, 'Enable', 'on', 'UserData', [], ...
            'FontName', 'helvetica', ...
            'FontUnits', 'pixels', 'FontSize', 14, ...
            'FontWeight', 'normal');

        % Centre line profiles

        %subplot(3, 2, [1 2]);
        h_ax1 = nexttile(1, [1 2]);
        zz = img.StraightenFilament(filament);
        xx = 1:size(zz, 2);
        yy = (1:size(zz, 1))'-ceil(size(zz, 1)/2);
        h_img = imagesc(xx, yy, zz);
        set(gca, 'yDir', 'normal', ...
            'box', 'on', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', 'FontSize', 12);
        colormap(img.cMap);
        set(gca, 'XLim', [0 l(end)]);
        hold('on');
        h_peaks = plot(peak_loc, zeros(size(peak_loc)), 'c|', 'MarkerSize', 60, 'LineWidth', 1);
        xlabel('Length / pixel', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', ...
            'FontUnits', 'points', 'FontSize', 14);
        ylabel('Width / pixel', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', ...
            'FontUnits', 'points', 'FontSize', 14);
        title('Straightened filament with peak list  –  Left click to add peak, Right click to remove peak(s)', ...
            'interpreter', 'none', ...
            'FontName', 'helvetica', 'FontWeight', 'bold', ...
            'FontUnits', 'points', 'FontSize', 14);

        h_img.ButtonDownFcn = @SelPeaks;
        h_peaks.ButtonDownFcn = @SelPeaks;


        %subplot(3, 2, [3 4]);
        h_ax3 = nexttile(3, [1 2]);
        findpeaks(z_smooth, l_fine, ...
            'MinPeakDistance', periodicity*0.5);
        hold('on');
        plot(l, z, 'r-');
        h_peaks2 = plot(peak_loc, pchip(l_fine, z_smooth, peak_loc), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
        h_data = plot(l_fine, z_smooth, 'b-');
        set(gca, 'box', 'on', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', 'FontSize', 12);
        xlabel('Peak locations / pixel', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', ...
            'FontUnits', 'points', 'FontSize', 14);
        ylabel('Height / nm', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', ...
            'FontUnits', 'points', 'FontSize', 14);
        legend({'Smoothed data' 'Initial estimated peaks' 'Image data' 'Final (edited) peaks'}, ...
            'Location', 'southeast', 'Orientation', 'horizontal', 'FontSize', 12);
        
        if ~isempty(x_lims)
            h_ax1.XLim = x_lims;
            h_ax3.XLim = x_lims;
        end
        linkaxes([h_ax1 h_ax3], 'x');

        h_data.ButtonDownFcn = @SelPeaks;
        h_peaks2.ButtonDownFcn = @SelPeaks;

        %{
        subplot(3, 2, [5 6]);
        findpeaks(-z_smooth, l, ...
            'MinPeakDistance', periodicity*0.5);
        hold('on');
        plot(l, -z, 'r-');
        xlabel('Minima locations / pixel', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', ...
            'FontUnits', 'pixels', 'FontSize', 14);
        ylabel('-Height / nm', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', ...
            'FontUnits', 'pixels', 'FontSize', 14);
        %}



        % 1D PSD and periodicity
        %subplot(3, 2, 5);
        h_ax5 = nexttile(5, [1 1]);
        h_freq = plot(ff, amp, 'b-');
        hold('on');
        plot(locs(current_peak), pks(current_peak), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
        set(gca, 'box', 'on', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', 'FontSize', 12);
        xlabel('Frequency / pixel^-^1', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', ...
            'FontUnits', 'points', 'FontSize', 14);
        ylabel('Power', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', ...
            'FontUnits', 'points', 'FontSize', 14);
        title('1D power spectrum', 'interpreter', 'none', ...
            'FontName', 'helvetica', 'FontWeight', 'bold', ...
            'FontUnits', 'points', 'FontSize', 14);
        if ~isempty(freq_lims)
            h_ax5.XLim = freq_lims;
        end

        h_freq.ButtonDownFcn = @SelFreq;

        %subplot(3, 2, 6);
        h_ax6 = nexttile(6, [1 1]);
        h_period = plot(1./ff, amp, 'b-');
        hold('on');
        plot(1./locs(current_peak), pks(current_peak), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
        set(gca, 'box', 'on', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', 'FontSize', 12);
        xlabel('Periodic length / pixel', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', ...
            'FontUnits', 'points', 'FontSize', 14);
        ylabel('Power', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', ...
            'FontUnits', 'points', 'FontSize', 14);
        title('Periodicity', 'interpreter', 'none', ...
            'FontName', 'helvetica', 'FontWeight', 'bold', ...
            'FontUnits', 'points', 'FontSize', 14);
        if ~isempty(period_lims)
            h_ax6.XLim = period_lims;
        end

        h_period.ButtonDownFcn = @SelPeriod;

        peak_dist = diff(peak_loc);
        h_tt = text(0.95, 0.9, sprintf('Periodicity: %g\nMean p-p dist: %g\nClick to select freq. peak', ...
            1/locs(current_peak), mean(peak_dist)), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', 'FontSize', 12, ...
            'FontName', 'helvetica', 'FontWeight', 'normal', 'Color', 'k');

    end


% This part with the old ginput method is now obsolete
%{
subplotClicked = 0;
while ~quitEval && subplotClicked ~= 1 && subplotClicked ~= 2
    [xx, ~, button] = ginput(1);

    axesClicked = gca;
    currentfig = gcf;
    allAxes = findobj(currentfig.Children, 'Type', 'axes');
    subplotClicked = find(allAxes == axesClicked);
    % The indexes are in reverse order of subplots added
    if isempty(button) || button == 27 || button == 113
        quitEval = 1;

        % If click on top filament image
    elseif (button == 1 || button == 3) && subplotClicked == 4
        % If left buton click, add peak
        if button == 1
            peak_loc = sort([peak_loc; xx]);
        end

        % If right button click, remove peak within 5 pixels
        if button == 3
            peak_loc = peak_loc(abs(peak_loc-xx) > 5);
        end

        subplot(3, 2, [1 2]);
        delete(h_peaks);
        h_peaks = plot(peak_loc, size(zz, 1)/2, 'r+', 'MarkerSize', 10, 'LineWidth', 2);

        peak_dist = diff(peak_loc);
        set(tt, 'String', ...
            sprintf('Periodicity: %g\nMean p-p dist: %g\nClick to select freq. peak\nPress return to exit', ...
            1/locs(current_peak), mean(peak_dist)))
    end

end

% Next query if needed
if ~quitEval && subplotClicked == 1
    loc_query = 1./xx;
    current_peak = find(abs(locs-loc_query) == min(abs(locs-loc_query)), 1);
elseif ~quitEval && subplotClicked == 2
    loc_query = xx;
    current_peak = find(abs(locs-loc_query) == min(abs(locs-loc_query)), 1);
end

%}


    function [] = SelPeaks(~, evt)
        xx = evt.IntersectionPoint(1);
        %yy = evt.IntersectionPoint(2);
        button = evt.Button;

        if button == 1
            peak_loc = sort([peak_loc; xx]);

        elseif button == 3 || button == 2
            % If right button click, remove peak within 5 pixels
            peak_loc = peak_loc(abs(peak_loc-xx) > 5);
        end

        delete(h_peaks);
        h_peaks = plot(h_ax1, peak_loc, zeros(size(peak_loc)), 'c|', 'MarkerSize', 60, 'LineWidth', 1);
        h_peaks.ButtonDownFcn = @SelPeaks;
        delete(h_peaks2);
        h_peaks2 = plot(h_ax3, peak_loc, pchip(l_fine, z_smooth, peak_loc), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
        h_peaks2.ButtonDownFcn = @SelPeaks;

        delete(h_tt);
        peak_dist = diff(peak_loc);
        h_tt = text(h_ax6, 0.95, 0.9, sprintf('Periodicity: %g\nMean p-p dist: %g\nClick to select freq. peak', ...
            1/locs(current_peak), mean(peak_dist)), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', 'FontSize', 12, ...
            'FontName', 'helvetica', 'FontWeight', 'normal', 'Color', 'k');
    end


    function [] = SelFreq(~, evt)
        xx = evt.IntersectionPoint(1);
        loc_query = xx;
        current_peak = find(abs(locs-loc_query) == min(abs(locs-loc_query)), 1);
        x_lims = h_ax3.XLim;
        freq_lims = h_ax5.XLim;
        period_lims = h_ax6.XLim;
        clf;
        UpdatePlots;
    end

    function [] = SelPeriod(~, evt)
        xx = evt.IntersectionPoint(1);
        loc_query = 1./xx;
        current_peak = find(abs(locs-loc_query) == min(abs(locs-loc_query)), 1);
        x_lims = h_ax3.XLim;
        freq_lims = h_ax5.XLim;
        period_lims = h_ax6.XLim;
        clf;
        UpdatePlots;
    end

    function [] = Done(~, ~)
        % Tripe check here for status, a bit of overkill
        h_fig.UserData = 1;
        delete(doneBtn);
        pause(0.1);
        close(h_fig);
    end



waitfor(h_fig, 'UserData', 1);

% Saving the average values of peaks and troughs
peak_heights = interp1(l_fine, z_smooth, peak_loc);
min_heights = abs(min_peaks);

z_max = mean(peak_heights);
z_min = mean(min_heights);



end