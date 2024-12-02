function [img, currentFilament] = MakeHelicalFilamentModel(img, ...
    filamentIndex, handedness, symmetry, tipData, refineTip, smoothness, refine)

%
% DESCRIPTION
% – Semi automatic helical 3D reconstruction of traced filaments. Assumes 
% the filaments have helical symmetry. Outputs a 3D model of the filament
% along with the cross-section point-cloud.
% – Algorithm first described in Lutter, L. et al. Three-dimensional 
% reconstruction of individual helical nano-filament structures from 
% atomic force microscopy topographs. Biomol Concepts 11, 102-115, (2020). 
% https://doi.org/10.1515/bmc-2020-0009
% – Kernel density bandwidth estimation for x/z cross-section point cloud 
% uses a non-parametric bandwidth selection method for ksdensity adapted 
% from Kernel Density Estimator kde.m on
% https://www.mathworks.com/matlabcentral/fileexchange/14034-kernel-density-estimator
% and Kernel density estimation via diffusion Z. I. Botev, 
% J. F. Grotowski, and D. P. Kroese (2010) Annals of Statistics, Volume 38,
% Number 5, pages 2916-2957, doi:10.1214/10-AOS799
% – Part of Trace_y
%
% USAGE
% – Standard method usage with user spplied twist handedness and 
% cross-section symmetry number.
% >> img = img.MakeHelicalFilamentModel(filamentIndex, handedness, symmetry);
% – Also user-supply initial tip information and repine tip geometry (slow)
% >> img = img.MakeHelicalFilamentModel(filamentIndex, ...
%   handedness, symmetry, tipData, refineTip);
% – Also set how smooth and refined the model surface appears
% >> img = img.MakeHelicalFilamentModel(filamentIndex, ...
%   handedness, symmetry, tipData, refineTip, smoothness, refine);
% – Get the filament object data separatly from the img object
% >> [img, currentFilament] = ...
%   img.MakeHelicalFilamentModel(filamentIndex, handedness, symmetry);
%
% INPUTS
% img  –  AFMimage object. Implicit with dot notation
% filamentIndex  –  The index number of filament to be modelled
% handedness  –  The twist handedness, 'left' or 'right'
% symmetry  –  Cross-sectional symetry number
% tipData  –  Tip geometry information, e.g. [tip_radius tip_angle]. The 
% default tip side angle is 18 deg, is optional
% refineTip  –  To refine tip radius parameter or not, boolean. Optional,
% with default to skip the slow itterative refine process.
% smoothness, refine  –  How smooth the model surface appears. Affects
% surface coordinate interpolation. The parameters are integers. The 
% smoothness parameter is how many angular points to sample for the 
% cross-section. It controls how smooth the spline interpolation for the 
% cross-section is, 1 is original angular density, 2 is twice as smooth 
% etc. The refine (or y_refine) parameter is how many subpixels to sample 
% in the filament y-axis, 1 is the original pixel density, 2 is twice as 
% many etc., for example 4 is four sub-divisions per pixel. These are 
% optional with default 4 and 2, respectively.
%
% OUTPUTS
% img  –  Updated img object with the 3D model data.
% currentFilament  –  The Filament object with the model separately (also
% contained in img.features.filaments property).
%
% DEPENDENCIES
% – Uses Statistical and Machine Learning, Image Processing, 
% Optimization, Signal Processing Matlab toolboxes.
% – Method for Trace_y's @AFMimage/ object.
% – Uses CalcMinMaxHeight.m, CalcFilamentFFT1D.m, 
% CalcTipFilamentConvolution.m and FilamentModel.m
%
% AUTHORS
% Wei-Feng Xue, Liisa Lutter, Liam D Aubrey
%
% HISTORY
% 2019.08  –  WFX drafted the first version for Aubrey and Blakeman et al
% (2020), https://doi.org/10.1038/s42004-020-00372-3.
% 2021.01  –  LL drafted changes for incorporating non-gridded (in the x/z
% cross-section) reconstruction. Was used for Lutter et al. (2020), 
% https://doi.org/10.1515/bmc-2020-0009.
% 2021.03  –  WFX edited LL draft changes and updated method. Also 
% added posibility to input handedness and symmetry. Display options are 
% modified to includ LDA surface modifications
% 2024.09  –  WFX Trace_y update including updates to the visualisation
% code and replaced old ginput method with callback functions for various 
% UI elements. Also fixed few bugs relating to vector initiation ranges.
%



%
% Constants and settings
%

% Test max symetry to max_sym, 2 or more, standard is 4
% This UI is no longer used as users will have to explicitly test
%max_sym = 4;

if ~exist('tipData', 'var') || isempty(tipData)
    a_tip = 18;
else
    a_tip = tipData(2);
end

% Smoothness
if ~exist('smoothness', 'var') || isempty(smoothness)
    smoothness = 4;
end

% Refine or not
if ~exist('refine', 'var') || isempty(refine)
    refine = 2;
end

% Refine tip_r or not
if ~exist('refineTip', 'var') || isempty(refineTip)
    refineTip = false;
elseif strcmpi(refineTip, 'no')
    refineTip = false;
elseif strcmpi(refineTip, 'yes')
    refineTip = true;
end


% Get the filament to work on
%{
if (~exist('filamentIndex', 'var') || strcmp(filamentIndex, 'last')) && ~isempty(img.features.lastFilament)
    %filamentIndex = 0;
    filamentIndex = numel(img.features.filaments);
elseif ~exist('filamentIndex', 'var') && isempty(img.features.lastFilament)
    filamentIndex = length(img.features.filaments);
end

if filamentIndex == 0
    f = img.features.lastFilament;
else
    f = img.features.filaments(filamentIndex);
end
%}
f = img.features.filaments(filamentIndex);


% Evaluate few basic parameters
fprintf('\nHelical 3D reconstruction of %s filament %g\n', img.dataFile, filamentIndex);
fprintf('Estimating peak locations and periodicity...\n');
[z_max, z_min, periodicity, peaks_loc] = img.CalcMinMaxHeight(filamentIndex);
r_max = z_max/2;

fprintf('Estimating tip-sample convolution and tip radius...\n\n');
%figure(3);
%set(3, 'OuterPosition', [501 601 500 500]);
%clf;
%[~, r_tip] = img.Calc_Filament_Corr(filamentIndex, a_tip);
% Initial convolution estimate code updated
[~, r_tip, ~, appwidth] = img.CalcTipFilamentConvolution(filamentIndex, a_tip);
img.features.filaments(filamentIndex).appWidth = appwidth;

fprintf('Average max height: z_max = %g nm\n', z_max);
fprintf('Average min height: z_min = %g nm\n', z_min);
fprintf('Average max filament radius: r_max = z_max/2 = %g nm\n', r_max);
fprintf('Periodicity: periodicity = %g pixels\n', periodicity);
fprintf('Tip radius initial estimate: r_tip = %g nm\n', r_tip);
fprintf('Tip angle: a_tip = %g deg\n', a_tip);

% Get image data for the filament
z = img.StraightenFilament(filamentIndex);
z = imrotate(z, -90);
x = 0:img.xResolution:img.xResolution*(size(z, 2)-1);
y = (0:img.xResolution:img.xResolution*(size(z, 1)-1))';
imgResolution = img.xResolution;
px = img.xResolution;


% Correct for tip-sample convolution
% Technically this should be done on raw image but this order on
% straightened filament is faster, but assumes relative straight filament
% and symetric tip. A point for future developements
%[xc, yc, zc, ~, ~, ~, unreliability] = AFM_Image_Correction(r_tip, a_tip, x, y, z, 0.1);
% As with convolution estimation, contact point reconstruction code also 
% updated. Everything is now faster without long loops etc.
img = img.SetTipModel('rounded_cone', [r_tip, a_tip]);
[xc, yc, zc, unreliability] = CalcContactPoints(x, y, z, img.zT{1}, img.zT{2}, 20);
unreliability_rotated = imrotate(unreliability, -90);
centre_line = ceil(size(xc, 2)/2);


% User input needed
% Deciding on pixel data to use, i.e. selecting pixel lines to use

% Old UI for selecting lines
%{
figure(4);
set(4, 'OuterPosition', [501 51 500 960]);
%set(4, 'OuterPosition', [501 51 500 500]);
clf;

imagesc(unreliability);
colormap(gca, jet);
colorbar;
set(gca, 'YDir', 'normal');
hold('on');
centre_line = ceil(size(xc, 2)/2);
ylims = ylim;
plot([centre_line centre_line], ylims, 'y-', 'LineWidth', 2);
xlabel('pixel', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'pixels', 'FontSize', 14);
ylabel('pixel', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'pixels', 'FontSize', 14);
title(sprintf('\nBased on the uncertainty map, how many lines to zero pad?\nMark either left or right outer bounds of the filament!\n'), ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontUnits', 'pixels', 'FontSize', 14, 'Color', 'r');
%}

% Continue to use the figure window from CalcTipFilamentConvolution.m
% instead
pause(0.2);
h_fig = figure(14);
h_img = findobj(h_fig, 'Type', 'image');
% Set background colour to red for overlay with the unreliability as Alpha 
% map
h_img.Parent.Color = [0 0 1];
h_img.AlphaData = 1-(unreliability_rotated./max(unreliability_rotated(:)));
h_tt = text(h_img.Parent, 0.025, 0.025, ...
    sprintf('Blue overlay shows the unreliability of the pixels.\nSelect image lines to use for 3D reconstruction'), ...
    'Units', 'normalized', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom', 'FontSize', 12, ...
    'FontName', 'helvetica', 'FontWeight', 'normal', 'Color', 'w');
h_img.ButtonDownFcn = @SelLines;
h_pa = [];
n_lines = [];

doneBtn = uicontrol(h_fig, 'Style', 'pushbutton', ...
    'String', 'Done', ...
    'Units', 'pixels', 'Position', [15 10 100 25], ...
    'Callback', @Done, 'Enable', 'off', 'UserData', [], ...
    'FontName', 'helvetica', ...
    'FontUnits', 'pixels', 'FontSize', 14, ...
    'FontWeight', 'normal');


% Call back functions for lines selection

    function [] = SelLines(~, evt)
        %xx = evt.IntersectionPoint(1);
        yy = evt.IntersectionPoint(2);
        %button = evt.Button;

        yy = round(yy/imgResolution);
        n_lines = 2*abs(yy)+1;
        xx = h_img.Parent.XLim+[0.002 -0.002]*(h_img.Parent.XLim(2)-h_img.Parent.XLim(1));
        if yy > 0
            yy = [-yy-0.5 yy+0.5];
        else
            yy = [yy-0.5 -yy+0.5];
        end

        delete(h_pa);
        h_pa = patch([xx(1) xx(2) xx(2) xx(1)], ...
            imgResolution.*[yy(1) yy(1) yy(2) yy(2)], 'c', ...
            'FaceColor', 'none', 'EdgeColor', 'c', 'LineWidth', 2);
        h_tt.String = ...
            sprintf('Blue overlay shows the unreliability of the pixels.\nSelected %g image lines to use for 3D reconstruction', ...
            n_lines);
        doneBtn.Enable = 'on';
    end


    function [] = Done(~, ~)
        % Tripe check here for status, a bit of overkill
        h_fig.UserData = 1;
        delete(doneBtn);
        pause(0.1);
        close(h_fig);
    end


waitfor(h_fig, 'UserData', 1);
fprintf('Number of pixel lines to use: n_lines = %g\n', n_lines);



% Old UI method no longer needed
%{
%fprintf('\nBased on the uncertainty map, how many lines to zero pad?\n');
%fprintf('Mark either left or right outer bounds of the filament!\n');
%pause
[xx, ~] = ginput(1);
xx = round(xx);
if xx > centre_line
    xx = xx-2*(xx-centre_line);
end
z(:, [1:xx end-xx+1:end]) = 0;
plot([xx xx], ylims, 'r-', 'LineWidth', 2);
plot([size(z, 2)-xx+1 size(z, 2)-xx+1], ylims, 'r-', 'LineWidth', 2);

title(sprintf('\nHow many lines to use for 3D model reconstruction?\nMark either left or right inner bounds of the filament!\n'), ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontUnits', 'pixels', 'FontSize', 14, 'Color', 'r');

%fprintf('\nBased on the uncertainty map, how many lines to use for 3D model reconstruction?\n');
%fprintf('Mark either left or right inner bounds of the filament!\n');
[xx, ~] = ginput(1);
xx = round(xx);
plot([xx xx], ylims, 'r-', 'LineWidth', 2);
plot([size(z, 2)-xx+1 size(z, 2)-xx+1], ylims, 'r-', 'LineWidth', 2);
%}


% Use apparent width from simulatios to do zero padding
z_orig = z;
bg_lines = (size(z, 2)-appwidth)/2;
z(:, [1:bg_lines end-bg_lines+1:end]) = 0;
%n_lines = 2*abs(xx-centre_line)+1;

% LL 2021 Jan edit for non-gridded reconstruction
side_lines = floor(n_lines/2);
selected_lines = centre_line-side_lines:centre_line+side_lines;
%centre_axis = [mean(xc(:)) z_max/2];

% Data points weighted by w which scales between (0,1)
% High uncertainty means low w, w-> 1 when uncertainty->0
w = 1./(1+unreliability.^2);

% Original gridded y
%{
% Interpolate slightly to get y on grid, i.e. the slices
yci = y*ones(1, n_lines);
xci = xc(:, centre_line-floor(n_lines/2):centre_line+floor(n_lines/2));
zci = griddata(xc, yc, zc, xci, yci, 'cubic');
centre_axis = [mean(xc(:)) z_max/2];

unreliability_i = griddata(xc, yc, unreliability, xci, yci, 'linear');
w = 1./(1+unreliability_i.^2);
%}

% Test code
%{
% Check points in cylindrical coordinates
yci = (min(y):img.yResolution/4:max(y)+0.9*img.yResolution/4)';
xci = xc(:, centre_line-floor(n_lines/2):centre_line+floor(n_lines/2));
yc = yc(:, centre_line-floor(n_lines/2):centre_line+floor(n_lines/2));
zci = zc(:, centre_line-floor(n_lines/2):centre_line+floor(n_lines/2));

%zci = griddata(xc, yc, zc, xci, yci, 'cubic');
centre_axis = [mean(xc(:)) z_max/2];

[theta, rho] = cart2pol(xci-centre_axis(1), zci-centre_axis(2));
theta = rad2deg(theta);
mean_arc = mean(theta, 1);

thetai = round(min(mean_arc)):1:round(max(mean_arc));
[thetai, yci] = meshgrid(thetai, yci);
rhoi = griddata(theta, yc, rho, thetai, yci, 'cubic');

plot(theta, rho, '+');
%}



% User input needed
% Deciding on handedness and symmetry based on user id of largest crosspeak
% in the 2d fft image

% Calculating 2D FFT map for later use
%{
%trim = ceil(periodicity*2);
% trim periodicity*x pixels on each sidemake it ready for symmery number of
% up to 2*x
% Trimming not needed anymore since the simlated image will be of the same
% size as the data image
trim = 0;
z_trim = z(trim+1:end-trim, :);
% z_trim must have odd number of rows because it will have odd number of
% columns, after padding both dimensions will be odd
if mod(size(z_trim, 1), 2) == 0
    z_trim = z_trim(1:end-1, :);
end
z_size = size(z_trim);
% Pad to square
zpad_size = [max(z_size), max(z_size)];
zpad = padarray(z_trim, (zpad_size-z_size)/2);


% 2d fft
psd2d = fftshift(fft2(zpad));
psd2d = log10(abs(psd2d).^2);

% LL added frequency scale display
% Sampling freq / Å^-1
px = img.xResolution;
freq = 1/(10*px);
% Frequnency scale
sz = max(size(zpad));
s_c = ceil(sz/2);
r = 1:1*s_c;
fx = freq.*r'/sz;
fxx = [-fliplr(fx') fx(2:end)'];

cmap = jet(512);
%}


%
% Handedness evaluation if needed
%
% No longer done this way. User should need to suply the value
if ~exist('handedness', 'var') || isempty(handedness)
    
    %{
    figure(5);
    set(5, 'OuterPosition', [501 151 500 960]);
    clf;
    subplot(1, 2, 1);
    imagesc(fxx, fxx, psd2d);
    set(gca, 'YDir', 'normal');
    colormap(cmap);
    xlabel('Freqency / Å^-^1', ...
        'FontName', 'helvetica', 'FontWeight', 'normal', ...
        'FontUnits', 'pixels', 'FontSize', 14);
    ylabel('Freqency / Å^-^1', ...
        'FontName', 'helvetica', 'FontWeight', 'normal', ...
        'FontUnits', 'pixels', 'FontSize', 14);
    title(sprintf('2D FFT\n'), ...
        'interpreter', 'none', ...
        'FontName', 'helvetica', 'FontWeight', 'bold', ...
        'FontUnits', 'pixels', 'FontSize', 14);
    
    text(0.05, 0.08, sprintf('Select twist handedness...\nFourier space main peaks:\n     / is left;  \\ is right\nDirect image pattern:\n     \\ is left;  / is right\n\nPress (l)eft or (r)ight'), ...
        'Units', 'normalized', 'HorizontalAlignment', 'left', 'Interpreter', 'none', ...
        'FontName', 'helvetica', 'FontWeight', 'normal', 'Color', 'r');
    
    subplot(1, 2, 2);
    imagesc(z);
    set(gca, 'YDir', 'normal');
    colormap(cmap);
    xlabel('x / pixel', ...
        'FontName', 'helvetica', 'FontWeight', 'normal', ...
        'FontUnits', 'pixels', 'FontSize', 14);
    ylabel('y / pixel', ...
        'FontName', 'helvetica', 'FontWeight', 'normal', ...
        'FontUnits', 'pixels', 'FontSize', 14);
    title(sprintf('Straightened filament\n'), ...
        'interpreter', 'none', ...
        'FontName', 'helvetica', 'FontWeight', 'bold', ...
        'FontUnits', 'pixels', 'FontSize', 14);
    
    handedness = [];
    while isempty(handedness)
        [~, ~, button] = ginput(1);
        
        
        if button == 108
            % (l)eft
            handedness = 'left';
        elseif button == 114
            % (r)ight
            handedness = 'right';
        end
    end
    fprintf('Filament twist handedness appears to be: %s-handed\n', handedness);
    %}

    handedness = 'left';
    fprintf('Filament twist handedness assumed to be: %s-handed\n', handedness);
else
    % Use user supplied handedness
    fprintf('Filament twist handedness is set to: %s-handed\n', handedness);
end

%{
fprintf('\nPlease look for the clearest cross-peak in the 2D-FFT image!\n');
fprintf('When ready, press return to mark the spot!\n');
pause
[~, yy] = ginput(1);
centre = ceil(size(psd2d, 1)/2);
clf;
hh1 = plot(psd2d(centre, :), 'r-');
hold('on');
hh2 = plot(psd2d(round(yy), :), 'k-');
legend([hh1 hh2], ...
    'Equatorial centre', ...
    'Selected horizontal line, please select first off-centre peak on this line', ...
    'Location', 'southeast');
[xx, ~] = ginput(1);
crosspeak1 = [xx, yy];
crosspeak2 = [centre-xx+centre, centre-yy+centre];
% Detect quadrant
if (crosspeak1(1) > centre && crosspeak1(2) > centre) || ...
        (crosspeak1(1) < centre && crosspeak1(2) < centre)
    handedness = 'left';
else
    handedness = 'right';
end
    
fprintf('Filament twist handedness appears to be: %s-handed\n', handedness);
%}



% User input needed
% Finally deciding on symmetry based on user comparison of 2D-FFT of
% simulated images vs the original data

%
% Symmetry evaluation if needed
%
% Also no longer done this way. User should need to suply the value
if ~exist('symmetry', 'var') || isempty(symmetry) || length(symmetry) > 1
    
    %{
    if ~exist('symmetry', 'var') || isempty(symmetry)
        symmetry_test = 1:max_sym;
    else
        % Use input range
        symmetry_test = symmetry;
        max_sym = max(symmetry_test);
    end
    
    % The 2D-FFT of the data
    % Reusing fig 5
    figure(5);
    set(5, 'OuterPosition', [251 151 960 960]);
    clf;
    colormap(cmap);
    ax = zeros(1, max_sym+1);
    ax2 = zeros(1, max_sym+1);
    ax(1) = subplot(2, max_sym+1, 1);
    imagesc(fxx, fxx, psd2d);
    set(gca, 'YDir', 'normal');
    xlabel('x / Å^-^1', ...
        'FontName', 'helvetica', 'FontWeight', 'normal', ...
        'FontUnits', 'pixels', 'FontSize', 14);
    ylabel('y / Å^-^1', ...
        'FontName', 'helvetica', 'FontWeight', 'normal', ...
        'FontUnits', 'pixels', 'FontSize', 14);
    %title('Data');
    tt = get(gca, 'Title');
    set(tt, 'HorizontalAlignment', 'left');
    %colormap(ax(1), cmap);
    %hold('on');
    %plot([crosspeak1(1) crosspeak2(1)], [crosspeak1(2) crosspeak2(2)], 'k+-', 'LineWidth', 2);
    %
    ax2(1) = subplot(2, max_sym+1, max_sym+2);
    imagesc(z);
    set(gca, 'YDir', 'normal');
    %colormap(ax2(1), img.cMap);
    title('Data');
    clim = get(gca, 'CLim');
    xlabel('x / pixel', ...
        'FontName', 'helvetica', 'FontWeight', 'normal', ...
        'FontUnits', 'pixels', 'FontSize', 14);
    ylabel('y / pixel', ...
        'FontName', 'helvetica', 'FontWeight', 'normal', ...
        'FontUnits', 'pixels', 'FontSize', 14);
    
    
    for symmetry = symmetry_test
        %{
    % Set one single spline setting for scouting, not tested and optimised
    % here to quickly test the different symmetry settings
    if strcmp(method, 'spap2')
        n_spline_pieces = round(periodicity/2);
    else
        % 'pchip'
        n_spline_pieces = -1;
    end
    %[xm, ym, zm] = Filament_Model(xci(:, ceil(n_lines/2)), yci(:, ceil(n_lines/2)), zci(:, ceil(n_lines/2)), centre_axis, ...
    %    symmetry, handedness, periodicity, n_spline_pieces, w(:, ceil(n_lines/2)));
    %[xm, ym, zm] = Filament_Model(xci, yci, zci, centre_axis, ...
    %    symmetry, handedness, periodicity, n_spline_pieces, w);
        %}
        
        % Make the 3D model, quick without additional y-refine
        [xm, ym, zm] = FilamentModel(xc(:, selected_lines), ...
            yc(:, selected_lines), zc(:, selected_lines), ...
            centre_axis, symmetry, handedness, periodicity, peaks_loc, ...
            w(:, selected_lines), px, smoothness, 1);
        
        % Get the top part of the model
        %[xm_top, ym_top, zm_top] = Filament_Model_Top(xm, ym, zm, z_max/2);
        [xm_top, ym_top, zm_top] = FilamentModelTop(xm, ym, zm, z_max/2);
        
        % Quick image simulation for symmetry check
        % precision set to -1 for speed, only do the dialation with no refining
        % of the coordinates
        [~, ~, zsim] = AFM_Image_Sim(r_tip, a_tip, [], -1, ...
            'data3Dmodel', {x; y}, 0, xm_top+centre_axis(1), ym_top, zm_top);
        
        % 2D-FFT of the simulated images
        zsim_trim = zsim(trim+1:end-trim, :);
        % zsim_trim must have odd number of rows because it will have odd number of
        % columns, after padding both dimensions will be odd, same as for data
        if mod(size(zsim_trim, 1), 2) == 0
            zsim_trim = zsim_trim(1:end-1, :);
        end
        z_size = size(zsim_trim);
        % Pad to square
        zpad_size = [max(z_size), max(z_size)];
        zsim_pad = padarray(zsim_trim, (zpad_size-z_size)/2);
        
        psd2d_test = fftshift(fft2(zsim_pad));
        psd2d_test = log10(abs(psd2d_test).^2);
        
        % Plotting all of the symmetries
        ax(symmetry+1) = subplot(2, max_sym+1, symmetry+1);
        imagesc(fxx, fxx, psd2d_test);
        set(gca, 'YDir', 'normal');
        xlabel('x / Å^-^1', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', ...
            'FontUnits', 'pixels', 'FontSize', 14);
        %ylabel('y / Å^-^1', ...
        %    'FontName', 'helvetica', 'FontWeight', 'normal', ...
        %    'FontUnits', 'pixels', 'FontSize', 14);
        %title(sprintf('sym = %g', symmetry));
        %colormap(ax(symmetry+1), cmap);
        %hold('on');
        %plot([crosspeak1(1) crosspeak2(1)], [crosspeak1(2) crosspeak2(2)], 'k+-', 'LineWidth', 2);
        %{
        ax2(symmetry+1) = subplot(2, max_sym+1, max_sym+1+symmetry+1);
        imagesc(zsim);
        set(gca, 'YDir', 'normal');
        colormap(ax2(symmetry+1), img.cMap);
        %}
        
        ax2(symmetry+1) = subplot(2, max_sym+1, max_sym+2+symmetry);
        imagesc(zsim);
        set(gca, 'YDir', 'normal');
        caxis(clim);
        %colormap(ax2(1), img.cMap);
        title(sprintf('Symmetry = %g', symmetry));
        xlabel('x / pixel', ...
            'FontName', 'helvetica', 'FontWeight', 'normal', ...
            'FontUnits', 'pixels', 'FontSize', 14);
        %ylabel('y / pixel', ...
        %    'FontName', 'helvetica', 'FontWeight', 'normal', ...
        %    'FontUnits', 'pixels', 'FontSize', 14);
        
    end
    %linkaxes(ax, 'xy');
    %linkaxes(ax2, 'xy');
    %fprintf('\nWhen ready, press return to select symmetry by choosing closest match to the data!\n');
    %pause;
    subplot(2, max_sym+1, 1);
    title({'{\color{red}Select symmetry by choosing closest match to the data}'; ''; ''}, 'FontSize', 16);
    
    
    currentfig = gcf;
    subplotClicked = 2*(max_sym+1);
    while ceil(subplotClicked/2) == max_sym+1
        % While clicked subplot is not the data plots
        [~, ~] = ginput(1);
        axesClicked = gca;
        allAxes = findobj(currentfig.Children, 'Type', 'axes');
        subplotClicked = find(allAxes == axesClicked);
        % The indexes are in reverse order of subplots added
    end
    %title('Selected');
    symmetry = max_sym+1-ceil(subplotClicked/2);
    fprintf('Filament screw axis symmetry estimate appears to be: %g\n', symmetry);
    %}

    symmetry = 2;
    fprintf('Filament cross-section symmetry assumed to be: %g\n', symmetry);
else
    % Use user supplied symmetry
    fprintf('Filament cross-section symmetry estimate is set to: %g\n', symmetry);
end



%
% Finalise 3D model
%
%close(1);
%close(2);
%close(3);
%close(4);
tic;

%{
figure;
if strcmp(method, 'spap2')
    n_spline_pieces = round(symmetry*periodicity/2);
else
    % 'pchip'
    n_spline_pieces = -1;
end
%}
%[xm, ym, zm, xarea] = Filament_Model(xci, yci, zci, centre_axis, ...
%    symmetry, handedness, periodicity, n_spline_pieces, w);
%xarea_mean = mean(xarea);
%[xm_top, ym_top, zm_top] = Filament_Model_Top(xm, ym, zm, z_max/2);



% Grid search for better r_tip estimate

if refineTip
    %n_r_tip_test = 10;
    n_r_tip_test = 5;
    i_r_tip_test = 1:n_r_tip_test+1;
    %r_tip_test(1:10) = linspace(1, round(2*r_tip), 10);
    r_tip_test = zeros(n_r_tip_test+1, 1);
    %r_tip_test(1:n_r_tip_test) = linspace(0.5*r_tip, 2*r_tip, n_r_tip_test);
    r_tip_test(1:n_r_tip_test) = [r_tip-2 r_tip-1 r_tip r_tip+1 r_tip+2];

    fprintf('\nRefining tip radius estimation...\n');
else
    n_r_tip_test = 1;
    i_r_tip_test = 1;
    r_tip_test = r_tip;
    fprintf('\nFinalising 3D model...\n');
end
    
d_rmsd = zeros(size(r_tip_test))+Inf;
d_corr = zeros(size(r_tip_test))+Inf;
%d_xcorr2 = zeros(size(r_tip_test))+Inf;

%{
figure(5);
set(5, 'OuterPosition', [51 51 1440 960]);
clf;
ax1 = subplot(6, 2, 1);
imagesc(imrotate(z, 90));
colormap(ax1, img.cMap);
set(gca, 'YDir', 'normal');
set(gca, 'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1]);
clim = get(gca, 'CLim');
title(sprintf('Image data'));
%}

if n_r_tip_test > 1
    hh_tip = waitbar(0, 'Refining tip radus estimate');
    set(hh_tip, 'OuterPosition', get(hh_tip, 'OuterPosition')+[0 100 0 0]);
end

for aa = i_r_tip_test

    if numel(i_r_tip_test) > 1
        waitbar(aa/(n_r_tip_test+1), hh_tip, 'Refining tip radus estimate');
    end

    % For the final model
    if aa == n_r_tip_test+1
        
        idx = find(d_rmsd == min(d_rmsd), 1);
        if idx == 1
            idx = 2;
        elseif idx == n_r_tip_test
            idx = n_r_tip_test-1;
        end
        rr = (r_tip_test(idx-1):0.01:r_tip_test(idx+1))';
        dd = spline(r_tip_test(1:end-1), d_rmsd(1:end-1), rr);
        
        r_tip_test(aa) = rr(find(dd == min(dd), 1));
        r_tip = r_tip_test(end);
        
        fprintf('Tip radius refined estimate: r_tip = %g nm\n', r_tip);
        fprintf('Finalising 3D model...\n');
        
        % Test code
        %{
        clf;
        subplot(2, 1, 1);
        plot(rr, dd, '+-');
        subplot(2, 1, 2);
        dd = spline(r_tip_test(1:10), d_corr(1:10), rr);
        plot(rr, dd, '+-');
        %}
    end
    
    
    img = img.SetTipModel('rounded_cone', [r_tip_test(aa), a_tip]);
    [xc, yc, zc, unreliability] = CalcContactPoints(x, y, z, img.zT{1}, img.zT{2}, 20);
    centre_axis = [mean(xc(:)) z_max/2];
    
    % Make the 3D model
    [xm, ym, zm, theta_m, rho_m, x_untwisted, z_untwisted, twist_angle_y] = ...
        FilamentModel(xc(:, selected_lines), ...
        yc(:, selected_lines), zc(:, selected_lines), ...
        centre_axis, symmetry, handedness, periodicity, peaks_loc, ...
        w(:, selected_lines), px, smoothness, refine);
    
    % Get the top part of the model
    %[xm_top, ym_top, zm_top] = Filament_Model_Top(xm, ym, zm, z_max/2);
    [xm_top, ym_top, zm_top] = FilamentModelTop(xm, ym, zm, z_max/2);
    
    % Simulate the AFM image for comparison
    [~, ~, zsim] = SimHeightImageData(r_tip_test(aa), a_tip, [], refine, ...
        'data3Dmodel', {x; y}, 2, xm_top+centre_axis(1), ym_top, zm_top);
    
    % Trim not needed anymore
    %trim = (size(z, 1)-size(zm, 1))/2;
    %zsim([1:trim end-trim+1:end], :) = 0;
    %zf_sim = zsim(trim+1:end-trim, :);
    
    d_rmsd(aa) = sqrt(sum((zsim(:)-z(:)).^2)./numel(z));
    d_corr(aa) = pdist([zsim(:)'; z(:)'], 'correlation');
    %d_xcorr2_all = xcorr2(zsim, z);
    %d_xcorr2(aa) = d_xcorr2_all(size(z, 1), size(z, 2));

    %{
    ax2 = subplot(6, 2, aa+1);
    imagesc(imrotate(zsim, 90));
    colormap(ax2, img.cMap);
    set(gca, 'YDir', 'normal');
    set(gca, 'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1]);
    caxis(clim);
    title(sprintf('r_t_i_p = %g, d_s_s_r = %g', r_tip_test(aa), d_rmsd(aa)));
    %}
end
if n_r_tip_test > 1
    close(hh_tip);
end

fprintf('Image-model RMSD: d_rmsd = %g nm\n', d_rmsd(end));
fprintf('Image-model correlation distance: d_corr = %g \n\n', d_corr(end));

% Cross section calculations
[xm_untwisted, zm_untwisted] = pol2cart(theta_m, rho_m);
% 95% bounds for xsection = [2.5 procentile, mean, 97.5 procentile]
[xm_xsection, zm_xsection] = pol2cart(theta_m(1:3, :), ...
    [prctile(rho_m, 2.5, 1); mean(rho_m, 1); prctile(rho_m, 97.5, 1)]);
xarea = 0.5*sum(diff(theta_m, [], 2).*rho_m(:, 2:end).^2, 2);
xarea_mean = mean(xarea);


% Kernel density bandwidth estimation for x/z cross-section point cloud 
% The default estimator for ksdensity is (probably) Scott’s or Silverman’s 
% rule and it seems to overestimate the bandwidth so use a Non-parametric 
% bandwidth selection for ksdensity instead.
% Using the method from Kernel Density Estimator kde.m
% https://www.mathworks.com/matlabcentral/fileexchange/14034-kernel-density-estimator
% Kernel density estimation via diffusion
% Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
% Annals of Statistics, Volume 38, Number 5, pages 2916-2957
% doi:10.1214/10-AOS799

bandwidth = [0 0];
for aa = [1 2]

    if aa == 1
        data = x_untwisted(:);
    elseif aa == 2
        data = z_untwisted(:);
    end

    % Compute the optimal bandwidth^2 using the referenced method

    % Recomended number of mesh points
    n = 2^14;
    range = max(data)-min(data);
    n_data = length(unique(data));
    ii = (1:n-1)'.^2;
    dx = range/(n-1);
    xmesh = min(data)+(0:dx:range);
    %initial_f = histc(data, xmesh)/n_data;
    initial_f = [histcounts(data, xmesh)'/n_data; 0];
    initial_f = initial_f/sum(initial_f);

    % computes the discrete cosine transform of the column vector data
    [nrows, ~]= size(initial_f);
    % Compute weights to multiply DFT coefficients
    weight = [1; 2*(exp(-1i*(1:nrows-1)*pi/(2*nrows))).'];
    % Re-order the elements of the columns of x
    a = [initial_f(1:2:end, :); initial_f(end:-2:2, :)];
    % Multiply FFT by weights:
    a = real(weight.*fft(a));

    a2=(a(2:end)/2).^2;

    % Use fzero to solve the equation t=zeta*gamma^[5](t)
    t_star = find_root(@(t)fixed_point(t,n_data,ii,a2),n_data);
    bandwidth(aa) = sqrt(t_star)*range;

end
% Use a uniform bandwidth across x/z
xz_bandwidth = min(bandwidth);


% Subroutines for bandwidth estimator
    function  out = fixed_point(t, N, I, a2)
        % this implements the function t-zeta*gamma^[l](t)
        l=7;
        ff=2*pi^(2*l)*sum(I.^l.*a2.*exp(-I*pi^2*t));
        for s=l-1:-1:2
            K0=prod(1:2:2*s-1)/sqrt(2*pi);  const=(1+(1/2)^(s+1/2))/3;
            time=(2*const*K0/N/ff)^(2/(3+2*s));
            ff=2*pi^(2*s)*sum(I.^s.*a2.*exp(-I*pi^2*time));
        end
        out=t-(2*N*sqrt(pi)*ff)^(-2/5);
    end

    function t = find_root(f, N)
        % try to find smallest root whenever there is more than one
        N=50*(N<=50)+1050*(N>=1050)+N*((N<1050)&(N>50));
        tol=10^-12+0.01*(N-50)/1000;
        flag=0;
        while flag==0
            try
                t=fzero(f,[0,tol]);
                flag=1;
            catch
                tol=min(tol*2,.1); % double search interval
            end
            if tol==.1 % if all else fails
                t=fminbnd(@(x)abs(f(x)),0,.1); flag=1;
            end
        end
    end



% Save result
% Straightened x, y, z values
f.xf = x;
f.yf = y;
%f.zf = z;
f.zf = z_orig;
% Corrected x, y, z values
f.xc = xc;
f.yc = yc;
f.zc = zc;
f.unreliability = unreliability;

% 3D reconstruction parameter settings
%{
if ~isprop(f, 'hf_3DModel')
    addprop(f, 'hf_3DModel');
end
%}
% Old structure
%{
f.hf3DModel = struct(...
    'r_max', r_max, ...
    'periodicity', periodicity, ...
    'r_tip', r_tip, ...
    'a_tip', a_tip, ...
    'n_lines', n_lines, ...
    'handedness', handedness, ...
    'symmetry', symmetry, ...
    'xm', xm, ...
    'ym', ym, ...
    'zm', zm, ...
    'n_spline_pieces', n_spline_pieces, ...
    'xarea', xarea, ...
    'xarea_mean', xarea_mean, ...
    'zf_sim', zf_sim ...
    );
%}
f.helicalFilament3DModel = struct(...
    'r_max', r_max, ...
    'z_max', z_max, ...
    'z_min', z_min, ...
    'periodicity', periodicity, ...
    'peaks_loc', peaks_loc, ...
    'r_tip', r_tip, ...
    'a_tip', a_tip, ...
    'n_lines', n_lines, ...
    'smoothness', smoothness, ...
    'refine', refine, ...
    'handedness', handedness, ...
    'symmetry', symmetry, ...
    'x_untwisted', x_untwisted, ...
    'z_untwisted', z_untwisted, ...
    'twist_angle_y', twist_angle_y, ...
    'xm', xm, ...
    'ym', ym, ...
    'zm', zm, ...
    'xm_untwisted', xm_untwisted, ...
    'zm_untwisted', zm_untwisted, ...
    'theta_m', theta_m, ...
    'rho_m', rho_m, ...
    'xarea', xarea, ...
    'xarea_mean', xarea_mean, ...
    'xm_xsection', xm_xsection, ...
    'zm_xsection', zm_xsection, ...
    'zsim', zsim, ...
    'r_tip_test', r_tip_test, ...
    'd_rmsd', d_rmsd, ...
    'd_corr', d_corr, ...
    'xz_bandwidth', xz_bandwidth ...
    );


currentFilament = f;
%if filamentIndex == 0
%    img.features.lastFilament = f;
%else
%    img.features.filaments(filamentIndex) = f;
%end
img.features.filaments(filamentIndex) = f;

% Display final result
%DispFilament3DModel(img, filamentIndex, [], [], 0, 500, 4);
DispFilament3DModel(img, filamentIndex);

toc;
fprintf('Done!\n\n');



end