function [hh] = DispFilament3DModel(img, filamentIndex, clim_img, cc_rho, y_start, y_length, refine)

%
% DESCRIPTION
% – Display the model made with MakeHelicalFilamentModel.m along with the 
% straightened filament image and visualisation of other model-depedent 
% info.
% – Also produces a back-simulated image from the model for comparison and
% a cross-section density visualisation
% – Part of Trace_y by WFX
%
% USAGE
% Standard method usage with filament number as input to get 500 nm long
% model
% >> img.DispFilament3DModel(filamentIndex);
%   Specify the colour map coding for image and for the model based on
% distance to the central axis
% >> img.DispFilament3DModel(filamentIndex, clim_img, cc_rho);
%   Specify the length of the model
% >> DispFilament3DModel(filamentIndex, [], [], y_start, y_length)
%
% INPUTS
% img  –  AFMimage object. Implicit with dot notation.
% filamentIndex  –  The index number of the filament.
% clim_img  –  The colour scale limits of the image, [min, max].
% cc_rho  –  The colour scale limits of the model, [min, max].
% y_start  –  The left end of the model shown, default is 0 (left end).
% y_length  –  The length of the model shown, default is 500 nm.
% refine  –  How smooth the surface of the model look, higher means
%   smoother, default is 4 (4 points per pixel)
%
% OUTPUTS
% Plot of the image and the 3D model as well as simulated image from the
% model etc.
%
% DEPENDENCIES
% – Statistics and Machine Learning Toolbox
% – Method for Trace_y's @AFMimage/ object, including StraightenFilament.m
%
% AUTHORS
% Wei-Feng Xue
%
% HISTORY
% 2021.02  –  First draft
% 2024.09  –  WFX Trace_y update including updates to the visualisation
% code.
%



%
% Defaults
%

if ~exist('refine', 'var') || refine == 0
    refine = 4;
end

if ~exist('clim_img', 'var') || isempty(clim_img)
    clim_img = img.cScale;
end

if ~exist('cc_rho', 'var') || isempty(cc_rho)
    rho_m = img.features.filaments(filamentIndex).helicalFilament3DModel.rho_m;
    %cc_rho = [floor(min(rho_m(:))) max(ceil(max(rho_m(:))), 2*floor(min(rho_m(:))))];
    %cc_rho = [0 clim_img(2)/2];
    cc_rho = [0 max(ceil(max(rho_m(:))), 2*floor(min(rho_m(:))))];
end

if ~exist('y_start', 'var') || isempty(y_start)
    y_start = 0;
end
if ~exist('y_length', 'var') || isempty(y_length)
    y_length = min(500, img.xResolution*(size(img.features.filaments(filamentIndex).zf, 1)-1));
end
y_min = y_start;
y_max = y_start+y_length;
if y_max > img.xResolution*(size(img.features.filaments(filamentIndex).zf, 1)-1)
    y_max = img.xResolution*(size(img.features.filaments(filamentIndex).zf, 1)-1);
end


% Update old file version
% This should now be handled by @AFMImage/UpdateVersion.m
%{
if isprop(img.features.filaments(filamentIndex), 'hf3DModel') && isempty(img.features.filaments(filamentIndex).helicalFilament3DModel)
    fprintf('\nOld format, re-make 3D model...')
    img.MakeHelicalFilamentModel(filamentIndex, ...
        img.features.filaments(filamentIndex).hf3DModel.handedness, ...
        img.features.filaments(filamentIndex).hf3DModel.symmetry, ...
        true, ...
        img.features.filaments(filamentIndex).hf3DModel.smoothness, ...
        img.features.filaments(filamentIndex).hf3DModel.refine);
end
%}


%
% Get data
%

% Get image data for the filament
% Straightened image
%px = img.xResolution;
z = img.StraightenFilament(filamentIndex);
z = imrotate(z, -90);
% Assumes as always square pixels everywhere
x = 0:img.xResolution:img.xResolution*(size(z, 2)-1);
x = x-img.xResolution*(size(z, 2)-1)/2;
y = (0:img.xResolution:img.xResolution*(size(z, 1)-1))';


zsim = img.features.filaments(filamentIndex).helicalFilament3DModel.zsim;
xsim = img.features.filaments(filamentIndex).xf-img.xResolution*(size(zsim, 2)-1)/2;
ysim = img.features.filaments(filamentIndex).yf;
xm = img.features.filaments(filamentIndex).helicalFilament3DModel.xm;
ym = img.features.filaments(filamentIndex).helicalFilament3DModel.ym;
zm = img.features.filaments(filamentIndex).helicalFilament3DModel.zm;
rho_m = img.features.filaments(filamentIndex).helicalFilament3DModel.rho_m;
theta_m = img.features.filaments(filamentIndex).helicalFilament3DModel.theta_m;
x_untwisted = img.features.filaments(filamentIndex).helicalFilament3DModel.x_untwisted;
z_untwisted = img.features.filaments(filamentIndex).helicalFilament3DModel.z_untwisted;
%xm_xsection = img.features.filaments(filamentIndex).helicalFilament3DModel.xm_xsection;
%zm_xsection = img.features.filaments(filamentIndex).helicalFilament3DModel.zm_xsection;

%zz = img.StraightenFilament(filament);
%yy = px*((1:size(zz, 2))-1);
%xx = px*((1:size(zz, 1))-1);



%
% Display final result
%

h_fig = figure(12);
%set(6, 'OuterPosition', [1 51 1440 960]);
%set(6, 'OuterPosition', [251 151 1000 800]);
set(h_fig, 'NumberTitle', 'off', 'Name', 'Trace_y – DispFilament3DModel');
h_fig.OuterPosition(3:4) = [1100 880];
h_fig.OuterPosition(1:2) = h_fig.OuterPosition(1:2)+[25 -25];

%drawnow;
screenSize = get(0, 'ScreenSize');
if h_fig.OuterPosition(2)+h_fig.OuterPosition(4) > screenSize(4)
    h_fig.OuterPosition(2) = screenSize(4)-h_fig.OuterPosition(4)-10;
end
clf;

tl = tiledlayout(15, 6, 'Padding', 'compact', 'TileSpacing', 'tight');
title(tl, '   ', ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontSize', 16);
tx = xlabel(tl, sprintf('%s, filament %d', img.dataFile, filamentIndex), ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontSize', 16);
set(tx, 'FontWeight', 'normal', 'FontAngle', 'italic', 'FontSize', 16);


% The image
h_ax1 = nexttile(1, [4 3]);
%ax1 = subplot(5, 4, 1:4);
imagesc(y, x, imrotate(z, 90));
hold('on');
colormap(h_ax1, img.cMap);
set(h_ax1, 'YDir', 'normal');
%axis('equal');
set(h_ax1, 'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1]);
set(h_ax1, 'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 12);
set(h_ax1, 'CLim', clim_img);
%set(h_ax1, 'XTickMode', 'manual');
ylims = h_ax1.YLim;
plot([y_min y_min], ylims, 'c-');
plot([y_max y_max], ylims, 'c-');
h_ax1.XLim = [y_min y_max];
%clim_img = get(gca, 'CLim');
%ylabel('x / nm', ...
%    'FontName', 'helvetica', 'FontWeight', 'normal', ...
%    'FontUnits', 'pixels', 'FontSize', 14);
xlabel('y / nm', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'pixels', 'FontSize', 14);
title(sprintf('Image Data'), ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontUnits', 'points', 'FontSize', 14);

h_ax2 = nexttile(4, [4 3]);
%ax2 = subplot(5, 4, 5:8);
imagesc(ysim, xsim, imrotate(zsim, 90));
hold('on');
colormap(h_ax2, img.cMap);
set(h_ax2, 'YDir', 'normal');
%axis('equal');
set(h_ax2, 'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1]);
set(h_ax2, 'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 12);
%set(h_ax2, 'YTickMode', 'manual');
set(h_ax2, 'YLim', h_ax1.YLim);
set(h_ax2, 'CLim', clim_img);
ylims = h_ax2.YLim;
plot([y_min y_min], ylims, 'c-');
plot([y_max y_max], ylims, 'c-');
h_ax2.XLim = [y_min y_max];
%ylabel('x / nm', ...
%    'FontName', 'helvetica', 'FontWeight', 'normal', ...
%    'FontUnits', 'pixels', 'FontSize', 14);
xlabel('y / nm', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 16);
title(sprintf('Simulation from model'), ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontUnits', 'pixels', 'FontSize', 14);

linkaxes([h_ax1 h_ax2], 'xy');



%
% The model
%

% Refine the 3D display model by interpolating with pchip
delta_y = ym(2, 1)/refine;
%yi = (min(ym(:)):delta_y:max(ym(:)))'*ones(1, size(ym, 2));
yi = (y_min:delta_y:y_max)'*ones(1, size(ym, 2));
xi = zeros(size(yi));
zi = zeros(size(yi));

for aa = 1:size(ym, 2)
    [theta, d] = cart2pol(xm(:, aa), zm(:, aa));
    di = pchip(ym(:, aa), d, yi(:, aa));

    % theta is periodical so need to check changes
    d_theta = [0; diff(theta)];
    idx = find(d_theta < -pi);
    d_theta(idx) = d_theta(idx)+2*pi;
    idx = find(d_theta > pi);
    d_theta(idx) = d_theta(idx)-2*pi;
    theta = theta(1)+cumsum(d_theta);

    thetai = pchip(ym(:, aa), theta, yi(:, aa));
    [xi(:, aa), zi(:, aa)] = pol2cart(thetai, di);
end

%ax_m = subplot(5, 4, 9:12);
h_ax_m = nexttile(25, [2 6]);
cc = sqrt(xi.^2+zi.^2);
surf(xi, yi, zi, cc, 'EdgeColor', 'none', 'AmbientStrength', 0.5);
%surf(xm, ym, zm, cc, 'EdgeColor', 'none', 'AmbientStrength', 0.4);
%cc_rho = [min(rho_m(:)) max(rho_m(:))];
%cc_rho = [4 10];
%cc_rho = [floor(min(rho_m(:))) max(ceil(max(rho_m(:))), 2*floor(min(rho_m(:))))];
%caxis(cc_rho);
set(h_ax_m, 'CLim', cc_rho);
axis('equal');
camlight('left');
%camlight('right');
colormap(h_ax_m, 'parula');
set(gca,'visible', 'off');
view(90, 90);
ylim([y_min y_max]);
%title(sprintf('3D Model'), ...
%    'interpreter', 'none', ...
%    'FontName', 'helvetica', 'FontWeight', 'bold', ...
%    'FontUnits', 'points', 'FontSize', 14);



%
% Helical polar coordinates
%

%ax5 = subplot(5, 4, [13 17]+2);
h_ax5 = nexttile(37, [3 6]);
%imagesc(rad2deg(theta_m(1, :)), ym(:, 1), rho_m, cc_rho);
imagesc(ym(:, 1)', rad2deg(theta_m(1, :))', imrotate(rho_m, 90), cc_rho);
set(h_ax5, 'YDir', 'normal');
set(h_ax5, 'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 12);
set(h_ax5, 'YTickMode', 'manual', 'YTick', [-180 -90 0 90 180]);
xlim([y_min y_max]);
colormap(h_ax5, 'parula');
%colorbar;
h_cb = colorbar(h_ax5, 'eastoutside');
set(h_cb, 'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontSize', 12)
set(h_cb.Label, 'String', 'rho / nm', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontSize', 16)
xlabel('y / nm', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 16);
ylabel('theta / deg', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 16);
title(sprintf('Model'), ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontUnits', 'points', 'FontSize', 14);

linkaxes([h_ax1 h_ax2 h_ax5], 'x');



%
% Fourier space
%
% Lifted out code originally from MakeHelicalFilamentModel.m

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
%cmap = jet(512);

%ax3 = subplot(5, 4, [13 17]);
h_ax3 = nexttile(55, [6 2]);
imagesc(fxx, fxx, psd2d);
set(h_ax3, 'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 12);
%set(h_ax3, 'YTickMode', 'manual');
%set(h_ax3, 'YDir', 'normal');
set(h_ax3, 'yDir', 'normal', ...
    'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1], ...
    'box', 'on');
clims = get(h_ax3, 'CLim');
colormap(h_ax3, 'jet');
xlabel('x / Å^-^1', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'pixels', 'FontSize', 16);
ylabel('y / Å^-^1', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'pixels', 'FontSize', 16);
title(sprintf('Data 2D-FFT'), ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontUnits', 'pixels', 'FontSize', 14);

%ax4 = subplot(5, 4, [13 17]+1);
h_ax4 = nexttile(57, [6 2]);
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

psd2d = fftshift(fft2(zsim_pad));
psd2d = log10(abs(psd2d).^2);

imagesc(fxx, fxx, psd2d);
%set(h_ax4, 'YDir', 'normal');
set(h_ax4, 'yDir', 'normal', ...
    'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1], ...
    'box', 'on');
colormap(h_ax4, 'jet');
%caxis(clims);
set(h_ax4, 'CLim', clims);
%set(h_ax4, 'YTickMode', 'manual');
set(h_ax4, 'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 12);
xlabel('x / Å^-^1', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'pixels', 'FontSize', 16);
ylabel('y / Å^-^1', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'pixels', 'FontSize', 16);
title(sprintf('Simulation 2D-FFT'), ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontUnits', 'pixels', 'FontSize', 14);

linkaxes([h_ax3 h_ax4], 'xy');



%
% Cross-section
%

%ax6 = subplot(5, 4, [13 17]+3);
h_ax6 = nexttile(59, [6 2]);

%[xm_untwisted, zm_untwisted] = pol2cart(theta_m, rho_m);
% xsection = [2.5 procentile, mean, 97.5 procentile], essentiall 95% bound
%[xm_xsection, zm_xsection] = pol2cart(theta_m(1:3, :), ...
%    [prctile(rho_m, 2.5, 1); mean(rho_m, 1); prctile(rho_m, 97.5, 1)]);

%{
pts = pointCloud([x_untwisted(:) z_untwisted(:) 0.*x_untwisted(:)]);
dist = zeros(numel(x_untwisted), 1);
for aa = 1:numel(x_untwisted)
    [~, dist2] = findNearestNeighbors(pts, pts.Location(aa, :), 2);
    dist(aa) = max(dist2);
end
%}


% Calc ksdensity
bandwidth = img.features.filaments(filamentIndex).helicalFilament3DModel.xz_bandwidth;
%xz_range = max([x_untwisted(:); z_untwisted(:)]);
%xz_range = 5;

% For visualisation
% Cross-section plot lims are set to 10, 15, 20... etc.
if max([x_untwisted(:); z_untwisted(:)]) < 10
    xz_range = 10;
else
    xz_range = ceil(max([x_untwisted(:); z_untwisted(:)])/5)*5;
end

xx = 0:0.1:1.5*xz_range;
xx = unique([xx -xx]);
%xy_range = max(max(x_untwisted(:)), max(z_untwisted(:)))-min(min(x_untwisted(:)), min(z_untwisted(:)));
%xx = min(min(x_untwisted(:)), min(z_untwisted(:)))-xy_range/2:0.1:max(max(x_untwisted(:)), max(z_untwisted(:)))+xy_range/2;
yy = xx;
[xx, yy] = meshgrid(xx, yy);
ff = ksdensity([x_untwisted(:) z_untwisted(:)], [xx(:) yy(:)], 'Bandwidth', bandwidth);
ff = reshape(ff, size(xx));
ff = ff./max(ff(:));
%contour(xx(1, :), yy(:, 1), ff, 0.3:0.1:1);
imagesc(xx(1, :), yy(:, 1), ff);
colormap(h_ax6, 'sky');
%colormap(h_ax6, [0 0 0; hot(511)]);
hold('on');

%plot(xm_xsection', zm_xsection', 'k-', 'LineWidth', 0.5);
%plot(xm_xsection(2, :)', zm_xsection(2, :)', 'k-', 'LineWidth', 2);
%plot(xm_xsection(2, :)', zm_xsection(2, :)', 'k-', 'LineWidth', 1);
plot(0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1);
%set(gca, 'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1]);
set(h_ax6, 'yDir', 'normal', ...
    'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1], ...
    'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1], ...
    'box', 'on');
set(h_ax6, 'XGrid', 'on', 'yGrid', 'on', 'XTick', -50:5:50, 'YTick', -50:5:50);
set(h_ax6, 'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 12);
xlim([-xz_range xz_range]);
ylim([-xz_range xz_range]);
%axis('equal');
h_cb = colorbar(h_ax6, 'eastoutside');
set(h_cb, 'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontSize', 12)
set(h_cb.Label, 'String', 'Relative density', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontSize', 16)
xlabel('x / nm', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'pixels', 'FontSize', 16);
ylabel('z / nm', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'pixels', 'FontSize', 16);
title(sprintf('Cross-section'), ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontUnits', 'pixels', 'FontSize', 14);



if nargout == 1
    hh = gcf;
end


% For zoom-in segments
% Need to lift out somewhere else for summary figure production
%{
if exist('y_start', 'var')

    if ~exist('y_length', 'var') || refine == 0
        y_length = 500;
    end

    % Refine the 3D display model by interpolating with pchip
    delta_y = ym(2, 1)/refine;
    yi = (min(ym(:)):delta_y:max(ym(:)))'*ones(1, size(ym, 2));
    xi = zeros(size(yi));
    zi = zeros(size(yi));

    for aa = 1:size(ym, 2)
        [theta, d] = cart2pol(xm(:, aa), zm(:, aa));
        di = pchip(ym(:, aa), d, yi(:, aa));

        % theta is periodical so need to check changes
        d_theta = [0; diff(theta)];
        idx = find(d_theta < -pi);
        d_theta(idx) = d_theta(idx)+2*pi;
        idx = find(d_theta > pi);
        d_theta(idx) = d_theta(idx)-2*pi;
        theta = theta(1)+cumsum(d_theta);

        thetai = pchip(ym(:, aa), theta, yi(:, aa));
        [xi(:, aa), zi(:, aa)] = pol2cart(thetai, di);
    end


    % Zooms
    y_min = y_start;
    y_max = y_start+y_length;
    z = imrotate(z, 90);

    figure(10);
    clf;
    set(10, 'OuterPosition', [501 301 1201 700]);
    imagesc(y, x, z);
    hold('on');
    set(gca, 'yDir', 'normal', 'YTickMode', 'manual', ...
        'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1], ...
        'box', 'on');
    % Colour map same as in the image
    colormap(gca, img.cMap);
    caxis(clim_img);
    set(gca,'visible','off')
    %xlabel('l / nm', ...
    %    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    %    'FontUnits', 'pixels', 'FontSize', 14);
    xlim([y_min y_max]);
    %colorbar('northoutside');

    hh = [hh gcf];


    figure(11);
    clf;
    set(11, 'OuterPosition', [501 51 1201 250]);
    c = sqrt(xi.^2+zi.^2);

    %{
    twist_angle_y = img.features.filaments(filamentIndex).helicalFilament3DModel.twist_angle_y;
    lines_idx = (size(img.features.filaments(filamentIndex).yc, 2)-...
        img.features.filaments(filamentIndex).helicalFilament3DModel.n_lines)/2+1:...
        (size(img.features.filaments(filamentIndex).yc, 2)-...
        img.features.filaments(filamentIndex).helicalFilament3DModel.n_lines)/2+...
        img.features.filaments(filamentIndex).helicalFilament3DModel.n_lines;
    yc = img.features.filaments(filamentIndex).yc(:, lines_idx);
    [yc, idx] = sort(yc(:));
    twist_angle_y = twist_angle_y(idx);
    d_twist_angle_y = diff(twist_angle_y)./diff(yc);
    d_twist_angle_y = [d_twist_angle_y(1); d_twist_angle_y];

    c = pchip(yc, d_twist_angle_y, yi);

    %c = pchip(ym(:, 1), img.features.filaments(filamentIndex).helicalFilament3DModel.xarea, yi);
    %}

    surf(xi, yi, zi, c, 'EdgeColor', 'none', 'AmbientStrength', 0.5);

    axis('equal');
    caxis(cc_rho);
    %colorbar;
    view(90, 90)
    camlight('left')
    colormap(gca, 'bone');
    set(gca,'visible','off')
    ylim([y_min y_max]);
    %colorbar('northoutside');

    hh = [hh gcf];
end

figure(14);
clf;
% Calc ksdensity (again) with fixed range
%xx = 0:0.1:7;
xx = 0:0.1:xz_range+2;
xx = unique([xx -xx]);
yy = xx;
[xx, yy] = meshgrid(xx, yy);
ff = ksdensity([x_untwisted(:) z_untwisted(:)], [xx(:) yy(:)], 'Bandwidth', bandwidth);
ff = reshape(ff, size(xx));
ff = ff./max(ff(:));
imagesc(xx(1, :), yy(:, 1), ff);
%colorbar;
%colormap(gca, 'parula');
colormap(gca, 'bone');
hold('on');
set(gca, 'yDir', 'normal', ...
    'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1], ...
    'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1], ...
    'box', 'on');
xlabel('x / nm', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'pixels', 'FontSize', 14);
ylabel('z / nm', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'pixels', 'FontSize', 14);
%set(gca,'visible','off')
%plot([0 0], [-7 7], '-w');
%plot([-7 7], [0 0], '-w');
plot(0, 0, '+w', 'MarkerSize', 24, 'LineWidth', 2);
%colorbar('northoutside');

hh = [hh gcf];
%}



end