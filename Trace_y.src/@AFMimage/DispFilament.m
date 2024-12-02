function [hh] = DispFilament(img, filament_number, unit)

%
% DESCRIPTION
% – Display more detailed info on filament traces
% – Part of Trace_y by WFX
%
% USAGE
% Standard method usage with filament number as input
% >> img.DispFilament(filament_number);
%
% INPUTS
% img  –  AFMimage object. Implicit with dot notation
% filament_number  –  The index of the filament to be ploted in the object,
% essentially the index of img.features.filaments.
% unit  –  'px' (default) or 'nm'
%
% OUTPUTS
% An image of strightened fibril is ploted in a figure with information 
% such as 1D-FFT and central line height plots
%
% DEPENDENCIES
% – Method for Trace_y's @AFMimage/ object and uses other methods including 
% StraightenFilament.m and CalcFilamentFFT1D.m
% – Used also by other routines in Trace_y
%
% AUTHORS
% Wei-Feng Xue
%
% HISTORY
% 2019.01  –  First draft
% 2024.09  –  WFX Trace_y update including updates to the visualisation
%



%
% Defaults
%
% Default to last filament if not specified
if ~exist('filament_number', 'var') || isempty(filament_number)
    filament_number = numel(img.features.filaments);
end
%{
if (~exist('filament_number', 'var') || strcmp(filament_number, 'last')) && ~isempty(img.features.lastFilament)
    filament_number = 0;
elseif ~exist('filament_number', 'var') && isempty(img.features.lastFilament)
    filament_number = length(img.features.filaments);
end
%}

% Default to direct pixel length unit
if ~exist('unit', 'var') || isempty(unit)
    unit = 'pixel';
end
% Handling units
switch unit
    case {'pixel' 'px'}
        % 1 pixel/pixel conversion
        px = 1;
        lableText1 = 'length / pixel';
        lableText2 = 'frequency pixel^-^1';
    case 'nm'
        % nm/pixel conversion factor
        px = img.xResolution;
        lableText1 = 'length / nm';
        lableText2 = 'frequency nm^-^1';
end

% Get data
filament = img.features.filaments(filament_number);
zz = img.StraightenFilament(filament_number);
xx = 1:size(zz, 2);
yy = 1:size(zz, 1);
ll = filament.l;
z = filament.z;




% This is not needed anymore with ExploreFeatures.m
% Plot all filaments in the file
%{
figure(1);
img.DispImage;
set(1, 'OuterPosition', [51 51 960 960]);
hold on

for aa = 1:length(img.features.filaments)
    filament = img.features.filaments(aa);
    plot(filament.x, filament.y, 'c-', 'LineWidth', 2);
end

%{
if ~isempty(img.features.lastFilament)
    filament = img.features.lastFilament;
    plot(filament.x, filament.y, 'c-', 'LineWidth', 1);
end
%}

% Plot the selected filament in detail

if filament_number == 0
    filament = img.features.lastFilament;
else
    filament = img.features.filaments(filament_number);
end


% Markers in the image

x = filament.x;
y = filament.y;
l = filament.l;
z = filament.z;

plot(x(1), y(1), 'b+', 'MarkerSize', 20);
plot(x(end), y(end), 'r+', 'MarkerSize', 20);
plot(x, y, 'r-', 'LineWidth', 2);
%}


%
% Data plots
%

h_fig = figure(11);
set(11, 'NumberTitle', 'off', 'Name', 'Trace_y – DispFilament');
%h_fig.OuterPosition(3:4) = [1440 960];
h_fig.OuterPosition(3:4) = [1100 880];

%drawnow;
screenSize = get(0, 'ScreenSize');
if h_fig.OuterPosition(2)+h_fig.OuterPosition(4) > screenSize(4)
    h_fig.OuterPosition(2) = screenSize(4)-h_fig.OuterPosition(4)-10;
end
clf;

tl = tiledlayout(7, 2, 'Padding', 'loose', 'TileSpacing', 'compact');
title(tl, '   ', ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontSize', 16);
tx = xlabel(tl, sprintf('%s, filament %d', img.dataFile, filament_number), ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontSize', 16);
set(tx, 'FontWeight', 'normal', 'FontAngle', 'italic', 'FontSize', 16);


% Straightened filament
nexttile([1 2]);
imagesc(px*xx, px*yy, zz);
hold('on');
set(gca, 'yDir', 'normal', 'YTickMode', 'manual', ...
    'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1], ...
    'box', 'on', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', 'FontSize', 12);

% Colour map same as in the image
colormap(img.cMap);
set(gca, 'CLim', img.cScale);
% Start/end markers
plot(px*1, px*(size(zz, 1)+1)/2, 'b+', 'MarkerSize', 10, 'LineWidth', 2);
plot(px*size(zz, 2), px*(size(zz, 1)+1)/2, 'r+', 'MarkerSize', 10, 'LineWidth', 2);

xlabel(lableText1, ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'pixels', 'FontSize', 14);
title('Straightened filament', 'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontUnits', 'pixels', 'FontSize', 14);
%xlim([0 200]);


% Height profile
nexttile([2 2]);
plot(px*ll, z, 'b-', 'LineWidth', 1);
set(gca, 'XLim', [0 px*ll(end)]);
set(gca, 'YLim', [min(z)-0.1.*abs(min(z)) max(z)+0.1.*abs(max(z))]);
set(gca, 'FontName', 'helvetica', 'FontWeight', 'normal', 'FontSize', 12);
xlabel(lableText1, ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'pixels', 'FontSize', 14);
ylabel(sprintf('z / %s', img.zUnit), ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'pixels', 'FontSize', 14);
title('Height profile', 'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontUnits', 'pixels', 'FontSize', 14);


% Length and EE distance
h_ax = nexttile([2 2]);
plot(px*ll(2:end), filament.lStep(2:end), 'b-', 'LineWidth', 1);
hold('on');
d_rEE = diff([0; filament.rEE]);
if length(d_rEE) > 1
    plot(px*ll(2:end), d_rEE(2:end), 'r-', 'LineWidth', 1);
end
plot(px*ll(2:end), filament.rEE(2:end)./ll(2:end), '-', 'LineWidth', 1);
set(gca, 'FontName', 'helvetica', 'FontWeight', 'normal', 'FontSize', 12);
set(gca, 'XLim', [0 px*ll(end)]);
h_ax.YLim(1) = min(h_ax.YLim(1), 0.98);
h_ax.YLim(2) = h_ax.YLim(2)+0.1*(h_ax.YLim(2)-h_ax.YLim(1));
xlabel(lableText1, ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'pixels', 'FontSize', 14);
%ylabel('Delta / pixels', ...
%    'FontName', 'helvetica', 'FontWeight', 'normal', ...
%    'FontUnits', 'pixels', 'FontSize', 14);
title('Contour length and E-E distance change between points', 'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontUnits', 'pixels', 'FontSize', 14);
legend('Length per pixel / px', 'EE distance per pixel / px', ...
    'EE distance - Length ratio', 'Location', 'southwest', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', 'FontSize', 10);


% 1D PSD and periodicity
[ff, psd] = img.CalcFilamentFFT1D(filament_number);
% Sampling freq
fs = 1/px;
% Nyquist freq
fn = fs/2;

nexttile([2 1]);
plot(ff/px, psd, 'b-', 'LineWidth', 1);
%hold('on');
%ylims = ylim;
%plot([fn fn], ylims, 'r-');
set(gca, 'XLim', [0 fn]);
set(gca, 'FontName', 'helvetica', 'FontWeight', 'normal', 'FontSize', 12);
xlabel(lableText2, ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'pixels', 'FontSize', 14);
ylabel('Power', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'pixels', 'FontSize', 14);
title('1D power spectrum', 'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontUnits', 'pixels', 'FontSize', 14);
nexttile([2 1]);
plot(px*1./ff, psd, 'b-', 'LineWidth', 1);
hold('on');
ylims = ylim;
plot([1/fn 1/fn], ylims, 'r-');
set(gca, 'FontName', 'helvetica', 'FontWeight', 'normal', 'FontSize', 12);
set(gca, 'XLim', [0 px*ll(end)/4]);
xlabel(lableText1, ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'pixels', 'FontSize', 14);
ylabel('Power', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'pixels', 'FontSize', 14);
title('Periodicity', 'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontUnits', 'pixels', 'FontSize', 14);



if nargout == 1
    hh = h_fig;
end



end