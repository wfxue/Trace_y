function [hh] = DispFilamentCS(img, filamentIndex, ncalc, precision)

%
% DESCRIPTION
% – Displays the filament cross-section contact point-cloud density map
% along with its FRC curve and resolution estimate.
% – The resolution estimate is based on the 1/2 bit information criterion,
% which is recommended by van Heel and Schatz, Fourier shell correlation 
% threshold criteria, Journal of Structural Biology 151 (2005) 250–262, 
% https://doi.org/10.1016/j.jsb.2005.05.009
% – Part of Trace_y
%
% USAGE
% Standard method usage
% >> img.DispFilamentCS(filamentIndex);
% Also specify the precision of the estimate
% >> img.DispFilamentCS(filamentIndex, ncalc, precision);
%
% INPUTS
% img  –  AFMimage object. Implicit with dot notation
% filamentIndex  –  The index number of filament to be deleted
% ncalc  –  How many Monty Carlo repeats to do for the FRC estiate. For
% each estimate, the point-cloud is indepedently split randomly in half 
% sets. Optional with default of 25. The FRC estimate seems to converge
% faily quickly
% precision  –  The pixel precision of the density map in nm unit. Optional 
% with default of 1Ä = 0.1 nm.
%
% OUTPUTS
% Figure showing the filament cross-section contact point-cloud density map
% and its FRC curve
%
% DEPENDENCIES
% – Uses Statistics and Machine Learning Matlab toolbox
% – Method for Trace_y's @AFMimage/ object.
%
% AUTHORS
% Wei-Feng Xue, Liisa Lutter
%
% HISTORY
% 2019.09  –  Original FRC_Res.m drafted by LL and edited by WFX.
% Originally the often used, simple but scientifically not valid 0.5
% threshold was used.
% 2023.06  –  WFX edit to the code and to swap for the 1/2 bit information
% criterion threshold recommended by van Heel and Schatz, Fourier shell 
% correlation threshold criteria, Journal of Structural Biology 151 (2005) 
% 250–262, https://doi.org/10.1016/j.jsb.2005.05.009
% 2024.09  –  WFX incorporated FRC and density map visualisation into a
% separate method as part of Trace_y update
%



%
% Defaults
%

if ~exist('ncalc', 'var') || isempty(ncalc)
    ncalc = 25;
end

if ~exist('precision', 'var') || isempty(precision)
    % precision in nm per pixel units, 0.1 is 1 Å/px
    precision = 0.1;
end



%
% Get data
%

% Helical axis aligned contact point cloud
cpc = [img.features.filaments(filamentIndex).helicalFilament3DModel.x_untwisted(:) ...
    img.features.filaments(filamentIndex).helicalFilament3DModel.z_untwisted(:)];
xy_range = max(abs(([cpc(:, 1); cpc(:, 2)])));

ncpc = size(cpc, 1);
bandwidth = img.features.filaments(filamentIndex).helicalFilament3DModel.xz_bandwidth;

% Calc ksdensity with intervals defined by precision
xx = 0:precision:1.5*xy_range;
xx = unique([xx -xx]);
yy = xx;
[xx, yy] = meshgrid(xx, yy);



%
% Set up figure window
%
h_fig = figure(15);
set(h_fig, 'NumberTitle', 'off', 'Name', 'Trace_y – DispFilamentCS');
h_fig.OuterPosition(3:4) = [1200 440];

%drawnow;
screenSize = get(0, 'ScreenSize');
if h_fig.OuterPosition(2)+h_fig.OuterPosition(4) > screenSize(4)
    h_fig.OuterPosition(2) = screenSize(4)-h_fig.OuterPosition(4)-10;
end
clf;


tl = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tl, '   ', ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontSize', 16);
tx = xlabel(tl, sprintf('%s, filament %d', img.dataFile, filamentIndex), ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontSize', 16);
set(tx, 'FontWeight', 'normal', 'FontAngle', 'italic', 'FontSize', 16);


% For visualisation
% Cross-section plot lims are set to 10, 15, 20... etc.
if max([cpc(:, 1) cpc(:, 2)]) < 10
    plot_range = 10;
else
    plot_range = ceil(max([cpc(:, 1); cpc(:, 1)])/5)*5;
end

% Kernel density estimation image for the filament
h_ax1 = nexttile(1, [1 1]);
ff = ksdensity([cpc(:, 1) cpc(:, 2)], [xx(:) yy(:)], 'Bandwidth', bandwidth);
ff = reshape(ff, size(xx));
ff = ff./max(ff(:));
%contour(xx(1, :), yy(:, 1), ff, 0.3:0.1:1);
imagesc(xx(1, :), yy(:, 1), ff);
colormap(h_ax1, 'sky');
hold('on');

%plot(xm_xsection', zm_xsection', 'k-', 'LineWidth', 0.5);
%plot(xm_xsection(2, :)', zm_xsection(2, :)', 'k-', 'LineWidth', 2);
%plot(xm_xsection(2, :)', zm_xsection(2, :)', 'k-', 'LineWidth', 1);
plot(0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1);
%set(gca, 'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1]);
hold('on');
set(h_ax1, 'yDir', 'normal', ...
    'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1], ...
    'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1], ...
    'box', 'on');
set(h_ax1, 'XGrid', 'on', 'yGrid', 'on', 'XTick', -50:5:50, 'YTick', -50:5:50);
set(h_ax1, 'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 12);
xlim([-plot_range plot_range]);
ylim([-plot_range plot_range]);
%axis('equal');
h_cb = colorbar(h_ax1, 'eastoutside');
set(h_cb, 'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontSize', 12)
set(h_cb.Label, 'String', 'Relative density', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontSize', 16)
xlabel('x / nm', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 16);
ylabel('z / nm', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 16);
title(sprintf('Cross-section'), ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontUnits', 'points', 'FontSize', 16);


%h_wb = waitbar(0);
frc_all = 0;
for bb = 1:ncalc
    %waitbar(bb/ncalc, h_wb);

    % Randomly 50/50 split cpc points
    idx = zeros(ncpc, 1);
    idx(randsample(ncpc, floor(ncpc/2))) = 1;
    idx2 = [idx ~idx];

    im_cpc = cell(2, 1);

    % Work independently on each half-set point-cloud
    for aa = 1:2

        x_untwisted = cpc(logical(idx2(:, aa)), 1);
        y_untwisted = cpc(logical(idx2(:, aa)), 2);

        % Create the kernel density estimation
        ff = ksdensity([x_untwisted(:) y_untwisted(:)], [xx(:) yy(:)], ...
            'Bandwidth', bandwidth);
        ff = reshape(ff, size(xx));
        ff = ff./max(ff(:));
        im_cpc{aa} = ff;

        %Test code
        %{
        figure(aa);
        imagesc(xx(1, :), yy(:, 1), ff);
        colorbar;
        %colormap(gca, 'parula');
        colormap(gca, 'bone');
        hold('on');
        set(gca, 'yDir', 'normal', ...
            'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1], ...
            'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1], ...
            'box', 'on');
        set(gca,'visible','off')
        %plot([0 0], [-7 7], '-w');
        %plot([-7 7], [0 0], '-w');
        plot(0, 0, '+w', 'MarkerSize', 24, 'LineWidth', 2);
        %}

    end

    % Old code for FRC, ddited version to use as sub-routine below
    %frc = CalcFRCcpr(cpc);
    %FRC_Res(im_cpc{1}, im_cpc{2}, 0.1, 'FRC2');
    [frc, fx, nx] = FRC_Res(im_cpc{1}, im_cpc{2}, precision);

    % Average
    frc_all = frc+frc_all;
    frc_mean = frc_all./bb;

    %
    % Plot FRC progress
    %
    if bb == 1
    h_ax2 = nexttile(2, [1 2]);
    %plot(fx(1:end-5), frc(1:end-5));
    h_frc = plot(fx, frc_mean, 'b-', 'LineWidth', 1);
    hold('on');
    plot([h_ax2.XLim], [0 0], 'k-', 'LineWidth', 0.5);
    set(h_ax2, 'FontName', 'helvetica', 'FontWeight', 'normal', ...
        'FontUnits', 'points', 'FontSize', 12);
    xlabel('Spatial frequency (1/Å)', ...
        'FontName', 'helvetica', 'FontWeight', 'normal', ...
        'FontUnits', 'points', 'FontSize', 16);
    ylabel('FRC', ...
        'FontName', 'helvetica', 'FontWeight', 'normal', ...
        'FontUnits', 'points', 'FontSize', 16);
    title(sprintf('FRC for cross-section contact-point density map'), ...
        'interpreter', 'none', ...
        'FontName', 'helvetica', 'FontWeight', 'bold', ...
        'FontUnits', 'points', 'FontSize', 16);
    % Nyquist limit is half the sampling rate so don't plot more beyond it
    h_ax2.XLim = [0 0.5/(precision*10)];
    %h_ax2.YLim(2) = 1.1;

    % Progress information
    h_tt = text(0.01, 0.025, sprintf('Estimating FRC... %g %%', 100*bb/ncalc), ...
        'Units', 'normalized', 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom', 'FontSize', 12, ...
        'FontName', 'helvetica', 'FontWeight', 'normal', 'Color', 'k');


    % 1/2-bit threshold curve
    % For FSC, nx should be square root but for FRC, no square root needed
    ndim = 2;
    t_halfbit = (0.2071+1.9102./nx.^(1/(ndim-1)))./(1.2071+0.9102./nx.^(1/(ndim-1)));

    h_th = plot(fx, t_halfbit, '-r', 'LineWidth', 0.5);
    legend([h_frc h_th], {'FRC' '1/2 bit threshold'}, ...
        'Location', 'northeast', 'Orientation', 'vertical', 'FontSize', 12);

    else
        delete(h_frc);
        h_frc = plot(fx, frc_mean, 'b-', 'LineWidth', 1);
        h_tt.String = sprintf('Estimating FRC... %g %%', 100*bb/ncalc);

        legend([h_frc h_th], {'FRC' '1/2 bit threshold'}, ...
            'Location', 'northeast', 'Orientation', 'vertical', 'FontSize', 12);

    end
    drawnow;

end
%delete(h_wb);
delete(h_tt);


% Find and display where FRC and threshold curve intersects
% Interpolate for a bit more presision
fx1 = (fx(1):(fx(2)-fx(1))/10:fx(end))';
th1 = interp1(fx, t_halfbit, fx1, 'linear');
%th1 = frc_mean-t_halfbit;
id1 = find(interp1(fx, frc_mean, fx1, 'linear')-th1 < 0, 1, 'first');
frc_halfbit = 1/fx1(id1);
h_th1 = plot(fx1(id1), th1(id1), '.');
datatip(h_th1);
dtRows = [...
    dataTipTextRow('Spatial frequency / Ä^-^1', fx1(id1)); ...
    dataTipTextRow('Resolution estimate / Å', frc_halfbit); ...
    ];
h_th1.DataTipTemplate.DataTipRows = dtRows;
set(h_th1.DataTipTemplate, 'FontName', 'helvetica', 'FontSize', 12);

legend([h_frc h_th], {'FRC' '1/2 bit threshold'}, ...
    'Location', 'northeast', 'Orientation', 'vertical', 'FontSize', 12);



% FRC calculations
    function [frc, fx, nx] = FRC_Res(z1, z2, precision)
        % WFX Aug 2019 2 image FRC calculation update

        sz = max(size(z1));

        % Make sure size is odd and get centre coordinate s_c
        if mod(sz, 2) == 0
            sz = sz-1;
        end
        s_c = ceil(sz/2);

        lookup_x = (1:sz)-s_c;
        lookup_y = ((1:sz)-s_c)';
        [lookup_x, lookup_y] = meshgrid(lookup_x, lookup_y);
        [~, lookup_d] = cart2pol(lookup_x, lookup_y);

        % 2D FFT and shift 0 frequency component to the centre
        % Size of fft not needing to change
        %sz = size(z1, 1);
        t1 = fftshift(fft2(z1, sz, sz));
        t2 = fftshift(fft2(z2, sz, sz));

        % Test code
        % View 2D FFTs
        %imagesc(abs(t1));
        %imagesc(log10(abs(t1).^2));

        % FRC numerator FFT multiplication
        f1f2 = t1.*conj(t2);

        % Store it in a way that we can use it with accumarray
        %lookup_M = uint32(lookup_M(:));
        %lookup_d_int = ceil(lookup_d(:));
        % This line is incorrect as index distance 0 (centre should be
        % indexed to 1
        %lookup_d_int(lookup_d_int == 0) = 1;
        % This in the correct way
        lookup_d_int = ceil(lookup_d(:))+1;

        % Number of pixels per ring
        nx = accumarray(lookup_d_int, ones(size(lookup_d_int)));

        % sum the values of data points of each concentric circle
        rsum_num = accumarray(lookup_d_int, f1f2(:));

        % Compute fourier ring correlation curve
        % FRC numerator
        % t1.*conj(t2) should be all real, WFX checked imag component is
        % very small so ok
        % Numerator
        frc_num = real(rsum_num);

        in1 = abs(t1).^2;
        in2 = abs(t2).^2;

        % Denominator
        rsum_d1 = accumarray(lookup_d_int, in1(:));
        rsum_d2 = accumarray(lookup_d_int, in2(:));
        frc_denom = sqrt(abs(rsum_d1.*rsum_d2));

        % FRC calculation
        frc = double(frc_num)./double(frc_denom);
        %r = (1:max(lookup_M))';


        %
        % Find resolution scale from curve
        %

        % Convert x axis to spatial frequency
        % Sampling rate (px/Angstrom)
        sampling = 1/(precision*10);
        % Each Fourier pixel represents wavelengths per total length in Angstroms
        %fx = (0:sz-1)/(sz/sampling);
        %fx = fx(1:length(frc))';

        % Only use full circles
        fx = ((0:floor(sz/2))/(sz/sampling))';
        %fx = fx(1:length(frc))';
        frc = frc(1:numel(fx));
        nx = nx(1:numel(fx));

    end



if nargout == 1
    hh = h_fig;
end


end