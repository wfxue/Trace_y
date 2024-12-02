function [r, r_tip, tip_angle, appWidth, resnorm, residual, exitflag, output] = ...
    CalcTipFilamentConvolution(img, filament, tip_angle, precision)

%
% DESCRIPTION
% – CalcTipFilamentConvolution estimates the filament-tip convolution 
% effect, thereby also the helical filament cross-sectional radius and the
% tip radius by least-squares fitting of a cylindrical cross-sectional 
% model to the filament data averaged across its length. 
% – This method provides initial estimates of the tip geometry and an be
% used to generate a initial rudimentry point-cloud of the filament as the
% magnitude of the convolution as function of position compared to the 
% filament axis (interaction between r and r_tip) is estimated. This 
% methods assumes cylindrical or helical symetry.
% – Used by MakeHelicalFilamentModel.m as one of the first step in 
% helical 3D reconstruction.
% – Part of Trace_y
%
% USAGE
% Standard method usage to estimate r (cross-section) vs. r_tip
% >> [r, r_tip] = img.CalcTipFilamentConvolution(filament, tip_angle);
% Also get an apparent width of the filament in the image based on
% simulated image:
% >> [r, r_tip, tip_angle, appWidth] = img.CalcTipFilamentConvolution(filament, tip_angle);
% Also get full least-squares stats:
% >> [r, r_tip, tip_angle, appWidth, resnorm, residual, exitflag, output] = ...
%   img.CalcTipFilamentConvolution(img, filament)
%
% INPUTS
% img  –  AFMimage object. Implicit with dot notation
% filament  –  The index number of the traced filament to be analysed.
% tip_angle  –  The side angles of the tip. Can be estimated by inputing 
% [] but this is highly not recommended as the sides do not interact much 
% with small stuff on the size range of the tip-radius.
% precision  –  Optional sub-pixel precision for visualisation only. 
%
% OUTPUTS
% r  –  Cross-sectional radius of the filament
% r_tip  –  AFM tip radius
% tip_angle  –  AFM tip side (average) angles
% appWidth  –  Apparent width of the filament in the image in pixels
% [resnorm, residual, exitflag, output]  –  Least squares estimation stats
% from lsqnonlin function (Optimization toolbox)
%
% DEPENDENCIES
% – Uses Optimization Matlab toolbox
% – Method for Trace_y's @AFMimage/ object.
% – Used by @AFMImage/MakeHelicalFilamentModel.m
%
% AUTHORS
% Wei-Feng Xue
%
% HISTORY
% 2019.05  –  Initial draft of Calc_Filament_Corr.m by WFX as the first 
% attempt to estimate tip radius and to reconstruct a point-cloud of a 
% cylindrical or helical filament.
% 2023.11  –  Updated the code for use in updated 3D reconstruction
% algorithm MakeHelicalFilamentModel.m
% 2024.09  –  WFX Trace_y update including updates to the visualisation
% code and UI for MakeHelicalFilamentModel.m
%



%
% Defaults
%

%{
if (~exist('filament', 'var') || strcmp(filament, 'last')) && ~isempty(img.features.lastFilament)
    filament = 0;
elseif ~exist('mode', 'var') && isempty(img.features.lastFilament)
    filament = length(img.features.filaments);
end
%}

if ~exist('filament', 'var') || isempty(filament)
    filament = 1;
end

% Tip angle estimation is not robust for small stuff and should not be used
if ~exist('tip_angle', 'var') || isempty(tip_angle)
    tip_angle0 = 18;
    fit_tip_angle = true;
else
    fit_tip_angle = false;
end

if ~exist('precision', 'var') || isempty(precision)
    precision = 20;
end

% For plotting
h_fig = figure(14);
clf;
set(h_fig, 'NumberTitle', 'off', 'Name', 'Trace_y – CalcTipFilamentConvolution');
h_fig.OuterPosition(3:4) = [1320 660];

%drawnow;
screenSize = get(0, 'ScreenSize');
if h_fig.OuterPosition(2)+h_fig.OuterPosition(4) > screenSize(4)
    h_fig.OuterPosition(2) = screenSize(4)-h_fig.OuterPosition(4)-10;
end
clf;

tl = tiledlayout(1, 5, 'Padding', 'loose', 'TileSpacing', 'compact');
title(tl, '   ', ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontSize', 16);
tx = xlabel(tl, sprintf('%s, filament %d', img.dataFile, filament), ...
    'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontSize', 16);
set(tx, 'FontWeight', 'normal', 'FontAngle', 'italic', 'FontSize', 16);



%
% Get data
%

% Get straightened filament pixels
zz = img.StraightenFilament(filament);
% Pre-rotation to vertical
yy = 0:img.xResolution:img.xResolution*(size(zz, 2)-1)';
xx = (0:img.xResolution:img.xResolution*(size(zz, 1)-1))-img.xResolution*ceil(size(zz, 1)/2-1);

imgResolution = img.xResolution;
%pixelDensity = 1./img.xResolution;

%subplot(2, 1, 1);
h_ax1 = nexttile(2, [1 4]);
imagesc(yy, xx, zz);
colormap(img.cMap);
set(gca, 'yDir', 'normal', ...
    'box', 'on', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', 'FontSize', 12);
set(gca, 'YTickLabelMode', 'manual');
xlabel('y / nm', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'points', 'FontSize', 14);
%ylabel('x / nm', ...
%    'FontName', 'helvetica', 'FontWeight', 'normal', ...
%    'FontUnits', 'points', 'FontSize', 14);
title('Straightened filament image data (aspect ratio not to scale)', 'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontUnits', 'points', 'FontSize', 14);
h_cb = colorbar;
set(h_cb, 'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontSize', 12)
set(h_cb.Label, 'String', sprintf('%s / %s', img.imgType, img.zUnit), ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontSize', 16)

drawnow;

% Rotate to vertical
zz = imrotate(zz, -90);
za = mean(zz, 1);

% Test code
%{
figure(1000);
subplot(2, 1, 1);
imagesc(xx, yy, zz);
set(gca, 'YDir', 'normal');
xlabel('y / nm', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'pixels', 'FontSize', 14);
ylabel('x / nm', ...
    'FontName', 'helvetica', 'FontWeight', 'normal', ...
    'FontUnits', 'pixels', 'FontSize', 14);
title('Straightened filament image data', 'interpreter', 'none', ...
    'FontName', 'helvetica', 'FontWeight', 'bold', ...
    'FontUnits', 'pixels', 'FontSize', 14);
subplot(2, 1, 2);
plot(xx, za, 'bo');
%}

% Trim zz a bit for plotting lims
w = floor(length(za)/4);
%wp = w*precision;
%za = za(w+1:end-w);
%xx = xx(w+1:end-w);
%xcentre = mean(xx);
xcentre = 0;

% The x coordinates in nm
% Height values -> z already in nm
% za is averaged values across all the y cross sections
%xx = xx.*img.xResolution;

% p = [radius, tip_radius]
% p0: initial estimate
if fit_tip_angle
    p0 = [max(za)/2; 2; tip_angle0];
    lb = [0; 0; 1];
    ub = [max(za); 100; 27.5];
else
    p0 = [max(za)/2; 2];
    lb = [0; 0];
    ub = [max(za); 100];
end

% Weight as only top pixels are meaningful for this calculation
% Too stringent if only considering top pixels
%{
lsq_weight = zeros(size(za));
lsq_weight(za > 0.5*max(za)) = 1;
% Incase thinn filament with only one pixel ridge.
if sum(lsq_weight) == 1
    ind = find(lsq_weight == 1, 1, 'first');
    lsq_weight(ind-1:ind+1) = 1;
end
%}
% Try this simple weight instead
lsq_weight = abs(za)./max(za);

lsqOpt = optimoptions('lsqnonlin', ...
    'Display', 'off', ...
    'MaxFunEvals', 1000000, 'MaxIter', 500, ...
    'TolX', 1e-6, 'TolFun', 1e-6, ...
    'OutputFcn', @LsqProgressPlot);
%'Display', 'iter-detailed', ...

% Set up various other variables
hh_data = [];
hh1 = [];
hh2 = [];
hh3 = [];
h_tt = [];
x_sim = 0;
z_sim = 0;
x_calc = 0;
z_calc = 0;

% Estimation by LSQ
[p, resnorm, residual, exitflag, output] = lsqnonlin(@LsqFcn, p0, lb, ub, lsqOpt);
r = p(1);
r_tip = p(2);

% Apparent width can now be estimated based on the simulated cylindar
% cross-section
appWidth = sum(z_sim.*z_calc ~= 0)+2;



% Sub-routines
    function [residual, zcalc, ssr] = LsqFcn(pfit)
        
        r = pfit(1);
        r_tip = pfit(2);
        if fit_tip_angle
            tip_angle = pfit(3);
        end
        %fprintf('r = %g, tip_r = %g\n', r, tip_r);
        
        % Old sim code
        %
        % Target function
        % add 0.9*img.xResolution/precision for safety
        xcalc = min(xx):img.xResolution/precision:max(xx)+0.9*img.xResolution/precision;
        %xcalc = [-fliplr(xcalc(2:end)) xcalc];
        zcalc = Cylinder_Upper_Surface(xcalc, r, xcentre);
        
        % Using old AFM_Image_Sim code (updated in SimHeightImageData.m)
        % for a single scan line 2D simulation. 
        % This function default to 0 based x values so need to adjust x
        xcentre2Dsim = xcalc(1);
        xcalc = xcalc-xcentre2Dsim;
        %[s_x, ~, s_z_sim, s_x_cor, ~, s_z_cor] = ...
        %    AFM_Image_Sim(r_tip, tip_angle, img.xResolution, precision, 'data2D', [xcalc; zcalc], 0);
        [x_sim, ~, z_sim, x_calc, ~, z_calc] = ...
            SimHeightImageData(r_tip, tip_angle, img.xResolution, ...
            precision, 'data2D', [xcalc; zcalc], 0);
        % Adjust x back to 0 centrerd
        x_sim = x_sim+xcentre2Dsim;
        x_calc = x_calc+xcentre2Dsim;
        % Above is an akward patch needing tyding up later
        %}

        % Using SimHeightImage works but is inefficient and potentially 
        % really slow since we don't need an image just one scanline.
        %{
        img_sim = SimHeightImage(numel(za), pixelDensity, ...
            'cylinder', r, 'rounded_cone', [r_tip tip_angle], ...
            [], [], [], false);
        x_sim = img_sim.x;
        z_sim = img_sim.z(1, :);
        x_calc = img_sim.xC(1, :);
        z_calc = img_sim.zC(1, :);
        %}

        % Check corrected data
        zcalc = z_calc;
        residual = z_sim-za;
        residual = residual.*lsq_weight;
        ssr = sum(residual(:).^2);
    end

    function z = Cylinder_Upper_Surface(x, r, c)
        % Cylinder centrerd at x = 0, across y
        x = x-c;
        z = sqrt(r.^2-x.^2)+r;
        z = real(z)-(imag(z) > 0).*r;
    end

    function stop = LsqProgressPlot(pfit, optimVal, state)
        % Standard info
        stop = false;
        switch state
            case 'iter'
                % Make updates to plot or guis as needed
            case 'interrupt'
                % Probably no action here. Check conditions to see
                % whether optimization should quit.
            case 'init'
                % Setup for plots or guis
                
                %subplot(2, 1, 2);
                h_ax2 = nexttile(1, [1 1]);
                %hh_data = plot(za(w+1:end-w), xx(w+1:end-w), 'bo');
                hh_data = plot(za, xx, 'bo', 'MarkerSize', 10);
                set(h_ax2, 'box', 'on', ...
                    'FontName', 'helvetica', 'FontWeight', 'normal', 'FontSize', 12);
                box('on');
                hold('on');
                %set(gca, 'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1]);
                %set(gca, 'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1]);
                xlabel('mean z / nm', ...
                    'FontName', 'helvetica', 'FontWeight', 'normal', ...
                    'FontUnits', 'points', 'FontSize', 14);
                ylabel('x / nm', ...
                    'FontName', 'helvetica', 'FontWeight', 'normal', ...
                    'FontUnits', 'points', 'FontSize', 14);
                title('Convolution estimation', 'interpreter', 'none', ...
                    'FontName', 'helvetica', 'FontWeight', 'bold', ...
                    'FontUnits', 'points', 'FontSize', 14);
                h_ax1.YLim = imgResolution*[-w w];
                h_ax2.YLim = imgResolution*[-w w];
                linkaxes([h_ax1 h_ax2], 'y');
                
            case 'done'
                % Cleanup of plots, guis, or final plot
                %title(sprintf('LSQ done: iteration = %g, r = %g, r_t_i_p = %g, angle_t_i_p = %g', ...
                %    optimVal.iteration, r, r_tip, tip_angle), ...
                %    'interpreter', 'tex', ...
                %    'FontName', 'helvetica', 'FontWeight', 'bold', ...
                %    'FontUnits', 'pixels', 'FontSize', 14);
                %drawnow;
                %return
            otherwise
        end
        
        % Plotting
        r = pfit(1);
        r_tip = pfit(2);
        % Target function
        %xcalc = 0:img.xResolution/precision:max(xx)+0.9*img.xResolution/precision;
        %xcalc = [-fliplr(xcalc(2:end)) xcalc];
        xcalc = (0:imgResolution/precision:max(xx)+imgResolution/precision);
        xcalc = unique([-xcalc xcalc]);
        zcalc = Cylinder_Upper_Surface(xcalc, r, xcentre);
        
        delete(hh1);
        delete(hh2);
        delete(hh3);
        
        %hh1 = plot(xcalc(wp+1:end-wp), zcalc(wp+1:end-wp), 'k-');
        %hh2 = plot([x_sim(w+1:end-w); x_calc(w+1:end-w)], [z_sim(w+1:end-w); z_calc(w+1:end-w)], 'r+-');
        %hh3 = plot(x_sim(w+1:end-w), z_sim(w+1:end-w), 'r+-');
        hh1 = plot(zcalc, xcalc, 'k-');
        hh2 = plot([z_sim; z_calc], [x_sim; x_calc], 'r+-', 'MarkerSize', 10, 'LineWidth', 1);
        hh3 = plot(z_sim, x_sim, 'r-', 'LineWidth', 1);

        legend([hh_data hh1 hh2(1)], ...
            'Image data', 'Estimated crossection', 'Estimated convolution', ...
            'Location', 'northeast', 'FontSize', 10);
        %title(sprintf('LSQ: iteration = %g, r = %g, r_t_i_p = %g, angle_t_i_p = %g', ...
        %    optimVal.iteration, r, r_tip, tip_angle), ...
        %    'interpreter', 'tex', ...
        %    'FontName', 'helvetica', 'FontWeight', 'bold', ...
        %    'FontUnits', 'pixels', 'FontSize', 14);

        delete(h_tt);
        h_tt = text(0.95, 0.025, sprintf('LSQ iteration %g\n\nr_c_s = %g\nr_t_i_p = %g nm\nangle_t_i_p = %g deg', ...
            optimVal.iteration, r, r_tip, tip_angle), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'bottom', 'FontSize', 12, ...
            'FontName', 'helvetica', 'FontWeight', 'normal', 'Color', 'k');
        drawnow;
        
    end


end