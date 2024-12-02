function [img] = TraceFilamentM(img, appWidth, zThreshold, maxSearchPhiDeg, xy_nodes, targetFig)

%
% DESCRIPTION
% – Semi automatic filament tracer, with image modelling method, to
% estimate filaments' central axis.
% – User can select node points on the image to select the filament(s) to 
% be traced.
% – Part of Trace_y
%
% USAGE
% Standard method usage
% >> img = img.TraceFilamentM(appWidth, zThreshold);
% Tracing the filament with aditional information on the angles and node
% points given apriori
% >> img = img.TraceFilamentM(appWidth, zThreshold, , maxSearchPhiDeg, xy_nodes);
%
% INPUTS
% img  –  AFMimage object. Implicit use with dot notation
% appWidth  –  Apparent width of the filament as appers in the image in 
% pixels
% zThreshold  –  z-threshold for background noise to be ignored
% maxSearchPhiDeg  –  Aaximum ±angle change in degrees per pixel
% nodes  –  [x y] coordinates that trace must go through, should have 2 
% columns and at least two rows (i.e. two points)
% targetFig  –  Target figure window handle to use
%
% OUTPUTS
% img  –  An instance of Trace_y @AFMimage/ AFM image object with the newly
% traced filament(s) added in img.features.filaments
%
% DEPENDENCIES
% – Optimization, Signal processing, Statistical and Machine Learning,
% Image Processing Matlab toolboxes
% – Method for Trace_y's @AFMimage/ object.
% – @AFMimage/MakeImgModelLSQ.m, TraceFilamentMC.m
% – Used also by other routines
%
% AUTHORS
% Wei-Feng Xue
%
% HISTORY
% 2009.06  –  FibAFM scripts by WFX used in Fibril fragmentation enhances 
% amyloid cytotoxicity, JBC, 284 (2009) 34272-34282, 
% http://doi.org/10.1074/jbc.m109.049809 and Amyloid fibril length 
% distribution quantified by atomic force microscopy single-particle 
% image analysis, Protein Engineering Design and Selection 22 (2009), 
% 489-496, http://doi.org/10.1093/protein/gzp026.
% 2019.02  –  Initial draft of TraceFilamentM using updated Gaussian wall
% model. Modernised the otherwise old 2009 WFX tracing scripts.
% 2022.04  –  WFX Updated with better GUI filament selection approach and 
% so no more end condition checks needed
% 2022.09  –  WFX updated method to increase run speed using a step free 
% approach to estimate rough contour as intitial input for applying the 
% Gausian wall model as refinement step. Should be much more efficient.
% 2024.09  –  WFX updated the UI to allow consecutive traces and to save
% the trace by default to streamline the workflow 
%



%
% Defaults
%

% Setting for plotting every plotstep or not
%htrace = [];
plotstep = 20;
% Step size in pixels
step_size = 0.5;


% Init default parameters
if ~exist('appWidth', 'var') || isempty(appWidth)
    appWidth = 5;
end
if ~exist('zThreshold', 'var') || isempty(zThreshold)
    zThreshold = 3*img.zNoiseStd;
end

% Convert to max angular change per step in radians
if ~exist('maxSearchPhiDeg', 'var') || isempty(maxSearchPhiDeg)
    maxSearchPhiDeg = 10;
    maxSearchPhi = deg2rad(maxSearchPhiDeg)*step_size;
else
    maxSearchPhi = deg2rad(maxSearchPhiDeg)*step_size;
end

% This part is no longer needed
%{
if ~exist('maxStep', 'var') || isempty(maxStep)
    maxStep = 2*img.nLinesImg;
end
if ~exist('stepSize', 'var') || isempty(minStepSize)
    % standard stepSize is 1 pixel
    minStepSize = 1;
end
%}

% Get image data
z = img.z;

% Set up main window and show image data
%figure(1);
%set(1, 'OuterPosition', [51 51 960 960]);
%clf;
%img.DispImage;
if ~exist('targetFig', 'var') || isempty(targetFig)
    h_fig = img.ExploreFeatures;
else
    % Save the current view
    h_img = findobj(targetFig, 'Type', 'Image');
    h_ax = h_img.Parent;
    current_xlims = h_ax.XLim;
    current_ylims = h_ax.YLim;

    h_fig = img.ExploreFeatures(targetFig);

    % Recall the current view
    h_img = findobj(h_fig, 'Type', 'Image');
    h_ax = h_img.Parent;
    h_ax.XLim = current_xlims;
    h_ax.YLim = current_ylims;
end
hold('on');
%h_fig = img.ExploreFeatures;

set(h_fig, 'NumberTitle', 'off', 'Name', 'Trace_y – TraceFilamentM');

% Modify callback functions
h_img = findobj(h_fig, 'Type', 'Image');
h_img.ButtonDownFcn = @SetXYnodes;

h_lines = findobj(h_fig, 'Type', 'Line');
for aa = 1:numel(h_lines)
    h_lines(aa).ButtonDownFcn = [];
end

% Add a UI button for start
traceBtn = uicontrol('Style', 'pushbutton', ...
    'String', 'Trace', ...
    'Units', 'pixels', 'Position', [75 70 150 25], ...
    'Callback', @TraceNow, 'Enable', 'off', 'UserData', [], ...
    'FontName', 'helvetica', ...
    'FontUnits', 'pixels', 'FontSize', 14, ...
    'FontWeight', 'normal');



% Plot all filaments in the file
% No longer needed with ExploreFeatures.m
%{
hold('on');

for aa = 1:length(img.features.filaments)
    filament = img.features.filaments(aa);
    plot(filament.x, filament.y, 'b-', 'LineWidth', 2);
end
%{
if ~isempty(img.features.lastFilament)
    filament = img.features.lastFilament;
    plot(filament.x, filament.y, 'c-', 'LineWidth', 2);
end
%}
%}


xy_nodes = [];
ll = [];
is_node = [];
hh = [];
hh2 = [];
hh3 = [];


%
% Callback function for node selection
%

    function [] = SetXYnodes(~, evt)

        xx = evt.IntersectionPoint(1);
        yy = evt.IntersectionPoint(2);
        button = evt.Button;


        % For centre and range for UI zooming [centre range]
        %xlims = xlim;
        %ylims = ylim;
        %zoom_range = (xlims(2)-xlims(1))/2;

        if button == 1 && (z(round(yy), round(xx)) > zThreshold) && (xx >= 1 && xx <= size(z, 2)) && (yy >= 1 && yy <= size(z, 1))
            % Normal left button click Check point within image bounds and
            % point above noise
            %xy_nodes = sortrows([xy_nodes; xx yy]);

            % Put the new nodes in to be in seqence
            dd_min = Inf;
            for bbb = 0:size(xy_nodes, 1)
                xy_nodes_test = [xy_nodes(1:bbb, :); xx yy; xy_nodes(bbb+1:end, :)];
                dd = pdist(xy_nodes_test);
                dd = squareform(dd);
                dd = sum(diag(dd, -1));
                if dd < dd_min
                    xy_nodes_new = xy_nodes_test;
                    dd_min = dd;
                end
            end
            xy_nodes = xy_nodes_new;


        elseif button == 3 || button == 2
            % If right button click, remove peak within 5 pixels
            dd = pdist([xy_nodes; xx yy]);
            dd = squareform(dd);
            xy_nodes = xy_nodes(dd(end, 1:end-1) > 5, :);

        %elseif zoom_range < size(z, 1)/2
        %    % If nothing is clicked then pan window to follow click
        %    xlim(round([xx-zoom_range xx+zoom_range]));
        %    ylim(round([yy-zoom_range yy+zoom_range]));
        %    % This is very disorienting so disable
        end

        % Plot the nodes
        delete(hh);
        if ~isempty(xy_nodes)
            hh = plot(xy_nodes(:, 1), xy_nodes(:, 2), 'c+', 'MarkerSize', 16, 'LineWidth', 2);
            hh.ButtonDownFcn = @SetXYnodes;
        end

        % Plot the initial contour
        delete(hh2);
        if size(xy_nodes, 1) > 1
            % compute the splines with added interpolation with approx 0.5
            % pixel distance
            dd = pdist(xy_nodes);
            dd = squareform(dd);
            dd = [0; cumsum(diag(dd, -1))];

            %ll = [0:step_size:max(dd) dd(end)]';
            %ll = [(0:step_size:max(dd))'; dd];
            ll = linspace(0, max(dd), round(max(dd)/step_size))';

            %is_node = logical([zeros(size((0:step_size:max(dd))')); ones(size(dd))]);
            %[ll, idx] = unique(ll, 'last');
            %is_node = is_node(idx);
            xx = spline(dd, xy_nodes(:, 1), ll);
            yy = spline(dd, xy_nodes(:, 2), ll);

            is_node = ones(size(dd));
            for bbb = 1:numel(dd)
                is_node(bbb) = find(abs(dd(bbb)-ll) == min(abs(dd(bbb)-ll)));
            end

            % draw final result
            hh2 = plot(xx, yy, 'c-', 'LineWidth', 1);
            hh2.ButtonDownFcn = @SetXYnodes;
            %plot(xx(is_node), yy(is_node), 'ro', 'MarkerSize', 16);
            %drawnow

            traceBtn.Enable = 'on';
        else
            traceBtn.Enable = 'off';
        end

    end


% The UL node selection using gin method is slow and obsolete now with
% callback functions
%{
% Centre and range for UI zooming [centre range]
x_zoom_full = xlim;
y_zoom_full = ylim;
zoom_range_full = x_zoom_full(1)+(x_zoom_full(2)-x_zoom_full(1))/2;
zoom_range = zoom_range_full;


% Graphical input node point if needed
if ~exist('xy_nodes', 'var') || isempty(xy_nodes)
    
    hold('on')
    xy_nodes = [];
    
    % Set up button check  flag
    quitEval = 0;
    hh = [];
    hh2 = [];
    
    while ~quitEval
        
        title(sprintf('Left click: mark, Right click: un-mark, (z)oom in, (f)ull zoom, Return: start\n'));
        [xx, yy, button] = ginput(1);
        
        % enter returns empty, 27 = esc, 113 = 'q'
        if (size(xy_nodes, 1) > 1) && (isempty(button) || button == 27 || button == 113)
            quitEval = 1;
            
            % toggle, 116 = 't'
        %elseif ~isempty(button) && button == 116
        %    zoom('on');
        %    pause;
            % pause does not work with uifigure for compiled version so try
            % ths following zooming instead when clicking on bg or 'f','z'
            % Revert to full zoom, 102 = 'f'
        elseif ~isempty(button) && button == 102
            xlim(x_zoom_full);
            ylim(y_zoom_full);
            zoom_range = zoom_range_full;

            % Zoom in, 122 = 'z'
        elseif ~isempty(button) && button == 122
            zoom_range = zoom_range/2;
            x_zoom = round([xx-zoom_range xx+zoom_range]);
            y_zoom = round([yy-zoom_range yy+zoom_range]);
            xlim(x_zoom);
            ylim(y_zoom);

            % Move
        elseif z(round(yy), round(xx)) < zThreshold
            if zoom_range == zoom_range_full
                zoom_range = zoom_range/4;
            end
            x_zoom = round([xx-zoom_range xx+zoom_range]);
            y_zoom = round([yy-zoom_range yy+zoom_range]);
            xlim(x_zoom);
            ylim(y_zoom);


            % Check point within image bounds
        elseif ~isempty(button) && button == 1 && (z(round(yy), round(xx)) > zThreshold) && (xx >= 1 && xx <= img.pixelPerLine) && (yy >= 1 && yy <= img.nLines)
            %xy_nodes = sortrows([xy_nodes; xx yy]);
            
            % Put the new nodes in to be in seqence
            dd_min = Inf;
            for aa = 0:size(xy_nodes, 1)
                xy_nodes_test = [xy_nodes(1:aa, :); xx yy; xy_nodes(aa+1:end, :)];
                dd = pdist(xy_nodes_test);
                dd = squareform(dd);
                dd = sum(diag(dd, -1));
                if dd < dd_min
                    xy_nodes_new = xy_nodes_test;
                    dd_min = dd;
                end
            end
            xy_nodes = xy_nodes_new;
            
            
            % If right button click, remove peak within 5 pixels
        elseif ~isempty(button) && button == 3
            dd = pdist([xy_nodes; xx yy]);
            dd = squareform(dd);
            xy_nodes = xy_nodes(dd(end, 1:end-1) > 5, :);
            
        end
        
        if ~isempty(xy_nodes)
            delete(hh);
            hh = plot(xy_nodes(:, 1), xy_nodes(:, 2), 'b+', 'MarkerSize', 16, 'LineWidth', 2);
        end
        
        if size(xy_nodes, 1) > 1
            % compute the splines with added interpolation with approx 0.5
            % pixel distance
            dd = pdist(xy_nodes);
            dd = squareform(dd);
            dd = [0; cumsum(diag(dd, -1))];
            
            %ll = [0:step_size:max(dd) dd(end)]';
            %ll = [(0:step_size:max(dd))'; dd];
            ll = linspace(0, max(dd), round(max(dd)/step_size))';

            %is_node = logical([zeros(size((0:step_size:max(dd))')); ones(size(dd))]);
            %[ll, idx] = unique(ll, 'last');
            %is_node = is_node(idx);
            xx = spline(dd, xy_nodes(:, 1), ll);
            yy = spline(dd, xy_nodes(:, 2), ll);

            is_node = ones(size(dd));
            for bb = 1:numel(dd)
                is_node(bb) = find(abs(dd(bb)-ll) == min(abs(dd(bb)-ll)));
            end
            
            % draw final result
            delete(hh2);
            hh2 = plot(xx, yy, 'c-', 'LineWidth', 2);
            %plot(xx(is_node), yy(is_node), 'ro', 'MarkerSize', 16);
            %drawnow
            
        end
    end
    
end
%}


%
% Callback for starting the trace procedure
%
    function [] = TraceNow(~, ~)

        traceBtn.Enable = 'inactive';

        % First pass: preparing rough contour
        [xx, yy] = TraceFilamentMC(h_fig, z, xx, yy, ll, is_node, appWidth, zThreshold);
        delete(hh2);
        hh2 = plot(xx, yy, 'c-', 'LineWidth', 2);

        % Second pass, refining
        traceBtn.String = 'Refining trace...';

    end


%
% Resuming tracing procedure
%

traceNext = true;
fprintf('\n');

while traceNext

    waitfor(traceBtn, 'String', 'Refining trace...');
    % Taking over from the callback of TraceNow button

    if isgraphics(traceBtn)

        % Second pass, refining
        % Tracing...
        %
        %hh_text = text(0.02, 0.04, 'Refine trace...', 'Units', 'normalized', 'FontSize', 16, 'Color', 'c');
        % Note variable name change to xf and yf, xx and yy for use as crop coords
        % This is due to chronological developement reasons
        % Last elements are loop dummies
        xf = [xx; xx(end)];
        yf = [yy; yy(end)];
        % Also reusing the normals soo rename theta variable.
        %theta_norm = theta(:, 1);

        maxStep = size(xx, 1);

        % Distans between points
        %xy_dist = sqrt((xf(2:end)-xf(1:end-1)).^2+(yf(2:end)-yf(1:end-1)).^2);

        % The differences from spline with manual nodes from first pass
        xd = zeros(maxStep+1, 1);
        yd = zeros(maxStep+1, 1);

        theta = zeros(maxStep+1, 1);
        hf = zeros(maxStep+1, 1);
        sigma = ones(maxStep+1, 1);
        model = zeros(maxStep+1, 1);
        exitflag = zeros(maxStep+1, 1);
        rmsd = zeros(maxStep+1, 1);

        % Some initial values guesses
        %theta(1) = cart2pol(xf(2)-xf(1), yf(2)-yf(1));
        theta(1:end-1) = cart2pol((xf(2:end)-xf(1:end-1)), (yf(2:end)-yf(1:end-1)));
        theta(end-1:end) = theta(end-2);
        %sigma(1) = appWidth/4;
        sigma = sigma.*appWidth/4;
        %sigma_std = appWidth/15;

        theta_norm = theta+pi/2;
        drho = zeros(maxStep+1, 1);
        % For checking normals
        %{
        [norm_x, norm_y] = pol2cart(theta_norm, appWidth);
        for aa = 1:size(theta_norm, 1)-1
            %plot([xx(aa) xx(aa)+normal(aa, 1)], [yy(aa) yy(aa)+normal(aa, 2)], 'r');
            %plot(xx(aa)+norm_x(:, aa), yy(aa)+norm_y(:, aa), 'r-');
            plot([xf(aa)-norm_x(aa) xf(aa) xf(aa)+norm_x(aa)], [yf(aa)-norm_y(aa) yf(aa) yf(aa)+norm_y(aa)], 'r-');
        end
        %}

        % Save original values from 1st pass
        xf1 = xx;
        yf1 = yy;

        % Masking not used at the moment
        z_mask = ones(size(z));

        for aa = 1:maxStep
            % Calc optimum model (LSQ)

            % initial guess
            %x0 = xf(aa)+xd(aa);
            %y0 = yf(aa)+yd(aa);
            x0 = xf(aa);
            y0 = yf(aa);
            theta0 = theta(aa);
            sigma0 = sigma(aa);

            % Crop the relevant part of the img
            % Decided to use input parameter appWidth insted for local crop to
            % allow useer some control of its size, wich  will effect the trace
            % quality slightly
            xx = (...
                max(round(x0)-ceil(appWidth), 1):...
                min(round(x0)+ceil(appWidth), size(z, 2))...
                )';
            yy = (...
                max(round(y0)-ceil(appWidth), 1):...
                min(round(y0)+ceil(appWidth), size(z, 1))...
                )';
            %{
    xx = (...
        max(round(x0)-2*ceil(sigma0), 1):...
        min(round(x0)+2*ceil(sigma0), size(z, 2))...
        )';
    yy = (...
        max(round(y0)-2*ceil(sigma0), 1):...
        min(round(y0)+2*ceil(sigma0), size(z, 1))...
        )';
            %}
            zz = z(yy, xx).*z_mask(yy, xx);

            % Weights not used, set to 1
            w = ones(size(zz));

            % Reminder: for MakeImgModelLSQ the parameters are:
            % p: [x, y, sigma, h, theta, l, c2]
            % only p(1) to p(5) are used for the models 1 and 2 (2 used here)
            % p(7) = c2 is the bending parameter, can be enabled but not nessasary
            % nessasary, not activated here anyway
            % Only use model 2 in MakeImgModelLSQ: @FibSegmentLSQ

            % p0 and bounds estimate
            %p0 = [x0 y0 sigma0 z(round(y0), round(x0)) theta0 appWidth 0];

            % Old bounds with masking idea, not used
            %{
    if aa == 1
        % Old bounds
        p0 = [x0,        y0,        sigma(aa),  z(round(y0), round(x0)), theta0,              appWidth,   0];
        ub = [x0+sigma0, y0+sigma0, appWidth/2, 2*max(zz(:)),            theta0+maxSearchPhi, 2*appWidth, 0.02];
        lb = [x0-sigma0, y0-sigma0, 0.25,       zThreshold/2,            theta0-maxSearchPhi, 0,          -0.02];
        sigma_mean = sigma(aa);
        hf_mean = z(round(y0), round(x0));
    elseif aa >=2 && aa <= 2*appWidth
        p0 = [x0,               y0,               sigma(aa-1),     z(round(y0), round(x0)), theta0,               appWidth,   0];
        ub = [x0+xy_dist(aa)/2, y0+xy_dist(aa)/2, 1.2*sigma(aa-1), 2*max(zz(:)),            theta0+maxSearchPhi, 2*appWidth, 0.02];
        lb = [x0-xy_dist(aa)/2, y0-xy_dist(aa)/2, 0.8*sigma(aa-1), zThreshold/2,            theta0-maxSearchPhi, 0,          -0.02];
        sigma_mean = mean(sigma(1:aa-1));
        sigma_std = std(sigma(1:aa-1));
        sigma(aa:end) = sigma_mean;
        hf_mean = mean(hf(1:aa-1));
    else
        p0 = [x0,               y0,               sigma0,             z(round(y0), round(x0)), theta0,              appWidth,   0];
        ub = [x0+xy_dist(aa)/5, y0+xy_dist(aa)/5, sigma0+2*sigma_std, 2*max(zz(:)),            theta0+maxSearchPhi, 2*appWidth, 0.02];
        lb = [x0-xy_dist(aa)/5, y0-xy_dist(aa)/5, sigma0-2*sigma_std, zThreshold/2,            theta0-maxSearchPhi, 0,          -0.02];
    end
            %}

            p0 = [x0, y0, sigma0,   z(round(y0), round(x0)), theta0,              appWidth,   0,   theta_norm(aa), 0];
            ub = [x0, y0, 2*sigma0, 2*max(zz(:)),            theta0+maxSearchPhi, appWidth, 0.02,  theta_norm(aa), appWidth];
            lb = [x0, y0, 0.25,     zThreshold/2,            theta0-maxSearchPhi, appWidth, -0.02, theta_norm(aa), -appWidth];

            [modelc, p, ~, zcalc, ~, exitflagc] = ...
                img.MakeImgModelLSQ(xx, yy, zz, w, p0, lb, ub, 5);

            % Finalise current point
            model(aa) = modelc;
            zd = zz-zcalc;
            rmsd(aa) = sqrt(sum(zd(:).^2)./numel(zd));

            xd(aa) = p(1)-xf(aa);
            yd(aa) = p(2)-yf(aa);

            xf(aa) = p(1);
            yf(aa) = p(2);
            sigma(aa) = p(3);
            hf(aa) = p(4);
            theta(aa) = p(5);

            drho(aa) = p(9);
            exitflag(aa) = exitflagc;
            % z is not calculated here, h will be used and z will calc later

            % Masking, idea not used at the moment
            %z_mask(yy, xx) = z_mask(yy, xx) & (zcalc > hf_mean*normpdf(2*sigma_mean, 0, sigma_mean));

            % Guessing p0 for next point
            %[x_next, y_next] = pol2cart(theta(aa), xy_dist(aa));
            %xf(aa+1) = xf(aa)+x_next;
            %yf(aa+1) = yf(aa)+y_next;
            %xd(aa+1) = xd(aa);
            %yd(aa+1) = yd(aa);
            %sigma(aa+1) = sigma(aa);
            %theta(aa+1) = theta(aa);

            %[dtheta, drho] = cart2pol(xd(aa), yd(aa));
            %drho_test = drho*normpdf(ll, ll(aa), 2*appWidth)/normpdf(0, 0, 2*appWidth);
            %[xx_test, yy_test] = pol2cart(dtheta, drho_test);

            %xx_test = xd(aa)*normpdf(ll, ll(aa), 2*appWidth)/normpdf(0, 0, 2*appWidth);
            %yy_test = yd(aa)*normpdf(ll, ll(aa), 2*appWidth)/normpdf(0, 0, 2*appWidth);

            %xd(aa+1:end-1) = xd(aa+1:end-1)+xx_test(aa+1:end);
            %yd(aa+1:end-1) = yd(aa+1:end-1)+yy_test(aa+1:end);
            %xf(aa+1:end-1) = xf(aa+1:end-1)+xx_test(aa+1:end);
            %yf(aa+1:end-1) = yf(aa+1:end-1)+yy_test(aa+1:end);
            %theta(aa+1) = theta(aa);


            % Get a better estimate of sigma
            if aa <= appWidth
                sigma_mean = mean(sigma(1:aa));
                sigma(aa+1:end) = sigma_mean;
            end

            %{
    if aa <= appWidth
        %isgoodfit(aa) =  1;
        %sigma_std = std(sigma(1:aa));
        sigma_mean = mean(sigma(1:aa));
        sigma(aa+1:end) = sigma_mean;
        %rmsd_mean = mean(rmsd(1:aa));
        %rmsd_std = std(rmsd(1:aa));

    elseif aa > appWidth
        
        %sigma_mean = mean(sigma(1:aa));
        %sigma(aa+1:end) = sigma_mean;

        rmsd_mean = mean(rmsd(1:aa-1));
        rmsd_std = std(rmsd(1:aa-1));
        if rmsd(aa) < rmsd_mean+3*rmsd_std
            isgoodfit(aa) =  1;
        end

    end
            %}

            % Plotting current trace
            if mod(aa, plotstep) == 0 || aa == maxStep-1
                %figure(1);
                delete(hh3);
                hh3 = plot(xf(1:aa+1), yf(1:aa+1), 'c-+', 'LineWidth', 1);
                traceBtn.String = sprintf('Refining... %g%%', round(100*aa/maxStep));
                drawnow;
            end

            %{
    % For checking and debugging
    figure(2);
    %set(2, 'OuterPosition', [1011 51 500 960]);
    set(2, 'OuterPosition', [811 51 850 960]);
    clf;
    subplot(3, 2, 1);
    [uu, vv] = pol2cart(theta(aa), 2);
    imagesc(xx, yy, zz);
    set(gca, 'yDir', 'normal', ...
        'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1], ...
        'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1], ...
        'box', 'on');
    colormap(img.cMap);
    set(gca, 'CLim', img.cScale);
    hold('on');
    plot(xf(aa), yf(aa), 'r+', 'MarkerSize', 16, 'LineWidth', 2);
    quiver(xf(aa), yf(aa), uu, vv, 'r', 'LineWidth', 2);
    plot(x0, y0, 'ro', 'MarkerSize', 16);
    plot(xf(aa+1)+xd(aa+1), yf(aa+1)+yd(aa+1), 'b+', 'MarkerSize', 16);

    subplot(3, 2, 2);
    imagesc(xx, yy, zcalc);
    set(gca, 'yDir', 'normal', ...
        'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1], ...
        'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1], ...
        'box', 'on');
    set(gca, 'CLim', img.cScale);
    hold('on');
    plot(xf(aa), yf(aa), 'r+', 'MarkerSize', 16, 'LineWidth', 2);
    quiver(xf(aa), yf(aa), uu, vv, 'r', 'LineWidth', 2);

    subplot(3, 2, 3);
    imagesc(xx, yy, zz-zcalc);
    set(gca, 'yDir', 'normal', ...
        'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1], ...
        'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1], ...
        'box', 'on');
    set(gca, 'CLim', img.cScale);
    hold('on');
    plot(xf(aa), yf(aa), 'r+', 'MarkerSize', 16, 'LineWidth', 2);
    quiver(xf(aa), yf(aa), uu, vv, 'r', 'LineWidth', 2);
    
    ax = subplot(3, 2, 4);
    imagesc(xx, yy, z_mask(yy, xx));
    set(gca, 'yDir', 'normal', ...
        'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1], ...
        'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1], ...
        'box', 'on');
    colormap(ax, gray);
    set(gca, 'CLim', [0 1]);
    hold('on');
    plot(xf(aa), yf(aa), 'r+', 'MarkerSize', 16, 'LineWidth', 2);
    quiver(xf(aa), yf(aa), uu, vv, 'r', 'LineWidth', 2);

    subplot(3, 2, 5);
    imagesc(z.*z_mask);
    set(gca, 'yDir', 'normal', ...
        'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1], ...
        'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1], ...
        'box', 'on');
    colormap(img.cMap);
    set(gca, 'CLim', img.cScale);
    hold('on');
    plot(xf(aa), yf(aa), 'r+', 'MarkerSize', 16, 'LineWidth', 2);
    xlim(x_zoom);
    ylim(y_zoom);
    
    ax = subplot(3, 2, 6);
    imagesc(z_mask);
    set(gca, 'yDir', 'normal', ...
        'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1], ...
        'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1], ...
        'box', 'on');
    colormap(ax, gray);
    set(gca, 'CLim', [0 1]);
    hold('on');
    plot(xf(aa), yf(aa), 'r+', 'MarkerSize', 16, 'LineWidth', 2);
    plot(xf(aa), yf(aa), 'r+', 'MarkerSize', 16, 'LineWidth', 2);
    xlim(x_zoom);
    ylim(y_zoom);

    %pause;
            %}

        end

        % Delete the last dummy value in the vector
        model = model(1:maxStep);
        xf = xf(1:maxStep);
        yf = yf(1:maxStep);
        sigma = sigma(1:maxStep);
        hf = hf(1:maxStep);
        theta = theta(1:maxStep);
        theta_norm = theta_norm(1:maxStep);
        drho = drho(1:maxStep);
        rmsd = rmsd(1:maxStep);
        exitflag = exitflag(1:maxStep);

        % Deciding if the refinement fit is good or not
        rmsd_distr = ksdensity(rmsd, 0:0.01:max(rmsd), 'Support', 'positive');
        [~, locs, ~, prominence] = findpeaks(rmsd_distr, 0:0.01:max(rmsd));
        rmsd_mode = locs(find(prominence == max(prominence), 1));
        rmsd_test = rmsd(rmsd <= rmsd_mode);
        rmsd_test2 = [rmsd_test; rmsd_mode+(rmsd_mode-rmsd_test)];
        rmsd_test = rmsd_mode+3*std(rmsd_test2);

        isgoodfit = zeros(maxStep, 1);
        isgoodfit(rmsd < rmsd_test) = 1;
        isgoodfit(1) = 1;
        isgoodfit(end) = 1;

        %{
delete(hh3);
hh3 = plot(xf, yf, 'r+-', 'LineWidth', 1);
plot(xf(logical(isgoodfit)), yf(logical(isgoodfit)), 'c+', 'LineWidth', 1);
        %}

        % Interpolate the bad segments
        % Inclusive index
        goodfitSegmentStart = find(isgoodfit(2:end)-isgoodfit(1:end-1) == -1)+1;
        goodfitSegmentEnd = find(isgoodfit(2:end)-isgoodfit(1:end-1) == 1);
        % Intrapolate rho diference
        for aa = 1:numel(goodfitSegmentStart)
            drho_final = linspace(0, 1, goodfitSegmentEnd(aa)-goodfitSegmentStart(aa)+1+2)';
            drho_final = drho(goodfitSegmentStart(aa)-1).*(1-drho_final)+...
                drho(goodfitSegmentEnd(aa)+1).*drho_final;
            drho(goodfitSegmentStart(aa):goodfitSegmentEnd(aa)) = drho_final(2:end-1);
            [xx_test, yy_test] = pol2cart(theta_norm(goodfitSegmentStart(aa):goodfitSegmentEnd(aa)), drho_final(2:end-1));
            xf(goodfitSegmentStart(aa):goodfitSegmentEnd(aa)) = xf1(goodfitSegmentStart(aa):goodfitSegmentEnd(aa))+xx_test;
            yf(goodfitSegmentStart(aa):goodfitSegmentEnd(aa)) = yf1(goodfitSegmentStart(aa):goodfitSegmentEnd(aa))+yy_test;
        end

        %plot(xf, yf, 'c-', 'LineWidth', 1);


        % Finalise centre contour line trace
        % compute the splines with added interpolation
        % for refinement, point density approx 1 pixel apart
        l = [0; sqrt((xf(2:end)-xf(1:end-1)).^2+(yf(2:end)-yf(1:end-1)).^2)];
        l = cumsum(l);
        % precision in 1/pixels
        precision = 100;
        ll = (0:1/precision:max(l))';
        xx = spline(l, xf, ll);
        yy = spline(l, yf, ll);
        % interpolate back to
        % Find points on refined curve that are closest to 1 px apart
        %l = (0:1:max(ll))';
        l = linspace(0, max(ll), ceil(max(ll)))';
        x = spline(ll, xx, l);
        y = spline(ll, yy, l);


        % Apparent width set to +/- 2*sigma = 4*sigma wide
        %currentFilament = Filament(img, x, y, [], [], 4*mean(sigma), ...
        %    {'TraceFilamentM', appWidth, zThreshold, maxSearchPhiDeg}');
        %img.features.lastFilament = currentFilament;

        tracingData = table(model, xf, yf, sigma, hf, theta, theta_norm, drho, rmsd, exitflag);
        %img.features.lastFilament = Filament(img, x, y, ...
        %    {'TraceFilamentM', appWidth, zThreshold}', ...
        %    tracingData);
        filamentNumber = numel(img.features.filaments)+1;
        img.features.filaments(filamentNumber, 1) = Filament(img, x, y, ...
            {'TraceFilamentM', appWidth, zThreshold}', ...
            tracingData);
        fprintf('Traced filament segment %g: %g pixels long, %g nm mean height.\n', ...
            filamentNumber, ...
            img.features.filaments(filamentNumber).lContour, ...
            mean(img.features.filaments(filamentNumber).z));

        % Save updated img in UserData for Trace_y UI if needed
        if exist('targetFig', 'var')
            targetFig.UserData = img;
        end

        %figure(2);
        %img.DispFilament(filamentNumber);
        %drawnow;
        %set(2, 'OuterPosition', [201 51 1260 960]);
        %currentfig = gcf;

        % This bit is not used any more. All filaments are automatically saved
        %{
%figure(1);
allAxes = findobj(currentfig.Children, 'Type', 'axes');
subplot(allAxes(end));
title(sprintf('Save filament? (y): yes, otherwise: no \n'), 'FontSize', 16, 'Color', 'r');
[~, ~, button] = ginput(1);
if button == 121 || button == 89
    % toggle, 89 = 'Y', 121 = 'y'
    img.features.filaments = ...
        [img.features.filaments; img.features.lastFilament];
    img.features.lastFilament = [];
    title(sprintf('Filament #%g saved! \n', length(img.features.filaments)), 'Color', 'k');
else
    title(sprintf('Last Filament \n'), 'Color', 'k');
end
        %}



        %
        % Get ready for the next trace
        %

        % Should not use img methods recursively
        %{
        % Save the current view
        h_ax = h_img.Parent;
        current_xlims = h_ax.XLim;
        current_ylims = h_ax.YLim;

        h_fig = img.ExploreFeatures;
        set(h_fig, 'NumberTitle', 'off', 'Name', 'Trace_y – TraceFilamentM');

        % Modify callback functions
        h_img = findobj(h_fig, 'Type', 'Image');
        h_img.ButtonDownFcn = @SetXYnodes;

        h_lines = findobj(h_fig, 'Type', 'Line');
        for aa = 1:numel(h_lines)
            h_lines(aa).ButtonDownFcn = [];
        end

        % Recall the current view
        h_ax = h_img.Parent;
        h_ax.XLim = current_xlims;
        h_ax.YLim = current_ylims;

        % Add a UI button for start
        traceBtn = uicontrol('Style', 'pushbutton', ...
            'String', 'Trace', ...
            'Units', 'pixels', 'Position', [75 70 150 25], ...
            'Callback', @TraceNow, 'Enable', 'off', 'UserData', [], ...
            'FontName', 'helvetica', ...
            'FontUnits', 'pixels', 'FontSize', 14, ...
            'FontWeight', 'normal');
        %}
        traceBtn.String = 'Trace';
        traceBtn.Enable = 'off';

        % draw final result of the previous trace
        %delete(hh_text);
        if exist('h_recentFilament', 'var') && isgraphics(h_recentFilament)
            h_recentFilament.LineWidth = 1.5;
        end
        h_recentFilament = plot(x, y, 'r-', 'LineWidth', 3);
        %drawnow

        delete(hh);
        delete(hh2);
        delete(hh3);

        xy_nodes = [];
        ll = [];
        is_node = [];
        hh = [];
        hh2 = [];
        hh3 = [];

        h_filament = img.DispFilament(filamentNumber, 'nm');
        uicontrol(h_filament, 'Style', 'pushbutton', ...
            'String', 'Delete filament', ...
            'Units', 'pixels', 'Position', [25 20 150 25], ...
            'Callback', @DeleteFilament, 'Enable', 'on', 'UserData', [], ...
            'FontName', 'helvetica', ...
            'FontUnits', 'pixels', 'FontSize', 14, ...
            'FontWeight', 'normal');

    else
        % Finishing
        traceNext = false;
        %fprintf('\n');

    end


end



% Call back for not saving the last trace. Saving is now the default
    function [] = DeleteFilament(~, ~)
        img = img.DeleteFilament(filamentNumber);
        fprintf('  Deleted!\n');
        close(h_filament);
        delete(h_recentFilament);
        %delete(h_lines(filamentNumber));
    end



end