function [xm, ym, zm, theta_m, rho_m, x_untwisted, z_untwisted, twist_angle_y] = ...
    FilamentModel(x_data, y_data, z_data, centre_axis, ...
    symmetry, handedness, periodicity, peaks_loc, w, px, smoothness, y_refine)

%
% DESCRIPTION
% – This function creates a 3D model of a fibril imaged with AFM.
% – One of the main routines used by MakeHelicalFilamentModel.m
% – Part of Trace_y
%
% USAGE
% Standard method usage with full inputs and outputs
% >> [xm, ym, zm, theta_m, rho_m, x_untwisted, z_untwisted, twist_angle_y] = ...
%   FilamentModel(x_data, y_data, z_data, centre_axis, ...
%   symmetry, handedness, periodicity, peaks_loc, w, px, smoothness, y_refine)
%
% INPUTS
% x_data, y_data, z_data  –  x, y and z values of a fibril line(s)
% centre_axis  –  [x, z] coordinates estimate of the centre_axis location
% symmetry  –   Cross-sectional symetry number, e.g. 2, 3 or more screw 
% axis, loosely but not always related to the nr of protofilaments in the 
% fibril
% handedness  –  Helical filament twist handedness ('left' or 'right')
% periodicity  –  Initial estimate of pitch length in pixels, this method 
% of course only works for helical filaments with twist
% peaks_loc  –  The coordinates of peaks in the central line.
% w  –  Optional weights for the pixel data
% px  –  Image resolution in nm/pixel
% smoothness, y-refine  –  Controls how smooth the spline interpolation is
% in the cross-section and in the filament axis.
%
% OUTPUTS
% xm, ym, zm  –  The model coordinates in nm in cartesian coordinates
% theta_m, rho_m  –  The model coordinates in rad and nm respectively in
% helical polar coordinates that rotates around the central filament axis 
% acording to peaks_loc data.
% x_untwisted, z_untwisted, twist_angle_y  –  The model cross-section 
% coordinates in cartesian coordinates with the x/z plane that rotates 
% around the central filament y-axis acording to peaks_loc data
%
% DEPENDENCIES
% – Used by MakeHelicalFilamentModel.m.
% – Uses Curve Fitting, Statistical and Machine Learning, Image Processing, 
% Optimization, Signal Processing Matlab toolboxes.
% – Method for Trace_y's @AFMimage/ object.
%
% AUTHORS
% Liisa Lutter, Wei-Feng Xue
%
% HISTORY
% 2019.05  –  LL initial draft 3D model construction solely based on 
% central line profile. With WFX edit to simplify the central line model, 
% using a moving window approach as described in Lutter et al 2020.
% 2019.07  –  WFX edit for using corrected x, y and z coordinates in 
% multiple pixel lines.
% 2019.08  –  WFX cleaned up the code, incorporated least squares spline
% fitting, updated the algorithm for calculations in helical polar
% coordinates that twists acording to peaks_loc data and incorporated in 
% updated MakeHelicalFilamentModel.m
% 2024.09  –  WFX Trace_y update and minor updated with 
% MakeHelicalFilamentModel.m changes.
%



% Optional weights
if ~exist('w', 'var') || isempty(w)
    w = ones(size(z_data));
end

% WFX July 2019: new way to find centre of filament, i.e the centre
% coordinate on corrected (not x gridded) coordinates is now input
x = x_data-centre_axis(1);
y = y_data;
z = z_data-centre_axis(2);
%centre_line = ceil(size(y, 2)/2);




% This is the old approach
%{
% WFX: calc the angle shift per pixel in deg depeding on
% periodicity in pixels and symetry
delta_angle = 360/(periodicity*symmetry);
angle = [90:delta_angle:270 90:-delta_angle:-90];
angle = unique(angle)';

% WFX: Not sure if I got the correct handedness
if strcmp(handedness, 'left')
    angle = flip(angle);
end
% delta_angle used for rotation around filment centre
delta_angle = (angle-90)*ones(1, size(z, 2));

% Vectors from filament centre at (0, 0) in the x, z, plane in polar
% coordinates
% y is the same in each row (slice)
[theta, rho] = cart2pol(x, z);

% For moving window approach, this is the final number of slices
y_ind = (ceil(length(angle)/2):length(z)-floor(length(angle)/2))';

% Pre-allocate variables
% Model coords xm, ym and zm have number of rows as y-slices and in each
% y-slice there are length(angle)*size(z, 2) numebr of coordinates from the
% data points used
xm = zeros(length(y_ind), 361);
ym = xm;
zm = xm;
theta_m = zeros(length(y_ind), 361);
rho_m = theta_m;
xarea = zeros(length(y_ind), 1);


% For moving window approach, this is the final number of slices
y_ind = (ceil(length(angle)/2):length(z)-floor(length(angle)/2))';

% Pre-allocate variables
% Model coords xm, ym and zm have number of rows as y-slices and in each
% y-slice there are length(angle)*size(z, 2) numebr of coordinates from the
% data points used
xm = zeros(length(y_ind), 361);
ym = xm;
zm = xm;
theta_m = zeros(length(y_ind), 361);
rho_m = theta_m;
xarea = zeros(length(y_ind), 1);



% WFX July 2019:
% Lots of data points, therefore smooth cross-section of the slice by
% Least square cubic spline estimation as an option
% Number of spline pieces is npieces
% Order: k = 4 for cubic
k = 4;

if exist('npieces', 'var') && npieces == 0
    %waitbar(0, hwb, 'Optimising least squares splines...');
    % For spline fitting, estimate best number of pieces to use using AICc
    % scoring approach
    % x3 comes from the workaround for periodic fits
    ndata = 3*length(delta_angle(:));
    npieces_max = min(ceil(size(delta_angle, 1)/3)+4, 30);
    aicc = zeros(npieces_max, 1);
    
    for npieces = 1:npieces_max
        for aa = 1:length(y_ind)
            % For each slice
            theta_slice = theta(aa:aa+length(angle)-1, :)+deg2rad(delta_angle);
            rho_slice = rho(aa:aa+length(angle)-1, :);
            w_slice = w(aa:aa+length(angle)-1, :);
            
            % Merge data for least square cubic spline estimation
            % Can't specify periodic edge condition so fit 3 repeating periods
            % as workaround.
            theta_slice_fit = [theta_slice(:)-2*pi; theta_slice(:); theta_slice(:)+2*pi];
            rho_slice_fit = [rho_slice(:); rho_slice(:); rho_slice(:)];
            w_slice_fit = [w_slice(:); w_slice(:); w_slice(:)];
            
            % Define knots, since periodic workaround above, npieces is x3
            knots = augknt(-2.5*pi:6*pi/(npieces*3):3.5*pi, k);
            % Make spline (B-form)
            sp = spap2(knots, k, theta_slice_fit, rho_slice_fit, w_slice_fit);
            % Get spline between -0.5pi to 1.5 pi, i.e. centred on pi/2 (top)
            
            theta_ind = find((theta_slice_fit >= -0.5*pi) & (theta_slice_fit < 1.5*pi));
            rho_slice_residual = fnval(sp, theta_slice_fit(theta_ind))-rho_slice_fit(theta_ind);
            
            % Calculating AICc score and sum all scores form the slices up
            ncoeff = sp.number;
            ssr = 3*sum(rho_slice_residual.^2);
            aicc(npieces) = aicc(npieces)+ndata.*log(ssr/ndata)+2*ncoeff+(2*ncoeff*(ncoeff+1))/(ndata-ncoeff-1);
            
        end
        %waitbar(npieces./npieces_max, hwb);
    end
    
    % Choose final npieces with minimum AICc score
    npieces = find(aicc == min(aicc), 1, 'last');
    
    % Test code
    %{
    plot(1:npieces_max, aicc, '+-');
    hold('on');
    ylims = ylim;
    plot([npieces npieces], [aicc(npieces) ylims(1)], '-or');
    %pause;
    %}
end



% Evaluating one slice at time
for aa = 1:length(y_ind)
    % Working in cylindrical coordinates, do for each slice
    theta_slice = theta(aa:aa+length(angle)-1, :)+deg2rad(delta_angle);
    rho_slice = rho(aa:aa+length(angle)-1, :);
    w_slice = w(aa:aa+length(angle)-1, :);
    
    if exist('npieces', 'var') && npieces > 0
        % Merge data for least square cubic spline estimation
        % Can't specify periodic edge condition so fit 3 repeating periods
        % as workaround.
        theta_slice_model = [theta_slice(:)-2*pi; theta_slice(:); theta_slice(:)+2*pi];
        rho_slice_model = [rho_slice(:); rho_slice(:); rho_slice(:)];
        w_slice_fit = [w_slice(:); w_slice(:); w_slice(:)];
        
        % Define knots, since periodic workaround above, npieces is x3
        knots = augknt(-2.5*pi:6*pi/(npieces*3):3.5*pi, k);
        % Make spline (B-form)
        sp = spap2(knots, k, theta_slice_model, rho_slice_model, w_slice_fit);
        
        % Interpolating with the resulting spline
        % Get spline between -0.5pi to 1.5 pi, i.e. centred on pi/2 (top)
        % Doing 1 deg spacing per point
        theta_slice_model = -0.5*pi:pi/180:1.5*pi;
        rho_slice_model = fnval(sp, theta_slice_model);
        [x_slice_model, z_slice_model] = pol2cart(theta_slice_model, rho_slice_model);
    else
        % Averaging over multiple lines
        % First treat each line
        theta_slice_model = -0.5*pi:pi/180:1.5*pi;
        rho_slice_model = 0.*theta_slice_model;
        
        for bb = 1:size(theta_slice, 2)
            theta_slice_pp = [theta_slice(:, bb)-2*pi; theta_slice(:, bb); theta_slice(:, bb)+2*pi];
            [theta_slice_pp, idx] = sort(theta_slice_pp);
            rho_slice_pp = [rho_slice(:, bb); rho_slice(:, bb); rho_slice(:, bb)];
            rho_slice_pp = rho_slice_pp(idx);
            %pp = csaps(theta_slice_pp, rho_slice_pp, 0.5);
            %rho_slice_model = fnval(pp, theta_slice_model);
            rho_slice_model = rho_slice_model+pchip(theta_slice_pp, rho_slice_pp, theta_slice_model);
            
            % Test code
            %{
            clf;
            plot(theta_slice(:, bb), rho_slice(:, bb), '+-');
            hold('on');
            plot(theta_slice_model, rho_slice_model, '-');
            pause;
            %}
            
        end
        rho_slice_model = rho_slice_model./size(theta_slice, 2);
        [x_slice_model, z_slice_model] = pol2cart(theta_slice_model, rho_slice_model);
        
    end
    
    % Cross sectional area by polar integration
    xarea(aa) = 0.5*sum((pi/180)*rho_slice_model(1:360).^2);
    
    % Some test code
    %{
    clf;
    subplot(1, 2, 1)
    plot(theta_slice, rho_slice, '+');
    hold('on');
    plot(theta_slice_model, rho_slice_model, '-');
    subplot(1, 2, 2)
    [x_slice, z_slice] = pol2cart(theta_slice, rho_slice);
    plot(x_slice, z_slice, '+');
    axis('equal');
    hold('on');
    plot(x_slice_model, z_slice_model, '-');
    pause;
    %}
    
    % Writes the coordinates of the slice
    xm(aa, :) = x_slice_model;
    ym(aa, :) = y(y_ind(aa), 1);
    zm(aa, :) = z_slice_model;
    theta_m(aa, :) = theta_slice_model;
    rho_m(aa, :) = rho_slice_model;
    
    %waitbar(aa./length(y_ind), hwb, 'Constructing 3D model...');
end

xm_old = xm;
ym_old = ym;
zm_old = zm;
%}

%{
surf(xm, ym, zm, 'EdgeColor', 'none', 'AmbientStrength', 0.6);
axis('equal');
camlight('right');
set(gca,'visible', 'off');
view(0, 10);
%}



% New version
% Try new approach, gridding in cylindrical rotating axis

% "Untwisting" first

% Periodicity in nm = periodicity in pixels * resolution in nm/pixel
periodicity_nm = periodicity*px;

nperiods = ceil(length(peaks_loc)/symmetry);

% Twist in deg per nm
% The nodes are the peak positions, set up for getting local twist angles
% Symmetry defines which peak to set for each full turn, every peak for
% symmetry = 1, every other peak for symmetry = 2 etc
twist_nodes = zeros(nperiods, 1);
twist_loc = zeros(nperiods, 1);
for aa = 1:nperiods
    twist_nodes(aa) = 360*(aa-1);
    % Location in nm = location in pixels * resolution in nm/pixel
    twist_loc(aa) = peaks_loc(1+(aa-1)*symmetry)*px;
end
% Add initial and end point using general mean periodicity value
twist_loc = [min(y(:)); twist_loc; max(y(:))];
twist_nodes = [...
    (twist_loc(1)-twist_loc(2))*360/(periodicity_nm*symmetry)+twist_nodes(1); ...
    twist_nodes; ...
    (twist_loc(end)-twist_loc(end-1))*360/(periodicity_nm*symmetry)+twist_nodes(end)];

%twist_angle = spline(twist_loc, twist_angle, 0:max(y(:)));
%twist_angle = 360/(periodicity_nm*symmetry);

% WFX: Again not sure if I got the correct handedness
% Left handed is counter-clockwise, counter-clockwise is angle increase
% Right handed is clockwise, clockwise is angle decrease
if strcmp(handedness, 'right')
    twist_nodes = -twist_nodes;
end

%untwist_angle = y.*(-twist_angle);

% Local twist angles
%twist_angle = spline(twist_loc, twist_nodes, y);
%sp_twist = csapi(twist_loc, twist_nodes);
%twist_angle = fnval(sp_twist, y);
% pchip works the best, no problematic oschilations
twist_angle = pchip(twist_loc, twist_nodes, y);
twist_angle_y = deg2rad(twist_angle);

% Test code
%{
figure;
periodicity_nm = periodicity*px;
plot(twist_loc, twist_nodes, '+');
hold('on');
plot(y, twist_angle, '-');
plot(y(:), y(:).*360/(periodicity_nm*symmetry)+min(twist_angle(:)), 'k-');
%}

% Vectors from filament centre at (0, 0) in the x, z, plane in polar
% coordinates
% All coordinate points are in vector format, x, y, z in cartesian and
% theta, rho, y in cylindrical coordinates
[theta, rho] = cart2pol(x, z);
%theta_untwisted = theta+deg2rad(untwist_angle);

% "Untwisting" by taking away the twist angle, note conversion to rad units
% as twist angles are in deg but theta are in rad
theta_untwisted = theta-deg2rad(twist_angle);
[x_untwisted, z_untwisted] = pol2cart(theta_untwisted, rho);

% Set up variables
% For y, any y coordinate to calc can be used now, here every nm
% Set up the data in triplicate for continuity for the middle during spline
% fitting
ym = (0:(px/y_refine):max(y_data(:)))'*ones(1, 361);
theta_m = deg2rad(ones(size(ym, 1), 1)*(-180:1:180));
%rho_m = zeros(size(ym));

% The corresponding twist angles in the rotating polar cylindrical coords
twist_angle_m = pchip(twist_loc, twist_nodes, ym);


% Setup for spline fitting using spap2

% LL introduced smoothness for smoothing filament model using spap2 to
% define nknots based on pixel size. 
% WFX simplifid to just 1 based smoothness. 1 is no additonal smoothing
%smoothness = 4;
% Order: k = 4 for cubic
k = 4;

% Knots set up still uses the global periodicity value, should be ok
delta_angle = 360/(periodicity*symmetry);
% Alternatively go with the following for mean pixel angular diffeerence in
% a line. This is biggeer so retain delta angle from line to line I think
%delta_angle = rad2deg(mean(diff(mean(theta, 1))));

s_angle = delta_angle*smoothness;
th_knt = [-3*pi:deg2rad(s_angle):3*pi 3*pi];

knotsx = optknt(th_knt, k);
coefsy = zeros(size(ym, 1), length(knotsx)-k);

hwb = waitbar(0, 'Constructing 3D model...');

for aa = 1:size(ym, 1)
    
    % Get the current angle in the cylindrical coordinates
    %y_current = ym(aa, 1);
    %twist_angle_current = fnval(sp_twist, y_current);
    %twist_angle_current = pchip(twist_loc, twist_nodes, y_current);
    twist_angle_current = twist_angle_m(aa, 1);
    
    % Find all data points in an angular distance to the y current location
    % Go with +/- full period to be able to estimate the edges indepedent
    % of the symmetry constraints
    idx = find(abs(twist_angle-twist_angle_current) < 360);
    x_current = x_untwisted(idx);
    z_current = z_untwisted(idx);
    
    % Weight data based on how far the points is from the current location
    % Weight w is already indication of the data point from uncertainty of
    % deconvolution
    delta_angle = abs(twist_angle(idx)-twist_angle_current);
    w_slice = w(idx);
    % (cos(x/4))^c gives 1 at x = 0 and 0 at x = 2*pi
    % c = 4 seems to give a good compromise between fall off from 180 to 
    % 360 but not overly discount around +/- 360
    % This is a bit arbitrary but seems to work practically, c can be
    % changed
    w_slice = w_slice.*(cosd(0.25*delta_angle)).^4;
    
    [theta_slice, rho_slice] = cart2pol(x_current, z_current);
    theta_current = [theta_slice-2*pi; theta_slice; theta_slice+2*pi];
    rho_current = [rho_slice; rho_slice; rho_slice];
    w_current = [w_slice; w_slice; w_slice];
    
    sp = spap2(knotsx, k, theta_current, rho_current, w_current);
    coefsy(aa, :) = fnbrk(sp, 'coefficients');
    
    
    % Test code
    %{
    clf;
    rho_m(aa, :) = fnval(sp, theta_m(aa, :));
    subplot(1, 2, 1)
    plot(theta_current, rho_current, '.');
    hold('on');
    plot(theta_m(aa, :), rho_m(aa, :), '-');
    
    subplot(1, 2, 2)
    [xm, zm] = pol2cart(theta_m(aa, :), rho_m(aa, :));
    plot(x_current, z_current, '.');
    axis('equal');
    hold('on');
    plot(xm, zm, '-');
    [xx, zz] = pol2cart(deg2rad(twist_angle_current), mean(rho(:)));
    plot([0 xx], [0 zz], 'k-'); 
    drawnow;
    %}
    
    % 
    if mod(aa, 100) == 0
        waitbar(aa./size(ym, 1), hwb, 'Constructing 3D model...')
    end
    
end
close (hwb);

% Spline fitting the second dimension (y)
% LL define y smoothness based on pixel size
%y_step = px*smoothness;
% WFX suggest on a quarter of periodicity 
y_step = smoothness*periodicity/4;
% In practice, y smoothness has very little effect seems due to the noise
% is mainly in the angular dimension
y_knt = [0:y_step:max(y_data(:)) max(y_data(:))];
knotsy = optknt(y_knt, k);
sp2 = spap2(knotsy, k, ym(:, 1), coefsy.');
coefs = fnbrk(sp2, 'c').';

rho_m = spcol(knotsy, k, ym(:, 1))*coefs*spcol(knotsx, k, theta_m(1, :)).';
% Harmonise -180° and 180° point
rho_m(:, 1) = mean(rho_m(:, [1 end]), 2);
rho_m(:, end) = rho_m(:, 1);

% Test code
%{
figure;
imagesc(theta_m(1, :), ym(:, 1), rho_m);
set(gca, 'YDir', 'normal');
colormap('jet');
colorbar;

figure;
[xm, zm] = pol2cart(theta_m, rho_m);
plot(xm, zm);
%}


% Not needed
%{
y_window = 0:periodicity_nm*symmetry/2;
w_window = y_window.*(0.5*pi/y_window(end));
y_window = unique([-y_window y_window]');
w_window = cos(unique([-w_window w_window]')).^2;

x_untwisted = repmat(x_untwisted, length(y_window), 1);
z_untwisted = repmat(z_untwisted, length(y_window), 1);
[theta_untwisted, rho_untwisted] = cart2pol(x_untwisted, z_untwisted);

y_untwisted = repmat(y, length(y_window), 1);
w_untwisted = repmat(w, length(y_window), 1);

row_start = 1;
row_length = size(y, 1);
for aa = 1:length(y_window)
    y_untwisted(row_start:row_start+row_length-1, :) = y+y_window(aa);
    w_untwisted(row_start:row_start+row_length-1, :) = w.*w_window(aa);
    row_start = row_start+row_length;
end


ym = (0:1:max(y_data(:)))'*ones(1, 3*361-2);
theta_m = deg2rad(ones(size(ym, 1), 1)*(-180-360:1:180+360));
%rho_m = zeros(size(theta_m));

%rho_m = griddata(theta_untwisted(:), y_untwisted(:), rho_untwisted(:), theta_m, ym, 'cubic');
interpolant = scatteredInterpolant(...
    [theta_untwisted(:)-2*pi; theta_untwisted(:); theta_untwisted(:)+2*pi], ...
    [y_untwisted(:); y_untwisted(:); y_untwisted(:)], ...
    [rho_untwisted(:); rho_untwisted(:); rho_untwisted(:)], ...
    'natural', 'linear');
rho_m = interpolant(theta_m, ym);

interpolant = scatteredInterpolant(...
    [theta_untwisted(:)-2*pi; theta_untwisted(:); theta_untwisted(:)+2*pi], ...
    [y_untwisted(:); y_untwisted(:); y_untwisted(:)], ...
    [w_untwisted(:); w_untwisted(:); w_untwisted(:)], ...
    'natural', 'linear');
wm = interpolant(theta_m, ym);
%wm(wm < 0) = min(wm(wm > 0));

% Smoothing
rho_m = SplineSurfaceSmoothing(theta_m, rho_m, ym, wm, px, periodicity, symmetry, 2);
%}


% Twisting back by adding back the twist angle
%retwist_angle = ym.*twist_angle;
%twist_angle = fnval(sp_twist, ym);
%twist_angle = pchip(twist_loc, twist_nodes, ym);
theta_retwisted = theta_m+deg2rad(twist_angle_m);
[xm, zm] = pol2cart(theta_retwisted, rho_m);


% Test code
% Can turn on and off the following segments but need the old section

% All of the slices
%{
for aa = 1:size(ym, 1)
    clf;
    plot(xm(aa, :), zm(aa, :), '-');
    axis('equal');
    drawnow;
end
%}

% Compared to the old method and data
%{
y_row_data = round(mean(y_data, 2));
y_row_old = round(mean(ym_old, 2));
for aa = 1:length(y_row_old)
    clf;
    plot(xm(y_row_old(aa)+1, :), zm(y_row_old(aa)+1, :), '-');
    hold('on');
    plot(xm_old(aa, :), zm_old(aa, :), '+');
    
    idx = find(abs(y_row_data-y_row_old(aa)) == min(abs(y_row_data-y_row_old(aa))), 1);
    plot(x_data(idx, :)-centre_axis(1), z_data(idx, :)-centre_axis(2), 'r+');
    
    axis('equal');
    drawnow;
end
%}

% Compared to the data
%{
for aa = 2:size(y_data, 1)
    clf;
    plot(xm(round(mean(y_data(aa, :)))+1, :), zm(round(mean(y_data(aa, :)))+1, :), '-');
    hold('on');
    plot(x_data(aa, :)-centre_axis(1), z_data(aa, :)-centre_axis(2), '+');
    axis('equal');
    pause;
end
%}



end