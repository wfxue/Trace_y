function [s_x_sim, s_y_sim, s_z_sim, s_x_cor, s_y_cor, s_z_cor] = ...
    SimHeightImageData(t_r, t_a, p, precision, surf_f_name, surf_f_p, doPlot, x, y, z)

%
% DESCRIPTION
% – Simulate AFM image data of a object or a regular geometry based on hard
% tip-sample contacts.
% – This is an old code bus is still used for single line simulations which
% is sometimes needed and faster than simulation of a whole image. For
% other applications, SimHeightImage.m is generally used and this function
% will need updating and incorporate into SimHeightImage.m gradually
% – Part of Trace_y
%
% USAGE
% – Full inputs for standard usage
% >> [s_x_sim, s_y_sim, s_z_sim, s_x_cor, s_y_cor, s_z_cor] = ...
%   SimHeightImageData(t_r, t_a, p, precision, surf_f_name, ...
%   surf_f_p, doPlot, x, y, z);
%
% INPUTS:
% t_r  –  Tip radius, nm
% t_a  –  Tip side angle, deg
% p  –  Pixel density, nm/pixel
% precision  –  precision for surface function evaluation, expressed as
% fraction of pixel, e.g. 2 is 1/2 of a pixel, 3 is 1/3 etc
% surf_f  –  The object type to test
% surf_f_p  –  Parameter for the object, e.g. radius
% The function parameter pairs are:
%   'sphere'         : radius
%   'cylinder'       : radius
%   'circle'         : 2D radius
%   'rhombus'        : side length
%   'data2D'         : {data_x, data_y}
%   'data3D'         : {data_x, data_y, data_z}
%   'data3Dmodel'    : {data_x, data_y, data_z} supplied separatly in x, y
%                      and z inputs
% 
% OUTPUTS
% s_x_sim, s_y_sim, s_y_sim  –  Simulated x, y, z data
% s_x_cor, s_x_cor, s_x_cor  –  The original surface x, y, z coordinates
%
% DEPENDENCIES
% – Uses Optimization Matlab toolbox
% – Uses Trace_y's @AFMimage/ object defs and methods.
%
% AUTHORS
% Liisa Lutter, Wei-Feng Xue
% 
% HISTORY
% 2019.04  –  LL initial draft edited by WFX. Each of the steps were very 
% slow, so WFX added code that speed things up. Original code in 
% old AFM_Image_Sim function.
% 2019.05  –  LL drafted the aditional option to use a 3D model, with WFX 
% edited the new 3D model option and tidied up code a bit.
% 2019.07  –  Moved all pixel coordinates to go from [1 1] to introduce 
% pixeldensity and size parameters etc., and moved to dilation based 
% method, together with other edits should make things go a bit faster.
%



%tic

% WFX 2019 April: I have moved tese parameters to input and suggested
% surface as a function handle to accomodate any object to test
%{
% Define tip and pixel parameters
% These are typical values for use in tip function RCone
t_r = 2;
t_a = 18;
p = 5;
%}

% LL code
%{
    % Create a sphere
    [x, z, y] = sphere(100);
    x = x*r + 2*r;
    y = y*r + r ;
    z = z*r + 2*r ;
    
    % View the sphere
    surf(x, y, z);
    axis equal;
%}

if nargin < 7
    doPlot = 1;
end
if precision <= 1 && doPlot == 1
    doPlot = 2;
end

% The default
if nargin <= 4
    % Function for spherical dome
    surf_f_name = 'sphere';
end

switch surf_f_name
    case 'sphere'
        surf_f = @Sphere_upper_surface;
        r = surf_f_p;
        is2D = false;
        s_x = 0:p:6*r;
        s_y = 0:p:6*r;
    case 'cylinder'
        surf_f = @Cylinder_upper_surface;
        r = surf_f_p;
        is2D = false;
        s_x = 0:p:6*r;
        s_y = 0:p:6*r;
    case 'circle'
        surf_f = @Cylinder_upper_surface;
        r = surf_f_p;
        is2D = true;
        s_x = 0:p:20*r;
        s_y = 0;
    case 'rhombus'
        surf_f = @Rhombus_upper_surface;
        r = surf_f_p(1);
        is2D = true;
        s_x = 0:p:6*r;
        s_y = 0;
    case 'data2D'
        surf_f = @UserData2D_upper_surface;
        data_x = surf_f_p(1, :);
        data_z = surf_f_p(2, :);
        is2D = true;
        if min(data_x) < 0
            s_x = 0:p:max(data_x);
            s_x = [-fliplr(p:p:-min(data_x)) s_x];
        else
            s_x = min(data_x):p:max(data_x);
        end
        s_y = 0;
    case 'data3D'
        surf_f = @UserData3D_upper_surface;
        data_x = surf_f_p{1};
        data_y = surf_f_p{2};
        data_z = surf_f_p{3};
        s_x = data_x;
        s_y = data_y;
        is2D = false;
        p = s_x(2)-s_x(1);
    case 'data3Dmodel'
        surf_f = @UserData3Dmodel_upper_surface;
        data_x = x(:);
        data_y = y(:);
        data_z = z(:);
        data_x = data_x(~isnan(data_z));
        data_y = data_y(~isnan(data_z));
        data_z = data_z(~isnan(data_z));
        is2D = false;
        % s_x and s_y here are vectors defining the grid for simulated
        % image
        s_x = surf_f_p{1};
        s_y = surf_f_p{2};
        p = s_x(2)-s_x(1);
end

% WFX 2019 April: sphere() gives coords but a function instead here would
% make things easier. Check this out:



% Surface functions
    function z = Sphere_upper_surface(x, y)
        % Dome centrerd at middle of x and y range
        x0 = max(x(:))./2;
        y0 = max(y(:))./2;
        z = sqrt(r.^2-(x-x0).^2-(y-y0).^2)+r;
        z = real(z)-(imag(z) > 0).*r;
    end

    function z = Cylinder_upper_surface(x, ~)
        % Cylinder centrerd at middle of x
        x0 = max(x(:))./2;
        z = sqrt(r.^2-(x-x0).^2)+r;
        z = real(z)-(imag(z) > 0).*r;
    end

    function z = Rhombus_upper_surface(x, ~)
        % Rhombus centrerd at middle of x
        x0 = max(x(:))./2;
        z = -abs(x-x0)+2*r;
        z(abs(x-x0) > r) = 0;
    end

    function z = UserData2D_upper_surface(x, ~)
        % User test line data
        z = interp1(data_x, data_z, x, 'pchip');
    end

    function z = UserData3D_upper_surface(x, y)
        % User test surface data
        [data_x_grid, data_y_grid] = meshgrid(data_x, data_y);
        z = interp2(data_x_grid, data_y_grid, data_z, x, y, 'cubic', 0);
    end

    function z = UserData3Dmodel_upper_surface(x, y)
        % User 3D model
        % x, y are the gridded query points
        % data_x, data_y, data_z not gridded so use griddata
        warning('off', 'MATLAB:griddata:DuplicateDataPoints');
        warning('off', 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
        %z = griddata(data_x, data_y, data_z, x, y, 'cubic');
        z = griddata(data_x, data_y, data_z, x, y, 'natural');
        warning('on', 'MATLAB:griddata:DuplicateDataPoints');
        warning('on', 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
        % NaN values outside griddata bounds, i.e. outside top area are 0
        z(isnan(z)) = 0;
    end


% WFX 2019 April: This section could be much simpler based on r

%{
% Find all the pixels to sample
s_x = min(x, [], 'all')-2*p:p:max(x, [], 'all')+2*p;
s_y = min(y, [], 'all')-p:p:max(y, [], 'all')+p;

% Get layers of sphere at each pixel "line"
for i = 1:length(s_y)
    ii = s_y(i);
    M = contour(x, z, y, [ii ii]);
    layer{i} = M(:, 2:end)';
end

% For the line at the 0 isoline, add it to the matrix separately
[~, idx_f] = ismember(y(1, 1), s_y);
[~, idx_l] = ismember(y(end, 1), s_y);

if isempty(layer{idx_f})
    M = [x(1, :) ; z(1, :)];
    layer{idx_f} = M';
end
%}



% WFX 2019 April: Tip function is symetric around z so no need to define,
% input x can be seen as distance r to centre
% Instead make the tip function in 3D centred at tip_x, tip_y

%{
% Create a 2D tip

tip_x = linspace(-5, 5, 501);
tip_z = RCone(tip_x, t_r, t_a);

% View the tip
% plot(tip_x, tip_z)

% Rotate the tip in the y axis
tip_y = zeros(size(tip_z(:)))';
V = [tip_x ; tip_y ; tip_z];
Rz=[cosd(90) -sind(90) 0; sind(90) cosd(90) 0; 0 0 1];
V_rot = [Rz*V]';
tip_x_90(:, 1) = V_rot(:, 1);
tip_y_90(:, 1) = V_rot(:, 2);
tip_z_90(:, 1) = V_rot(:, 3);
    
%}


% WFX 2019 April: Stepwise simulation is slow. Also do not need to evaluate
% distance as you go, which also adds unessasary operations. Check out my
% custome solver to accomplish the same thing below


%{
% Define the height translations

s_z = max(z, [], 'all'):-0.1:0;
if s_z(end) ~= 0
    s_z(end+1) = 0;
end

% Loop through each height translation i.e. lower the tip to the surface
% until the tip makes contact with the surface. Contact is here
% defined as separation of 0.02 nm as there aren't infinite nr of defined
% points on the tip or sphere surface

contact = [];
image = [];

for i = 1:length(s_y)
    % If the pixel is before or after the object layer, then use the first
    % or last layer to find the contact point
        if s_y(i) < min(x, [], 'all')
            s1 = layer{idx_f};
            s1(:, end+1) = repmat(s_y(idx_f), size(s1, 1), 1);
        elseif s_y(i) > max(x, [], 'all')
            s1 = layer{idx_l};
            s1(:, end+1) = repmat(s_y(idx_l), size(s1, 1), 1);
        else
            s1 = layer{i};
            % Add the y-values to the sphere i.e. which layer it is
            s1(:, end+1) = repmat(s_y(i), size(s1, 1), 1);
        end
    
    for ii = s_x
        % Translate the tip in the x plane to the pixel
        if s_y(i) < min(x, [], 'all')
            t_x = tip_x_90 + ii;
            % Add y axis values to the tip
            t_y = tip_y_90 + s_y(i);
        elseif s_y(i) > max(x, [], 'all')
            t_x = tip_x_90 + ii;
            t_y = tip_y_90 + s_y(i);
        else
            t_x = tip_x + ii;
            t_y = zeros(size(t_x, 2), 1) + s_y(i);
        end
        
        for iii = s_z
            % Translate the tip in the z plane
            if s_y(i) < min(x, [], 'all') || s_y(i) > max(x, [], 'all')
                t_z = tip_z_90 + iii;
                s2 = [t_x t_z t_y];
            else
                t_z = tip_z + iii;
                s2 = [t_x' t_z' t_y];
            end
            % Look for a contact point between the layer and the tip
            [~, distance] = knnsearch(s1, s2);
            % Determine minimum distance value and the index
            [min_dist, min_idx] = min(distance);
                if min_dist <= 0.02
                    % Save the coordinates of the contact point and the point at the bottom of
                    % the tip, which would be recorded by the AFM
                    contact = [contact; s2(min_idx, :)];
                    image = [image; s2(ceil(size(tip_x, 2)/2), :)];
                    break
                elseif iii == s_z(end)
                    % If there is no contact point save the bottom of the tip at
                    % the last point
                    contact = [contact; s2(ceil(size(tip_x, 2)/2), :)];
                    image = [image; s2(ceil(size(tip_x, 2)/2), :)];
                end
        end
    end
end
%}



% Here we solve the equation Tip_lower_surface - zf = 0 at z with the
% condition that Tip_lower_surface - zf is >= 0 at z by looking at first
% collision based on predefined precision

% WFX Jul 2009
% Incorporated dialation first to generate simulated image, much faster


% Tip function, to be called for any current positions tip_x, tip_y
%{
    function z = Tip_lower_surface(x, y, tip_x, tip_y, t_r, t_a)
        % Tip centrerd at x = tip_x, y = tip_y
        % distance to tip_x, tip_y
        [~, dist] = cart2pol(x-tip_x, y-tip_y);
        z = RCone(dist, t_r, t_a);
    end
%}


if doPlot > 0
    hwb = waitbar(0, 'Simulating image, dialation...');
end


% Dialation
[s_x_sim, s_y_sim] = meshgrid(s_x, s_y);
s_z_sim = 0.*s_x_sim;
s_z = surf_f(s_x_sim, s_y_sim);

[s_tip_x, s_tip_y] = meshgrid(unique([-s_x s_x]), unique([-s_y; s_y]));
[~, tip_lower_surface] = cart2pol(s_tip_x, s_tip_y);
%tip_lower_surface = RCone(tip_lower_surface, t_r, t_a);
tip_lower_surface = TipModel_RoundedCone(tip_lower_surface, [], [t_r t_a]);
s_tip_xcentre = length(s_x);
s_tip_ycentre = length(s_y);

for aa = 1:size(s_x_sim, 2)
    for bb = 1:size(s_x_sim, 1)
        
        % Translate the tip to pixel coordinates
        tip_z_r = -tip_lower_surface(s_tip_ycentre-bb+1:end-bb+1, s_tip_xcentre-aa+1:end-aa+1);
        
        % Test code
        % From old code that re-calcs the tip surface
        %{
        tip_z_r = -Tip_lower_surface(s_x_sim, s_y_sim, s_x_sim(bb, aa), s_y_sim(bb, aa), t_r, t_a);
        %
        tip_z_r_test = -Tip_lower_surface(s_x_sim, s_y_sim, s_x_sim(bb, aa), s_y_sim(bb, aa), t_r, t_a);
        imagesc(tip_z_r_test-tip_z_r);
        colorbar;
        set(gca, 'YDir', 'normal');
        %}
        
        % No need to rotate for symetric tip
        s_z_sim = max(s_z_sim, s_z(bb, aa)+tip_z_r);
        
        %waitbar(bb./size(xf, 1), hwb);
    end
    if doPlot > 0
        waitbar(aa./size(s_x_sim, 2), hwb);
    end
end



% Refine and find contact point with finer grids
% If presision parameter is set to positive value, otherwise do not need to
% refine

if doPlot > 0
    waitbar(0, hwb, 'Simulating image, refine...');
end

if precision > 1
    % Surface function, fine grid
    xf = min(s_x(:)):p/precision:max(s_x(:));
    yf = min(s_y(:)):p/precision:max(s_y(:));
    % Tip function fine grid
    [s_tip_xf, s_tip_yf] = meshgrid(unique([-xf xf]), unique([-yf; yf]));
    [~, tip_lower_surface] = cart2pol(s_tip_xf, s_tip_yf);
    %tip_lower_surface = RCone(tip_lower_surface, t_r, t_a);
    tip_lower_surface = TipModel_RoundedCone(tip_lower_surface, [], [t_r t_a]);
    s_tip_xfcentre = length(xf);
    s_tip_yfcentre = length(yf);
    
    % Discretize surface function, fine grid: zf
    [xf, yf] = meshgrid(xf, yf);
    zf = surf_f(xf, yf);
    
    % Set up corrected image variables
    [s_x_cor, s_y_cor] = meshgrid(s_x, s_y);
    s_z_cor = s_z_sim;
    
    % Solver stats, just for checking and debugging if needed
    solver_stats = ones(size(s_z_sim));
    
    for aa = 1:size(s_x_sim, 2)
        for bb = 1:size(s_x_sim, 1)
            
            tip_zf = tip_lower_surface(s_tip_yfcentre-precision*(bb-1):end-precision*(bb-1), s_tip_xfcentre-precision*(aa-1):end-precision*(aa-1));
            
            % Test code
            % From old code that re-calcs the tip surface
            %{
            tip_zf = Tip_lower_surface(xf, yf, s_x_sim(bb, aa), s_y_sim(bb, aa), t_r, t_a);
            %
            tip_zf_test = Tip_lower_surface(xf, yf, s_x_sim(bb, aa), s_y_sim(bb, aa), t_r, t_a);
            imagesc(tip_zf_test-tip_zf);
            colorbar;
            set(gca, 'YDir', 'normal');
            %}
            
            % Fine adjustments
            diff = inf;
            while abs(diff) > 1e-4
                diff = (tip_zf+s_z_sim(bb, aa))-zf;
                diff = min(diff(:));
                s_z_sim(bb, aa) = s_z_sim(bb, aa)-diff;
            end
            
            % Find contact point
            diff = (tip_zf+s_z_sim(bb, aa))-zf;
            [y_ind, x_ind] = find(diff == min(diff(:)), 1);
            s_x_cor(bb, aa) = xf(y_ind, x_ind);
            s_y_cor(bb, aa) = yf(y_ind, x_ind);
            s_z_cor(bb, aa) = zf(y_ind, x_ind);
            solver_stats(bb, aa) = min(diff(:));
        end
        if doPlot > 0
            waitbar(aa./size(s_x_sim, 2), hwb);
        end
    end
    
end



% WFX Jul 2009: This solver below is no longer needed
%{

s_z_sim = 0.*s_x_cor;
s_z_cor = s_z_sim;

% Solver stats
solver_stats = ones(size(s_z_sim));

hwb = waitbar(0, 'Progress');

for aa = 1:length(s_x)
    for bb = 1:length(s_y)
        %counter = 0;
        % Discretize tip function: tip_zf at coordinate s_x(aa), s_y(bb)
        tip_zf = Tip_lower_surface(xf, yf, s_x(aa), s_y(bb), t_r, t_a);
        
        % initial tip location
        zz1 = min(s_zf(:))-precision;
        zz2 = max(s_zf(:))+precision;
        
        % Setup custom solver
        contact = false;
        while ~contact
            %counter = counter+1;
            %if counter == 100
            %    counter
            %end
            % This is the target function
            % find zz1 when ncontact1 is 1 (or maybe 2)
            ncontact1 = sum(sum((tip_zf+zz1) <= s_zf));
            ncontact2 = sum(sum((tip_zf+zz2) <= s_zf));
            
            if ncontact1 == 1
                % contact point found
                contact = true;
            elseif ncontact1 == ncontact2 && ncontact1 == 0
                zz1 = zz1-precision;
            elseif ncontact2 == 0 && zz2 > zz1
                % test point 2 above surface, lower towards test point 1
                zz2 = (zz2+zz1)/2;
            elseif ncontact1 == 0 && zz1 > zz2
                % test point 1 above surface, lower towards test point 2
                zz1 = (zz2+zz1)/2;
            elseif ncontact2 == 0 && zz2 < zz1
                % above surface below test point 1, lower test point 1
                zz1 = (zz2+zz1)/2;
            elseif ncontact1 == 0 && zz2 < zz1
                % above surface below test point 2, lower test point 2
                zz2 = (zz2+zz1)/2;
            elseif ncontact1 > 0 && ncontact1 < ncontact2
                % below surface above test point 2, improve test point 2
                zz2 = interp1([ncontact1 ncontact2], [zz1 zz2], 1, 'linear', 'extrap');
            elseif ncontact1 > 0 && ncontact1 > ncontact2
                % below surface below test point 2, improve test point 1
                zz1 = interp1([ncontact1 ncontact2], [zz1 zz2], 1, 'linear', 'extrap');
                %elseif ncontact1 == ncontact2 && ncontact1 > 0
                %    zz1 = zz1+precision;
                %elseif ncontact1 == ncontact2 && ncontact1 == 2
                %    % rare occurance due to discretization
                %    contact = true;
            else
                % Something went wrong
                %zz1 = (min(s_zf(:))+zz1)/2;
                %zz2 = (max(s_zf(:))+zz2)/2;
                %close(hwb);
                %error('Solver stuck in an infinite loop, try increase precision');
                % Now gives a warning at the end instead
                %warning('Solver found %g clash points', ncontact1);
                solver_stats(bb, aa) = ncontact1;
                contact = true;
            end
            
        end
        
        % Find corrected x and y and z
        [y_ind, x_ind] = find((tip_zf+zz1) <= s_zf, 1);
        s_x_cor(bb, aa) = xf(y_ind, x_ind);
        s_y_cor(bb, aa) = yf(y_ind, x_ind);
        s_z_cor(bb, aa) = s_zf(y_ind, x_ind);
        
        % Slightly correct and save z
        % s_z(bb, aa) = zz1-(tip_zf(y_ind, x_ind)+zz1-s_zf(y_ind, x_ind));
        % =>
        s_z_sim(bb, aa) = -tip_zf(y_ind, x_ind)+s_zf(y_ind, x_ind);
        %waitbar(bb./length(s_y), hwb);
        if is2D
            break
        end
    end
    waitbar(aa./length(s_x), hwb);
end
%}

% WFX July 2019: Not needed
%{

if is2D
    solver_stats = solver_stats(1, :);
end

if ~noPlot
    % Gives a warning if there are potential solver issues
    if sum(solver_stats(:)) > length(solver_stats(:))
        warning('Solver found 2 or more solutions at %g pixels (%g %%)', ...
            sum(solver_stats(:) > 1), ...
            100*sum(solver_stats(:) > 1)/length(solver_stats(:)));
        warning('Solver found 3 or more solutions at %g pixels (%g %%)', ...
            sum(solver_stats(:) > 2), ...
            100*sum(solver_stats(:) > 2)/length(solver_stats(:)));
        warning('Solver found max %g solutions at %g pixels (%g %%), try increasing precision', ...
            max(solver_stats(:)), ...
            sum(solver_stats(:) == max(solver_stats(:))), ...
            100*sum(solver_stats(:) == max(solver_stats(:)))/length(solver_stats(:)));
    end
end
%}

% WFX 2019 April: I suggest using imagesc for simulated image plots


%{
% Arrange the data points in separate x, y, z variable for surface plotting
contact_x = [];
for i = 0:length(s_y)-1
    ldp = [1:length(s_x)] + length(s_x) * i;
    contact_x = [contact_x  contact(ldp, 1)];
end

contact_z = [];
for i = 0:length(s_y)-1
    ldp = [1:length(s_x)] + length(s_x) * i;
    contact_z = [contact_z  contact(ldp, 2)];
end

contact_y = [];
for i = 0:length(s_y)-1
    ldp = [1:length(s_x)] + length(s_x) * i;
    contact_y = [contact_y  contact(ldp, 3)];
end


image_x = [];
for i = 0:length(s_y)-1
    ldp = [1:length(s_x)] + length(s_x) * i;
    image_x = [image_x  image(ldp, 1)];
end

image_z = [];
for i = 0:length(s_y)-1
    ldp = [1:length(s_x)] + length(s_x) * i;
    image_z = [image_z  image(ldp, 2)];
end

image_y = [];
for i = 0:length(s_y)-1
    ldp = [1:length(s_x)] + length(s_x) * i;
    image_y = [image_y  image(ldp, 3)];
end

% Plot the surface of the simulated AFM image and the contact point image
surf(contact_x, contact_y, contact_z); axis equal;
figure(2);
surf(image_x, image_y, image_z); axis equal;



toc
%}

if is2D
    % Get the single row from the results
    s_z_sim = s_z_sim(1, :);
    s_x_cor = s_x_cor(1, :);
    s_z_cor = s_z_cor(1, :);
    xf = xf(1, :);
    zf = zf(1, :);
end



if doPlot == 0
    return
elseif doPlot == 2
    close(hwb);
    return
else
   % doPlot == 1
    close(hwb);
end



figure;
if ~is2D
    subplot(2, 3, 1);
    surf(xf, yf, zf, 'LineStyle', 'none');
    axis('equal');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    colorbar;
    
    ax2 = subplot(2, 3, 2);
    imagesc(s_x, s_y, s_z_sim);
    set(gca, 'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1]);
    set(gca, 'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1]);
    set(gca, 'yDir', 'normal');
    %set(gca, 'CLim', [0 1]);
    title('Simulated image');
    
    ax3 = subplot(2, 3, 3);
    % interpolate nn
    s_zf_cor = griddata(s_x_cor, s_y_cor, s_z_cor, xf, yf, 'nearest');
    imagesc(xf(1, :), yf(:, 1), s_zf_cor);
    set(gca, 'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1]);
    set(gca, 'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1]);
    set(gca, 'yDir', 'normal');
    %set(gca, 'CLim', [0 1]);
    title('Corrected image');
    
    
    ax4 = subplot(2, 3, 4);
    imagesc(xf(1, :), yf(:, 1), zf);
    set(gca, 'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1]);
    set(gca, 'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1]);
    set(gca, 'yDir', 'normal');
    title('The specimen surface');
    
    
    ax5 = subplot(2, 3, 5);
    hold('on');
    % Corrected grid lines
    for aa = 1:length(s_x)
        plot(s_x_cor(:, aa), s_y_cor(:, aa), '-k');
    end
    for bb = 1:length(s_y)
        plot(s_x_cor(bb, :), s_y_cor(bb, :), '-k');
    end
    set(gca, 'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1]);
    set(gca, 'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1]);
    box('on');
    title('Corrected grid lines');
    
    %{
    imagesc(xf(1, :), yf(:, 1), solver_stats);
    set(gca, 'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1]);
    set(gca, 'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1]);
    set(gca, 'yDir', 'normal');
    title('Solver Stats');
    %}
    
    ax6 = subplot(2, 3, 6);
    % fine interpolate
    s_zf_cor = griddata(s_x_cor, s_y_cor, s_z_cor, xf, yf, 'v4');
    imagesc(xf(1, :), yf(:, 1), s_zf_cor);
    set(gca, 'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1]);
    set(gca, 'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1]);
    set(gca, 'yDir', 'normal');
    %set(gca, 'CLim', [0 1]);
    title('Corrected interpolated image');
    
    linkaxes([ax2, ax3, ax4, ax5, ax6], 'xy');
    colormap(ax4, 'gray');
    
else
    % Is 2D
    subplot(1, 2, 1);
    h1 = plot(xf, zf, 'k-');
    hold('on');
    for aa = 1:length(s_z_sim)
        %tip_zf = tip_lower_surface+s_z_sim(aa);
        tip_zf = tip_lower_surface(s_tip_yfcentre:end, s_tip_xfcentre-precision*(aa-1):end-precision*(aa-1))+s_z_sim(aa);
        h2 = plot(xf, tip_zf, '-', 'Color', [0.8 0.8 0.8]);
    end
    h3 = plot(s_x, s_z_sim, 'o');
    h4 = plot(s_x_cor, s_z_cor, '+');
    plot([s_x; s_x_cor], [s_z_sim; s_z_cor], 'r-');
    xlim([min(s_x) max(s_x)]);
    
    ylim([min(s_z_sim)-0.25*(max(s_x)-min(s_x)) min(s_z_sim)+0.75*(max(s_x)-min(s_x))]);
    set(gca, 'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1]);
    set(gca, 'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1]);
    legend([h1 h2 h3 h4], ...
        'surface', 'tip', 'pixel value', 'corrected pixel value');
    
    
    subplot(1, 2, 2);
    plot(xf, zf, 'k-');
    hold('on');
    plot(s_x, s_z_sim, '-o');
    plot(s_x_cor, s_z_cor, '-+');
    xlim([min(s_x) max(s_x)]);
    ylim([min(s_z_sim)-0.25*(max(s_x)-min(s_x)) min(s_z_sim)+0.75*(max(s_x)-min(s_x))]);
    %ylim(xlim+3*r);
    set(gca, 'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1]);
    set(gca, 'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1]);
    legend('surface', 'pixel values', 'corrected pixel values');
end



end