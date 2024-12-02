function z = SurfaceModel_EMmap(x, y, p)

%
% DESCRIPTION
% – Part of SimHeightImage.m for calc top of EM map models
% – Part of Trace_y
%
% AUTHORS
% Wei-Feng Xue and Liisa Lutter
% 
% HISTORY
% Updated with SimHeightImage.m and see notes in SimHeightImage.m
%



xyz = p{1};
voxSize = p{2};

% Construct alpha shape
% 2*voxel size to save a bit of computation as tip radius likely to
% be a lot bigger
% Also use alphashape to find largest region
% Basically suppresing noise by using region 1 (the biggest region)
alpha_val = 2*voxSize;

%em_extended = alphaShape(p(:, 1), p(:, 2), p(:, 3), tip_radius, ...
%    'HoleThreshold', max(mapfv_aligned.UserData.CELLA./10).^2', 'RegionThreshold', 0.01);
p = alphaShape(xyz(:, 1), xyz(:, 2), xyz(:, 3), alpha_val, ...
    'HoleThreshold', 0, 'RegionThreshold', 10*voxSize.^3);


% Find upper boundaries for region 1
[~, p_top] = boundaryFacets(p, 1);
z_step = p.Alpha/2;


z = zeros(size(x));
z_lb = zeros(size(x));
z_ub = z_lb+max(p_top(:, 3));
%z_b = logical(z);
z_b = (z > 0);

for zz = fliplr([0:z_step:max(p_top(:, 3)) max(p_top(:, 3))])
    z_test = zeros(size(x))+zz;
    tf = inShape(p, x, y, z_test, 1);
    z_ub(~tf & ~z_b) = z_test(~tf & ~z_b);
    z_lb(tf & ~z_b) = z_test(tf & ~z_b);
    
    z_b(tf & ~z_b) = true;
end

% Refine
while sum(z_ub(:)-z_lb(:) > z_step/100) > 0
    % voxSize/100 is the tolerance
    z_test = z_lb+(z_ub-z_lb)./2;
    tf = inShape(p, x, y, z_test);
    z_ub(~tf) = z_test(~tf);
    z_lb(tf) = z_test(tf);
end

% Use uper bound value on top and lower bound values (incl 0) at bottom
z(z_test > max(p_top(:, 3))/2) = z_ub(z_test > max(p_top(:, 3))/2);
z(z_test <= max(p_top(:, 3))/2) = z_lb(z_test <= max(p_top(:, 3))/2);

end
