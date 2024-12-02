function [zz, xx, yy] = StraightenFilament(img, filamentNumber)

%
% DESCRIPTION
% – Digitally straighten a filament based on its traced central axis
% – Based the method in Egelman EH. An algorithm for straightening images 
% of curved filamentous structures. Ultramicroscopy. 1986 19(4), 367–73 
%
% USAGE
% Standard method usage with filament number as input:
% >> zz = img.StraightenFilament(filamentNumber)
% Also get the full strightened xyz coordinates
% >> [zz, xx, yy] = img.StraightenFilament(filamentNumber)
%
% INPUTS
% img  –  AFMimage object. Implicit with dot notation
% filamentNumber  –  The index number of the traced filament.
%
% OUTPUTS
% [zz, xx, yy]  –  The xyz coordinates of the straightend filament image
%
% DEPENDENCIES
% – Method for Trace_y's @AFMimage/ object
%
% AUTHORS
% Liisa Lutter and Wei-Feng Xue
%
% HISTORY
% 2018.12  –  First draft by LL.
% 2019.01  –  WFX edited and incorporated method into AFMimage class.
% 2024.09  –  Minor update to make image width estimate more effective.
%



if ~exist('filamentNumber', 'var')
    filamentNumber = 1;
end

%figure(1);
%img.DispImage;
%hold on

% Importing x, y values of the fibril.
%{
if filamentNumber == 0
    x = img.features.lastFilament.x;
    y = img.features.lastFilament.y;
    l = img.features.lastFilament.l;
    %z = img.features.lastFilament.z;
    
    % Width of the straightened filament img
    w = ceil(img.features.lastFilament.appWidth);
else
    x = img.features.filaments(filamentNumber).x;
    y = img.features.filaments(filamentNumber).y;
    l = img.features.filaments(filamentNumber).l;
    w = ceil(img.features.filaments(filamentNumber).appWidth);
end
%}
x = img.features.filaments(filamentNumber).x;
y = img.features.filaments(filamentNumber).y;
l = img.features.filaments(filamentNumber).l;
% Get the sigma value from tracing data, set apparent width as 2*2sigma
if ~isempty(img.features.filaments(filamentNumber).tracingData)
    w = round(4*max(img.features.filaments(filamentNumber).tracingData{:, 4}));
else
    % Old method of using user supplied apparent width
    w = round(img.features.filaments(filamentNumber).appWidth);
end


% Calc total arclength and segment length of interpolated gridpoints

% compute the splines with added interpolation

% precision in 1/pixels
precision = 100;
ll = (0:1/precision:max(l))';
xx = spline(l, x, ll);
yy = spline(l, y, ll);

% interpolate back to
% Find points on refined curve that are closest to 1 px apart

%ll = sqrt(sum(diff(xx).^2+diff(yy).^2 ,2));
%ll = [0; cumsum(ll)];
%filamentLength = ll(end);
%idx = find(diff([-1; floor(round(ll*precision)/precision)]) == 1);
l = (0:1:max(ll))';
x = spline(ll, xx, l);
y = spline(ll, yy, l);



% Next step is to find the normals to these points
% Vectors to the relevant points
v = [x y];
normal = zeros(size(v));
% Calculate tangent vectors
tangent = v(1:end-1, :)-v(2:end, :);

%{
% 1/Distance weighting for weighted central differences 
% Points which are closer give a more accurate estimates
% may not be needed since all points are the same 1 px apart
distance = sqrt(tangent(:, 1).^2+tangent(:, 2).^2);
tangent(:, 1) = tangent(:, 1)./max(distance.^2, eps);
tangent(:, 2) = tangent(:, 2)./max(distance.^2, eps);
%}

% Calculate central differences
tangent = [0 0; tangent]+[tangent; 0 0];

% Find normalised normal vectors
%magnitude = sqrt(tangent(:, 1).^2+tangent(:, 2).^2);
%normal(:, 1)= -tangent(:, 2)./magnitude;
%normal(:, 2)= tangent(:, 1)./magnitude;
% No need to normalise
normal(:, 1)= -tangent(:, 2);
normal(:, 2)= tangent(:, 1);

% Find normal angles and points using polar tranformations
theta = ones(4*w+1, 1)*cart2pol(normal(:, 1), normal(:, 2))';
rho = flip(-2*w:2*w)'*ones(1, size(theta, 2));

% Set up new pixel coordinates
[dx, dy] = pol2cart(theta, rho);
xx = ones(4*w+1, 1)*x'+dx;
yy = ones(4*w+1, 1)*y'+dy;

% interpolate
% Crop out a sub-image z to speed up interpolation
xc = floor(min(x))-2*w:ceil(max(x)+2*w);
xc = xc(xc > 0 & xc <= img.pixelPerLineImg);
yc = floor(min(y))-2*w:ceil(max(y)+2*w);
yc = yc(yc > 0 & yc <= img.nLinesImg);
zc = img.z(yc, xc);
% calc new z values
zz = interp2(xc, yc, zc, xx, yy, 'cubic', 0);



end