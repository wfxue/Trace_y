function y = TipModel_RoundedCone(x, y, tip_par)

%
% DESCRIPTION
% – A tip geometry function.
% – The tip geometry is a truncated cone with a matching spherical cap. The
% Rounded cone is centred on x = 0 and lowest point at y = 0. Used for 3D 
% symmetric tip in polar coordinates with x = r as only coordinate input.
% – Part of Trace_y
%
% USAGE
% Standard function usage
% >> y = TipModel_RoundedCone(x, y, tip_par);
% Symmetric tip in polar coordinates with x is the radius from the centre
% >> y = TipModel_RoundedCone(x, [], tip_par);
%
% INPUTS
% x, y  –  x and y values
% tip_par  –  The parameters for the tip function [tip_r tip_a]. tip_r is 
% the tip radius, scalar in nm, and tip_a is the tip side angle to the tip 
% centre axis, scalar in degrees
%
% OUTPUTS
% y  –  Tip surface (technically z)
%
% DEPENDENCIES
% – Used for image simulations in Trace_y
%
% AUTHORS
% Wei-Feng Xue
% 
% HISTORY
% 2020.05  –  Separated out from other code files and added for 
% convenience as used simulations
%



tip_r = tip_par(1);
tip_a = tip_par(2);

if ~isempty(y)
    % Convert to polar coordinates
    [x, y] = meshgrid(x, y);
    [~, x] = cart2pol(x, y);
end
% Otherwise, x is the distance to centre in polar coordinates


tip_a = deg2rad(tip_a);

% the cone
yc = abs(x.*(cos(tip_a)./sin(tip_a)));
% the rounded tip
yr = -sqrt(tip_r.^2-x.^2)+tip_r./sin(tip_a);

% the tangent point
xx = tip_r.*cos(tip_a);

% complete function and move down so y = 0 at x = 0
y = yc.*(x <= -xx | x >= xx)+yr.*(x > -xx & x < xx)-tip_r./sin(tip_a)+tip_r;

% Test code
%{
% Check result on plot
plot(x, yc, '-b');
hold on
plot(x, y, '+r');

set(gca, ...
    'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1],...
    'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1]);
%}
end