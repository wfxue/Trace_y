function [img] = Dilation(img, precision, do_wb)

%
% DESCRIPTION
% – Simulation of an AFM image with known surface S and tip T using approx
% the standars dilation approach described in Villarrubia (1997) 
% https://doi.org/10.6028/jres.102.030
% – Assumes S and T coordinates or functions are already set in Trace_y 
% img object
% – Part of Trace_y
%
% USAGE
% – To use method, set S and T in an instance of @AFMimage/ object, e.g. 
% in img, then run with defaults:
% img.Dilation;
% – Can set pixel precision with the following:
% img.Dilation(precision);
% – Also turns off the progress bar:
% img.Dilation(precision, false);
%
% INPUTS
% img  –  An instance of Trace_y @AFMimage/ AFM image object.
% precision  –  Subdivisions of pixel for better dilation calc etc. Default
% is 4.
% do_wb  –  Flag for turning on/off the visual wait bar.
%
% OUTPUTS
% img  –  An instance of Trace_y @AFMimage/ AFM image object. The dilated
% image coordinates are stored in xC, yC and zC in img.
%
% DEPENDENCIES
% – Uses Trace_y's @AFMimage/ object defs and methods.
%
% AUTHORS
% Wei-Feng Xue
%
% HISTORY
% 2019.06  –  Initial code draft AFM_Image_Correction containing dilation 
% and reconstruction algorithms
% 2020.05  –  Code edit and incorporated into Trace_y @AFMimage/ object, 
% i.e. the img object
% 2021.03  –  Separated dilation and the reconstruuction algorithms, this 
% is for better organisation of data in the img object
% 2023.11  –  Code edited with the aim of some added visualisation method 
% for img object
% 2024.08  –  Edited to reduce memory and increase speed by removing the
% need for many expensive AlphaShape evaluations needed for if zS contain 
% EM map coordinates. Also fixed rare bugs regarding the range of x/y 
% vectors upon initiation
%

if ~exist('precision', 'var') || isempty(precision)
    precision = 4;
end

if ~exist('do_wb', 'var') || isempty(do_wb)
    do_wb = true;
end


% Generate simulated surface S if needed
[xS, yS] = meshgrid(img.x, img.y);

% Surface S, fine grid
finer_grid_spacing = (img.x(2)-img.x(1))/precision;
%xSf = img.x(1):finer_grid_spacing:img.x(end)+eps;
xSf = img.x(1):finer_grid_spacing:img.x(end)+0.5*(finer_grid_spacing);
%ySf = (img.y(1):finer_grid_spacing:img.y(end)+eps)';
ySf = (img.y(1):finer_grid_spacing:img.y(end)+0.5*(finer_grid_spacing))';

% Discretize surface function, fine grid: zf
[xSf, ySf] = meshgrid(xSf, ySf);

% If S is stored as function
if iscell(img.zS)
    zSf = img.zS{1}(xSf, ySf, img.zS{2});
else
    zSf = interp2(xS, yS, img.zS, xSf, ySf, 'spline');
end


% If S is stored as function
%{
if iscell(img.zS)
    zS = img.zS{1}(xS, yS, img.zS{2});
else
    zS = img.zS;
end
%}
% Just want to evaluate zS = img.zS{1}(xS, yS, img.zS{2}) once, as there is
% a cost of big alpha shape in there. So just do the finer grid and get the
% normal grid from there.
if iscell(img.zS)
    zS = interp2(xSf, ySf, zSf, xS, yS, 'linear');
else
    zS = img.zS;
end



% Init image
z = zeros(size(zS))-Inf;


% If T is stored as function
if iscell(img.zT)
    
    % Make a tip array around twice as big as the image so it can always 
    % cover, apex needs to be in the centre of the matrix
    %xT = 0:img.x(2)-img.x(1):img.x(end)-img.x(1)
    xT = 0:img.x(2)-img.x(1):img.x(end)-img.x(1)+0.5*(img.x(2)-img.x(1));
    xT = unique([-xT xT]);
    %yT = 0:img.y(2)-img.y(1):img.y(end)-img.y(1);
    yT = 0:img.y(2)-img.y(1):img.y(end)-img.y(1)+0.5*(img.y(2)-img.y(1));
    yT = unique([-yT yT])';
    
    % tip_lower_surface
    zT = img.zT{1}(xT, yT, img.zT{2});
    tip_xcentre = length(xS);
    tip_ycentre = length(yS);
    
    % else
    % Need to sort out later if tip stored as coordinates or other functions
    
end

% Reflect the tip (P = -T)
%zP = fliplr(zT);
%zP = flipud(zP);
%zP = -zP;
zP = -rot90(zT, 2);

if do_wb
    hwb = waitbar(0, 'Simulating image, dilate...');
end

for aa = 1:size(zS, 2)
    for bb = 1:size(zS, 1)
        
        % Translate the tip to pixel coordinates
        zP_translate = zP(tip_ycentre-bb+1:end-bb+1, tip_xcentre-aa+1:end-aa+1);
        
        % Dilate
        z = max(z, zS(bb, aa)+zP_translate);
        
    end
    
    % Progressbar
    if do_wb && mod(aa-1, 10) == 0
        waitbar(aa./size(zS, 2), hwb);
    end
end



% Refine and find contact point with finer grids if presision parameter
% is set to positive value, otherwise do not need to refine. 
% This means if precision is 1, everything is very fast since this below is
% skipped.

% If T is stored as function
if iscell(img.zT)
    
    % Tip function T fine grid
    %xTf = 0:finer_grid_spacing:img.x(end)-img.x(1);
    xTf = 0:finer_grid_spacing:img.x(end)-img.x(1)+0.5*(finer_grid_spacing);
    xTf = unique([-xTf xTf]);
    %yTf = 0:finer_grid_spacing:img.y(end)-img.y(1);
    yTf = 0:finer_grid_spacing:img.y(end)-img.y(1)+0.5*(finer_grid_spacing);
    yTf = unique([-yTf yTf])';
    % tip_lower_surface
    zTf = img.zT{1}(xTf, yTf, img.zT{2});
    tip_xcentre = length(xSf);
    tip_ycentre = length(ySf);
end
%zPf = -rot90(zTf, 2);

% Set up reconstituted contact points variables
[xC, yC] = meshgrid(img.x, img.y);
zC = z;

% Solver stats, just for checking and debugging if needed
solver_stats = ones(size(z));

if do_wb
    waitbar(0, hwb, 'Simulating image, dilate and refine...');
end

for aa = 1:size(zS, 2)
    for bb = 1:size(zS, 1)
        
        % Translate the tip (T with finer sampling) to pixel coordinates
        zT_translate = zTf(tip_ycentre-precision*(bb-1):end-precision*(bb-1), tip_xcentre-precision*(aa-1):end-precision*(aa-1));
        
        % Fine adjustments
        diff = Inf;
        while abs(diff) > 1e-4
            diff = (zT_translate+z(bb, aa))-zSf;
            diff = min(diff(:));
            z(bb, aa) = z(bb, aa)-diff;
        end
        
        % Find contact point
        diff = (zT_translate+z(bb, aa))-zSf;
        [y_ind, x_ind] = find(diff == min(diff(:)), 1);
        xC(bb, aa) = xSf(y_ind, x_ind);
        yC(bb, aa) = ySf(y_ind, x_ind);
        zC(bb, aa) = zSf(y_ind, x_ind);
        solver_stats(bb, aa) = min(diff(:));
    end
    
    if do_wb && mod(aa-1, 10) == 0
        waitbar(aa./size(zS, 2), hwb);
    end
end

if do_wb
    close(hwb);
end


img.z = z;

img.xC = xC;
img.yC = yC;
img.zC = zC;


end