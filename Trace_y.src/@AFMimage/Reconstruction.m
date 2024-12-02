function [img] = Reconstruction(img, tipModel, tipData, precision)

%
% DESCRIPTION
% – Contact point reconstruction algorithm for correcting tip-sample 
% convolution effect in AFM topology images. This is the method for
% @AFMimage class object and uses the new CPR algorithm.
% – The algorithm used is the new direct estimate and replaces the 
% previous inverse problem and slow itterative approach, which also rely on
% helical or cylindrical symetry.
% – Part of Trace_y
%
% USAGE
% Standard method usage
% >> img = img.Reconstruction(tipModel, tipData, precision)
%
% INPUTS
% img  –  AFMimage object. Implicit use with dot notation
% tipModel  –  Function handle to the tip function to be used, e.g. 
% @TipModel_RoundedCone;
% tipData  –  Parameters for the tipModel function, e.g. 
% [tip_radius tip_angle];
% precision  –  Sup-pixel precision for the calculations. Optional and set
% to default of 10. Actual contact points are interpolated between these
% sub-pixel points so the actual precision is better.
%
% OUTPUTS
% img  –  AFMimage object updated with the reconstructed point cloud
% coordinates
%
% DEPENDENCIES
% – Used by other methods for Trace_y's @AFMimage/ object
% – Uses CalcContactPoints.m
%
% AUTHORS
% Wei-Feng Xue
%
% HISTORY
% 2019.02  –  Initial code draft (AFM_Image_Correction) wich is by old inverse
% itterative method similar to old Calc_Filament_Corr.m and
% CalcTipFilamentConvolution.m and is very slow and rely on the helical
% symetry.
% 2020.05  –  Initial test for the new fast direct estimation algorithm and 
% incorporation into @AFMimage object method. This and CalcContactPoints.m
% are updated together.
% 2021.03  –  Further edits to make algorithm avaialable to any xyz image
% in an @AFMimage class object.
%



% Defaults
if ~exist('precision', 'var') || isempty(precision)
    precision = 10;
end

% Setting up the tip model if not already in the img object
% The default tip geometry if no inputs
if isempty(img.zT) && (~exist('tipModel', 'var') || isempty(tipModel))
    % Function for cone with spherical dome
    tipModel = 'rounded_cone';
    % Parameters are tip_r and tip_a (deg)
    tipData = [2, 18];
    img = img.SetTipModel(tipModel, tipData);
elseif exist('tipModel', 'var') && ~isempty(tipModel)
    % Use supplied tip model as priority over existing in img
    img = img.SetTipModel(tipModel, tipData);
end



% Get image coordinates
x = img.x;
y = img.y;
z = img.z;

% Contact point reconstruction
[xr, yr, zr, uncertainty] = CalcContactPoints(x, y, z, img.zT{1}, img.zT{2}, precision);

img.xR = xr;
img.yR = yr;
img.zR = zr;
img.zRUncertainty = uncertainty;

%figure(1);
%img.DispImage;



end