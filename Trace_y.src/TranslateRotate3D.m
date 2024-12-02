function xyz_tr = TranslateRotate3D(xyz0, p)

%
% DESCRIPTION
% – Standard translation and rotation transformations used by algorithms 
% for shape manipulations such as CalcHelicalAxisEM.m
% – Part of Trace_y
%
% USAGE
% Standard function usage
% >> xyz_tr = TranslateRotate3D(xyz0, [tx ty tz rx ry rz]);
%
% INPUTS
% xyz0  –  Pre-transformation coordinates in tri format
% p  –  Translation and rotation parameters [tx ty tz rx ry rz]. Angles in 
% degrees
%
% OUTPUTS
% xyz_tr  –  Transformed coordinates
%
% DEPENDENCIES
% – Used for shape manipulations in Trace_y
%
% AUTHORS
% Wei-Feng Xue
% 
% HISTORY
% 2021.02  –  Added for convenience as used for 3D shape transformations
%



% p = [tx ty tz rx ry rz];
xyz0(:, 1) = xyz0(:, 1)+p(1);
xyz0(:, 2) = xyz0(:, 2)+p(2);
xyz0(:, 3) = xyz0(:, 3)+p(3);
xyz_tr = Rotate3D(xyz0, [p(4) p(5) p(6)]);

end
