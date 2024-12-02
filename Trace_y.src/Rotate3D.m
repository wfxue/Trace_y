function xyz = Rotate3D(xyz, r_angles)

%
% DESCRIPTION
% – Standard 3D rotation using rotation matrixes
% – Part of Trace_y
%
% USAGE
% Standard function usage
% >> xyz_tr = Rotate3D(xyz0, [rx ry rz]);
%
% INPUTS
% xyz  –  coordinates in tri format
% r_angles  –  [rx ry rz] in degrees
%
% OUTPUTS
% xyz  –  Transformed coordinates
%
% DEPENDENCIES
% – Used for shape manipulations in Trace_y
%
% AUTHORS
% Wei-Feng Xue
% 
% HISTORY
% 2020.06  –  Added for convenience as used for 3D shape transformations
%


% Rotation matrixes, angles in degrees
rx = [1 0 0; 0 cosd(r_angles(1)) -sind(r_angles(1)); 0 sind(r_angles(1)) cosd(r_angles(1))];
ry = [cosd(r_angles(2)) 0 sind(r_angles(2)); 0 1 0; -sind(r_angles(2)) 0 cosd(r_angles(2))];
rz = [cosd(r_angles(3)) -sind(r_angles(3)) 0; sind(r_angles(3)) cosd(r_angles(3)) 0; 0 0 1];

xyz = (rz*ry*rx*xyz')';

end

