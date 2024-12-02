function  [x_top, y_top, z_top] = FilamentModelTop(x, y, z, r)

%
% DESCRIPTION
% – Indicate the topside of a filament model for simulation AFM
% – Part of Trace_y by WFX
%
% USAGE
% Standard method usage with filament number as input
% >> [x_top, y_top, z_top] = FilamentModelTop(x, y, z, r)
%
% INPUTS
% x, y, z  –  Model x, y and z coodinates
% r  –  The estimated max radius of the model, for moving the model to
% surface for simulation
%
% OUTPUTS
% x_top, y_top, z_top  –  The top coordinates marked (others are NaN)
%
% DEPENDENCIES
% – Is used by and updated with MakeHelicalFilamentModel.m
%
% AUTHORS
% Liisa Lutter, Wei-Feng Xue
%
% HISTORY
% 2019.07  –  Original draft by LL Filament_Model_Top
% 2019.08  –  WFX edit of the draft and used to make filament model
% 2021.01  –  LL edit for non-gridded model improvements as WFX assumed 
% fixed structure in x, z, which is not fully correct anymore.
% 2021.03  –  WFX edit for non-gridded x/z reconstructions
%



% Assume height is along the z-axis and model has been interpolated onto a
% y-grid

% Init variables
x_top = x;
y_top = y;
z_top = z;

for aa = 1:size(y, 1)
    
    % Get widest points
    [~, mn_id] = min(x(aa, :));
    [~, mx_id] = max(x(aa, :));
    mn_mx = [mn_id, mx_id];
    
    % Check which part is the top
    mn_mx = sort(mn_mx);
    half_id = mn_mx(1):mn_mx(2);
    other_id = [1:mn_mx(1)-1 mn_mx(2)+1:size(x, 2)];
    sum1 = sum(z(aa, half_id));
    sum2 = sum(z(aa, other_id));
    
    % Remove the other part, actually just mark them with NaN
    if sum1 > sum2
        z_top(aa, other_id) = NaN;
    else
        z_top(aa, half_id) = NaN;
    end
end

% Move up filament so that it's bottom is at z = 0
z_top = z_top + r;
end