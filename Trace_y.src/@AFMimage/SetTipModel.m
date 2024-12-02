function [img] = SetTipModel(img, tipType, tipData)

%
% DESCRIPTION
% – Method to handle the tip function, can be expanded in the future with 
% other tip models.
% – Part of Trace_y
%
% USAGE
% Standard method usage
% >> img = img.SetTipModel(tipType, tipData);
%
% INPUTS
% tipType  –  E.g. 'rounded_cone'
% tipData  –  The parameters for the tip function. E.g. [tip_r tip_a] for 
% 'rounded_cone' where tip_r is the tip radius, scalar in nm, and tip_a is 
% the tip side angle to the tip centre axis, scalar in degrees
%
% OUTPUTS
% img  –  AFM image object.
%
% DEPENDENCIES
% – Used for tip manipulations and simulations in Trace_y
%
% AUTHORS
% Wei-Feng Xue
% 
% HISTORY
% 2020.05  –  Separated out from other code files and added for 
% convenience as used simulations
%



switch tipType
    case 'rounded_cone'
        img.zT = {@TipModel_RoundedCone tipData};
    case 'coordinates'
        % Future update
end


end
