function [img] = DeleteFilament(img, filamentIndex)

%
% DESCRIPTION
% – Simple function to delete the indicated filament
% – Part of Trace_y
%
% USAGE
% Standard method usage
% >> img = img.DeleteFilament(filamentIndex);
%
% INPUTS
% img  –  AFMimage object. Implicit with dot notation
% filamentIndex  –  The index number of filament to be deleted
%
% OUTPUTS
% img  –  Updated object
%
% DEPENDENCIES
% – Method for Trace_y's @AFMimage/ object.
%
% AUTHORS
% Wei-Feng Xue
%
% HISTORY
% 2024.09  –  Utility for Trace_y update
%



if numel(img.features.filaments) > 1
    img.features.filaments = img.features.filaments([1:filamentIndex-1 filamentIndex+1:end]);
else
    img.features.filaments = Filament.empty;
end


end