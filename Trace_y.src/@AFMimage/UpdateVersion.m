function [img_updated] = UpdateVersion(img, filamentIndex)

%
% DESCRIPTION
% – A self-update method, mainly useful incase major changes are added in 
% the future on the data structure and back compatibility is needed
% – Part of Trace_y by WFX
%
% USAGE
% Standard method usage
% >> img = img.UpdateVersion;
%
% INPUTS
% img  –  AFMimage object. Implicit with dot notation
%
% OUTPUTS
% img_updated  –  Updated object
%
% DEPENDENCIES
% – Method for Trace_y's @AFMimage/ object.
%
% AUTHORS
% Wei-Feng Xue
%
% HISTORY
% See individual update comments



% Current version of @AFMimage is 7
currentVersion = 7;
img_updated = img;


% 2024.08
% Simplify the features property to a struct. Referencing should be the
% same with the dot notation but more efficient and got rid of legacy stuff
% no longer needed
if img.version <= 6
    img_updated.features = struct;
    img_updated.features.filaments = img.features.filaments;
end


% 2024.08
% Minor change to move format specific meta info to one structure. Mostly
% organisational and have no practical implications. Upon oppening, some
% previous scaling properties will have been lost though due to object dep
% update.
if img.version <= 6
    [~, ~, fileExt] = fileparts(img.dataFile);
    img_updated.dataFileProp = struct;
    img_updated.dataFileProp.format = fileExt;
end


% 2023.12
% Major update to how filament 3D model is generated and saved. This update
% requires the 3D model to be re-generated so some additional user action
% is needed
if img.version <= 6 && ~isempty(img.features.filaments) ...
        && (isprop(img.features.filaments(1), 'hf3DModel') || ...
        ~isprop(img.features.filaments(1, 1).helicalFilament3DModel, 'xz_bandwidth'))
    fprintf('Pre-v7 .trcy format, helical filament 3D models need to be re-reconstructed...\n')
    %fprintf('Use UpdateVersion with filament index input to update specific filament 3D model...\n')

    %{
    for filamentIndex = 1:numel(img.features.filaments)

        if isprop(img.features.filaments(filamentIndex), 'hf3DModel') && ~isempty(img.features.filaments(filamentIndex).hf3DModel)
            img_updated = img_updated.MakeHelicalFilamentModel(filamentIndex, ...
                img.features.filaments(filamentIndex).hf3DModel.handedness, ...
                img.features.filaments(filamentIndex).hf3DModel.symmetry, ...
                [], false);
        elseif ~isempty(img.features.filaments(filamentIndex).helicalFilament3DModel)
            img_updated = img_updated.MakeHelicalFilamentModel(filamentIndex, ...
                img.features.filaments(filamentIndex).helicalFilament3DModel.handedness, ...
                img.features.filaments(filamentIndex).helicalFilament3DModel.symmetry, ...
                [], false, ...
                img.features.filaments(filamentIndex).helicalFilament3DModel.smoothness, ...
                img.features.filaments(filamentIndex).helicalFilament3DModel.refine);
        end
    end
    %}
    
end

if exist('filamentIndex', 'var') && ~isempty(filamentIndex)
    filamentIndex = floor(filamentIndex(1));
    fprintf('Updating filament 3D model for filament %g...\n', filamentIndex);
    if isprop(img.features.filaments(filamentIndex), 'hf3DModel') && ~isempty(img.features.filaments(filamentIndex).hf3DModel)
        img_updated = img_updated.MakeHelicalFilamentModel(filamentIndex, ...
            img.features.filaments(filamentIndex).hf3DModel.handedness, ...
            img.features.filaments(filamentIndex).hf3DModel.symmetry, ...
            [], false);
    elseif ~isempty(img.features.filaments(filamentIndex).helicalFilament3DModel)
        img_updated = img_updated.MakeHelicalFilamentModel(filamentIndex, ...
            img.features.filaments(filamentIndex).helicalFilament3DModel.handedness, ...
            img.features.filaments(filamentIndex).helicalFilament3DModel.symmetry, ...
            [], false, ...
            img.features.filaments(filamentIndex).helicalFilament3DModel.smoothness, ...
            img.features.filaments(filamentIndex).helicalFilament3DModel.refine);
    end
    %img.features.filaments(filamentIndex) = rmfield(img.features.filaments(filamentIndex), 'hf3DModel');
end



% Update version number
img_updated.version = currentVersion;



end