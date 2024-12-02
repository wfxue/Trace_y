function [img] = ImportAFMImageData(varargin)

%
% DESCRIPTION
% – Imports AFM data files. Format can be auto-detected or manually
% selected in the open file dialogue 
% – Part of Trace_y by WFX
%
% USAGE
% Standard usage with open file dialogue
% >> img = ImportAFMImageData;
% With additional inputs to be passed on, e.g. all channels for a Bruker 
% Nanoscope file
% >> img = ImportAFMImageData('all');
%
% INPUTS
% varargin  –  Optional variable inputs for the individual import
% functions that will be passed on.
%
% OUTPUTS
% img  –  An instance of Trace_y @AFMimage/ AFM image object. For multiple
% channles, img is a stack in the form of a object array
%
% DEPENDENCIES
% – Uses Trace_y's @AFMimage/ object defs and methods.
% – Uses separate import function for different formats
%
% AUTHORS
% Wei-Feng Xue
%
% HISTORY
% 2024.08  –  Initial version. Move to this single function for AFM image 
% data imports. The format is detected and appropriate import function is 
% then invoked. Make it more straight forward to add import formats in the 
% future.
%



% Open file dialogue
[fileName, pathName] = uigetfile({...
    '*.*' 'Auto detect format'
    '*.spm' 'Bruker Nanoscope (.spm)'; ...
    '*.*' 'Bruker Nanoscope legacy (.[number])'; ...
    '*.csv;*.txt' 'Delimited text'});

if fileName == 0
    fprintf('No data imported! \n');
    return
end

filePath = strcat(pathName, fileName);



% Determining the format based on file extension. Can add more advanced
% header checks if nessasary
[~, ~, fileExt] = fileparts(fileName);
fileExt = lower(fileExt);

% Any special treatment
% For Bruker Nanoscope legacy files, the file extention is a number
if ~isnan(str2double(fileExt))
    fileExt = '.spm';
end



% Action the appropriate imports
switch fileExt
    case '.spm'
        if ~isempty(varargin)
            % Optional channel input
            img = ImportBrukerSPM(filePath, varargin{1});
            if strcmpi(varargin{1}, 'all')
                fprintf('Bruker Nanoscope file (all channels) imported: \n');
            else
                fprintf('Bruker Nanoscope file (%g selected channels) imported: \n', numel(varargin{1}));
            end
        else
            img = ImportBrukerSPM(filePath);
            fprintf('Bruker Nanoscope file (first channel) imported: \n');
        end
        
    case {'.csv', '.txt'}
        if ~isempty(varargin)
            % Optional but importaint and recomended inputs
            % varargin = {imgType zUnit scanSizex scanSizey scanSizeUnit}
            img = ImportCSV(filePath, varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5});
        else
            img = ImportCSV(filePath);
        end
        fprintf('Delimited text file imported: \n');

    otherwise
        fprintf('No data imported! \n');
        return
end
fprintf('%s\n\n', filePath);
%img(1).DispImage;



end