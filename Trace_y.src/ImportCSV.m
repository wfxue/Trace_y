function [img, filePath] = ImportCSV(varargin)

%
% DESCRIPTION
% – Imports AFM image data in csv or other delimited text files
% – Part of Trace_y by WFX
%
% INPUTS
% varargin  –  {filePath imgType zUnit scanSizex scanSizey scanSizeUnit} or
%   {imgType zUnit scanSizex scanSizey scanSizeUnit}
% filePath  –  Optional full file path.
% imgType  –  z-data type, e.g. 'Height'
% zUnit  –  z-data unit, e.g. 'nm'
% scanSizex  –  x scan size
% scanSizey  –  y scan size
% scanSizeUnit  –  x/y scan size unit, e.g. 'nm'
%
% OUTPUTS
% img  –  An instance of Trace_y @AFMimage/ AFM image object. 
% filePath  –  The full path to the file selected
%
% DEPENDENCIES
% – Uses Trace_y's @AFMimage/ object defs and methods.
% – Used by ImportAFMImageData.m
%
% AUTHORS
% Wei-Feng Xue
%
% HISTORY
% 2022.10  –  Initial draft
% 2024.08  –  Slight update to make compatible with ImportAFMImageData.m
%



% Open file dialogue if needed
if numel(varargin) == 5
    [fileName, pathName] = uigetfile('*.csv;*.txt');
    dataFile = fileName;
    filePath = strcat(pathName, fileName);

    imgType = varargin{1};
    zUnit = varargin{2};
    scanSizex = varargin{3};
    scanSizey = varargin{4};
    scanSizeUnit = varargin{5};

    % Import data
    z = readmatrix(filePath);

elseif numel(varargin) == 6
    filePath = varargin{1};

    imgType = varargin{2};
    zUnit = varargin{3};
    scanSizex = varargin{4};
    scanSizey = varargin{5};
    scanSizeUnit = varargin{6};
    if ispc
        dataFile = ...
            filePath(find(filePath == '\', 1, 'last')+1:end);
    else
        dataFile = ...
            filePath(find(filePath == '/', 1, 'last')+1:end);
    end

    % Import data
    z = readmatrix(filePath);

elseif isscalar(varargin)
    % No other parameters input case
    filePath = varargin{1};
    if ispc
        dataFile = ...
            filePath(find(filePath == '\', 1, 'last')+1:end);
    else
        dataFile = ...
            filePath(find(filePath == '/', 1, 'last')+1:end);
    end

    % Import data
    z = readmatrix(filePath);

    imgType = 'Text data';
    zUnit = 'AU';
    scanSizex = size(z, 2)-1;
    scanSizey = size(z, 1)-1;
    scanSizeUnit = 'px';

else
    % Nothing is passed in
    [fileName, pathName] = uigetfile('*.csv;*.txt');
    dataFile = fileName;
    filePath = strcat(pathName, fileName);

    % Import data
    z = readmatrix(filePath);

    imgType = 'Text data';
    zUnit = 'AU';
    scanSizex = size(z, 2)-1;
    scanSizey = size(z, 1)-1;
    scanSizeUnit = 'px';

end


%
% Create new AFMimage object
%
img = AFMimage(z, zUnit, scanSizex, scanSizey, scanSizeUnit);
img.dataFileFormat = 'Delimited text (.csv or .txt)';

[~, ~, fileExt] = fileparts(filePath);
img.dataFileProp.format = fileExt;
img.dataFile = dataFile;
img.imgType = imgType;



end