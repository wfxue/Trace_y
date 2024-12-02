function [img] = ImportBrukerSPM(filePath, getChannel)

%
% DESCRIPTION
% – Imports Bruker / Nanoscope data files (.spm)
% – Also opens legacy nanoscope files (.[number])
% – Part of Trace_y by WFX
%
% INPUTS
% filePath  –  Optional full file path.
% getChannel  –  The image channel to be imported, default is 1. 'all' for
%   importing all channels
%
% OUTPUTS
% img  –  An instance of Trace_y @AFMimage/ AFM image object. For multiple
% channles, img is a stack in the form of a object array
%
% DEPENDENCIES
% – Uses Trace_y's @AFMimage/ object defs and methods.
% – Used by ImportAFMImageData.m
%
% AUTHORS
% Wei-Feng Xue
%
% HISTORY
% 2015.08  –  Imitial draft of ImportBruker.m
% 2020.06  –  Updated to match v4 of AFMimage object
% 2022.10  –  Updated to match v6 of AFMimage constructor changes
% 2024.08  –  Full rewite to tidy up all code, to store most header 
% information and to make it compatible with newer .spm files (as WFX got
% sent newer Nanoscope files that the old ImportBruker.m cannot cope with).
% The new ImportBrukerSPM.m also include updats to match v7 of AFMimage 
% constructor mods.
%



%
% Defaults
%
if ~exist('getChannel', 'var') || isempty(getChannel)
    % The get all channels
    getChannel = 1;
end

% Open file dialogue if needed
if ~exist('filePath', 'var') || isempty(filePath)
    [fileName, pathName] = uigetfile('*.*');
    dataFile = fileName;
    filePath = strcat(pathName, fileName);
else
    if ispc
        dataFile = ...
            filePath(find(filePath == '\', 1, 'last')+1:end);
    else
        dataFile = ...
            filePath(find(filePath == '/', 1, 'last')+1:end);
    end
end

% This may need to be expanded with proper version information in the
% future if needed
[~, ~, fileExt] = fileparts(filePath);
if strcmpi(fileExt, '.spm')
    legacyFormat = false;
else
    legacyFormat = true;
end

fid1 = fopen(filePath, 'r');



%
% Get the all the header metadata
%

% .spm files' headers start with "\*File list" and ends with 
% "\*File list end". Find how many lines of headers we have
nl = 0;

% Count number of props and mage channels
% Number of image channels, starts with '\*Ciao image list'
nchannel = 0;
headerFinished = false;

while ~headerFinished
    nextLine = fgetl(fid1);
    nl = nl+1;
    if strcmp(nextLine, '\*Ciao image list')
        nchannel = nchannel+1;
    end
    if strcmp(nextLine, '\*File list end')
        headerFinished = true;
    end
end

% Get the header metadata for the mage channels
imageHeader = cell(nchannel, 1);
fileHeader = struct;
equipmentHeader = struct;
scannerHeader = struct;
scanHeader = struct;
fastscanHeader = struct;

frewind(fid1);
% Set a state flag, -1 is not relevant, 0 is image header line, 1 is image
% channel etc for pairs of headers and header properties
recordProp = -1;
channel = 0;
for aa = 1:nl

    nextLine = fgetl(fid1);
    if strcmp(nextLine, '\*Ciao image list')
        channel = channel+1;
        imageHeader{channel} = struct;
        recordProp = 0;
    elseif strcmp(nextLine, '\*Equipment list')
        recordProp = 2;
    elseif strcmp(nextLine, '\*Scanner list')
        recordProp = 4;
    elseif strcmp(nextLine, '\*Ciao scan list')
        recordProp = 6;
    elseif strcmp(nextLine, '\*File list')
        recordProp = 8;
    elseif strcmp(nextLine, '\*Fast Scan list')
        recordProp = 10;
    elseif regexp(nextLine, '\\\*')
        % Other stuff or the end of header line
        recordProp = -1;
    end

    prop = strsplit(nextLine, ': ');
    propName = regexprep(prop{1}, '@', 'a');
    propName = regexprep(propName, '\', '');
    propName = regexprep(propName, '[^a-zA-Z0-9]', '_');
    propName = regexprep(propName, '\s', '_');
    if isscalar(prop)
        propVal = [];
    elseif numel(prop) == 3
        propVal = prop{2}+prop{3};
    else
        propVal = prop{2};
    end

    switch recordProp
        case 0
            imageHeader{channel}.channel = channel;
            recordProp = 1;
        case 1
            imageHeader{channel}.(propName) = propVal;
        case 2
            recordProp = 3;
        case 3
            equipmentHeader.(propName) = propVal;
        case 4
            recordProp = 5;
        case 5
            scannerHeader.(propName) = propVal;
        case 6
            recordProp = 7;
        case 7
            scanHeader.(propName) = propVal;
        case 8
            recordProp = 9;
        case 9
            fileHeader.(propName) = propVal;
        case 10
            recordProp = 11;
        case 11
            fastscanHeader.(propName) = propVal;
            
    end
end



%
% Get the image data in LSB integer format
%

% For setting image size 
pixelPerLine = zeros(nchannel, 1);
nLines = zeros(nchannel, 1);
z = cell(nchannel, 1);
imgType = cell(nchannel, 1);
zUnit = cell(nchannel, 1);

%zSens = textscan(scanHeader.aSens__ZsensSens, '%s %f %s');
%zSens = zSens{2};

if (ischar(getChannel) || isstring(getChannel)) && strcmp(getChannel, 'all')
    % Then get all channels
    getChannel = 1:nchannel;
end

for aa = getChannel

    dataOffset = str2double(imageHeader{aa}.Data_offset);
    pixelPerLine(aa) = str2double(imageHeader{aa}.Samps_line);
    nLines(aa) = str2double(imageHeader{aa}.Number_of_lines);
    bytePerPixel = str2double(imageHeader{aa}.Bytes_pixel);
    bitsPerPixel = bytePerPixel*8;

    % Read data in LSB
    fseek(fid1, dataOffset, 'bof');
    if bytePerPixel == 2
        zz = fread(fid1, pixelPerLine(aa).*nLines(aa), 'int16');
    elseif bytePerPixel == 4
        zz = fread(fid1, pixelPerLine(aa).*nLines(aa), 'int32');
    end
    z{aa} = zeros(pixelPerLine(aa), nLines(aa));

    % Get scaling and sensitivity values
    zScale = regexprep(imageHeader{aa}.a2_Z_scale, '[[]()]', '"');
    zScale = textscan(zScale, '%s %q %q %f %s');
    if sum(cellfun(@isempty, zScale)) == 0
        zSensProp = zScale{2};
        zSensProp = regexprep(zSensProp, '[^a-zA-Z0-9]', '_');
        zSensProp = "a"+zSensProp{1};
        % Legacy format issue: check field names 
        if isfield(scanHeader, zSensProp)
            zSens = textscan(scanHeader.(zSensProp), '%s %f %s');
        else
            % This is in legacy format
            zSens = textscan(scannerHeader.(zSensProp), '%s %f %s');
        end

        if sum(cellfun(@isempty, zSens)) == 0
            zUnit{aa} = zSens{3}{1};
            zUnit{aa} = strsplit(zUnit{aa}, '/');
            zUnit{aa} = zUnit{aa}{1};
        else
            zUnit{aa} = zScale{5}{1};
        end
        zSens = zSens{2};
        zScale = zScale{4}/(2.^bitsPerPixel);
    else
        % For some reason some channels have less information
        zSens = 1;
        zUnit{aa} = '';
        zScale = str2double(zScale{3})/(2.^bitsPerPixel); 
    end

    % Calculate scaled z values, e.g. nm for height and V or nV for other
    % uncalibrated measures
    %z{aa}(:) = zz;
    z{aa}(:) = zz.*zScale.*zSens;
    z{aa} = z{aa}';

    imgType{aa} = regexprep(imageHeader{aa}.a2_Image_Data, '[[]()]', '"');
    imgType{aa} = textscan(imgType{aa}, '%s %q %q');
    imgType{aa} = imgType{aa}{3}{1};
end

fclose(fid1);



%
% Construct image object array, one for each channel
%

% Init the img objects
img(nchannel, 1) = AFMimage;

% Get needed scan properties
scanSizex = scanHeader.Scan_Size;
scanSizex = textscan(scanSizex, '%f %s');
% This is assuming same unit / aspect ratio of 1
scanSizeUnit = scanSizex{2}{1};
scanSizex = scanSizex{1};

% Legacy format issue: check field names
if isfield(scanHeader, 'Slow_Axis_Size')
    scanSizey = scanHeader.Slow_Axis_Size;
    scanSizey = textscan(scanSizey, '%f %s');
    scanSizey = scanSizey{1};
else
    scanSizey = scanSizex;
end

for aa = getChannel

    % Create new AFMimage object
    img(aa) = AFMimage(z{aa}, zUnit{aa}, scanSizex, scanSizey, scanSizeUnit);

    % Set some needed image properties
    if legacyFormat
        img(aa).dataFileFormat = 'Bruker Nanoscope legacy (.[number])';
    else
        img(aa).dataFileFormat = 'Bruker Nanoscope (.spm)';
    end
    img(aa).dataFile = dataFile;

    img(aa).imgType = imgType{aa};
    % Special
    if strcmpi(img(aa).imgType, 'Height Sensor')
        img(aa).imgType = 'Height';
    end

    img(aa).dataFileProp.format = fileExt;
    img(aa).dataFileProp.fileHeader = fileHeader;
    img(aa).dataFileProp.imageHeader = imageHeader{nchannel};
    img(aa).dataFileProp.equipmentHeader = equipmentHeader;
    img(aa).dataFileProp.scannerHeader = scannerHeader;
    img(aa).dataFileProp.scanHeader = scanHeader;
    img(aa).dataFileProp.fastscanHeader = fastscanHeader;

end
img = img(getChannel);



end