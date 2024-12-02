function [mapfv, map_x, map_y, map_z, map_v] = ImportMRC(contourLevel, filePath, visualise)

%
% DESCRIPTION
% – Imports data from MRC/CCP4 2014 file format. Also renders the 
% iso-surface with supplied contour level
% – Based on information from https://www.ccpem.ac.uk/mrc_format/mrc2014.php
% and ReadMRC.m. 
% – Part of Trace_y
%
% USAGE
% Choose map file and just outputting the iso surface and meta data:
% mapfv = ImportMRC(contourLevel);
% Include the raw map x/y/z/v data:
% [mapfv, map_x, map_y, map_z, map_v] = ImportMRC(contourLevel);
% Supply a map file without open file dialogue and visualise iso-surface
% [mapfv, map_x, map_y, map_z, map_v] = ImportMRC(contourLevel, filePath, true);
%
% INPUTS
% contourLevel: the isovalue for isosurface, default on EMDB entry
% filePath: Optional path to the map file with .mrc or .map filetypes
% visualise: Optional flag for whether do the iso-surface patch plot
%
% OUTPUTS
% mapfv  –  The data is-surface in faces and vertices format to be used as patch object.
% The main header and other meta info is in the data structure under UserData
% map_x, map_y, map_z, map_v  –  Volumetric data with voxels at x, y, z and
% value v
%
% DEPENDENCIES
% –
%
% AUTHORS
% Liisa Lutter and Wei-Feng Xue
%
% HISTORY
% 2019.11  –  Initial draft by LL as part of EM_import.m
% 2020.05  –  WFX edit to tidy up and to integrate with other code
%



mh = struct;

% open file
if nargin < 2 || isempty(filePath)
    [fileName, pathName, nfile] = uigetfile({'*.mrc;*.map' '.mrc or .map'});
    if nfile == 0
        return
    end
    filePath = strcat(pathName, fileName);
else
    %textscan(filePath, '%s', 'Delimiter', '/');
    [~, fileName, ext] = fileparts(filePath);
    fileName = strcat(fileName, ext);
end

% Defaults
if ~exist('visualise', 'var') || isempty(visualise)
    visualise = false;
end



% Information about the data file
mh.mapFileName = fileName;
mh.contourLevel = contourLevel;



% Try for little-endian data first
fh = fopen(filePath, 'r', 'ieee-le');

if fh == -1
    error(['The file could not be opened: ' filePath]);
end

% Check the file size
% Go to the end of the file
fseek(fh, 0, 'eof');
nBytes=ftell(fh);
% Back to the start of the file
fseek(fh, 0, 'bof');
if nBytes < 1024
    % Not big enaugh for the header
    warning('File does not appear to be in mrc format');
    return
end

% Get the first 10 values, which are integers:
% NX is the first variable
mh.NX = fread(fh, 1, 'int32');

% Check the nx value.
if abs(mh.NX) > 1e5
    % Possible the wrong endian format.
    fclose(fh);
    fh = fopen(filePath, 'r', 'ieee-be');
    mh.NX = fread(fh, 1, 'int32');
end

mh.NY = fread(fh, 1, 'int32');
mh.NZ = fread(fh, 1, 'int32');
mh.MODE = fread(fh, 1, 'int32');
mh.NXSTART = fread(fh, 1, 'int32');
mh.NYSTART = fread(fh, 1, 'int32');
mh.NZSTART = fread(fh, 1, 'int32');
mh.MX = fread(fh, 1, 'int32');
mh.MY = fread(fh, 1, 'int32');
mh.MZ = fread(fh, 1, 'int32');

% Words 11-16
mh.CELLA = fread(fh, 3, 'float32');
mh.CELLB = fread(fh, 3, 'float32');

% Words 17-22
mh.MAPC = fread(fh, 1, 'int32');
mh.MAPR = fread(fh, 1, 'int32');
mh.MAPS = fread(fh, 1, 'int32');

mh.DMIN = fread(fh, 1, 'float32');
mh.DMAX = fread(fh, 1, 'float32');
mh.DMEAN = fread(fh, 1, 'float32');

% Words 23-24
mh.ISPG = fread(fh, 1, 'int32');
mh.NSYMBT = fread(fh, 1, 'int32');

% Extra info in words 25-49
mh.EXTRA = fread(fh, 25, 'int32');
mh.EXTTYP = mh.EXTRA(3);
mh.NVERSION = mh.EXTRA(4);

% Words 50-52
mh.ORIGIN = fread(fh, 3, 'int32');

% Words 53 and 54 as single byte character strings
mh.MAP = char(zeros(1, 4));
mh.MAP(:) = fread(fh, 4, 'uchar');
mh.MACHST = fread(fh, 4, 'uchar');

% Words 55 and 56
mh.RMS = fread(fh, 1, 'float32');
mh.NLABL = fread(fh, 1, 'int32');

% up to 10 strings....
%ns = min(mh.NLABL, 10);
mh.LABEL = char(zeros(10, 80));
for aa = 1:10
    mh.LABEL(aa, :) = fread(fh, 80, 'uchar');
end

% Check data mode
switch mh.MODE
    case 0
        mode_string='*int8';
        
    case 1
        mode_string='*int16';
        
    case 2
        mode_string='*float32';
        
    case 6
        mode_string='*uint16';
        
    otherwise
        error(['Unknown data mode: ' num2str(mh.MODE)]);
end

% Go to the end of the main header.
fseek(fh, 1024, 'bof');

% File size sanity check
if nBytes < 1024+mh.NSYMBT
    % Not enough space for main and extended header
    warning('File is too short for the mrc extended header');
    return
end

% Skip extra header
if mh.NSYMBT > 0
    % extraHeader = fread(fh, mh.NSYMBT, 'uchar');
    fseek(fh, 1024+mh.NSYMBT, 'bof');
end

data_length = mh.NX*mh.NY*mh.NZ;
map = fread(fh, data_length, mode_string);
% close file
fclose(fh);

% Data size sanity check
if length(map(:)) ~= data_length
    warning('Not enough data points in the file: %d voxels, %d expected.', ...
        length(map(:)), data_length);
end

% Zero-fill the data array to expected size
map = [map zeros(1, data_length-length(map(:)))];

% Finally reshape data
map_v = reshape(map, mh.NY, mh.NX, mh.NZ);

% Matrix axis permuation needed as Matlab's row, column, page indexing is
% equivalent to y, x, z indexing
map_v = permute(map_v, [2 1 3]);
%map_v = rot90(permute(rot90(map_v, 3), [3 2 1]), 2);


% Get voxel size for x (assume aspect ratio of 1), in Å.
voxSize = mh.CELLA(1)/mh.MX;
% Convert to nm
mh.xyzUnit = 'nm';
voxSize = voxSize/10;
mh.voxSize = voxSize;

map_x = voxSize.*(0:mh.NX-1);
map_y = voxSize.*(0:mh.NY-1);
map_z = voxSize.*(0:mh.NZ-1);
[map_x, map_y, map_z] = meshgrid(map_x, map_y, map_z);



% Generate iso-surface
mapfv = isosurface(map_x, map_y, map_z, map_v, contourLevel);
% Estimate volume by voxel counting
mh.isoVolume = sum(map_v(:) > contourLevel).*voxSize.^3;
% Store data
mapfv.UserData = mh;



% Visualise
if visualise
    patch(mapfv, 'FaceColor', [0 0.7 0.7], 'FaceLighting', 'gouraud', ...
        'EdgeColor', 'none');
    view(3);
    camlight;
    axis('equal');
    xlabel('x / nm');
    ylabel('y / nm');
    zlabel('z / nm')
    set(gca, 'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on');
end

end