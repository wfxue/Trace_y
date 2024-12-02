function [v] = GetenvTrcy(prop)

%
% DESCRIPTION
% – Standarised function to set up and get various Trace_y properties and 
% paths
% – Part of Trace_y
%
% USAGE
% Standard method usage
% >> value = GetenvTrcy(property)
%
% INPUTS
% property  –  Property name
%
% OUTPUTS
% value  –  Property value
%
% DEPENDENCIES
% – Standard app utility used by various Trace_y code
%
% AUTHORS
% Wei-Feng Xue
%
% HISTORY
% 2024.10  –  Utility for Trace_y
%



% Set up dictionary
env = dictionary;

% Check OS and set up main documents folder
if ispc
    if ~isfolder([getenv('USERPROFILE') '\Documents\Trace_y'])
        mkdir([getenv('USERPROFILE') '\Documents\'], 'Trace_y');
    end

    mainPath = [getenv('USERPROFILE') '\Documents\Trace_y\'];

elseif ismac
    if ~isfolder('~/Documents/Trace_y')
        mkdir('~/Documents/', 'Trace_y');
    end

    mainPath = '~/Documents/Trace_y/';

end

% Main document folder path
env('mainPath') = mainPath;

% Logs
if ~isfolder([mainPath 'logs' filesep])
    mkdir(mainPath, 'logs');
end
env('logPath') = [mainPath 'logs' filesep];

% EMDB downloads
if ~isfolder([mainPath 'EMDB' filesep])
    mkdir(mainPath, 'EMDB');
end
env('emdbDataPath') = [mainPath 'EMDB' filesep];

% Get the property
v = env(prop);

end