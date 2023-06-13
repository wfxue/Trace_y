function [] = Trace_y(mode)

%
% Trace_y
% Tracing and individual filament 3D structural analysis
% and an app for AFM image analysis
%
% Author: WFX (Wei-Feng Xue)
% Contributors and testers: Liisa Lutter, Ben Blakeman, Liam Aubrey
%
% Version history
% 2014 March, v0.x: Project start, initially based on various old
%   WFX scripts (xAFMTools) for filament tracing on AFM images.
% 2014 August, v1: Initial version with rudimentary visualsations
% 2016 June, v2: Overhaul of tracing method to a image model based approach
%   which is faster and more acurate.
% 2019 August, v3: Refined tracing method for individual filament analysis.
%   Also first implementation of model building and image simulations
% 2020 June, v4: Overhaul of the simulation code as well as allowing
%   generic simulations with EM maps etc. Also visualisation updates
% 2021 March, v5: Updated MakeHelicalFilament method and other requisite
%   functions
% 2022 October, v6: Major updates for making and displaying
%   HelicalFilament and the Trace_y data format. Also added few more GUI 
%   elements for user sharing friendliness. 
% 
% Usage:
% Trace_y('-m'); For matlab command line usage
% Trace_y; Using the simple 'Functions' GUI (compiled or in matlab
%   subseqent to Trace_y('-m');)


% Setup directories and display description
if nargin == 0
    % Save log for running compiled version with UI
    if ispc
        mkdir(pwd, 'Trace_y logs');
        logPath = [pwd '\Trace_y logs\' char(datetime(clock, 'Format', 'yyyyMMdd-HHmmss')) '_log.txt'];
    elseif ismac
        mkdir('~/Documents/', 'Trace_y');
        logPath = ['~/Documents/Trace_y/Trace_y_log_' char(datetime(clock, 'Format', 'yyyyMMdd-HHmmss')) '_log.txt'];
    end
    diary(logPath);
end

dirName = 'Trace_y.src';
appName = 'Trace_y';
appAuthor = 'Wei-Feng Xue';
appVersion = 'v6.2023.0606';
%subDir = char('common');
%subDir = [];
%subDir = cellstr(subDir);
description = ...
    ['Tracing filament polymers and individual filament 3D structural analysis' ...
    newline 'Algorithms for AFM image analysis'];

fprintf('\n%s\n\n%s\n%s\n%s\n', appName, description, appAuthor, appVersion);
fprintf('\n');

% Comment/bypass addpath line before compiling, only need once per session
if nargin > 0
    addpath(sprintf('%s/%s', pwd, dirName));
end

%{
%fprintf('addpath(''%s/%s.src'')', pwd, dirName);
for aa = 1:size(subDir, 1)
    addpath(sprintf('%s/%s.src/%s.src', ...
        pwd, dirName, char(subDir(aa))));
    fprintf('addpath(''%s/%s.src/%s.src'')', ...
        pwd, dirName, char(subDir(aa)))
end
%}



% Initiate AFMimage objects
img = AFMimage(zeros(256), 'nm', 1000, 1000, 'nm');



% Check if doing UI (compiled) mode or runing in Matlab command line
if nargin == 0
    % Using simple UI for compiled app
    fprintf('Log file: %s\n', logPath);
    fprintf('Please start by loading an AFM height image file.\n\n')


elseif strcmp(mode, '-m')
    % Running in Matlab commandline
    clear('img');
    return
    
end



%
% The UI bit for simple function list
% Init main UI window


% Definitions

% Define curent positions from bottom right and width
screenSize = get(0, 'ScreenSize');
posUIwindow = [1 screenSize(4)-575 550 500];
% UI elements pos relative to the main window
posFcnListHeader = [25 posUIwindow(4)-40 200 20];
posFcnList = [25 225 225 posUIwindow(4)-25-25-225];

posInputTableHeader = [275 posUIwindow(4)-40 200 20];
posInputTable = [275 255 250 posUIwindow(4)-25-25-255];

posRunButton = [275 223 250 25];

posTextBox = [25 25 500 posUIwindow(4)-25-25-260];

% Functions list
fcnList = {...
    'Select functions:'; ...
    ''; ...
    '___ Suggested Workflow ___'; ...
    'Import CSV text file (.csv)'; ...
    'Trace filament'; ...
    'Display filament'; ...
    'Reconstruct helical 3D model'; ...
    'Display filament 3D model'; ...
    'Save session Trace_y file (.trcy)'; ...
    'Export filament results (.csv)'; ...
    ''; ...
    '___ Loading and Saving ___'; ...
    'Export filament results (.csv)'; ...
    'Import Bruker file (.000, .001 ...)'; ...
    'Import CSV text file (.csv)'; ...
    'Load Trace_y file (.trcy)'; ...
    'Save session Trace_y file (.trcy)'; ...
    ''; ...
    '___ Visualisation ___'; ...
    'Display filament'; ...
    'Display filament 3D model'; ...
    'Display image'; ...
    'Display scan line'; ...
    ''; ...
    '___ Structural Analysis ___'; ...
    'Reconstruct helical 3D model'; ...
    'Trace filament'; ...
    ''; ...
    '___ Exit ___'; ...
    'Quit'; ...
    ''...
    };

% Default values
inputImportCSV = ...
    {'Height (z) unit' 'nm'; ...
    'Scan Size (xy) unit' 'nm'; ...
    'Scan Size (x)' 2000; ...
    'Scan Size (y)' 2000};
inputDispFilament = ...
    {'Filament number' 1; ...
    'Unit' 'nm'};
inputDispFilament3DModel = ...
    {'Filament number' 1; ...
    'Segment start / nm' 0; ...
    'Segment legth / nm' 500};
inputTraceFilamentM = ...
    {'Apparent width / pixels' 7; ...
    'Height threshold / nm' 1};
inputMakeHelicalFilamentModel = ...
    {'Filament number' 1; ...
    'Handedess' 'left'; ...
    'Symmetry estimate' 2; ...
    'Refine tip radius estimate' 'no'; ...
    'Smoothness' 4; ...
    'Refine' 2};



% Main Window

hUIwindow = uifigure('Name', ['*** ' appName ' ' appVersion ' ***'], ...
    'WindowStyle', 'normal', ...
    'Color', [1 1 1], ...
    'NumberTitle', 'off', ...
    'Position', posUIwindow, ...
    'Resize', 'off', ...
    'Visible', 'on', ...
    'DeleteFcn', @quitTrace_y);

hZoom = zoom(hUIwindow);
set(hZoom, 'Enable', 'off');


% Listbox of functions
uicontrol(hUIwindow, 'Style', 'text', ...
    'String', '   Functions', ...
    'Units', 'pixels', 'Position', posFcnListHeader, ...
    'Callback', [], 'UserData', [], ...
    'FontName', 'helvetica', ...
    'FontUnits', 'pixels', 'FontSize', 16, ...
    'FontWeight', 'normal', ...
    'ForegroundColor', 'k', ...
    'BackgroundColor', 'w', ...
    'HorizontalAlignment', 'left');

%hFcnList = 
uicontrol(hUIwindow, 'Style', 'listbox', ...
    'String', fcnList, ...
    'Units', 'pixels', 'Position', posFcnList, ...
    'Callback', @fcnSelect, 'UserData', [], ...
    'FontName', 'helvetica', ...
    'FontUnits', 'pixels', 'FontSize', 14, ...
    'FontWeight', 'normal');



% Run button
hRunButton = uicontrol(hUIwindow, 'Style', 'pushbutton', ...
    'String', 'Run', 'Enable', 'off', ...
    'Units', 'pixels', 'Position', posRunButton, ...
    'Callback', [], 'UserData', [], ...
    'FontName', 'helvetica', ...
    'FontUnits', 'pixels', 'FontSize', 14, ...
    'FontWeight', 'normal');


% Text box
hEchoText = uitextarea(hUIwindow, ...
    'Value', ' ', ...
    'Position', posTextBox, ...
    'UserData', [], ...
    'FontName', 'helvetica', ...
    'FontSize', 12, ...
    'FontWeight', 'normal', ...
    'HorizontalAlignment', 'left');


% Input parameter table
uicontrol(hUIwindow, 'Style', 'text', ...
    'String', '   Input parameters', ...
    'Units', 'pixels', 'Position', posInputTableHeader, ...
    'Callback', [], 'UserData', [], ...
    'FontName', 'helvetica', ...
    'FontUnits', 'pixels', 'FontSize', 16, ...
    'FontWeight', 'normal', ...
    'ForegroundColor', 'k', ...
    'BackgroundColor', 'w', ...
    'HorizontalAlignment', 'left');

hInputTable = uitable(hUIwindow, 'Data', [], ...
    'Units', 'pixels', 'Position', posInputTable, ...
    'ColumnName', [], ...
    'ColumnEditable', logical([0 1]), ...
    'ColumnFormat', [], ...
    'ColumnWidth', {145 50}, ...
    'CellSelectionCallback', [], ...
    'CellEditCallback', [], ...
    'UserData', [], ...
    'FontName', 'helvetica', ...
    'FontUnits', 'pixel', 'FontSize', 12, ...
    'FontWeight', 'normal');



    function [] = fcnSelect(~, action)
        fcnIndx = action.Source.Value;
        action.Source.Value = 1;

        switch fcnList{fcnIndx}
            case 'Load Trace_y file (.trcy)'
                fprintf('|| %s\n', fcnList{fcnIndx});
                UpdateLog;
                set(hInputTable, 'ColumnName', [], 'Data', []);
                set(hRunButton, 'Enable', 'off', ...
                    'Callback', []);

                [trcyFile, trcyFilePath] = uigetfile('.trcy', 'Please select a .trcy file to open');
                load([trcyFilePath trcyFile], '-mat', 'img');

                fprintf('%s loaded\n\n', [trcyFilePath trcyFile]);

            case 'Save session Trace_y file (.trcy)'
                fprintf('-> %s\n', fcnList{fcnIndx});
                UpdateLog;
                set(hInputTable, 'ColumnName', [], 'Data', []);
                set(hRunButton, 'Enable', 'off', ...
                    'Callback', []);

                [trcyFile, trcyFilePath] = uiputfile(strcat(img.dataFile, '.trcy'));
                save([trcyFilePath trcyFile], '-mat', 'img');

                fprintf('%s saved\n\n', [trcyFilePath trcyFile]);

            case 'Import Bruker file (.000, .001 ...)'
                fprintf('|| %s\n', fcnList{fcnIndx});
                UpdateLog;
                set(hInputTable, 'ColumnName', [], 'Data', []);
                set(hRunButton, 'Enable', 'off', ...
                    'Callback', []);

                [fileName, filePath] = uigetfile('*.*', 'Please select a Bruker file to open');
                img = ImportBruker([filePath fileName]);

                fprintf('%s loaded\n\n', [filePath fileName]);


            case 'Import CSV text file (.csv)'
                fprintf('|| %s\n', fcnList{fcnIndx});
                UpdateLog;
                set(hInputTable, 'Data', inputImportCSV, ...
                    'ColumnName', {'Input parameter' 'Value'});
                set(hRunButton, 'Enable', 'on', ...
                    'Callback', @ImportCSVCallback);



             case 'Export filament results (.csv)'
                fprintf('-> %s\n', fcnList{fcnIndx});
                UpdateLog;
                set(hInputTable, 'ColumnName', [], 'Data', []);
                set(hRunButton, 'Enable', 'off', ...
                    'Callback', []);

                [exportFileName, filePath] = uiputfile(strcat(img.dataFile, '_filaments.csv'));
                img.CompileFilamentResults([filePath exportFileName]);

                fprintf('%s exported\n\n', [filePath exportFileName]);



            case 'Display image'
                fprintf('-> %s\n', fcnList{fcnIndx});
                UpdateLog;

                set(hInputTable, 'ColumnName', [], 'Data', []);
                set(hRunButton, 'Enable', 'off', ...
                    'Callback', []);
                figure(1);
                clf;
                img.DispImage;
                set(gcf, 'OuterPosition', [1 300 700 700]);

                fprintf('\n');

            case 'Display scan line'
                fprintf('-> %s\n', fcnList{fcnIndx});
                UpdateLog;

                set(hInputTable, 'ColumnName', [], 'Data', []);
                set(hRunButton, 'Enable', 'off', ...
                    'Callback', []);
                img.DispLine;

                fprintf('\n');

            case 'Display filament'
                fprintf('-> %s\n', fcnList{fcnIndx});
                UpdateLog;

                set(hInputTable, 'Data', inputDispFilament, ...
                    'ColumnName', {'Input parameter' 'Value'});
                set(hRunButton, 'Enable', 'on', ...
                    'Callback', @DispFilamentCallback);

                fprintf('\n');

            case 'Display filament 3D model'
                fprintf('-> %s\n', fcnList{fcnIndx});
                UpdateLog;

                set(hInputTable, 'Data', inputDispFilament3DModel, ...
                    'ColumnName', {'Input parameter' 'Value'});
                set(hRunButton, 'Enable', 'on', ...
                    'Callback', @DispFilament3DModelCallback);

                fprintf('\n');



            case 'Trace filament'
                fprintf('-> %s\n', fcnList{fcnIndx});
                UpdateLog;

                set(hInputTable, 'Data', inputTraceFilamentM, ...
                    'ColumnName', {'Input parameter' 'Value'});
                set(hRunButton, 'Enable', 'on', ...
                    'Callback', @TraceFilamentMCallback);
                
                fprintf('\n');

            case 'Reconstruct helical 3D model'
                fprintf('-> %s\n', fcnList{fcnIndx});
                UpdateLog;

                set(hInputTable, 'Data', inputMakeHelicalFilamentModel, ...
                    'ColumnName', {'Input parameter' 'Value'});
                set(hRunButton, 'Enable', 'on', ...
                    'Callback', @MakeHelicalFilamentModelCallback);
                
                fprintf('\n');


            case 'Quit'
                fprintf('|| %s\n', fcnList{fcnIndx});
                UpdateLog;

                close(hUIwindow);
                
                fprintf('\n');
                return
        end

        UpdateLog;
    end


    function [] = ImportCSVCallback(~, ~)
        inputImportCSV = hInputTable.Data;
        
        [fileName, pathName] = uigetfile('*.csv;*.txt');
        
        img = ImportCSV(strcat(pathName, fileName), ...
            inputImportCSV{1, 2}, inputImportCSV{2, 2}, ...
            inputImportCSV{3, 2}, inputImportCSV{4, 2});
        fprintf('%s loaded\n\n', strcat(pathName, fileName));


    end

    function [] = DispFilamentCallback(~, ~)
        inputDispFilament = hInputTable.Data;
        img.DispFilament(inputDispFilament{1, 2}, inputDispFilament{2, 2});
    end

    function [] = DispFilament3DModelCallback(~, ~)
        inputDispFilament3DModel = hInputTable.Data;
        img.DispFilament3DModel(inputDispFilament3DModel{1, 2}, [], [], ...
            inputDispFilament3DModel{2, 2}, inputDispFilament3DModel{3, 2});
    end

    function [] = TraceFilamentMCallback(~, ~)
        inputTraceFilamentM = hInputTable.Data;
        img.TraceFilamentM(inputTraceFilamentM{1, 2}, inputTraceFilamentM{2, 2});
    end

    function [] = MakeHelicalFilamentModelCallback(~, ~)
        inputMakeHelicalFilamentModel = hInputTable.Data;
        img.MakeHelicalFilamentModel(inputMakeHelicalFilamentModel{1, 2}, inputMakeHelicalFilamentModel{2, 2}, ...
            inputMakeHelicalFilamentModel{3, 2}, inputMakeHelicalFilamentModel{4, 2}, ...
            inputMakeHelicalFilamentModel{5, 2}, inputMakeHelicalFilamentModel{6, 2});
        UpdateLog;
    end


UpdateLog;
    function [] = UpdateLog()
        logText = readlines(logPath);
        set(hEchoText, 'Value', logText);
        scroll(hEchoText, 'bottom');
    end

    function [] = quitTrace_y(~, ~)
        fprintf('%s\nSession finished\n\n', char(datetime(clock, 'Format', 'yyyyMMdd-HHmmss')))
        diary('off');
        %close(hUIwindow);
    end


end