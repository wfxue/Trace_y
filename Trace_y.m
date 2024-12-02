function [img_out] = Trace_y(mode)

%
% DESCRIPTION
% – Trace_y: AFM image analysis with emphasis on filamentous features
% – Tracing and individual filament level 3D structural analysis
% – Integrative structural analysis with EMDB structural map data
% – Useful general algorithms for AFM image analysis
%
% USAGE
% Standard app with some "Functions" user interface
% >> Trace_y;
% Running commandline in Matlab
% >> Trace_y('-m');
%
% INPUTS
% mode  –  '-m' for running with commandlines in Matlab
%
% OUTPUTS
% img  –  @AFMImage object
% – AFMimage object can be saved in app as .trcy data file (.mat format).
% – Various visualisations and tabel (.csv) outputs of data avaiable.
%
% DEPENDENCIES
% – Uses several Matlab toolboxes: Statistical and Machine Learning, 
% Optimization, Image Processing, Signal Processing, Computer Vision 
% – Main data format is @AFMimage/ object includng various methods in the
% object def.
%
% AUTHORS
% Wei-Feng Xue
% Contributors and testers: Liisa Lutter, Ben Blakeman, Liam Aubrey,
% Caludia Chitty
%
% HISTORY
% 2014.03  –  v0.x. Project start, initially based on various old WFX 
%   scripts (xAFMTools) for filament tracing on AFM images.
% 2014.08  –  v1. Initial version with rudimentary visualsations
% 2016.06  –  v2. Overhaul of tracing method to a image model based approach
%   which is faster and more acurate.
% 2019  –  v3. Refined tracing method for individual filament analysis.
%   Also first implementation of model building and image simulations
% 2020  –  v4. Overhaul of the simulation code as well as allowing
%   generic simulations with EM maps etc. Also visualisation updates
% 2021.03  –  v5. Updated MakeHelicalFilament method and other requisite
%   functions.
% 2022.10  –  v6. Major updates for making and displaying HelicalFilament 
%   and the Trace_y data format. Also added few more GUI elements for 
%   user friendliness and for eventual open-source page on GitHub
% 2023.06  –  v6.2023.0606. Minor including UI updates, started Github 
%   page and uploaded this version onto GitHub
% 2024.09  –  v7. Major updates to various visualisation and UI elements.
%   Changed the main functions window from a uifigure to a normal figure as
%   the uifigure is really slow to load (known issue in Matlab). Also 
%   other updates and bug fixes for helical 3D reconstruction and related
%   data structures. Also cleaned up and added large amounts of comments
%   in the code for the Trace_y poject on Github.
%



%
% Setup directories and display description
%

% Check all of this before compiling

% Log file
if nargin == 0
    % Save log for running compiled version with UI
    if ispc
        %if ~isfolder(strcat(pwd, '\Trace_y logs'))
        %    mkdir(pwd, 'Trace_y logs');
        %end
        %userpath;
        if ~isfolder([getenv('USERPROFILE') '\Documents\Trace_y'])
            mkdir([getenv('USERPROFILE') '\Documents\'], 'Trace_y');
            mkdir([getenv('USERPROFILE') '\Documents\Trace_y\'], 'logs');
        end
        logPath = [getenv('USERPROFILE') '\Documents\Trace_y\logs\Trace_y_log_' char(datetime('now', 'Format', 'yyyyMMdd')) '.txt'];
    elseif ismac
        if ~isfolder('~/Documents/Trace_y')
            mkdir('~/Documents/', 'Trace_y');
            mkdir('~/Documents/Trace_y/', 'logs');
        end
        logPath = ['~/Documents/Trace_y/logs/Trace_y_log_' char(datetime('now', 'Format', 'yyyyMMdd')) '.txt'];
    end
    diary(logPath);
end

% Descriptions, read from the README.md file
if ~isdeployed
    descrTextPath = sprintf('%s%c%s', pwd, filesep, 'README.md');
else
    descrTextPath = 'https://raw.githubusercontent.com/wfxue/Trace_y/refs/heads/main/README.md';
end
descrText = readlines(descrTextPath);

% General information
dirName = 'Trace_y.src';
appName = 'Trace_y';
% WFX is in the 8th line of the README.md file
%appAuthor = 'Wei-Feng Xue';
appAuthor = char(descrText(8));
% Version is in the 9th line of the README.md file
appVersion = char(descrText(9));
%subDir = char('common');
%subDir = [];
%subDir = cellstr(subDir);
%description = ...
%    ['Tracing filament polymers and individual filament 3D structural analysis' ...
%    newline 'Algorithms for AFM image analysis'];
% Description is in the 3rd-6th line of the README.md file
description = descrText(3)+newline+descrText(4)+newline+descrText(5)+newline+descrText(6)+newline;

fprintf('\n%s\n\n%s\n%s\n%s\n', appName, description, appAuthor, appVersion);
fprintf('\n');

%if nargin > 0
%    addpath(sprintf('%s/%s', pwd, dirName));
%end
% Comment/bypass addpath line before compiling, only need once per session
if ~isdeployed
    addpath(sprintf('%s/%s', pwd, dirName));
    if isfolder(sprintf('%s/%s/_Testing_/', pwd, dirName))
        addpath(sprintf('%s/%s/_Testing_/', pwd, dirName));
    end
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
%img = AFMimage(zeros(256), 'nm', 1000, 1000, 'nm');
img = AFMimage;

% Check if doing UI (compiled) mode or runing in Matlab command line
if nargin == 0
    % Using simple UI for compiled app
    fprintf('\nSession Started %s\n', char(datetime('now', 'Format', 'yyyyMMdd-HHmmss')))
    fprintf('Log file: %s\n', logPath);
    fprintf('Trace_y recommends starting by importing an AFM height image file.\n\n');

    h_fig = img.DispImage;
    %posMainImage = [552 190 880 880];
    %h_fig.OuterPosition = posMainImage;
    %set(gcf, 'OuterPosition', [552 190 880 880]);
    h_fig.CloseRequestFcn = @CloseImage;
    h_fig.Visible = 'off';

elseif strcmp(mode, '-m')
    % Running in Matlab commandline
    clear('img');
    return
    
end



%
% The UI bit for simple function list
%

% Init main UI window

% Definitions

% Define curent positions from bottom right and width
% This portion is very messy, needing future edits
screenSize = get(0, 'ScreenSize');
posUIwindowBottom = screenSize(4)-880-10;
posUIwindow = [1 posUIwindowBottom 550 880];

% UI elements pos relative to the main window
posFcnListHeader = [25 850-40 200 20];
posFcnListBottom = 275;
posFcnList = [25 posFcnListBottom 260 850-25-25-posFcnListBottom];

posInputTableHeader = [275 850-40 200 20];
posInputTableBottom = 580;
posInputTable = [290 posInputTableBottom 235 850-25-25-posInputTableBottom];

posDescrHeader = [290 850-posInputTableBottom+275 200 20];
posDescrBoxBottom = 273+15+12;
posDescrBox = [290 posDescrBoxBottom 235 850-25-25-265-posDescrBoxBottom];

posRunButtonBottom = 273;
%posRunButton = [275 posRunButtonBottom 250 25];
posRunButton = [290 posRunButtonBottom 235 25];

posTextBox = [25 25 500 850-25-25-560];



% Functions list
fcnList = {...
    'Select function:'; ...
    ''; ...
    '___ Suggested Workflow ___'; ...
    'Import AFM image data'; ...
    'Trace filament'; ...
    'Reconstruct helical 3D model'; ...
    'Explore features'; ...
    'Display filament CS density map'; ...
    'Save session Trace_y file (.trcy)'; ...
    'Export filament results (.csv)'; ...
    ''; ...
    '___ Loading and Saving ___'; ...
    'Export filament results (.csv)'; ...
    'Import AFM image data'; ...
    'Import Bruker Nanoscope file (.spm)'; ...
    'Import CSV text file (.csv)'; ...
    'Load Trace_y file (.trcy)'; ...
    'Save session Trace_y file (.trcy)'; ...
    ''; ...
    '___ Visualisation ___'; ...
    'Display filament'; ...
    'Display filament 3D model'; ...
    'Display filament CS density map'; ...
    'Display image'; ...
    'Display scan line'; ...
    'Explore features'; ...
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
    {'Image Type' 'Height'
    'z Unit' 'nm'; ...
    'Scan Size (x)' 2000; ...
    'Scan Size (y)' 2000; ...
    'Scan Size (xy) unit' 'nm'};
inputDispFilament = ...
    {'Filament number' 1; ...
    'Unit' 'nm'};
inputDispFilament3DModel = ...
    {'Filament number' 1; ...
    'Segment start / nm' 0; ...
    'Segment legth / nm' 500};
inputDispFilamentCS = ...
    {'Filament number' 1};
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

hUIwindow = figure(4);
set(hUIwindow, 'Name', ['*** ' appName ' ' appVersion ' ***'], ...
    'WindowStyle', 'normal', ...
    'MenuBar', 'none', ...
    'ToolBar', 'none', ...
    'Color', [1 1 1], ...
    'NumberTitle', 'off', ...
    'Resize', 'off', ...
    'Visible', 'on', ...
    'DeleteFcn', @quitTrace_y);

%hUIwindow = uifigure('Name', ['*** ' appName ' ' appVersion ' ***'], ...
%    'WindowStyle', 'normal', ...
%    'Color', [1 1 1], ...
%    'NumberTitle', 'off', ...
%    'Resize', 'off', ...
%    'Visible', 'on', ...
%    'DeleteFcn', @quitTrace_y);
    %'Position', posUIwindow, ...
    
hUIwindow.OuterPosition = posUIwindow;
drawnow;
posUIwindow(2) = hUIwindow.OuterPosition(2)-10;
hUIwindow.OuterPosition = posUIwindow;
posMainImage = posUIwindow;
posMainImage(1) = posMainImage(1)+posMainImage(3)+1;
posMainImage(3) = posMainImage(4);
h_fig.OuterPosition = posMainImage;
h_fig.Visible = 'on';

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
%{
hEchoText = uitextarea(hUIwindow, ...
    'Value', ' ', ...
    'Position', posTextBox, ...
    'UserData', [], ...
    'FontName', 'helvetica', ...
    'FontSize', 12, ...
    'FontWeight', 'normal', ...
    'HorizontalAlignment', 'left');
%}
hEchoText = uicontrol(hUIwindow, ...
    'Style', 'listbox', ...
    'Enable', 'inactive', ...
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



% Function description
uicontrol(hUIwindow, 'Style', 'text', ...
    'String', '   Description', ...
    'Units', 'pixels', 'Position', posDescrHeader, ...
    'Callback', [], 'UserData', [], ...
    'FontName', 'helvetica', ...
    'FontUnits', 'pixels', 'FontSize', 16, ...
    'FontWeight', 'normal', ...
    'ForegroundColor', 'k', ...
    'BackgroundColor', 'w', ...
    'HorizontalAlignment', 'left');

%{
hDescrBox = uitextarea(hUIwindow, ...
    'Value', ' ', ...
    'Position', posDescrBox, ...
    'UserData', [], ...
    'FontName', 'helvetica', ...
    'FontSize', 12, ...
    'FontWeight', 'normal', ...
    'HorizontalAlignment', 'left');
%}
hDescrBox = uicontrol(hUIwindow, ...
    'Style', 'edit', ...
    'Max', 2, ...
    'Min', 0, ...
    'Enable', 'inactive', ...
    'Position', posDescrBox, ...
    'UserData', [], ...
    'FontName', 'helvetica', ...
    'FontSize', 12, ...
    'FontWeight', 'normal', ...
    'HorizontalAlignment', 'left');



    function [] = fcnSelect(~, action)

        % Retrieve latest img if stored in UserData. Is the case if a
        % Callback has a waitfor, e.g. TraceFilamentM.m
        if ~isempty(h_fig.UserData)
            img = h_fig.UserData;
            h_fig.UserData = [];
            h_fig = img.ExploreFeatures;
            h_fig.OuterPosition = posMainImage;
        end

        fcnIndx = action.Source.Value;
        %action.Source.Value = 1;
        UpdateDescrBox(fcnList{fcnIndx});

        switch fcnList{fcnIndx}

            %___ Loading and Saving ___

            case 'Load Trace_y file (.trcy)'
                set(hInputTable, 'ColumnName', [], 'Data', []);
                set(hRunButton, 'Enable', 'off', ...
                    'Callback', []);

                fprintf('\n# %s\n', fcnList{fcnIndx});

                [trcyFile, trcyFilePath] = uigetfile('*.trcy;*.mat', 'Please select a .trcy file');
                load([trcyFilePath trcyFile], '-mat', 'img');
                posMainImage = h_fig.OuterPosition;
                h_fig = img.ExploreFeatures;
                h_fig.OuterPosition = posMainImage;

                fprintf('Loaded: %s\n', [trcyFilePath trcyFile]);
                UpdateLog;

            case 'Save session Trace_y file (.trcy)'
                set(hInputTable, 'ColumnName', [], 'Data', []);
                set(hRunButton, 'Enable', 'off', ...
                    'Callback', []);

                fprintf('## %s\n', fcnList{fcnIndx});

                [trcyFile, trcyFilePath] = uiputfile(strcat(img.dataFile, '.trcy'));
                save([trcyFilePath trcyFile], '-mat', 'img');

                fprintf('Saved: %s\n', [trcyFilePath trcyFile]);
                UpdateLog;

            case 'Import AFM image data'
                set(hInputTable, 'ColumnName', [], 'Data', []);
                set(hRunButton, 'Enable', 'off', ...
                    'Callback', []);

                fprintf('\n# %s\n', fcnList{fcnIndx});
                
                img = ImportAFMImageData;
                posMainImage = h_fig.OuterPosition;
                h_fig = img.ExploreFeatures;
                h_fig.OuterPosition = posMainImage;

                UpdateLog;

            case 'Import Bruker Nanoscope file (.spm)'
                set(hInputTable, 'ColumnName', [], 'Data', []);
                set(hRunButton, 'Enable', 'off', ...
                    'Callback', []);

                fprintf('\n# %s\n', fcnList{fcnIndx});
                
                [fileName, filePath] = uigetfile('*.*', 'Please select a Bruker Nanoscope file');
                img = ImportBrukerSPM([filePath fileName]);
                posMainImage = h_fig.OuterPosition;
                h_fig = img.ExploreFeatures;
                h_fig.OuterPosition = posMainImage;

                fprintf('Imported: %s\n', [filePath fileName]);
                UpdateLog;

            case 'Import CSV text file (.csv)'
                set(hInputTable, 'Data', inputImportCSV, ...
                    'ColumnName', {'Input parameter' 'Value'});
                set(hRunButton, 'Enable', 'on', ...
                    'Callback', @ImportCSVCallback);

                fprintf('\n# %s\n', fcnList{fcnIndx});
                UpdateLog;


            case 'Export filament results (.csv)'
                set(hInputTable, 'ColumnName', [], 'Data', []);
                set(hRunButton, 'Enable', 'off', ...
                    'Callback', []);

                fprintf('## %s\n', fcnList{fcnIndx});

                [exportFileName, filePath] = uiputfile(strcat(img.dataFile, '_filaments.csv'));
                img.CompileFilamentResults([filePath exportFileName]);

                fprintf('Exported: %s\n', [filePath exportFileName]);
                UpdateLog;


            %___ Visualisation ___


            case 'Display image'
                set(hInputTable, 'ColumnName', [], 'Data', []);
                set(hRunButton, 'Enable', 'off', ...
                    'Callback', []);

                fprintf('## %s\n', fcnList{fcnIndx});

                posMainImage = h_fig.OuterPosition;
                h_fig = img.DispImage;
                h_fig.OuterPosition = posMainImage;

                %fprintf('\n');
                UpdateLog;

            case 'Explore features'
                set(hInputTable, 'ColumnName', [], 'Data', []);
                set(hRunButton, 'Enable', 'off', ...
                    'Callback', []);

                fprintf('## %s\n', fcnList{fcnIndx});

                posMainImage = h_fig.OuterPosition;
                h_fig = img.ExploreFeatures;
                h_fig.OuterPosition = posMainImage;

                %fprintf('\n');
                UpdateLog;

            case 'Display scan line'
                set(hInputTable, 'ColumnName', [], 'Data', []);
                set(hRunButton, 'Enable', 'off', ...
                    'Callback', []);

                fprintf('## %s\n', fcnList{fcnIndx});
                
                img.DispLine;

                %fprintf('\n');
                UpdateLog;
                
            case 'Display filament'
                set(hInputTable, 'Data', inputDispFilament, ...
                    'ColumnName', {'Input parameter' 'Value'});
                set(hRunButton, 'Enable', 'on', ...
                    'Callback', @DispFilamentCallback);

                fprintf('## %s\n', fcnList{fcnIndx});
                UpdateLog;

            case 'Display filament 3D model'
                set(hInputTable, 'Data', inputDispFilament3DModel, ...
                    'ColumnName', {'Input parameter' 'Value'});
                set(hRunButton, 'Enable', 'on', ...
                    'Callback', @DispFilament3DModelCallback);

                fprintf('## %s\n', fcnList{fcnIndx});
                UpdateLog;

            case 'Display filament CS density map'
                set(hInputTable, 'Data', inputDispFilamentCS, ...
                    'ColumnName', {'Input parameter' 'Value'});
                set(hRunButton, 'Enable', 'on', ...
                    'Callback', @DispFilamentCSCallback);

                fprintf('## %s\n', fcnList{fcnIndx});
                UpdateLog;


            %___ Structural Analysis ___

            
            case 'Trace filament'
                set(hInputTable, 'Data', inputTraceFilamentM, ...
                    'ColumnName', {'Input parameter' 'Value'});
                set(hRunButton, 'Enable', 'on', ...
                    'Callback', @TraceFilamentMCallback);
                
                fprintf('## %s\n', fcnList{fcnIndx});
                UpdateLog;


            case 'Reconstruct helical 3D model'
                set(hInputTable, 'Data', inputMakeHelicalFilamentModel, ...
                    'ColumnName', {'Input parameter' 'Value'});
                set(hRunButton, 'Enable', 'on', ...
                    'Callback', @MakeHelicalFilamentModelCallback);
                
                fprintf('## %s\n', fcnList{fcnIndx});
                UpdateLog;



            case 'Quit'
                fprintf('## %s\n', fcnList{fcnIndx});
                UpdateLog;

                close(hUIwindow);
                
                fprintf('\n');
                return
        end

        UpdateLog;
    end


    function [] = ImportCSVCallback(~, ~)
        inputImportCSV = hInputTable.Data;

        %[fileName, pathName] = uigetfile('*.csv;*.txt');
        %img = ImportCSV(strcat(pathName, fileName), ...
        [img, filePath] = ImportCSV(...
            inputImportCSV{1, 2}, inputImportCSV{2, 2}, ...
            inputImportCSV{3, 2}, inputImportCSV{4, 2}, ...
            inputImportCSV{5, 2});
        posMainImage = h_fig.OuterPosition;
        h_fig = img.ExploreFeatures;
        h_fig.OuterPosition = posMainImage;

        %fprintf('%s loaded\n\n', strcat(pathName, fileName));
        fprintf('Imported: %s\n', filePath);
        UpdateLog;
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

    function [] = DispFilamentCSCallback(~, ~)
        inputDispFilamentCS = hInputTable.Data;
        img.DispFilamentCS(inputDispFilamentCS{1, 2});
    end

    function [] = TraceFilamentMCallback(~, ~)
        inputTraceFilamentM = hInputTable.Data;

        % Retrieve latest img if stored in UserData. Is the case if a
        % Callback has a waitfor, e.g. TraceFilamentM.m
        if ~isempty(h_fig.UserData)
            img = h_fig.UserData;
            h_fig.UserData = [];
        end

        % Callback with a waitfor line, making sure not to double save
        % results as img is passed in the UserData figure property
        img.TraceFilamentM(inputTraceFilamentM{1, 2}, inputTraceFilamentM{2, 2}, [], [], h_fig);
    end

    function [] = MakeHelicalFilamentModelCallback(~, ~)
        inputMakeHelicalFilamentModel = hInputTable.Data;
        img.MakeHelicalFilamentModel(inputMakeHelicalFilamentModel{1, 2}, inputMakeHelicalFilamentModel{2, 2}, ...
            inputMakeHelicalFilamentModel{3, 2}, [], inputMakeHelicalFilamentModel{4, 2}, ...
            inputMakeHelicalFilamentModel{5, 2}, inputMakeHelicalFilamentModel{6, 2});
        UpdateLog;
    end



% Log text
UpdateLog;
    function [] = UpdateLog()
        logText = readlines(logPath);
        %set(hEchoText, 'Value', logText);
        %scroll(hEchoText, 'bottom');

        set(hEchoText, 'String', logText, 'Value', numel(logText));
    end



% Description text
UpdateDescrBox;
    function [] = UpdateDescrBox(fcnLine)

        % For finding ext in source code
        % Not used, use README.md instead
        %{
        lineIdx = [1 1];
        descrFound = -1;
        while descrFound ~= 1
            nextline = char(descrText(lineIdx(2)));

            if descrFound == -1 && ~isempty(nextline) && nextline(1) == '%'
                descrFound = 0;
                lineIdx(2) = lineIdx(2)+1;
            elseif descrFound == -1 && (isempty(nextline) || nextline(1) ~= '%')
                lineIdx = lineIdx+1;
            elseif descrFound == 0 && ~isempty(nextline) && nextline(1) == '%'
                lineIdx(2) = lineIdx(2)+1;
            else
                descrFound = 1;
            end

        end

        descrText = descrText(lineIdx(1):lineIdx(2));
        %}

        if ~exist('fcnLine', 'var') || isempty(fcnLine)
            % On app launch
            %set(hDescrBox, 'Value', descrText);
            %scroll(hDescrBox, 'top');
            set(hDescrBox, 'String', descrText);
        else
            % Match the function line text
            lineIdx = [1 1];
            descrFound = -1;
            while descrFound ~= 1
                nextline = char(descrText(lineIdx(2)));
                if strcmp(nextline, sprintf('### %s', fcnLine))
                    descrFound = 0;
                    lineIdx(2) = lineIdx(2)+1;
                elseif descrFound == 0 && ~isempty(regexp(nextline, '#', 'once'))
                    descrFound = 1;
                    lineIdx(2) = lineIdx(2)-1;
                elseif descrFound == 0
                    lineIdx(2) = lineIdx(2)+1;
                elseif lineIdx(2) >= numel(descrText)
                    lineIdx = [1 numel(descrText)];
                    descrFound = 1;
                else
                    lineIdx = lineIdx+1;
                end
            end
            %set(hDescrBox, 'Value', descrText(lineIdx(1):lineIdx(2)));
            set(hDescrBox, 'String', descrText(lineIdx(1):lineIdx(2)));
            %scroll(hDescrBox, 'top');
        end
    end


% On closing image window

    function [] = CloseImage(~, ~)

        % Open an empty new image window
        img = AFMimage;
        posMainImage = h_fig.OuterPosition;
        h_fig = img.DispImage;
        h_fig.OuterPosition = posMainImage;
        h_fig.CloseRequestFcn = @CloseImage;
    end



% On quit
    function [] = quitTrace_y(~, ~)
        fprintf('\n\nSession finished %s\n\n', char(datetime('now', 'Format', 'yyyyMMdd-HHmmss')))
        diary('off');
        %close(hUIwindow);
        delete(h_fig);
        if isdeployed
            close('all');
        end
    end



% Finishing output
waitfor(hUIwindow);
if nargout == 1
    img_out = img;
end



end