classdef AFMimage
    
%
% DESCRIPTION
% – Class def and constructor for @AFMimage class object centrally used in
% Trace_y
% – Part of Trace_y
%
% USAGE
% – Constructor of a blank img object
% >> img = AFMimage;
% – Constructor of an img object containing some data
% >> img = AFMimage(z, zUnit, scanSizex, scanSizey, scanSizeUnit)
%
% INPUTS
% z  –  The z (image) value, for example height in nm
% zUnit  –  The unit of z, for example 'nm'
% scanSizex  –  The scansize of the x-axis (the scan lines), for example 
% 2000 nm
% scanSizey  –  The scansize of the y-axis (slow scan axis), for example 
% 2000 nm
% scanSizeUnit  –  The unit of the scan sizes, for example 'nm'
%
% OUTPUTS
% obj  –  The AFMimage object created
%
% DEPENDENCIES
% Used thoughout in Trace_y and by its methods.
%
% AUTHORS
% Wei-Feng Xue
%
% HISTORY
% 2014.03  –  First version drafted with xAFMTools
% 2014.04  –  Used as part of Trace_y version 1.0
% 2018.11  –  WFX edit, adding various old functions in xAFMTools/Trace_y 
% as class methods
% 2019.02  –  WFX edit for v2 of Trace_y
% 2020.03  –  WFX edit for v3 of Trace_y, added zS, zR and zT for Sample,
% Reconstruction and and Tip data. This structure will be used for future 
% simulations, deconvolutions and tip espimates in a better organised way
% 2022.10  –  WFX edit for v6, revamped to make adding new formats easier
% 2024.09  –  Slight tweak to remove filetype specific parameters. Updates
% and comments added for v7 of Trace_y for Github page update
%


    
    properties
        % General
        version;
        
        % Data file
        dataFile;
        dataFileFormat;
        imgType;
        dataFileProp;
        pixelPerLine;
        nLines;
        scanSizex;
        scanSizey;
        scanSizeUnit;
        zUnit;
        %zSens
        %zScale
        
        % Master image size
        pixelPerLineImg;
        nLinesImg;
        
        % coordinate values I
        x;
        y;
        z;
        zBg;
        xMaster;
        yMaster;
        zMaster;
        
        % Sample surface S
        zS;
        
        % Tip T
        zT;
        
        % Reconstructed R
        xR;
        yR;
        zR;
        zRUncertainty;
        
        % Calculated contact-points C
        xC;
        yC;
        zC;
        
        % Other properties
        zNoiseStd;
        xResolution;
        yResolution;
        cMap;
        cScale;
        
        % Objects on the image
        features;
    end
    methods
        
        % Constructor
        function obj = AFMimage(z, zUnit, scanSizex, scanSizey, scanSizeUnit)


            %
            % Defaults
            %
            if ~exist('z', 'var') || isempty(z)
                z = zeros(256);
                obj.imgType = 'Height';
            end
            if ~exist('zUnit', 'var')
                zUnit = 'nm';
            end
            if ~exist('scanSizex', 'var') || isempty(scanSizex)
                scanSizex = 256;
            end
            if ~exist('scanSizey', 'var') || isempty(scanSizey)
                scanSizey = 256;
            end
            if ~exist('scanSizeUnit', 'var')
                scanSizeUnit = 'nm';
            end


            % Version
            obj.version = 7;

            % Image data and dimensions
            obj.z = z;
            obj.zUnit = zUnit;
            obj.scanSizex = scanSizex;
            obj.scanSizey = scanSizey;
            obj.scanSizeUnit = scanSizeUnit;

            obj.pixelPerLine = size(z, 2);
            obj.nLines = size(z, 1);

            % Estinate noise std based on 1 sigma procentile
            obj.zNoiseStd = prctile(abs(obj.z(:)), 100*erf(1/sqrt(2)));

            % Resolution
            obj.xResolution = obj.scanSizex/(obj.pixelPerLine-1);
            obj.yResolution = obj.scanSizey/(obj.nLines-1);

            % x and y
            obj.x = 0:obj.xResolution:obj.scanSizex+0.1*obj.xResolution;
            obj.y = (0:obj.xResolution:obj.scanSizey+0.1*obj.xResolution)';

            % Masters
            obj.zMaster = obj.z;
            obj.xMaster = obj.x;
            obj.yMaster = obj.y;
            obj.pixelPerLineImg = size(z, 2);
            obj.nLinesImg = size(z, 1);
            
            % Colour scale
            %{
            signal = sort(abs(obj.z(:)));
            idx = ceil(0.999*length(signal));
            obj.cScale = [-ceil(signal(idx)) ceil(signal(idx))];
            if obj.cScale(1) == obj.cScale(2)
                obj.cScale = [-1 1];
            end
            %}
            signal = max(abs(obj.z(:)));
            if signal > 1
                obj.cScale = [-ceil(signal) ceil(signal)];
            elseif signal == 0
                obj.cScale = [-1 1];
            else
                obj.cScale = [-signal signal];
            end

            % Default placeholders
            %obj.cScale = [-1 1];
            obj.cMap = AFMimage.CmapAFMstd(2^8);
            %obj.features = AFMimageObjects;
            obj.features = struct;
            obj.features.filaments = Filament.empty;
            obj.zBg = 0;
            %obj.zMaster = 0;
            obj.dataFileProp = struct;
            obj.dataFileProp.format = [];
        end


        % Other methods
        % All methods have been moved out to separate files

        %{
        % Ths section needs revamp

        % import a image in Bruker format
        obj = ImportImageBruker(obj, filePath);
        
        % write back (modify) an image in Bruker format
        SaveImageBruker(obj, filePath)
        
        % Display the image
        DispImage(obj);
        
        % Display a scan line, l is optional input line number
        DispLine(obj, l)
        
        % Calc and display a section, x1, y1, x2, y2 are the optional 
        % input for start and end points
        s = CalcSection(obj, x1, y1, x2, y2, doplot)
        %}
    end
    
    methods (Static)
        % Static methods
        % The default AFM colours, use n = 2^8 as default;
        c = CmapAFMstd(n);
        % Surface geometry model LSQ
        [model, p, residual, zcalc, ssr, exitflag] = ...
            MakeImgModelLSQ(x, y, z, w, p0, lb, ub, model)
    end
    
end