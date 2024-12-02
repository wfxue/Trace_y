classdef Filament < dynamicprops
    
%
% DESCRIPTION
% – Class def and constructor for @Filament class object centrally used in
% Trace_y
% – Is of a dynamicprops obj type so that new properties can be easily
% added in the future if needed
% – Part of Trace_y
%
% USAGE
% – Constructor of an Filament object containing some data
% >> obj = Filament(img, xf, yf, tracingMethod, tracingData);
%
% INPUTS
% img  –  @AFMimage object.
% xf, yf  –  The x and y subpixel coordinates of the filament central axis
% tracingMethod  –  Function handle to the method used to trace and
% estimate the central axis along with parameters used for reference
% tracingData  –  The any extra data saved from the trace of the central
% axis, for example sigma and angles etc, if a wall model is used
%
% OUTPUTS
% obj  –  The Filament object created
%
% DEPENDENCIES
% Used thoughout in Trace_y and by its methods. Also used as part of 
% @AFMimage object for defining segmented features.
%
% AUTHORS
% Wei-Feng Xue
%
% HISTORY
% 2019.01  –  Initial draft too add to the @AFMimage class def for v2 of 
% Trace_y. Updated with @AFMimage class def
% 2024.09  –  Slight tweaks. Updates and comments added for v7 of Trace_y 
% for Github page update.
%



    % Semi-flexible filament
    
    properties
        x;
        y;
        z;
        l;
        lStep;
        % phi is the angular change
        phi;
        % theta is the angle in the image plane
        theta;
        lContour
        rEE
        xResolution
        yResolution
        xyUnit
        zUnit
        tracingMethod;
        tracingData;
        appWidth;
        traceR;
        traceRMSD;
        isAtImageBoundary;
        isSegment;
        % Straightened filament data
        xf;
        yf;
        zf;
        % Tip-sample corrected pixel data
        xc;
        yc;
        zc;
        unreliability;
        % 3D reconstruction
        helicalFilament3DModel;
    end
    
    methods
        function obj = Filament(img, xf, yf, tracingMethod, tracingData)
        %function obj = Filament(img, xf, yf, xr, yr, appWidth, tracingMethod)
            
            obj.tracingMethod = tracingMethod;
            appWidth = tracingMethod{2};
            obj.appWidth = appWidth;
            obj.isSegment = true;
            
            % Image info
            % resolution in nm/pixel
            obj.xResolution = img.xResolution;
            obj.yResolution = img.yResolution;
            obj.xyUnit = 'pixel';
            obj.zUnit = img.zUnit;
            
            switch tracingMethod{1}
                case 'TraceFilament'
                    
                    % This is obsolete
                    %{
                    % Check trace and retrace length
                    if abs(length(xf)-length(xr)) > appWidth
                        obj.isSegment = true;
                    else
                        obj.isSegment = false;
                    end
                    
                    xf = flipud(xf);
                    yf = flipud(yf);
                    
                    steps = min([length(xf) length(xr)]);
                    
                    % Check trace and retrace differences
                    obj.traceR = sqrt((xf(1:steps)-xr(1:steps)).^2+...
                        (yf(1:steps)-yr(1:steps)).^2);
                    obj.traceRMSD = sqrt(mean(obj.traceR.^2));
                    traceRlarge = find(obj.traceR > appWidth, 1, 'first');
                    if ~isempty(traceRlarge)
                        steps = traceRlarge-1;
                        obj.isSegment = true;
                    end
                    
                    obj.x = mean([xf(1:steps) xr(1:steps)], 2);
                    obj.y = mean([yf(1:steps) yr(1:steps)], 2);
                    
                    % Angles
                    obj.phi = [atan2d((obj.y(2:end)-obj.y(1:end-1)), ...
                        (obj.x(2:end)-obj.x(1:end-1))); NaN];
                    obj.theta = [obj.phi(1); ...
                        MeanAngled2(obj.phi(2:end-1), obj.phi(1:end-2)); ...
                        obj.phi(end-1)];
                    %}
                    
                case 'TraceFilamentM'
                    obj.x = xf;
                    obj.y = yf;
                    obj.tracingData = tracingData;
                    
            end
            
            % Interpolate z values
            xrange = floor(min(obj.x))-4:ceil(max(obj.x))+4;
            xrange = xrange(xrange > 0 & xrange <= img.pixelPerLineImg);
            yrange = floor(min(obj.y))-4:ceil(max(obj.y))+4;
            yrange = yrange(yrange > 0 & yrange <= img.nLinesImg);
            obj.z = interp2(xrange, yrange, img.z(yrange, xrange), ...
                obj.x, obj.y, 'cubic', 0);
            
            % Distances
            obj.lStep = [0; sqrt((obj.x(2:end)-obj.x(1:end-1)).^2+...
                (obj.y(2:end)-obj.y(1:end-1)).^2)];
            obj.l = cumsum(obj.lStep);
            
            obj.lContour = obj.l(end);
            obj.rEE = sqrt((obj.x-obj.x(1)).^2+...
                (obj.y-obj.y(1)).^2);
            
            % Angles
            % Need to update
            %{
            % Vectors to the relevant points
            v = [obj.x obj.y];
            % Calculate tangent vectors
            tangent = v(1:end-1, :)-v(2:end, :);
            % Calculate central differences
            tangent = [0 0; tangent]+[tangent; 0 0];
            % Find tangent angles and points using polar tranformations
            obj.theta = cart2pol(tangent(:, 1), tangent(:, 2))';
            %}

            % Check if termination at boundary
            if (min(obj.x) < appWidth+1) || ...
                    (max(obj.x) > size(img.z, 2)-appWidth-1) || ...
                    (min(obj.y) < appWidth+1) || ...
                    (max(obj.y) > size(img.z, 1)-appWidth-1)
                obj.isAtImageBoundary = true;
                
            else
                obj.isAtImageBoundary = false;
            end
            
            
        end
    end
end