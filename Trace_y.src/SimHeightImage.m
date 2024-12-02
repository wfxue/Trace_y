function [img] = SimHeightImage(imgSize, pixelDensity, ...
    objectType, objectData, tipType, tipData, noiseType, noiseData, precision, imgShow)

%
% DESCRIPTION
% – Simulate an AFM image of a object or a regular geometry based on hard
% tip-sample contacts.
% – Can use iso-surfaces as input, e.g. generated from cryo-EM density maps
% – Part of Trace_y
%
% USAGE
% – Simulate a sphere with radius 5 on flat surface with spherical tip with
% radius 2 nm and Gaussian noise with stadard dev sigma of 0.5 nm. The 
% simulated image is 64x64 px and 2 pixels per nm and is shown.
% img = SimHeightImage(64, 2, 'sphere', 5, 'rounded_cone', [2 18], 'gaussian_noise', 0.5, 4, true);
% – Same as above but default no noise and with a random suface with 
% 'frequency' around per 10 and peaks of max 20 nm.
% img = SimHeightImage(64, 2, 'random', [10 20], 'rounded_cone', [2 18]);
% – Same as above with defaults with an EM map iso-surface in fv format.
% The shape is rotated 90 around y-axis before simulation. The map in fv
% format can be obtained using ImportEMDB.m
% img = SimHeightImage(64, 2, 'em_map', {mapfv [0 90 0]}, 'rounded_cone', [2 18]);
% – Same as above with defaults with an helical filament EM map 
% iso-surface in fv format. The map should in fv format using ImportEMDB.m
% should be axis aligned using CalcHelicalAxisEM.m and the maximum radius
% also calculated using CompileHelicalFilamentInfo.m
% img = SimHeightImage(64, 2, 'em_map_helical', {mapfv_aligned csa_max}, 'rounded_cone', [2 18]);
%
% INPUTS:
% imgSize  –  Number of pixels in each dimension, image will be square
% pixelDensity  –  Pixel density, pixel/nm
% objectType  –  The surface geometry. Current options are listed
% below.
% objectPar  –  Parameters describing the object, see list below
%   objectType       : objectPar
%   'sphere'         : radius
%   'cylinder'       : radius
%   'rhombic_prism'  : side length
%   'random'         : {rep_frequency, amplitude}
%   'random_x'       : {rep_frequency, amplitude}
%   'em_map'         : {mapfv, rotation_angles}
%   'em_map_helical' : {mapfv_aligned, maximum_radius}
% tipType  –  Tip function or data type, e.g. 'rounded_cone'.
% tipData  –  Parameters describing the tip geometry, for example [2 18] for 
% a rounded cone with 2 nm tip-radius and 18 deg tip angle from the 
% central axis.
% noiseType  –  Type of noise e.g. 'gaussian', or 'no-noise'
% precision  –  Subdivisions of pixel for better dilation calc etc. Default
% is 4
% imgShow  –  Flag for displaying the sAFM image or not.
% 
% OUTPUTS
% img  –  AFM image object containing a simulated image.
%
% DEPENDENCIES
% – Uses Trace_y's @AFMimage/ object defs and methods.
% – For EM maps, maps in fv format can be generated using: 
%     ImportEMDB.m
%     CalcHelicalAxisEM.m
%     CompileHelicalFilamentInfo.m
% – Uses:
%     SurfaceModel_EMmap.m
%     Rotate3D.m
%
% AUTHORS
% Wei-Feng Xue and Liisa Lutter
% 
% HISTORY
% 2019.04  –  LL initial draft edited by WFX. Each of the steps were very 
% slow, so WFX added code that speed things up. Original code in 
% old AFM_Image_Sim function.
% 2019.05  –  LL drafted the aditional option to use a 3D model, with WFX 
% edited the new 3D model option and tidied up code a bit.
% 2019.07  –  Moved all pixel coordinates to go from [1 1] to introduce 
% pixeldensity and size parameters etc., and moved to dilation based 
% method, together with other edits should make things go a bit faster.
% 2020.05  –  LL drafted the import and sim code for using PDB and MRC
% structures, with wfx edited the algorithms. Tidied up and organised the 
% code and harmonised the code with Trace_y in creating AFMimage object
% 2024.08  –  WFX updated code to be consistent with updates with updates 
% with helical EM map simulation code.
%



%
% Set up progress bar
%

tic;

% To show simulated image and progress bar or not
if ~exist('imgShow', 'var') || isempty(imgShow)
    imgShow = true;
end

if imgShow
    hwb = waitbar(0, sprintf('Setting up the image... Elapsed time: %.2f s', toc));
    pos = get(hwb, 'Position');
    set(hwb, 'Position', pos+[0 pos(4)+36 0 0]);
end

%
% Create new AFMimage object
%
img = AFMimage(zeros(imgSize), 'nm', (imgSize-1)/pixelDensity, (imgSize-1)/pixelDensity, 'nm');

% Set up image grid coordinates
% set x and y = 0 at centre pixel
x = (0:imgSize-1)./pixelDensity-floor(imgSize/2)./pixelDensity;
y = ((0:imgSize-1)./pixelDensity-floor(imgSize/2)./pixelDensity)';
img.x = x;
img.y = y;



%
% Defaults
%

% The default surface object to simulate is a sphere (for testing mainly)
if ~exist('objectType', 'var') || isempty(objectType)
    % Function for spherical dome
    objectType = 'sphere';
    objectData = 5;
end

% The default tip geometry is rounded cone
if ~exist('tipType', 'var') || isempty(tipType)
    % Function for cone with spherical dome
    tipType = 'rounded_cone';
    % Parameters are tip_r (2 nm default) and tip_a (18 deg default)
    % Tip angle is very insensitive as the sphere at the tip does all the
    % contact usually
    tipData = [2, 18];
end

% The default noise is no noise
if ~exist('noiseType', 'var') || isempty(noiseType)
    % Function for no noise
    noiseType = 'none';
    noiseData = 0;
end

% Default precition for dialation
if ~exist('precision', 'var') || isempty(precision)
    precision = 4;
end



%
% Handles the tip model
%

% Currently only usues the rounded cone
% But has future potential for custom/reconstructed tips
img = img.SetTipModel(tipType, tipData);



%
% Handles the surface object
%

% List can be extended
% Default input objectData should be supplied (if simple)
switch objectType
    case 'sphere'
        surf_f = @Sphere_upper_surface;
        % objectData is sphere radius, default 5 units
        if isempty(objectData)
            objectData = 5;
        end

    case 'cylinder'
        surf_f = @Cylinder_upper_surface;
         % objectData is cross-section radius
        if isempty(objectData)
            objectData = 5;
        end

    case 'rhombic_prism'
        % A symetric Right angled rhombic prism standing on one corner
        surf_f = @Rhombus_upper_surface;
        % objectData is corner to coner distance (same for both)
        if isempty(objectData)
            objectData = 5;
        end
    case 'random'
        % Random  surface
        % Node 'frequency' (number of pixels) and amplitude / nm as inputs
        surf_f = @Random_surface;

        % Input objectData is rough 'frequency' and amplitude
        if isempty(objectData)
            freq = 10;
            amp = 10;
        else
            freq = objectData(1);
            amp = objectData(2);
        end

        x_rand = [min(x(:)):round(freq):max(x(:)) max(x(:))];
        x_rand = unique(x_rand);
        n = length(x_rand);
        y_rand = x_rand;
        [x_rand, y_rand] = meshgrid(x_rand, y_rand);
        z_rand = amp.*rand(n, n);

        objectData = cat(3, x_rand, y_rand, z_rand);
        % The random numbers are then saved in objectData
        
    case 'random_x'
        % Random cross-section in x, basically random in 1D
        % Node 'frequency' (number of pixels) and amplitude / nm as inputs
        surf_f = @Random_x;

        % Input objectData is rough 'frequency' and amplitude
        if isempty(objectData)
            freq = 5;
            amp = 10;
        else
            freq = objectData(1);
            amp = objectData(2);
        end

        x_rand = [min(x(:)):round(freq):max(x(:)) max(x(:))]';
        z_rand = amp.*rand(size(x_rand));

        objectData = [x_rand z_rand];
        % The random numbers used to generate cross-section are then saved
        % in objectData
        
    case 'em_map'
        % Object from EM map, use ImportMRC or ImportEMDB
        % objectData input is the iso-surface (mapfv) and rotation angles 
        % [x y z]
        surf_f = @SurfaceModel_EMmap;
        
        % Progress bar for this one, will take a bit time
        if imgShow
            waitbar(0.2, hwb, sprintf('Calculating EM map object alpha shape... Elapsed time: %.2f s', toc));
        end
        
        % em_map is isosurface data from MRC file using ImportMRC function
        em_map = objectData{1};
        r_angles = objectData{2};
        % Rotate map, angles in degrees
        em_map.vertices = Rotate3D(em_map.vertices, r_angles);
        
        % Centre molecule for simulations
        % Also put min z value to 0 (at the surface)
        em_map.vertices(:, 1) = em_map.vertices(:, 1)-mean(em_map.vertices(:, 1));
        em_map.vertices(:, 2) = em_map.vertices(:, 2)-mean(em_map.vertices(:, 2));
        em_map.vertices(:, 3) = em_map.vertices(:, 3)-min(em_map.vertices(:, 3));
        
        % Assume aspect ratio of 1, it should always be 1 across xyz
        % Try alpha shape
        % 2*voxel size to save a bit of computation as tip radius likely to
        % be a lot bigger
        voxSize = em_map.UserData.voxSize;
        %alpha_val = 2*voxSize;
        %p = em_map.vertices;
        
        % First use alphashape to find largest region
        % Basically suppresing noise by using region 1 (the biggest region)
        %em_map = alphaShape(...
        %    em_map.vertices(:, 1), em_map.vertices(:, 2), em_map.vertices(:, 3), ...
        %    alpha_val, 'HoleThreshold', 0, 'RegionThreshold', 10*voxSize.^3);
        %em_map = alphaShape(p(:, 1), p(:, 2), p(:, 3), alpha_val, ...
        %    'HoleThreshold', max(em_map.UserData.CELLA./10).^3);
        
        %objectData = em_map;
        % The Alpha shape from the iso-surface is now saved in objectData
        % which later is saved in img.zS
        % Now saving the xyz tri-formated coordinates instead as Alpha
        % shape is evaluated on every load and access even if nothing is
        % changed.
        xyz = em_map.vertices;
        objectData = {xyz voxSize};


     case 'em_map_helical'
        % Object is axis aligned helical structure from EM map
        % Assume helical screw axis already aligned to z-axis using 
        % CalcHelicalAxisEM.m
        % objectData input is the axis aligned iso-surface (mapfv_aligned) 
        % with helical info in UserData and csa_max is also needed

        % The surface function is the same as 'em_map' just the input map
        % is helical
        surf_f = @SurfaceModel_EMmap;
        
        % The extension code below has bugs but if input is alphaShape 
        % already then the simulation works. 
        % Below needs fixing. 
        % Is now fixed and in a separate file:
        % MakeExtendedHelicalFilamentEM.m
        %{
        if numel(objectData) > 1
            % Progress bar for this one, will take a bit time
            waitbar(0.2, hwb, sprintf('Calculating EM map object alpha shape... Elapsed time: %.2f s', toc));

            em_map = objectData{1};
            twist = objectData{2};
            rise = objectData{3};
            radius = objectData{4};

            voxSize = em_map.UserData.voxSize;

            % Check symmetry
            % Not completly good for higher sym
            sym = abs(round(360/twist));
            if sym >= 2 && sym <= 12
                twist = twist*sym-360;
                rise = rise*sym;
            end
            % Convert rise to nm
            rise = rise/10;

            p = em_map.vertices;

            % Estimate fibril width if needed
            if isempty(radius)
                [~, rho] = cart2pol(p(:, 1), p(:, 2));
                radius = max(rho);
            end

            %
            % Calc max length in nm required
            l_img = imgSize/pixelDensity;

            % Extend the filament along the z axis
            % Overlap with 1/4
            % Making sure to extend integer number of rises
            %p = em_map.vertices;
            l_map = max(p(:, 3))-min(p(:, 3));
            n_rise = l_map/rise;
            n_rise = floor(3*n_rise/4);

            p2 = p;
            while l_map < 1.1*l_img
                p2(:, 3) = p2(:, 3)+n_rise*rise;
                p2 = Rotate3D(p2, [0 0 -twist*n_rise]);
                p = cat(1, p, p2);
                l_map = max(p(:, 3))-min(p(:, 3));
            end
            %

            % Reduce sampling
            % Use pointCloud class, requires Computer Vision toolbox
            p = pointCloud(p);
            p = pcdownsample(p, 'gridAverage', 2*voxSize);
            %p = pcdenoise(p);
            p = p.Location;


            % Rotate axis to y
            p = Rotate3D(p, [90 0 0]);

            % Centre molecule for simulations
            %p(:, 1) = p(:, 1)-mean(p(:, 1));
            p(:, 2) = p(:, 2)-mean(p(:, 2));
            % Also put min z value to 0 (at the surface)
            %p(:, 3) = p(:, 3)-min(p(:, 3));
            p(:, 3) = p(:, 3)+radius;

            % Assume aspect ratio of 1
            % Try alpha shape
            % 2*voxel size to save a bit of computation
            alpha_val = 2*voxSize;


            % First use alphashape to find largest region
            % Basically suppresing noise by not using small regions
            em_map = alphaShape(p(:, 1), p(:, 2), p(:, 3), alpha_val, ...
                'HoleThreshold', max(em_map.UserData.CELLA./10).^3);

            objectData = em_map;
        end
        %}

        % Progress bar for this one, will take a bit time
        if imgShow
            waitbar(0.2, hwb, sprintf('Calculating EM map object alpha shape... Elapsed time: %.2f s', toc));
        end
        
        % objectData input is the axis aligned iso-surface
        mapfv_aligned = objectData{1};
        csa_max = objectData{2};

        % Calc max length in nm required
        l_img = imgSize/pixelDensity;
        % Extend the aligned map and making an alpha shape
        xyz = MakeExtendedHelicalFilamentEM(mapfv_aligned, l_img);

        % Estimating max cross-sectional radius
        %[~, rho] = cart2pol(xyz(:, 1), xyz(:, 2));
        %rho = max(rho(:));
        % Now needs to be supplied as this should be calculated from
        % downsampled and denoised point set in cross-section analysis in
        % CompileHelicalFilamentInfo.m

        % Rotate axis from z to y as for AFM simulations the filament is laying
        % down and z axis is researvd for height
        xyz = Rotate3D(xyz, [90 0 0]);

        % Move the filament to the surface (the z axis now in AFM axes)
        %xyz(:, 3) = xyz(:, 3)+rho;
        xyz(:, 3) = xyz(:, 3)+csa_max;

        % Construct alpha shape
        % 2*voxel size to save a bit of computation as tip radius likely to
        % be a lot bigger
        % Also use alphashape to find largest region
        % Basically suppresing noise by using region 1 (the biggest region)
        voxSize = mapfv_aligned.UserData.voxSize;
        %alpha_val = 2*voxSize;

        %em_extended = alphaShape(p(:, 1), p(:, 2), p(:, 3), tip_radius, ...
        %    'HoleThreshold', max(mapfv_aligned.UserData.CELLA./10).^2', 'RegionThreshold', 0.01);
        %em_map = alphaShape(xyz(:, 1), xyz(:, 2), xyz(:, 3), alpha_val, ...
        %    'HoleThreshold', 0, 'RegionThreshold', 10*voxSize.^3);
        %objectData = em_map;

        objectData = {xyz voxSize};
end


%
% Surface functions
%
    function z = Sphere_upper_surface(x, y, r)
        % Sphere centrerd at middle of x and y range
        % Is a dome as lower surface not sensed
        z = sqrt(r.^2-(x).^2-(y).^2)+r;
        z = real(z)-(imag(z) > 0).*r;
    end

    function z = Cylinder_upper_surface(x, ~, r)
        % Cylinder centrerd at middle of x
        z = sqrt(r.^2-(x).^2)+r;
        z = real(z)-(imag(z) > 0).*r;
    end

    function z = Rhombus_upper_surface(x, ~, r)
        % Right angled symetric rhombus centrerd at middle of x
        z = -abs(x)+2*r;
        z(abs(x) > r) = 0;
    end

    function z = Random_surface(x, y, p)
        z = interp2(p(:, :, 1), p(:, :, 2), p(:, :, 3), x, y, 'spline');
    end

    function z = Random_x(x, ~, p)
        z = interp1(p(:, 1), p(:, 2), x, 'spline');
    end

% for EM maps, use:
% z = SurfaceModel_EMmap(x, y, p)
% External file

% zS stored in function format
img.zS = {surf_f objectData};
% Use this syntax to recal zS: zS = img.zS{1}(xS, yS, img.zS{2});
% This is done is Dilation.m

% Test code
%imagesc(img.zS{1}(xS, yS, img.zS{2}));
%set(gca, 'YDir', 'normal');



%
% Perform Dilation algorithm
%

% It will also calculate contact points
if imgShow
    waitbar(0.6, hwb, sprintf('Dilating surface function... Elapsed time: %.2f s', toc));
end
img = img.Dilation(precision, imgShow);



%
% Add noise if optioned
%

switch noiseType
    case 'no_noise'
        % Default
        img.zBg = zeros(size(img.z));
    case 'gaussian_noise'
        %rng('shuffle');
        sigma = noiseData;
        img.zBg = normrnd(0, sigma, size(img.z));
        img.zNoiseStd = sigma;
        img.z = img.z+img.zBg;
end



%
% Contact point reconstruction
%
% For comparison
% This is fast now
if imgShow
    waitbar(0.8, hwb, sprintf('Contact point reconstruction... Elapsed time: %.2f s', toc));
end
img = img.Reconstruction;



%
% Add extra info
%

% Save master image
img.xMaster = img.x;
img.yMaster = img.y;
img.zMaster = img.z;

img.imgType = 'Simulated height image';

img.pixelPerLine = imgSize;
img.nLines = imgSize;
img.pixelPerLineImg = imgSize;
img.nLinesImg = imgSize;
img.scanSizex = max(x)-min(x);
img.scanSizey = max(y)-min(y);

img.zUnit = 'nm';
img.scanSizeUnit = 'nm';

% Estinate noise std based on 1 sigma procentile
img.zNoiseStd = 0;

% Resolution
img.xResolution = 1./pixelDensity;
img.yResolution = 1./pixelDensity;

% Colour scale
signal = sort(abs(img.z(:)));
idx = ceil(0.999*length(signal));
img.cScale = [-ceil(signal(idx)) ceil(signal(idx))];



%
% Finishing
%
if imgShow
    
    fprintf('sAFM image simulation ');
    toc;

    delete(hwb);
    img.DispImage;
end



end