function [output] = CompileFilamentResults(img, exportFileName)

%
% DESCRIPTION
% – Compiles and exports all available filament parameters into a csv
% table.
% – Part of Trace_y
%
% USAGE
% Standard method usage with filament number as input:
% >> img.CompileFilamentResults;
% Specify the full path for the saved file instead of opening a save
% dialogue window:
% >> img.CompileFilamentResults(exportFileName);
%
% INPUTS
% img  –  AFMimage object. Implicit with dot notation
% exportFileName –  Full path to the save-file location
%
% OUTPUTS
% resultTable  –  The result table that is saved. This table is also saved
% onto a file.
%
% DEPENDENCIES
% – Method for Trace_y's @AFMimage/ object
%
% AUTHORS
% Wei-Feng Xue
%
% HISTORY
% 2022.08  –  First draft
%


% Number of filaments in the img object
nFilaments = numel(img.features.filaments);
% Image resolution nm/pixel
xyRes = img.xResolution;

% Prepare parameter variables
fileName = cell(nFilaments, 1);
filamentIdx = (1:nFilaments)';

l_contour = zeros(nFilaments, 1);
l_ee = zeros(nFilaments, 1);
h_mean = zeros(nFilaments, 1);
h_max = zeros(nFilaments, 1);
h_min = zeros(nFilaments, 1);
periodicity = zeros(nFilaments, 1);
cod_mean = zeros(nFilaments, 1);
r_tip = zeros(nFilaments, 1);
handedness = cell(nFilaments, 1);
dpf = zeros(nFilaments, 1);
symmetry = zeros(nFilaments, 1);
csa_mean = zeros(nFilaments, 1);
csr_mean = zeros(nFilaments, 1);
csr_max = zeros(nFilaments, 1);
csr_min = zeros(nFilaments, 1);
csjz = zeros(nFilaments, 1);
d_rmsd = zeros(nFilaments, 1);
d_corr = zeros(nFilaments, 1);

for aa = 1:nFilaments
    % Image file info
    fileName{aa} = img.dataFile;

    % Directly measured physical parameters
    l_contour(aa) = img.features.filaments(aa).lContour.*xyRes;
    l_ee(aa) = img.features.filaments(aa).rEE(end).*xyRes;
    h_mean(aa) = mean(img.features.filaments(aa).z);

    % Parrameters from 3D reconstruction
    if ~isempty(img.features.filaments(aa).helicalFilament3DModel)
        % From centre line fft analysis
        h_max(aa) = img.features.filaments(aa).helicalFilament3DModel.z_max;
        h_min(aa) = img.features.filaments(aa).helicalFilament3DModel.z_min;
        periodicity(aa) = img.features.filaments(aa).helicalFilament3DModel.periodicity.*xyRes;
        cod_mean(aa) = mean(diff(img.features.filaments(aa).helicalFilament3DModel.peaks_loc)).*xyRes;

        % Tip radius estimate
        r_tip(aa) = img.features.filaments(aa).helicalFilament3DModel.r_tip;

        % Helcal symetry parameters
        symmetry(aa) = img.features.filaments(aa).helicalFilament3DModel.symmetry;
        handedness{aa} = img.features.filaments(aa).helicalFilament3DModel.handedness;
        if strcmpi(handedness{aa}, 'left')
            %dpf(aa) = -1./periodicity(aa);
            dpf(aa) = -1./(cod_mean(aa)*symmetry(aa));
        elseif strcmpi(handedness{aa}, 'right')
            %dpf(aa) = 1/periodicity(aa);
            dpf(aa) = 1./(cod_mean(aa)*symmetry(aa));
        end

        % Average cross-section parameters
        csa_mean(aa) = img.features.filaments(aa).helicalFilament3DModel.xarea_mean;

        xsection = mean(img.features.filaments(aa).helicalFilament3DModel.rho_m, 1);
        csr_mean(aa) = mean(xsection);
        csr_max(aa) = max(xsection);
        csr_min(aa) = min(xsection);
        
        % Second polar moment of area Jz (or Iz) = Ix+Iy
        x_xs = img.features.filaments(aa).helicalFilament3DModel.xm_xsection(2, :)';
        y_xs = img.features.filaments(aa).helicalFilament3DModel.zm_xsection(2, :)';
        iy = (1/12)*sum(...
            (x_xs(1:360).*y_xs(2:361)-x_xs(2:361).*y_xs(1:360)).*...
            (x_xs(1:360).^2+x_xs(1:360).*x_xs(2:361)+x_xs(2:361).^2)...
            );
        ix = (1/12)*sum(...
            (x_xs(1:360).*y_xs(2:361)-x_xs(2:361).*y_xs(1:360)).*...
            (y_xs(1:360).^2+y_xs(1:360).*y_xs(2:361)+y_xs(2:361).^2)...
            );
        csjz(aa) = abs(ix)+abs(iy);

        % Recostruction quality parameters
        d_rmsd(aa) = img.features.filaments(aa).helicalFilament3DModel.d_rmsd(end);
        d_corr(aa) = img.features.filaments(aa).helicalFilament3DModel.d_corr(end);
    end

end

resultTable = table(fileName, filamentIdx, ...
    l_contour, l_ee, h_mean, h_max, h_min, periodicity, cod_mean, ...
    r_tip, handedness, dpf, symmetry, ...
    csa_mean, csr_mean, csr_max, csr_min, csjz, ...
    d_rmsd, d_corr);

resultTable.Properties.VariableDescriptions = resultTable.Properties.VariableNames;

descriptions = {...
    'Image file name', ...
    'Filament index number', ...
    'Contour length / nm', ...
    'End-end distace / nm', ...
    'Mean height / nm', ...
    'Mean peak height / nm', ...
    'Mean trogh height / nm', ...
    'Main FFT periodicity / nm', ...
    'Mean cross over distance / nm', ...
    'Tip radius estimate / nm', ...
    'Twist handedness', ...
    'Directional periodic freqency / nm^-1', ...
    'Estimated approximate cross-section symmetry', ...
    'Mean cross-sectional area / nm^2', ...
    'Mean cross-section radius / nm', ...
    'Max cross-section radius / nm', ...
    'Minimum cross-section radius / nm', ...
    '2nd polar moment of cross-sectional area / nm^4', ...
    'Model-data RMSD / nm', ...
    'Model-data Correlation distance', ...
    };

resultTable.Properties.VariableNames = descriptions;

if ~exist('exportFileName', 'var') || isempty(exportFileName)
    [exportFileName, pathName] = uiputfile(strcat(img.dataFile, '_filaments.csv'));
    exportFileName = strcat(pathName, exportFileName);
end

if exportFileName ~= 0
    writetable(resultTable, exportFileName);
end

% Flip so that output has variable names as headings
resultTable.Properties.VariableNames = resultTable.Properties.VariableDescriptions;
resultTable.Properties.VariableDescriptions = descriptions;

if nargout == 1
    output = resultTable;
end

end
