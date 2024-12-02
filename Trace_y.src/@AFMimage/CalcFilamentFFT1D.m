function [ff, psd, p1] = CalcFilamentFFT1D(img, filament)

%
% DESCRIPTION
% – Simple function to calculate the central line spatial freqencies with 
% FFT. Function used to calculate the periodicity of helical fibrils from 
% filaments traced from AFM image data.
% – Part of Trace_y
%
% USAGE
% Standard method usage
% >> [ff, psd] = img.CalcFilamentFFT1D(filament)
%
% INPUTS
% img  –  AFMimage object. Implicit with dot notation
% filament  –  The index number of filament to be analysed
%
% OUTPUTS
% ff  –  Freqencies
% psd  –  Power spectral density
% p1  –  Single sided amplitude
%
% DEPENDENCIES
% – Method for Trace_y's @AFMimage/ object.
% – Used by several @AFMimage/ methods
%
% AUTHORS
% Liisa Lutter, Wei-Feng Xue
%
% HISTORY
% 2018.10  –  Initial draft by LL
% 2019.01  –  WFX edited and incorporated function into @AFMimage class 
% methods
%



if ~exist('filament', 'var')
    %filament = 0;
    filament = 1;
end


% Importing l, z values of the fibril.

if filament == 0
    l = img.features.lastFilament.l;
    z = img.features.lastFilament.z;
else
    l = img.features.filaments(filament).l;
    z = img.features.filaments(filament).z;
end

% resampling data with spline interpolation
% current number of points 
nz = length(l);
% number of points to sample (sampling length)
nz = 2.^(nextpow2(nz));
% sampleing freq
freq = 1/(l(end)/(nz-1));
% new evenly spaced points
ll = linspace(0, l(end), nz)';

% re-sample
zz = spline(l, z, ll);

% Detrending for amplitude values
zdetrend = detrend(zz, 'constant');
% FFT of detrended zvalues
zfft = fft(zdetrend);



% Making amplitude single sided
% 2-sided
p2 = abs(zfft/nz);
% 1-sided
p1 = p2(1:nz/2+1);
% correct the amplitudes (i.e. x2)
p1(2:end-1) = 2*p1(2:end-1);
ff = (freq*(0:(nz/2))/nz)';

% Power spectral density
psd = (1/nz)*abs(zfft).^2;
psd = psd(1:nz/2+1);
psd(2:end-1) = 2*psd(2:end-1);



end
