function [c] = CmapAFMstd(n)

%
% DESCRIPTION
% – The default colour scheme for AFM images
% – Part of Trace_y
%
% USAGE
% Standard method usage
% >> cmap = CmapAFMstd;
%
% INPUTS
% n  –  Optional number of colours
%
% OUTPUTS
% c  –  AFM image standard colour map.
%
% DEPENDENCIES
% – Used for displaying images in Trace_y
%
% AUTHORS
% Wei-Feng Xue
% 
% HISTORY
% 2008.04  –  Code separated out and properly defined once and for all
%



if nargin == 0
    % Standard 256 colours
    n = 2^8;
end

n1 = floor(n/2);
n2 = n-n1;

c1 = zeros(n1, 3);
c1(:, 1) = 0:(0.07/(n1-1)):0.07;
c1(:, 2) = 1;
c1(:, 3) = 0.1:(0.4/(n1-1)):0.5;

c2 = zeros(n2, 3);
c2(:, 1) = 0.07:(0.10/(n2-1)):0.17;
c2(:, 2) = fliplr(0.1:((0.9/(n2-1))):1);
c2(:, 3) = 0.5:(0.5/(n2-1)):1;

c = hsv2rgb([c1; c2]);

end