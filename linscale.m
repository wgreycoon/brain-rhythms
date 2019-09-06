function [ xscaled ] = linscale( x, minmax )
%% linscale.m:
%
% Summary
% -------------------------------------------------------------------------
% This function rescales an input vector 'x' to span the range defined by
% input vector 'minmax'
%
%  - 'x'         is the original vector
%  - 'minmax'    is a pair of values specifying the new min and max (range)
%                on which to linearly scale the original input vector
%
% Example:
%
%                xscaled = ( [3 5 7 9], [1 2] )
%
% ------------------------------------------------------------------------- 
%
% Author(s):   William Coon   wgreycoon@gmail.com
% Last edited: March 26, 2017 1140hrs
%
% -------------------------------------------------------------------------

xscaled = ( ((minmax(2)-minmax(1))*(x-min(x))) / range(x) ) + minmax(1);


end