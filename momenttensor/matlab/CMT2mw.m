function [Mw,M0] = CMT2mw(M)
%CMT2MW convert from moment tensor to moment magnitude
%
% INPUT
%   M       6 x n set of moment tensors, M = [Mrr Mtt Mpp Mrt Mrp Mtp]
%
% The units of M0 are the same as the elements of Mij, which should be
% Newton-meter (N-m). See Latex notes for details.
%

% seismic moemnt
M0 = CMT2m0(1,M);   % use default formula

% moment magnitude
Mw = m02mw(1,M0);   % use default formula

%==========================================================================
