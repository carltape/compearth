function [M0,Mw,hdur,gamma,epsilon,trM] = CMT2all(M)
%CMT2ALL converts moment tensor to numerous scalar quantities
%
% INPUT
%   M   6 x n moment tensors in up-south-east (CMT) convention
%       M = [Mrr Mtt Mpp Mrt Mrp Mtp]; r=up, theta=south, phi=east
%
% EXAMPLE:
%   M = 1e15*[ 0.1230   -3.4360    3.3120   -0.3960    0.0120   -1.8860]';
%   [M0,Mw,hdur,gamma,epsilon,trM] = CMT2all(M)
%
% NOTE: This does not include quantities associated with the lune, such as
%       those presented in Tape and Tape (2013) -- see TT2013AppA.m
%
% Carl Tape, 01-April-2011
%

% seismic moment (row vector)
im0 = 1;
M0 = CMT2m0(im0,M);

% moment magnitude
imag = 1;
Mw = m02mw(imag,M0);

% half duration
hdur = m02hdur(M0);

% CLVD amount
[gamma, epsilon, trM] = CMT2epsilon(M);
