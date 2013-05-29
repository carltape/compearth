function [M0,Mw,hdur,gamma,epsilon,eCLVD,eDC,trM,eISO,eDEV] = CMT2all(M)
% CMT2ALL converts moment tensor to numerous scalar quantities
%
% INPUT
%   M1,M2   6 x n moment tensors in CMT convention
%           M = [Mrr Mtt Mpp Mrt Mrp Mtp]; r=up, theta=south, phi=east
%
% calls CMT2m0.m, m02mw.m, m02hdur.m, CMT2epsilon.m
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

% CLVD fraction
[gamma, epsilon, eCLVD, eDC, trM, eISO, eDEV] = CMT2epsilon(M);

%-------------------

if 0==1
    clear, clc
    dir1 = '/home/carltape/results/SOURCES/socal_04/CMT_files_post_inverted/';
    filename = [dir1 'CMTSOLUTION_9818433'];
    [date,tshift,hdur,lat,lon,dep,M,eid,elabel] = readCMT(filename,13,0);
    [M0, Mw, hdur, epsilon, eCLVD, eDC, trM] = CMT2all(M)
end

%=====================================================