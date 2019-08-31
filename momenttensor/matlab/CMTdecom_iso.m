function [Miso,Mdev,trM] = CMTdecom_iso(M)
%CMTDECOM_ISO decompose moment tensor into deviatoric and isotropic parts
%
% INPUT
%   M     6 x n moment tensors in CMT convention
%         M = [Mrr Mtt Mpp Mrt Mrp Mtp]; r=up, theta=south, phi=east
%
% OUTPUT
%   Miso  6 x n moment tensors representing isotropic (spherical) part
%   Mdev  6 x n moment tensors representing deviatoric (nonisotropic) part
%   trM   1 x n vector of trace = Mrr+Mtt+Mpp
%
% NOTE: This is sometimes needed for numerical stability in plotting
% moment tensors with GMT's psmeca.
%
% calls Mdim.m
%
% Carl Tape, 2010-11-01
%

% make sure M is 6 x n
[M,n] = Mdim(M);

%Mrr = M(1,:); Mtt = M(2,:); Mpp = M(3,:);
%Mrt = M(4,:); Mrp = M(5,:); Mtp = M(6,:);

% trace (row vector)
%trM = zeros(1,n);
trM = M(1,:) + M(2,:) + M(3,:);

% decompose into spherical and deviatoric parts
Miso = zeros(6,n);
Miso(1,:) = 1/3*trM';
Miso(2,:) = 1/3*trM';
Miso(3,:) = 1/3*trM';
Mdev = zeros(6,n);
Mdev(1,:) = M(1,:) - Miso(1,:);
Mdev(2,:) = M(2,:) - Miso(2,:);
Mdev(3,:) = M(3,:) - Miso(3,:);
Mdev(4,:) = M(4,:);
Mdev(5,:) = M(5,:);
Mdev(6,:) = M(6,:);

%==========================================================================