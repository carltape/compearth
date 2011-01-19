%
% function 
% CARL TAPE, 22-Sept-2006
% printed xxxx
%
% This function returns three index vectors that describe how to access the
% n*ndim datapoint.
%
% Convention is taken as:
%    3D: r-th-phi, that is, up-south-east
%    2D: th-phi, that is, south-east
%
% calls xxx
% called by test_platemodel2strain.m
%

function [inds1, inds2, inds3] = subindexing(n,ndim,opts)

iorder = opts{1};

inds1 = []; inds2 = []; inds3 = [];

iall = [1:n*ndim]';

if iorder==1            % vAx, vBx, ..., vAy, vBy, ..., vAz, vBz, ...
    if ndim==1
        inds1 = [1:n];
        
    elseif ndim == 2
        inds1 = [1:n];
        inds2 = [n+1:n*ndim];
        
    elseif ndim == 3   
        inds1 = [1:n];
        inds2 = [n+1:n*2];
        inds3 = [n*2+1:n*ndim];
    end
    
elseif iorder==2        % vAx, vAy, vAz, vBx, vBy, vBz, ...
    if ndim==1
        inds1 = [1:n];

    elseif ndim == 2
        inds2 = find(mod(iall,2)==1);
        inds3 = find(mod(iall,2)==0);

    elseif ndim == 3   
        inds1 = find(mod(iall,3)==1);
        inds2 = find(mod(iall,3)==2);
        inds3 = find(mod(iall,3)==0);
    end
end

if length(inds1) + length(inds2) + length(inds3) ~= n*ndim
    error('dimension mis-match');
end

inds1 = inds1(:);
inds2 = inds2(:);
inds3 = inds3(:);

%=====================================================================
