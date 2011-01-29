%
% function [th, ph] = xyz2tp(xyz)
%
% This function takes a matrix of 3-vectors in xyz
% and returns a vector of theta values and phi values.
%
% xyz  = 3 x n matrix of xyz vectors
% th   = col vector of polar angles in radians
% ph   = col vector of azimuthal angles in radians
%
% calls xxx
% called by xxx
%

function [th, ph, rho] = xyz2tp(xyz)

% ensure that xyz is 3 x n
[n, m] = size(xyz); if n~=3, xyz = xyz'; end

[azi, ele, rho] = cart2sph(xyz(1,:), xyz(2,:), xyz(3,:));

% ensure that th and ph are column vectors
th = pi/2 - ele(:);
ph = azi(:);

%===========================================================
