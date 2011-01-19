%
% function [lat_rot, lon_rot, R] = euler_rot_tec(lat,lon,evec)
% Carl Tape, 01-Sept-2005
% printed 01-Sept-2005
%
% This function inputs a set of latlon points and an euler rotation vector,
% and outputs the new positions of the points, as well as the rotation
% matrix R.
%
% INPUT:
%       lat  = latitude (deg) vector of input points
%       lon  = longitude (deg) vector of input points
%    evec(1) = latitude (deg) of euler pole
%    evec(2) = longitude (deg) of euler pole
%    evec(3) = rotation angle (deg)
%
% OUTPUT:
%    lat_rot = latitude (deg) vector of rotated points
%    lon_rot = longitude (deg) vector of rotated points
%          R = rotation matrix
%
% See test_euler_rot_tec.m for examples.
%
% See eulerREADME for related programs.
%
% calls euler2rotmat.m, latlon2xyz.m, xyz2latlon.m
% called by test_euler_rot_tec.m
%

function [lat_rot, lon_rot, R] = euler_rot_tec(lat,lon,evec)

deg = 180/pi;

lat = lat(:);
lon = lon(:);

R = euler2rotmat(evec);                     % get the rotation matrix

Pxyz = latlon2xyz(lat,lon,1);               % convert (lat,lon) to (x,y,z)

xyz_rot = R * Pxyz;                         % apply rotation

[lat_rot, lon_rot] = xyz2latlon(xyz_rot);   % convert (x,y,z) to (lat,lon)

%==============================================================
