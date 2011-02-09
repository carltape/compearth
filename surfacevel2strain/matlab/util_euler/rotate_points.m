%
% function [lat_rot2,lon_rot2,R] = rotate_points(lat,lon,lat1,lon1,lat2,lon2,espin)
% Carl Tape, 07-Feb-2011
%
% This function applies a finite rotation to a set of input points without
% directly specifying an euler vector. For example, if we want to
% superimpose an outline of California next to Alaska -- the rotated points
% will preserve the shape, whereas a simple translation in lon-lat
% coordinates will not.
% 
% INPUT
%   lat,lon     points to rotate
%   lat1,lon1   reference starting point
%   lat2,lon2   reference starting point, after translation
%   espin       spin about the translated reference point
%
% OUTPUT
%   lat_rot2,lon_rot2   final points
%   R                   net rotation matrix
%
% calls euler_rot_tec.m, latlons2pole.m, xyz2lonlat.m
% called by xxx
%

function [lat_rot2,lon_rot2,R] = rotate_points(lat,lon,lat1,lon1,lat2,lon2,espin)

deg = 180/pi;

lon = lon(:);
lat = lat(:);

% compute the pole to translate body
[Pxyz,Plat,Plon] = latlons2pole(lat1,lon1,lat2,lon2);
Pdist = distance(lat1,lon1,lat2,lon2);

% apply translation
evec = [Plat Plon Pdist];
[lat_rot, lon_rot, R1] = euler_rot_tec(lat,lon,evec);

% apply rotation about final point
evec2 = [lat2 lon2 espin];
[lat_rot2, lon_rot2, R2] = euler_rot_tec(lat_rot,lon_rot,evec2);

% note these will look distorted on a Cartesian grid
figure; hold on;
plot(lon,lat,'k.');
plot(lon1,lat1,'kp','markersize',14);
plot(lon_rot,lat_rot,'b.');
plot(lon2,lat2,'bp','markersize',14);
plot(lon_rot2,lat_rot2,'r.');
title('rotate_points.m','interpreter','none')

% check by applying a single rotation to the starting points
Bxyz = latlon2xyz(lat,lon);
R = R2*R1;
Rxyz = R*Bxyz;
[latc,lonc] = xyz2latlon(Rxyz);
plot(lonc,latc,'ro');

%==============================================================
