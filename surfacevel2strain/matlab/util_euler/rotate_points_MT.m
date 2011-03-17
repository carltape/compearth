%
% function [lat_rot2,lon_rot2,R,M] = rotate_points_MT(lat,lon,lat1,lon1,lat2,lon2,espin,M)
% Carl Tape, 07-March-2011
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
%   lat2,lon2   reference starting point, after translation/rotation
%   espin       spin about the translated reference point
%   Min         OPTIONAL: 6 x n matrix of moment tensors
%
% OUTPUT
%   lat_rot2,lon_rot2   final points
%   R                   net rotation matrix
%   Mout        OPTIONAL: 6 x n matrix of moment tensors, after translation/rotation
%
% calls
%   euler_rot_tec.m
%   latlons2pole.m
%   global2local_mat
%   CMTtransform.m
%   xyz2lonlat.m, latlon2xyz.m
%   
% called by xxx
%

function [lat_rot2,lon_rot2,R,Mout] = rotate_points_MT(lat,lon,lat1,lon1,lat2,lon2,espin,Min)

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

% FUTURE: Replace function calls with one analytical expression containing the
%         input variables (lat1,lon1,lat2,lon2,espin).

% note these will look distorted on a Cartesian grid
figure; hold on;
plot(lon,lat,'k.');
plot(lon1,lat1,'kp','markersize',14);
plot(lon_rot,lat_rot,'b.');
plot(lon2,lat2,'bp','markersize',14);
plot(lon_rot2,lat_rot2,'r.');
title('rotate_points.m','interpreter','none')

% check by applying a single rotation to the starting points
% (red circles should be centered on red dots)
Bxyz = latlon2xyz(lat,lon);
R = R2*R1;
Rxyz = R*Bxyz;
[latc,lonc] = xyz2latlon(Rxyz);
plot(lonc,latc,'ro');

if nargin==8
    disp('rotate_points_MT.m: rotating beach balls along with the points');
    
    % transform moment tensors from local to global basis at START POINTS
    Pxyz1 = latlon2xyz(lat,lon);
    Mglobal = global2local_mat(Min,Pxyz1,0);

    % transform moment tensors using R, which is defined for x-y-z basis
    Mglobal_rot = CMTtransform(R,Mglobal);

    % transform moment tensors global to local basis at FINAL POINTS
    Pxyz2 = latlon2xyz(lat_rot2,lon_rot2);
    Mlocal_rot = global2local_mat(Mglobal_rot,Pxyz2,1);

    Mout = Mlocal_rot;
end

%==============================================================
