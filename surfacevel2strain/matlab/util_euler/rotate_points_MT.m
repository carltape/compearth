function [latr,lonr,R,Mout] = rotate_points_MT(lat,lon,lat1,lon1,lat2,lon2,espin,arg8)
%ROTATE_POINTS_MT apply finite rotation to a set of points AND moment tensors
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
%   espin       spin about the translated reference point, degrees
%   arg8        OPTIONAL 8th argument
%                 OPTION 1: Min, 6 x n matrix of moment tensors
%                 OPTION 2: otag, string for writing output file
%                 note: at present you can't rotate MTs AND write output files
%
% OUTPUT
%   latr,lonr   final points
%   R           net rotation matrix
%   Mout        OPTIONAL: 6 x n matrix of moment tensors, after translation/rotation
%
% This program can easily be generalized to allow for non-symmetric 3 x 3
% matrices M (use global2local_mat.m and transform_mat.m).
%
% calls
%   euler_rot_tec.m
%   latlons2pole.m
%   global2local_MT.m
%   transform_MT.m
%   xyz2lonlat.m
%
% Carl Tape, 07-March-2011
%

ilon360 = 1;

if ilon360==1
    lon = wrapTo360(lon);
    lon1 = wrapTo360(lon1);
    lon2 = wrapTo360(lon2);
else
    lon = wrapTo180(lon);
    lon1 = wrapTo180(lon1);
    lon2 = wrapTo180(lon2);
end

lon = lon(:);
lat = lat(:);
n = length(lon);

% compute the pole to translate body
[Pxyz,Plat,Plon] = latlons2pole(lat1,lon1,lat2,lon2);
Pdist = sdistance(lat1,lon1,lat2,lon2);
%Pdist = distance(lat1,lon1,lat2,lon2);

% apply translation
evec = [Plat Plon Pdist];
[lat_rot, lon_rot, R1] = euler_rot_tec(lat,lon,evec);

% apply rotation about final point
evec2 = [lat2 lon2 espin];
[latr, lonr, R2] = euler_rot_tec(lat_rot,lon_rot,evec2);

if ilon360==1
   lon_rot = wrapTo360(lon_rot);
   lonr = wrapTo360(lonr);
else
   lon_rot = wrapTo180(lon_rot);
   lonr = wrapTo180(lonr);
end

% FUTURE: Replace function calls with one analytical expression containing the
%         input variables (lat1,lon1,lat2,lon2,espin). This will reduce the
%         possibility of numerical round-off errors.

% note these will look distorted on a Cartesian grid
figure; hold on;
plot(lon,lat,'k.');
plot(lon1,lat1,'kp','markersize',14);
plot(lon_rot,lat_rot,'b.');
plot(lon2,lat2,'bp','markersize',14);
plot(lonr,latr,'r.');
title(sprintf('rotate_points.m (%i points)',n),'interpreter','none');

% check by applying a single rotation to the starting points
% (red circles should be centered on red dots)
Bxyz = latlon2xyz(lat,lon);
R = R2*R1;
Rxyz = R*Bxyz;
[latc,lonc] = xyz2latlon(Rxyz);
if ilon360==1, lonc=wrapTo360(lonc); else lonc=wrapTo180(lonc); end
plot(lonc,latc,'ro');

if nargin==8
    if ischar(arg8)
        otag = arg8;
        disp('rotate_points_MT.m: write rotated points to file');
        write_xy_points(otag,lonr,latr);
        
    else
        Min = arg8;
        disp('rotate_points_MT.m: rotating beach balls along with the points');

        % transform moment tensors from local to global basis at START POINTS
        Mglobal = global2local_MT(Min,lat,lon,0);

        % transform moment tensors using R, which is defined for x-y-z basis
        Mglobal_rot = transform_MT(R,Mglobal);

        % transform moment tensors global to local basis at FINAL POINTS
        Mlocal_rot = global2local_MT(Mglobal_rot,latr,lonr,1);

        Mout = Mlocal_rot;
    end
end

%==========================================================================
