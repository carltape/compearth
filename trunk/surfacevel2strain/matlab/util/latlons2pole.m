%
% function Pxyz = latlons2pole(lat1,lon1,lat2,lon2)
% Carl Tape, 20-June-2008
%
% This function takes two latitude-longitude points in degress,
% and computes the pole of the great circle containing the two points.
%
% Example (see below).
%
% calls xyz2latlon.m, unit.m
% called by xxx
%

function [Pxyz,Plat,Plon] = latlons2pole(lat1,lon1,lat2,lon2)

Pvecs = latlon2xyz([lat1 lat2], [lon1 lon2], 1);

Pxyz = cross( Pvecs(:,1), Pvecs(:,2) );
[Plat,Plon] = xyz2latlon(Pxyz);
Pxyz = unit(Pxyz);

%norm(Pvecs(:,1))
%norm(Pvecs(:,2))
%norm(Pxyz)

%--------------

if 0==1
    lat1 = 20; lon1 = 0; lat2 = 80; lon2 = -40;
	[Pxyz,Plat,Plon] = latlons2pole(lat1,lon1,lat2,lon2)
    
    lat1 = 0; lon1 = 0; lat2 = 0; lon2 = 80;
    [Pxyz,Plat,Plon] = latlons2pole(lat1,lon1,lat2,lon2)
end

%===========================================================
