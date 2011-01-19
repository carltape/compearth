%
% function Pxyz = latlons2pole(lat1,lon1,lat2,lon2)
% Carl Tape, 20-June-2008
% printed xxxx
%
% This function takes two latitude-longitude points in degress,
% and computes the pole of rotation.
%
% Example (see below).
%
% calls xxx
% called by xxx
%

function [Pxyz,Plat,Plon] = latlons2pole(lat1,lon1,lat2,lon2)

deg = 180/pi;

Pvecs = latlon2xyz([lat1 lat2],[lon1 lon2],1);

Pxyz = cross( Pvecs(:,1), Pvecs(:,2) );
[Plat,Plon] = xyz2latlon(Pxyz);

%--------------

if 0==1
    lat1 = 20; lon1 = 0; lat2 = 80; lon2 = -40;
	[Pxyz,Plat,Plon] = latlons2pole(lat1,lon1,lat2,lon2)
    
    lat1 = 0; lon1 = 0; lat2 = 0; lon2 = 80;
    [Pxyz,Plat,Plon] = latlons2pole(lat1,lon1,lat2,lon2)
end

%===========================================================
