%
% function mat = latlon2xyz(lat,lon,r)
% CARL TAPE, 30-Aug-2005
% printed xxxx
%
% This function takes the latitude and longitude of a point on a sphere
% with radius r and outputs the (x,y,z) coordinates.
%
% lat  = col vector of latitudes (deg)
% lon  = col vector of longitudes (deg)
% r    = radial value for the sphere
%
% calls tp2xyz.m
% called by xxx
%

function mat = latlon2xyz(lat, lon, r)

deg = 180/pi;

lat = lat(:);
lon = lon(:);
th = (90-lat)/deg;
ph = lon/deg;

mat = tp2xyz(th, ph, r);

%===========================================================
