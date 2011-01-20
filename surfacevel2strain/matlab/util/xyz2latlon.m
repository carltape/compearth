%
% function [lat,lon] = xyz2latlon(xyz)
%
% xyz  = 3 x n matrix of xyz vectors
% lat  = col vector of latitudes (deg)
% lon  = col vector of longitudes (deg)
%
% calls xyz2tp.m
% called by xxx
%

function [lat,lon, rho] = xyz2latlon(xyz)

[th, ph, rho] = xyz2tp(xyz);

deg = 180/pi;
lat = (pi/2 - th)*deg;
lon = ph*deg;

%===========================================================
