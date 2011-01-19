%
% function mat = tp2xyz(th, ph, r)
% CARL TAPE, 01-Sept-2005
% printed xxxx
%
% This function takes the theta and phi of a set of vectors
% on a sphere with radius r and and returns the xyz coordinates
% as a 3 x n matrix.
%
% th    = col vector of polar angles in radians
% ph    = col vector of azimuthal angles in radians
% r     = radial value for the sphere
%
% calls xxx
% called by cities.m, geteuler.m, trigridN.m, spline_wang.m
%

function mat = tp2xyz(th, ph, r)

% column vectors
th = th(:);
ph = ph(:);

azi = ph;
ele = pi/2 - th;

if length(r) == 1, r = r*ones(length(th),1); end

[xx, yy, zz] = sph2cart(azi, ele, r );
mat = [xx yy zz]';

%===========================================================
