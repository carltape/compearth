%
% read_gps_japan.m
% CARL TAPE, 08-March-2009
% printed xxx
%
% This file reads in a GPS field from Japan.
%
% calls xxx
% called by japan_gps_dat.m
%

function [lon,lat,ve,vn,vu,se,sn,su,ren,reu,rnu,start_date,finish_date,name] = read_gps_japan

dir0 = '/home/carltape/gmt/gps_data/ASIA/japan/';

file = 'OUTFILE';        % original datafile provided by Takeo Ito
filename = [dir0 file];

% site name  lat. lon. N (m/y) E(m/y) U(m/y) EN(m/y) EE(m/y) EU(m/y)
[name,lat,lon,vn,ve,vu,sn,se,su] = textread(filename,'%s%f%f%f%f%f%f%f%f');
nobs = length(name);

% no covariances included
ren = zeros(nobs,1); reu = zeros(nobs,1); rnu = zeros(nobs,1);

% start dates and end dates
% "Observation period is from 1/4/1996  4 to 1/4/2000 that about 4 years."
start_date  = datenum(1996,1,4)*ones(nobs,1);
finish_date = datenum(2000,1,4)*ones(nobs,1);

% convert to mm/yr
vn = vn*1e3; ve = ve*1e3; vu = vu*1e3;
sn = sn*1e3; se = se*1e3; su = su*1e3;

%======================================================
