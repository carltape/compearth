%
% function [dvec,azvec,lon1,lat1,lon2,lat2] = lonlat2distaz(lon,lat)
% CARL TAPE, 23-Sept-2007
% printed xxx
%
% This function reads in a discretized lat-lon curve on the sphere, and
% ouputs a discretized vector of the azimuths and distances for each
% segment.
% 
% calls xxx
% called by surfacevel2strain.m
%

function [dvec,azvec,lon1,lat1,lon2,lat2] = lonlat2distaz(lon,lat)

n = length(lon);
nout = n - 1;

lon1 = lon(1:end-1);
lat1 = lat(1:end-1);
lon2 = lon(2:end);
lat2 = lat(2:end);
[dvec, azvec] = distance(lat1,lon1,lat2,lon2);

%-------------------------

if 0==1
    load('safdata');
    lon = lonsaf; lat = latsaf;
    figure; plot(lon,lat,'.');
    
    [dvec, azvec, lon1,lat1,lon2,lat2] = lonlat2distaz(lon,lat);
    disp('segment     point 1               point 2        -->    distance   azimuth');
    for ii=1:length(lon)-1
       disp(sprintf('%4i (%9.4f, %8.4f) (%9.4f, %8.4f) --> %10.4f %10.4f',...
           ii,lon1(ii),lat1(ii),lon2(ii),lat2(ii),dvec(ii),azvec(ii)));
    end
end

%===================================================================
