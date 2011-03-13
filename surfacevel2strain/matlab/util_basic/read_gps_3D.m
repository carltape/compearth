%
% function [lon,lat,ve,vn,vu,se,sn,su,ren,reu,rnu,start_date,finish_date,name] = read_gps_3D(filename)
%
% This program reads in 3D velocity files that were generated using write_gps_3D.m
%
% INPUT
%   filename   path to input file
%
% OUTPUT
%   lon,lat         location of observation (degrees)
%   ve,vn,vu        velocity field (east, north, up), MM/YR
%   se,sn,su        standard deviation for each component, MM/YR
%   ren,reu,rnu     covariances
%   start_date      start date (serial date)
%   finish_date     finish date (serial date)
%   name            station name
%
% calls xxx
% called by get_gps_dataset.m
%

function [lon,lat,ve,vn,vu,se,sn,su,ren,reu,rnu,start_date,finish_date,name] = read_gps_3D(filename)

if ~exist(filename,'file')
    error([filename ' does not exist']);
else
    disp(['read_gps_3D.m: ' filename]);
    [lon,lat,ve,vn,vu,se,sn,su,ren,reu,rnu,start_date,finish_date,name] ...
        = textread(filename,'%f%f%f%f%f%f%f%f%f%f%f%f%f%s','headerlines',1);
end    
    
%=======================================================================