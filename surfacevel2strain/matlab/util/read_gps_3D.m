%
% function read_gps_3D.m
% Carl Tape, 02-Aug-2007
%
% This program reads in 3D velocity files that were generated using write_gps_3D.m
%
% calls xxx
% called by xxx
%

function [lon,lat,ve,vn,vu,se,sn,su,ren,reu,rnu,start_date,finish_date,name] = read_gps_3D(filename)

if ~exist(filename), error([filename ' does not exist']); end    

[lon,lat,ve,vn,vu,se,sn,su,ren,reu,rnu,start_date,finish_date,name] ...
    = textread(filename,'%f%f%f%f%f%f%f%f%f%f%f%f%f%s','headerlines',1);
    
%=======================================================================