%
% load_mt.m
%
% Set paths to additional directories in compearth
% (These are needed for utilities such as rotation scripts.)

% path to the compearth head directory 
%compearth = strcat(getenv('HOME'),'/compearth/');  % absolute path (safest)
compearth = '../../';                               % relative path

addpath(strcat(compearth,'surfacevel2strain/matlab/util_basic'));
addpath(strcat(compearth,'surfacevel2strain/matlab/util_euler'));
addpath(strcat(compearth,'momenttensor/matlab'));
