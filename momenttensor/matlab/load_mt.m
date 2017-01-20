%
% load_mt.m
%
% Set paths to additional directories in compearth
% (These are needed for utilities such as rotation scripts.)
%
% load_mt.m can be called from a matlab startup.m script,
% so that you never need to execute it.

% path to the compearth head directory 
% it is safest to specify an absolute path to your compearth directory
%compearth = strcat(getenv('HOME'),'/REPOSITORIES/compearth/');
% this relative path will only work if you run Matlab from compearth/momenttensor/matlab
compearth = '../../';

addpath(strcat(compearth,'surfacevel2strain/matlab/util_basic'));
addpath(strcat(compearth,'surfacevel2strain/matlab/util_euler'));
addpath(strcat(compearth,'momenttensor/matlab'));
