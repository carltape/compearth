%
% user_path.m
%
% File to set Matlab paths into subdirectories.
%
% Probably it would be cleaner to all user_path.n only ONCE, since this
% will keep adding the same directories into the path each time you run
% surfacevel2strain.m.
%

% THIS ASSUMES THAT YOU HAVE compearth IN YOUR HOME DIRECTORY
% IF IT IS NOT, THEN CHANGE IT.
dir_home = getenv('HOME');
dir_compearth = strcat(dir_home,'/compearth/');
if ~exist(dir_compearth,'dir')
    dir_compearth
    error('user_path.m: compearth directory does not exist: ');
end

bdir = strcat(dir_compearth,'surfacevel2strain/');
bdir_matlab = strcat(bdir,'matlab/');

addpath(strcat(bdir_matlab));
addpath(strcat(bdir_matlab,'util_basic'));
addpath(strcat(bdir_matlab,'util_est'));
addpath(strcat(bdir_matlab,'util_euler'));
