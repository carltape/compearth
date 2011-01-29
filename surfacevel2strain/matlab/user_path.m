%
% user_path.m
%
% File to set Matlab paths into subdirectories.
%
%

% USER change
bdir = '/home/carltape/compearth/surfacevel2strain/';
bdir1 = [bdir 'matlab/'];

if ~exist(bdir,'dir')
    bdir
    error('user_path.m: dir does not exist: ');
end

% add path to additional matlab scripts
path(path,bdir1);
path(path,[bdir1 'util_basic']);
path(path,[bdir1 'util_est']);
path(path,[bdir1 'util_euler']);

%------------------------------------------------------------------------
