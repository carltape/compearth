%
% function Mout = global2local_MT(Min,lat,lon,iglob2loc)
% Carl Tape, 15-March-2011
%
% This function converts a set of symmetric matrices between
% a global basis (x-y-z) and a local basis (r-theta-phi).
%
% INPUT
%   Min         6 x n set of symmetric matrices
%   lat,lon     input points on the sphere
%   iglob2loc   =1 for (x,y,z) --> (r,th,ph)
%               =0 for (r,th,ph) --> (x,y,z)
% OUTPUT
%   Mout        6 x n set of symmetric matrices in new basis
%
% calls global2local_rotmat.m, Mdim.m, Mvec2Mmat.m
% called by rotate_points_MT.m
%

function Mout = global2local_MT(Min,lat,lon,iglob2loc)

% check that M is 6 x n
[Min,n1] = Mdim(Min);

% transform to 3 x 3 x n
Min = Mvec2Mmat(Min,1);

% transform M between r-theta-phi and x-y-z
Mout = global2local_mat(Min,lat,lon,iglob2loc);

% transform to 6 x n
Mout = Mvec2Mmat(Mout,0);

%==============================================================