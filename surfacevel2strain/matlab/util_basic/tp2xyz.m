%
% function mat = tp2xyz(th, ph, r)
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
% called by xxx
%

function mat = tp2xyz(th, ph, r)

% column vectors
th = th(:);
ph = ph(:);

azi = ph;
ele = pi/2 - th;

if length(r) == 1
    r = r*ones(length(th),1);
    disp('tp2xyz.m: uniform radial value');
end

% azi, ele, r are vectors
[xx, yy, zz] = sph2cart(azi, ele, r );
mat = [xx yy zz]';

%===========================================================
