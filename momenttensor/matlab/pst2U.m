function U = pst2U(phi,sigma,theta)
%PST2U phi-sigma-theta to rotation matrix U 
%
% INPUT
%   phi     n x 1 vector of azimuth spherical coordinate angle, degrees
%   sigma   n x 1 vector of longitudinal angle, degrees
%   theta   n x 1 vector of polar spherical angle, degrees
%
% OUTPUT
%   U       3 x 3 x n set of orientation matrices
%
% See examples in U2pst.m
%
% See WTape and CTape (2017) "Volume in moment tensor space in terms of distance"
%
% Carl Tape, 2017/11
%

% check input
np = length(phi);
ns = length(sigma);
nt = length(theta);
if length(unique([np ns nt])) > 2
    whos phi sigma theta
    error('input must have same dimension');
end
n = np;

% TT2017
U = NaN(3,3,n);
for ii=1:n
    % Tape and Tape (2017), Eq 40a
    Zp = rotmat(phi(ii),3);
    Yt = rotmat(theta(ii),2);
    Zs = rotmat(sigma(ii),3);
    U(:,:,ii) = Zp * Yt * Zs;
end

%==========================================================================
