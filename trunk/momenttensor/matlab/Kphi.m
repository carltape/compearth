function Klam = Kphi(phi)
%KPHI given phi angle on the lune, return eigenvalues of crack tensor K
%
% A crack tensor is a moment tensor with normal vector and slip vector
% parallel (opening/tensional crack) or opposite (closing/compressional crack).
% This function  returns the eigenvalues of a crack tensor
% given an input parameter phi. This implements Eq 39 of Tape and Tape (2013),
% "The Classical Model for Moment Tensors". (See Figure 6.)
% 
% INPUT
%   phi     polar angles on the lune (phi = 0 points toward +ISO), degrees
%
% OUTPUT 
%   Klam    3 x n set of normalized eignevalues of the crack tensors
%
% See example below.
%
% Carl Tape, 12/18/2013
%

% ensure row vector
phi = phi(:)';

sinphi = sin(phi*pi/180);
cosphi = cos(phi*pi/180);

% Eq 14 of Tape and Tape (2013)
V = 1/sqrt(6) * [ sqrt(3) 0 -sqrt(3) ;
                  -1 2 -1 ;
                  sqrt(2) sqrt(2) sqrt(2) ];

% Eq 39 of Tape and Tape (2013)
% note: V-transpose = V-inverse for this V
rphi = (4*sinphi.^2 + cosphi.^2).^(-1/2);
ivec = [sqrt(3)*abs(sinphi) ; -sinphi ; cosphi];
Klam =  repmat(rphi,3,1) .* (V'*ivec);

%==========================================================================
% EXAMPLE

if 0==1
    % values from Fig. 6 of TT2013
    phi = [15:30:165 -165:30:-15];
    Klam = Kphi(phi);

    % check
    lam1 = Klam(1,:);
    lam2 = Klam(2,:);
    lam3 = Klam(3,:);
    % check that all points are on the lune boundary (gamma = +/- 30) (Eq 20b, TT2013)
    gamma = 180/pi*atan( (-lam1 + 2*lam2 - lam3) ./ (sqrt(3)*(lam1-lam3)) );
    % recover phi value using Eq 24 of TT2013
    phicheck = 180/pi*atan2( (lam1 - 2*lam2 + lam3) , (sqrt(2)*(lam1 + lam2 + lam3)) );
    % display results
    [phi' Klam' gamma' phicheck']
end

%==========================================================================
