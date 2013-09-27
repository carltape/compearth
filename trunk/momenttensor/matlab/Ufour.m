function Ufour = Ufour(U)
%U2FOUR find four equivalent principal axis triples
%
% INPUT
%   U       3 x 3 (assumed) rotation matrix and output 
%
% OUTPUT
%   Ufour   3 x 3 x 4 set of four equivalent principal axes
%
% See TapeTape2012 "Angle between principal axis triples"
% 
% Carl Tape, 8/12/2012
%

Xpi = diag([1 -1 -1]);
Ypi = diag([-1 1 -1]);
Zpi = diag([-1 -1 1]);

Ufour = NaN(3,3,4);
Ufour(:,:,1) = U;
Ufour(:,:,2) = U * Xpi;
Ufour(:,:,3) = U * Ypi;
Ufour(:,:,4) = U * Zpi;
