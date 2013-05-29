function Ufour = Ufour(U)
%U2FOUR find four equivalent principal axis triples
%
% See TapeTape2013 "Angle between principal axis triples"
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

%==========================================================================
