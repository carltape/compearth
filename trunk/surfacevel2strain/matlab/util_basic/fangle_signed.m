%
% function stheta = fangle_signed(va,vb,vnor)
% Carl Tape, 12-April-2011
%
% This function returns the signed angle of rotation between two vectors.
%
% INPUT
%   va      3 x n set of initial vectors
%   vb      3 x n set of rotated vectors
%   vnor    3 x n set of "normal of rotation" vectors
%
% OUTPUT
%   theta   n x 1 set of rotation angles
% 
%
% calls fangle.m
% called by xxx
%

function stheta = fangle_signed(va,vb,vnor)

% get rotation angle (always positive)
theta = fangle(va,vb);
n = length(theta);

EPSVAL = 0;
stheta = theta;     % initialize to rotation angle
for ii=1:n
   if abs(theta - 180) <= EPSVAL
       stheta(ii) = 180;
   else
       Dmat = [va(:,ii) vb(:,ii) vnor(:,ii)];
       if det(Dmat) < 0
           stheta(ii) = -theta(ii);
       end
   end
end

%--------------
% EXAMPLES

if 0==1
    va = [1 0 0]'; vb = [0 0 3]'; vnor = -[0 2 0]'; stheta = signedfangle(va,vb,vnor)
    va = [1 0 0]'; vb = [0 0 3]'; vnor = [0 2 0]'; stheta = signedfangle(va,vb,vnor)
end

%===========================================================
