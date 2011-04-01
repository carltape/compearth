%
% function R = rotmat(x,ixyz)
% Carl Tape, 01-April-2011
%
% This returns a rotation matrix, given an angle and axis.
%
% INPUT
%   x       input angle, degrees
%   ixyz    index designating axis of rotation (=1,2,3)
% OUTPUT
%   R       rotation matrix
%
% calls xxx
% called xxx
%

function R = rotmat(x,ixyz)

cosx = cos(x * pi/180);
sinx = sin(x * pi/180);

if ixyz==1
    R = [1 0 0 ; 0 cosx -sinx ; 0 sinx cosx ];
elseif ixyz==2
    R = [cosx 0 sinx ; 0 1 0 ; -sinx 0 cosx];
elseif ixyz==3
    R = [cosx -sinx 0 ; sinx cosx 0 ; 0 0 1];
else
   error('ixyz = 1,2,3 only'); 
end

%==============================================================
