function R = rotmat(xdeg,ixyz)
%ROTMAT given index of axis (1,2,3) and angle, return rotation matrix
%
% INPUT
%   x       input angle, degrees
%   ixyz    index designating axis of rotation (=1,2,3)
% OUTPUT
%   R       rotation matrix
%
% Carl Tape, 01-April-2011
%

n = length(xdeg);
cosx = cos(xdeg * pi/180);
sinx = sin(xdeg * pi/180);

if n==1
    if ixyz==1
        R = [1 0 0 ; 0 cosx -sinx ; 0 sinx cosx ];
    elseif ixyz==2
        R = [cosx 0 sinx ; 0 1 0 ; -sinx 0 cosx];
    elseif ixyz==3
        R = [cosx -sinx 0 ; sinx cosx 0 ; 0 0 1];
    else
       error('ixyz = 1,2,3 only'); 
    end
    
else
    R = zeros(3,3,n);
    if ixyz==1
        R(1,1,:) = 1; R(1,2,:) = 0; R(1,3,:) = 0;
        R(2,1,:) = 0; R(2,2,:) = cosx(:)'; R(2,3,:) = -sinx(:)';
        R(3,1,:) = 0; R(3,2,:) = sinx(:)'; R(3,3,:) = cosx(:)';
    elseif ixyz==2
        R(1,1,:) = cosx(:)'; R(1,2,:) = 0; R(1,3,:) = sinx(:)';
        R(2,1,:) = 0; R(2,2,:) = 1; R(2,3,:) = 0;
        R(3,1,:) = -sinx(:)'; R(3,2,:) = 0; R(3,3,:) = cosx(:)';
    elseif ixyz==3
        R(1,1,:) = cosx(:)'; R(1,2,:) = -sinx(:)'; R(1,3,:) = 0;
        R(2,1,:) = sinx(:)'; R(2,2,:) = cosx(:)'; R(2,3,:) = 0;
        R(3,1,:) = 0; R(3,2,:) = 0; R(3,3,:) = 1;
    else
       error('ixyz = 1,2,3 only'); 
    end   
end

%==========================================================================
