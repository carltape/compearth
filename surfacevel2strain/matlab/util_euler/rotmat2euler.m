%
% function R = rotmat2euler(evec)
% Carl Tape, 01-Nov-2005
%
% Given a finite rotation matrix, this program returns the three components
% that describe the rotation (lat,lon,omega).
% 
% INPUT: matrix for finite rotation
%
% OUTPUT: euler pole vector
%    evec(1) = latitude (deg) of euler pole,   lat  = [-90,90]
%    evec(2) = longitude (deg) of euler pole,  lon  = [-180,180]
%    evec(3) = rotation angle (deg),           omeg = [0,180]
%
% Reverse program of euler2rotmat.m 
%
% See Cox & Hart (1986), p. 234
%
% calls arctan.m
% called by test_euler_rot_tec.m, c161C.m
%

function evec = rotmat2euler(R)

deg = 180/pi;

if R == eye(3)
    disp(' R is identity matrix, so no rotation is applied.');
    evec = [0 0 0]';
    
else

    R11 = R(1,1);
    R12 = R(1,2);
    R13 = R(1,3);
    R21 = R(2,1);
    R22 = R(2,2);
    R23 = R(2,3);
    R31 = R(3,1);
    R32 = R(3,2);
    R33 = R(3,3);

    % factor used in computing latitude and rotation angle
    fac = sqrt( (R32 - R23)^2 + (R13 - R31)^2 + (R21 - R12)^2 );

    elat = deg * asin( (R21 - R12) / fac );

    % BE CAREFUL: determining the longitude is tricky (see test_arctan.m)
    atop = R13 - R31;
    abot = R32 - R23;
    if atop >=0         % want lon to be [0,180]
        elon = deg * arctan( atop / abot );
    elseif atop <= 0    % want lon to be [-180,0]
        elon = -180 + deg * arctan( atop / abot );
    end

    % arctan.m ensures the output is [0,180]
    omeg = deg * arctan( fac / (R11 + R22 + R33 - 1) );

    % output euler vector
    evec = [elat elon omeg]';

end

%==============================================================
