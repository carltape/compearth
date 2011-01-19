%
% function evec_out = euler_mapping(evec)
% Carl Tape, 23-April-2007
% printed xxx
%
% This reads in an euler vector in "non-standard" form and outputs an euler
% vector in "standard" form.
%
% For example, an euler vector specified by
%    (latitude = -120, longitude = 465, and omega = -700) is
% retured as
%    (latitude = -60, longitude = -75, and omega = 700) 
% euler_mapping([-120 465 340])
%
% calls latlon2xyz.m, xyz2latlon.m
% called by euler2rotmat.m
%

function evec_out = euler_mapping(evec)

deg = 180/pi;

% input euler vector
elat = evec(1);
elon = evec(2);
omega = evec(3);

%--------------------------------------
% modification of rotation vector
% (this does not change the rotation matrix -- check it yourself)

% unit vector for rotation pole
exyz = latlon2xyz(elat,elon,1);

% updated rotation pole
[elat,elon] = xyz2latlon(exyz);

% % shift longitude to (-180,180] and latitude to [-90,90]
% elon = lonshift(elon,[1 1]);
% elat = latshift(elat,1);

% if omega is negative, then make it positive and take the antipode
if omega < 0
   [elat,elon] = antipode(elat,elon);
   omega = -omega;
end

% map the POSITIVE rotation angle onto [0,360)
omega = mod(omega,360);

% output euler vector
evec_out = [elat elon omega]';

%==============================================================
