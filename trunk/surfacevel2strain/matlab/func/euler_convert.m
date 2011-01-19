%
% function outvecs = euler_convert(invecs,opts)
% Carl Tape, 21-Oct-2005
% printed 21-Oct-2005
%
% This function converts between two types of Euler pole conventions:
% (1) (wx,wy,wz) --> (lat-deg,lon-deg,omega-deg/Myr)
% (2) (lat-deg,lon-deg,omega-deg/Myr) --> (wx,wy,wz)
%
% See example below (and test_euler_rot_tec.m, plate_model.m).
%
% See eulerREADME for related programs.
%
% calls euler_rot_tec.m, xyz2latlon.m
% called by xxx
%

function outvecs = euler_convert(invecs,opts)

deg = 180/pi;

% Specify the direction you want to convert.
% ctype = 1 : (wx,wy,wz) --> (lat-deg,lon-deg,omega-deg)
% ctype = 0 : (lat-deg,lon-deg,omega-deg) --> (wx,wy,wz)
ctype = opts(1);

% ensure that invecs is 3 x n
[n,m] = size(invecs); if n~=3, invecs = invecs'; end

if ctype == 1   % (wx,wy,wz) --> (lat-deg,lon-deg,omega-deg)
    
    wx = invecs(1,:);
    wy = invecs(2,:);
    wz = invecs(3,:);
    
    % vector lengths in deg/Myr
    omg = sqrt( wx.^2 + wy.^2 + wz.^2 );
    [lat,lon] = xyz2latlon([wx ; wy ; wz]);
    
    outvecs = [lat(:) lon(:) omg(:)]';  % 3 x n
    
%------------------------------------------------------------------
else            % (lat-deg,lon-deg,omega-deg) --> (wx,wy,wz)
    
    lat = invecs(1,:);
    lon = invecs(2,:);
    omg = invecs(3,:);      % deg/Myr
    
    % convert negative rates to positive rates using antipode
    ineg = find(omg < 0);
    [lat_new, lon_new] = antipode(lat(ineg), lon(ineg));
    lat(ineg) = lat_new;
    lon(ineg) = lon_new;
    omg(ineg) = -omg(ineg);

    th = (90-lat)/deg;      % radians
    ph = lon/deg;           % radians
    
    % convert to spherical (x,y,z) points with length omega
    azi = ph;
    ele = pi/2 - th;
    [xx, yy, zz] = sph2cart(azi, ele, omg);
    
    outvecs = [xx(:) yy(:) zz(:)]';     % 3 x n
end

%-----------------------

% EXAMPLE
if 0==1
    % load rotation vectors (deg/Myr)
    stdir = '/home/carltape/splines/eh_plate/Rotation_from_Craig_ONeill/';
    ww = 'ind2002.a1';
    [wx,wy,wz,names] = textread([stdir ww],'%f%f%f%s','headerlines',2);
    
    invec = [wx wy wz]';
    outvec1 = euler_convert(invec,1);
    outvec2 = euler_convert(outvec1,0);
    
    % check
    mean(mean(abs(invec-outvec2)))
end

%==============================================================
