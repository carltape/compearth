function Uout = U2pa(Uin,itype,iorthoU)
%U2PA convert between basis U and plunge/azimuth of three basis vectors
%
% INPUT
%   Uin     either 3 x 3 x n array U OR n x 6 set of plunge/azimuth angles
%   itype   =1 for U to plunge/azimuth
%           =0 for plunge/azimuth to U
%   iorthoU OPTIONAL: type of orthogonalization to apply to U (see Uorth.m)
%   
% OUTPUT
%   Uout    either 3 x 3 x n array U OR n x 6 set of plunge/azimuth angles
%
%   U       eigenbasis in SOUTH-EAST-UP convention
%           (but U does NOT have to be in the green zone)
%
%   pl1,az1 plunge and azimuth for 1st eigenvector
%   pl2,az2 plunge and azimuth for 2nd eigenvector
%   pl3,az3 plunge and azimuth for 3rd eigenvector
%   eigenvectors ordered as lam1 >= lam2 >= lam3
% 
% Carl Tape, 12-August-2011
%

% by default we will orthogonalize the input U and the output U
if nargin==2
   iorthoU = 1; 
end

if itype==1
    Uin = Uorth(Uin,iorthoU,0);
    Uin = Udetcheck(Uin);

    p1 = squeeze(Uin(:,1,:));
    p2 = squeeze(Uin(:,2,:));
    p3 = squeeze(Uin(:,3,:));

    [lat1,lon1] = xyz2latlon(p1);
    [lat2,lon2] = xyz2latlon(p2);
    [lat3,lon3] = xyz2latlon(p3);

    % plunge must be 0 to 90
    ipos1 = find(p1(3,:) < 0);
    ipos2 = find(p2(3,:) < 0);
    ipos3 = find(p3(3,:) < 0);
    [~,~,n] = size(Uin);
    disp(sprintf('%i/%i p1 vectors need to be flipped',length(ipos1),n));
    disp(sprintf('%i/%i p2 vectors need to be flipped',length(ipos2),n));
    disp(sprintf('%i/%i p3 vectors need to be flipped',length(ipos3),n));

    % note these angles are with respect to the LOCAL basis
    [lat1(ipos1),lon1(ipos1)] = antipode(lat1(ipos1),lon1(ipos1));
    [lat2(ipos2),lon2(ipos2)] = antipode(lat2(ipos2),lon2(ipos2));
    [lat3(ipos3),lon3(ipos3)] = antipode(lat3(ipos3),lon3(ipos3));

    pl1 = lat1;
    pl2 = lat2;
    pl3 = lat3;
    az1 = wrap360(-lon1);
    az2 = wrap360(-lon2);
    az3 = wrap360(-lon3);

    Uout = [pl1 az1 pl2 az2 pl3 az3];
    
else
    deg = 180/pi;
    [n,m] = size(Uin);
    if m==6
        % all three eigenvectors provided
        pl1 = Uin(:,1);
        az1 = Uin(:,2);
        pl2 = Uin(:,3);
        az2 = Uin(:,4);
        pl3 = Uin(:,5);
        az3 = Uin(:,6);
    else
        % neutral axis not provided
        pl1 = Uin(:,1);
        az1 = Uin(:,2);
        pl3 = Uin(:,3);
        az3 = Uin(:,4);
    end

    % convert to normal spherical angles
    ph1 = ph2az(az1);
    ph3 = ph2az(az3);
    th1 = 90 + pl1;
    th3 = 90 + pl3;
    u1 = tp2xyz(th1/deg,ph1/deg,1);
    u3 = tp2xyz(th3/deg,ph3/deg,1);
    
    if m==4
        u2 = cross(u3,u1);
    else
        ph2 = ph2az(az2);
        th2 = 90 + pl2;
        u2 = tp2xyz(th2/deg,ph2/deg,1);
    end
    
    Uout = zeros(3,3,n);
    Uout(:,1,:) = u1;
    Uout(:,2,:) = u2;
    Uout(:,3,:) = u3;
    
    % transform from east-north-up to south-east-up
    P = [0 1 0 ; -1 0 0 ; 0 0 1];
    for ii=1:n
        Uout(:,:,ii) = P'*squeeze(Uout(:,:,ii));
    end
    
    % ensure orthogonal U by using SVD
    Uout = Uorth(Uout,iorthoU,0);
    
    % ensure that these are rotation matrices
    Uout = Udetcheck(Uout);
end

%==========================================================================

function [lat,lon] = antipode(lat,lon)

lat = -lat;
lon = 180 - mod(-lon, 360);

%==========================================================================

function v_out = ph2az(v_in)
%PH2AZ convert between math and map conventions for azimuthal angle, in degrees
%
%   phi       azimuth
%   0           90
%   90          0
%   180         270
%   270         180
%
% NOTE: The formula is the same, whether you are going from ph2az or az2ph.
% NOTE: does matlab have a built-in function for this?

v_in  = wrap360(v_in);
v_out = -v_in + 90;
v_out = wrap360(v_out);

%==========================================================================
% EXAMPLE

if 0==1
    teid = 'M010176A';  % first event in GCMT catalog
    [otime,tshift,hdur,lat,lon,dep,M,~,~,eid,elabel,str1,dip1,rk1,str2,dip2,rk2,...
        lams,pl1,az1,pl2,az2,pl3,az3,icmt,dataused,itri] = readCMT(teid);
    [MDC,kap1,theta1,sig1,kap2,theta2,sig2,k1,d1,n1,k2,d2,n2,U0,lamsx] = CMT2dcfaultpar(M,0);
    
    % check
    Mt = U0 * diag(lamsx) * U0';
    disp([M convert_MT(5,1,Mvec2Mmat(Mt,0))])
    
    U = U0;
% check conversion to GCMT convention
%     U(1,:) = U0(3,:);
%     U(2,:) = U0(1,:);
%     U(3,:) = U0(2,:);
%     U * diag(lamsx) * U'
    
    Uout = U2pa(U,1);
    % compare Uout (top row) with angles listed in GCMT catalog (bottom row)
    [Uout ; pl1 az1 pl2 az2 pl3 az3]
    
    %----------------------------
    % convert back to U
    
    Ucheck = U2pa(Uout,0);
    Ucheck, U                       % note: no problem that the U do not match
    Mt, Ucheck*diag(lamsx)*Ucheck'  %       because the moment tensors do!
end

%==========================================================================
