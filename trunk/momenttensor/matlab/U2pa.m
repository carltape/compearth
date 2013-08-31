function Uout = U2pa(Uin,itype)
%U2PA convert between basis U and plunge/azimuth of three basis vectors
%
% INPUT
%   Uin     either 3 x 3 x n array U OR n x 6 set of plunge/azimuth angles
%   itype   =1 for U to plunge/azimuth; =0 for plunge/azimuth to U
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

if itype==1
    % TO BE SAFE, WE COULD CALL Uorth.m TO ENSURE THAT U IS ORTHOGONAL,
    % THEN CALL Udetcheck.m TO ENSURE THAT U IS A ROTATION MATRIX.
    [~,~,n] = size(Uin);

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
    az1 = wrapTo360(-lon1);
    az2 = wrapTo360(-lon2);
    az3 = wrapTo360(-lon3);

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
    Uout = Uorth(Uout,1,0);
    
    % ensure that these are rotation matrices
    Uout = Udetcheck(Uout);
end

%==========================================================================
% EXAMPLE

if 0==1
    clear, clc, close all
    teid = 'M010176A';  % first event in GCMT catalog
    [otime,tshift,hdur,lat,lon,dep,M,~,~,eid,elabel,str1,dip1,rk1,str2,dip2,rk2,...
        lams,pl1,az1,pl2,az2,pl3,az3,icmt,dataused,itri] = readCMT(teid);
    [MDC,kap1,theta1,sig1,kap2,theta2,sig2,k1,d1,n1,k2,d2,n2,U0,lamsx] = CMT2faultpar(M,0);
    
    % check
    Mt = U0 * diag(lamsx) * U0';
    disp([M CMTconvert(Mvec2Mmat(Mt,0),5,1)])
    
    U = U0;
% check conversion to GCMT convention
%     U(1,:) = U0(3,:);
%     U(2,:) = U0(1,:);
%     U(3,:) = U0(2,:);
%     U * diag(lamsx) * U'
    
    Uout = U2pa(U,1);
    % compare with angles listed in GCMT catalog
    [Uout ; pl1 az1 pl2 az2 pl3 az3]
    
    %----------------------------
    % convert back to U
    
    Ucheck = U2pa(Uout,0);
    Ucheck, U                       % note: no problem that the U do not match
    Mt, Ucheck*diag(lamsx)*Ucheck'  %       because the moment tensors do!
    
    %----------------------------
    % GCMT catalog
    
    clear, clc, close all
    % load CMT catalog
    isub = [1:10]';
    %isub = [1:1000]';
    %isub = [1:33507]';
    [otime,tshift,hdur,lat,lon,dep,M,~,~,eid,elabel,str1,dip1,rk1,str2,dip2,rk2,...
        lams,pl1,az1,pl2,az2,pl3,az3,icmt,dataused,itri] = readCMT(isub);
    [MDC,kap1,theta1,sig1,kap2,theta2,sig2,k1,d1,n1,k2,d2,n2,U,lamsx] = CMT2faultpar(M,0);
    n = length(isub);
    
    Uout = U2pa(U,1);
    if 0==1
        Ucheck = U2pa(Uout,0);
        Min    = CMTrecom(lamsx,U)              % south-east-up
        Mcheck = CMTrecom(lamsx,Ucheck)         % south-east-up
        (CMTconvert(Mcheck,5,1) - M) ./ M       % up-south-east check
    end
    pl1x = Uout(:,1);
    az1x = Uout(:,2);
    pl2x = Uout(:,3);
    az2x = Uout(:,4);
    pl3x = Uout(:,5);
    az3x = Uout(:,6);

    % compare eig angles from GCMT with eig angles computed by me from
    % GCMT M to U to eig angles
    % NOTE 1: There are 24 plunge angle discrepancies and 1051 azimuth
    % angle discrepancies.
    % NOTE 2: There is definitely a time dependence associated with the
    % discrepancies; they drop off after 2003.
    ATRSH = 2.0;
    ipbad1 = find(abs(pl1x-pl1) >= ATRSH);
    ipbad2 = find(abs(pl2x-pl2) >= ATRSH);
    ipbad3 = find(abs(pl3x-pl3) >= ATRSH);
    iabad1 = find(abs(az1x-az1) >= ATRSH);
    iabad2 = find(abs(az2x-az2) >= ATRSH);
    iabad3 = find(abs(az3x-az3) >= ATRSH);
    ipbad = unique([ipbad1 ; ipbad2 ; ipbad3]);
    iabad = unique([iabad1 ; iabad2 ; iabad3]);
    length(ipbad), length(iabad)
    
    figure; nr=3; nc=2;
    subplot(nr,nc,1); hold on; plot(isub,pl1x-pl1,'.'); xlim([0 n]);
    if ~isempty(ipbad1), plot(ipbad1,pl1x(ipbad1)-pl1(ipbad1),'ro'); end
    title(sprintf('U2pa.m, eig1 plunge: %i / %i with diff >= %.1f',length(ipbad1),n,ATRSH));
    subplot(nr,nc,3); hold on; plot(isub,pl2x-pl2,'.'); xlim([0 n]);
    if ~isempty(ipbad2), plot(ipbad2,pl2x(ipbad2)-pl2(ipbad2),'ro'); end
    title(sprintf('U2pa.m, eig2 plunge: %i / %i with diff >= %.1f',length(ipbad2),n,ATRSH));
    subplot(nr,nc,5); hold on; plot(isub,pl3x-pl3,'.'); xlim([0 n]);
    if ~isempty(ipbad3), plot(ipbad3,pl3x(ipbad3)-pl3(ipbad3),'ro'); end
    title(sprintf('U2pa.m, eig3 plunge: %i / %i with diff >= %.1f',length(ipbad3),n,ATRSH));
    
    subplot(nr,nc,2); hold on; plot(isub,az1x-az1,'.'); xlim([0 n]);
    if ~isempty(iabad1), plot(iabad1,az1x(iabad1)-az1(iabad1),'ro'); end
    title(sprintf('U2pa.m, eig1 azimuth: %i / %i with diff >= %.1f',length(iabad1),n,ATRSH));
    subplot(nr,nc,4); hold on; plot(isub,az2x-az2,'.'); xlim([0 n]);
    if ~isempty(iabad2), plot(iabad2,az2x(iabad2)-az2(iabad2),'ro'); end
    title(sprintf('U2pa.m, eig2 azimuth: %i / %i with diff >= %.1f',length(iabad2),n,ATRSH));
    subplot(nr,nc,6); hold on; plot(isub,az3x-az3,'.'); xlim([0 n]);
    if ~isempty(iabad3), plot(iabad3,az3x(iabad3)-az3(iabad3),'ro'); end
    title(sprintf('U2pa.m, eig3 azimuth: %i / %i with diff >= %.1f',length(iabad3),n,ATRSH));
    orient tall, wysiwyg
    
    figure; plot_histo(year(otime(iabad)),[1976:2011]);
end

%==========================================================================
