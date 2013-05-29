%
% function U = pa2U(pl1,az1,pl2,az2,pl3,az3)
% Carl Tape, 12-August-2011
%
% This function computes the eigenbasis rotation matrix from a set of
% (plunge, azimuth) angles describing three eigenvectors.
% This is the reverse program to U2pa.m.
%
% INPUT
%   pl1,az1     plunge and trend for 1st eigenvector
%   pl2,az2     plunge and trend for 2nd eigenvector
%   pl3,az3     plunge and trend for 3rd eigenvector
%
% OUTPUT
%   U       eigenbasis in SOUTH-EAST-UP convention and in the green zone
% 
% calls Ueigvec.m, latlon2xyz.m
% called by xxx
%

function U = pa2U(pl1,az1,pl2,az2,pl3,az3,EPSVAL)

n = length(pl1);

p1 = latlon2xyz(pl1,-az1);
p2 = latlon2xyz(pl2,-az2);
p3 = latlon2xyz(pl3,-az3);

U = zeros(3,3,n);
U(:,1,:) = p1;
U(:,2,:) = p2;
U(:,3,:) = p3;

% see EPSVAL in rotmat2rotvec.m
if nargin==6, EPSVAL = 1e-6; end
[U,igreen] = Ueigvec(U,EPSVAL);

%==========================================================================
% EXAMPLE

if 0==1
    clear, clc, close all
    % load CMT catalog
    %isub = [1:10]';
    isub = [1:33507]';
    %isub = [559 702 730]';
    %isub = [7059 7914 16810]';  % U = I causes problems
    n = length(isub);
    [otime,tshift,hdur,lat,lon,dep,M,~,~,eid,elabel,str1,dip1,rk1,str2,dip2,rk2,...
        lams,pl1,az1,pl2,az2,pl3,az3,icmt,dataused,itri] = readCMT_all(isub);
    [MDC,kap1,del1,lam1,kap2,del2,lam2,k1,d1,n1,k2,d2,n2,U0,lamsx] = CMT2faultpar(M,0);
    
    % take GCMT eigvec angles to compute U
    EPSVAL = 2e-2;
    U1 = pa2U(pl1,az1,pl2,az2,pl3,az3,EPSVAL);
    
    % take my U derived from GCMT M, then compute eigvec angles, then new U
    [pl1x,az1x,pl2x,az2x,pl3x,az3x] = U2pa(U0);
    U2 = pa2U(pl1x,az1x,pl2x,az2x,pl3x,az3x,1e-6);
    
    % compare U:
    % (1) U computed from the GCMT azimuth and plunge angles of eigenvectors;
    % (2) U computed from the GCMT moment tensor
    %
    % NOTE 1: 143 out of 33507 do NOT match, due to sign differences that
    % probably(?) arise because the angles computed from M are tied
    % to the precision of M, whereas GCMT provides angles with single-digit
    % precision.
    %
    % NOTE 2: If you use EPSVAL=2e-2, then 3 out of 33507 MTs will not
    % match the direct forward-reverse operation; this is because the 3 MTs
    % are so close to the base MT that rotmat2rotvec.m operation assumes U = I.
    %
    udiff1 = zeros(n,1);
    for ii=1:n
        u0 = squeeze(U0(:,:,ii));
        u1 = squeeze(U1(:,:,ii));
        u2 = squeeze(U2(:,:,ii));
        %disp(u0(:)'); disp(u1(:)'); %disp(u2(:)');
        udiff1(ii) = norm(u1(:) - u0(:));
        udiff2(ii) = norm(u2(:) - u0(:));
    end
    
    NTRSH = 0.5;
    ibad1 = find(abs(udiff1 >= NTRSH));
    figure; plot(isub,udiff1,'b.',ibad1,udiff1(ibad1),'ro');
    xlim([0 n]); ylabel('norm(U1 - U0)');
    title(sprintf('%i / %i have diff >= %.2f',length(ibad1),n,NTRSH));
    
    ibad2 = find(abs(udiff2 >= NTRSH));
    figure; plot(isub,udiff2,'b.',ibad2,udiff2(ibad2),'ro');
    xlim([0 n]); ylabel('norm(U1 - U0)');
    title(sprintf('%i / %i have diff >= %.2f',length(ibad2),n,NTRSH));
end

%==========================================================================
