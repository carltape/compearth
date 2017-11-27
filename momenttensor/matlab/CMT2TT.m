function [gamma,delta,M0,kappa,theta,sigma,K,N,S,thetadc,lam,U] = CMT2TT(M,bdisplay)
%CMT2TT converts a moment tensor to six parameters of TapeTape2012
%
% INPUT
%   M           6 x n moment tensors in CMT convention (UP-SOUTH-EAST)
%               M = [Mrr Mtt Mpp Mrt Mrp Mtp]; r=up, theta=south, phi=east
%   bdisplay    OPTIONAL (if present, will display details)
%               
% OUTPUT
%   gamma       angle from DC meridian to MT point (-30 <= gamma <= 30)
%   delta       angle from deviatoric plane to MT point (-90 <= delta <= 90)
%   M0          seismic moment, N-m
%   kappa       strike angle, degrees: [0,360]
%   theta       dip angle, degrees: [0,90]
%   sigma       slip (or rake) angle, degrees: [-90,90]
% optional:
%   K           strike vector (SOUTH-EAST-UP)
%   N           normal vector (SOUTH-EAST-UP)
%   S           slip vector (SOUTH-EAST-UP)
%   thetadc     angle to DC
%   lam         eigenvalues
%   U           basis (SOUTH-EAST-UP)
%
% Reverse program for TT2CMT.m
% See WTape and CTape (2012) "A geometric setting for moment tensors" (TT2012).
%
% Carl Tape, 2012/12
%

if nargin==1, bdisplay=false; end

BIGN = 1e5;

% make sure M is 6 x n
[M,n] = Mdim(M);

disp(sprintf('CMT2TT.m: %i moment tensors to convert into lune + strike/dip/rake',n));

% KEY: convert M into another basis
% YOU MUST ALSO CHANGE north AND zenith IN faultvec2ang BELOW
% --> U will be with respect to this basis (from CMTdecom.m)
%M = convert_MT(1,3,M);
M = convert_MT(1,5,M);  % moment tensor in south-east-up basis

%---------------------
% PART 1: moment tensor source type (or pattern)

% decomposition of moment tensor into eigenvalues + basis (M = U*lam*U')
% NOTE: ordering of eigenvalues is important
if n >= BIGN, disp('CMT2TT: M to lam + U...'); end
isort = 1;
[lam,U] = CMTdecom(M,isort);

% compute lune coordinates and magnitude from eigenvalues
if n >= BIGN, disp('CMT2TT: lam to lune...'); end
[gamma,delta,M0,thetadc] = lam2lune(lam);

%---------------------
% PART 2: moment tensor orientation
% TT2012, Section 6.3

[kappa,theta,sigma,K,N,S] = U2sdr(U,bdisplay);

%==========================================================================
% EXAMPLE

if 0==1
    clear, close all, clc
    n = 10;
    M0 = 1*ones(n,1);
    % random source types
    gamma = randomvec(-30,30,n); b = randomvec(-1,1,n); delta = asin(b)*180/pi;
    % random orientations
    kappa = randomvec(0,360,n); h = randomvec(0,1,n); sigma = randomvec(-90,90,n);
    theta = acos(h)*180/pi;
    %kappa = -10; theta = 30; sigma = 45; gamma = -20; delta = 70; M0 = 1;
    %kappa = -10; theta = 30; sigma = 45; gamma = -20; delta = 70; M0 = 1;
    M = TT2CMT(gamma,delta,M0,kappa,theta,sigma);
    
    [gammac,deltac,M0c,kappac,thetac,sigmac,K,N,S,thetadc,lam,U] = CMT2TT(M,1);
    disp([gamma gammac delta deltac M0 M0c kappa kappac theta thetac sigma sigmac]);
    
    % horizontal fault
    kappa = 30; theta = 0; sigma = 310;
    M = TT2CMT(0,0,1,kappa,theta,sigma)
    M = [0 0 0 -sqrt(3)/2 1/2 0]'
    [gammac,deltac,M0c,kappac,thetac,sigmac,K,N,S,thetadc,lam,U] = CMT2TT(M,1);
    disp([kappa kappac theta thetac sigma sigmac]);
    
    % GCMT catalog
    [otime,tshift,hdur,lat,lon,dep,M,M0,Mw,eid] = readCMT;
    display_eq_summary(otime,lon,lat,dep,Mw);
    [gammac,deltac,M0c,kappac,thetac,sigmac] = CMT2TT(M);
    
end

%==========================================================================
