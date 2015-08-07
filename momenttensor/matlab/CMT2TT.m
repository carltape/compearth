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
% Carl Tape, 12/2012
%

if nargin==1, bdisplay=false; end

global BIGN EPSVAL
BIGN = 1e5;
EPSVAL = 1e-6;

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

Yrot = rotmat(45,2);

% compute candidate fault vectors
S = zeros(3,n);
N = zeros(3,n);
for ii=1:n
   V = U(:,:,ii) * Yrot;    % V = U * Yrot (TT2012, p. 487)
   S(:,ii) = V(:,1);        % slip vector
   N(:,ii) = V(:,3);        % fault normal
end

% reassign ~0 elements to =0; ~1 to =1, ~-1 to =1
N = setzero(N); S = setzero(S);

% compute fault angles for four possible combinations (TT2012 Figure 15)
S1 =  S; N1 =  N;
S2 = -S; N2 = -N;
S3 =  N; N3 =  S;
S4 = -N; N4 = -S;
[theta1,sigma1,kappa1,K1] = faultvec2ang(S1,N1);
[theta2,sigma2,kappa2,K2] = faultvec2ang(S2,N2);
[theta3,sigma3,kappa3,K3] = faultvec2ang(S3,N3);
[theta4,sigma4,kappa4,K4] = faultvec2ang(S4,N4);

% % reassign ~0 elements to =0; ~1 to =1, ~-1 to =1
% K1 = setzero(K1); N1 = setzero(N1); S1 = setzero(S1);
% K2 = setzero(K2); N2 = setzero(N2); S2 = setzero(S2);
% K3 = setzero(K3); N3 = setzero(N3); S3 = setzero(S3);
% K4 = setzero(K4); N4 = setzero(N4); S4 = setzero(S4);

% display all four options
if bdisplay
    %xlab = '(kappa, theta, sigma)';
    xlab = '(strike, dip, slip)';
    stfmt = '(%7.1f, %7.1f, %7.1f)';
    %stfmt = '(%.16e, %.16e, %.16e)';
    for ii=1:n
        displayCMTshort(M(:,ii))
        disp(sprintf(['%4i  S, N %s = ' stfmt],...
            ii,xlab,kappa1(ii),theta1(ii),sigma1(ii)));
        disp('         K         N         S');
        disp([K1(:,ii) N1(:,ii) S1(:,ii)]);
        disp(sprintf(['     -S,-N %s = ' stfmt],...
            xlab,kappa2(ii),theta2(ii),sigma2(ii)));
        disp('         K         N         S');
        disp([K2(:,ii) N2(:,ii) S2(:,ii)]);
        disp(sprintf(['      N, S %s = ' stfmt],...
            xlab,kappa3(ii),theta3(ii),sigma3(ii)));
        disp('         K         N         S');
        disp([K3(:,ii) N3(:,ii) S3(:,ii)]);
        disp(sprintf(['     -N,-S %s = ' stfmt],...
            xlab,kappa4(ii),theta4(ii),sigma4(ii)));
        disp('         K         N         S');
        disp([K4(:,ii) N4(:,ii) S4(:,ii)]);
    end
end

% There are four combinations of N and S that represent a double couple
% moment tensor, as shown in Figure 15 of TT2012.
% From these four combinations, there are two possible fault planes.
% We want to isolate the combination that is within the bounding
% region shown in Figures 16 and B1.
thetaall = [theta1 theta2 theta3 theta4];
sigmaall = [sigma1 sigma2 sigma3 sigma4];
kappaall = [kappa1 kappa2 kappa3 kappa4];
btheta = thetaall <= 90+EPSVAL;       % dip angles
bsigma = abs(sigmaall) <= 90+EPSVAL;  % slip angle
bmatch = and(btheta,bsigma);
imatch = NaN(n,1);
for ii=1:n
    itemp = find(bmatch(ii,:)==1);
    switch length(itemp)
        case 0
            error('no match');
        case 1
            imatch(ii) = itemp;
        case 2
            % choose one of the two
            i1 = itemp(1);
            i2 = itemp(2);
            ipick = pickP1(thetaall(ii,i1),sigmaall(ii,i1),kappaall(ii,i1),thetaall(ii,i2),sigmaall(ii,i2),kappaall(ii,i2));
            imatch(ii) = itemp(ipick);
            if 0==1     % display output
                warning('moment tensor on boundary of orientation domain (%i candidates)',length(itemp));
                display_vals(thetaall(ii,:),sigmaall(ii,:),kappaall(ii,:),M(:,ii),ii,n);
                itemp, imatch(ii)
                disp(sprintf(' theta: %5.1f',thetaall(ii,imatch(ii))));
                disp(sprintf(' sigma: %5.1f',sigmaall(ii,imatch(ii))));
                disp(sprintf(' kappa: %5.1f',kappaall(ii,imatch(ii))));
            end
        case 3
            % this is a more unusual case, like for horizontal faults
            warning('moment tensor on boundary of orientation domain (%i candidates)',length(itemp));
            display_vals(thetaall(ii,:),sigmaall(ii,:),kappaall(ii,:),M(:,ii),ii,n);
            itemp
            % just take the first one in the list (for now)
            imatch(ii) = itemp(1);
        case 4
            error('not yet encountered');
    end
end

% get fault vectors
K = NaN(3,n); N = NaN(3,n); S = NaN(3,n);
kappa = NaN(n,1); theta = NaN(n,1); sigma = NaN(n,1); 
for ii=1:n
   kk = imatch(ii);
   switch kk
   case 1, K(:,ii) = K1(:,ii); N(:,ii) = N1(:,ii); S(:,ii) = S1(:,ii);
       kappa(ii) = kappa1(ii); theta(ii) = theta1(ii); sigma(ii) = sigma1(ii);
   case 2, K(:,ii) = K2(:,ii); N(:,ii) = N2(:,ii); S(:,ii) = S2(:,ii);
       kappa(ii) = kappa2(ii); theta(ii) = theta2(ii); sigma(ii) = sigma2(ii);
   case 3, K(:,ii) = K3(:,ii); N(:,ii) = N3(:,ii); S(:,ii) = S3(:,ii);
       kappa(ii) = kappa3(ii); theta(ii) = theta3(ii); sigma(ii) = sigma3(ii);
   case 4, K(:,ii) = K4(:,ii); N(:,ii) = N4(:,ii); S(:,ii) = S4(:,ii);
       kappa(ii) = kappa4(ii); theta(ii) = theta4(ii); sigma(ii) = sigma4(ii);
   end
end

%--------------------------------------------------------------------------

function [theta,sigma,kappa,K] = faultvec2ang(S,N)
% returns fault angles in degrees, assumes input vectors in south-east-up basis
global BIGN

deg = 180/pi;

[~,n] = size(S);
if n >= BIGN, disp('CMT2TT: running faultvec2ang...'); end

% for north-west-up basis (as in TT2012)
%zenith = [0 0 1]'; north  = [1 0 0]';

% for up-south-east basis (GCMT)
%zenith = [1 0 0]'; north  = [0 -1 0]';

% for south-east-up basis (as in TT2012)
zenith = [0 0 1]'; north  = [-1 0 0]';

kappa = NaN(n,1);
theta = NaN(n,1);
sigma = NaN(n,1);
K = NaN(3,n);
for ii=1:n
    % strike vector from TT2012, Eq. 29
    v = cross(zenith,N(:,ii));
    if norm(v)==0
        % TT2012 Appendix B
        disp('horizontal fault -- strike vector is same as slip vector');
        K(:,ii) = S(:,ii);
    else
        K(:,ii) = v / norm(v);
    end

    % Figure 14
    kappa(ii) = fangle_signed(north,K(:,ii),-zenith);

    % Figure 14
    costh = dot(N(:,ii),zenith);
    theta(ii) = acos(costh)*deg;

    % Figure 14
    sigma(ii) = fangle_signed(K(:,ii),S(:,ii),N(:,ii));
end

kappa = wrap360(kappa);

%--------------------------------------------------------------------------

function X = setzero(X)
% try to eliminate numerical round-off errors
global EPSVAL

% elements near zero
XN = X ./ max(abs(X(:)));
X(abs(XN) < EPSVAL) = 0;

% elements near +/- 1
X(abs(X - 1) < EPSVAL) = -1;
X(abs(X + 1) < EPSVAL) =  1;

%--------------------------------------------------------------------------

function display_vals(thetas,sigmas,kappas,M,ii,n)

disp(sprintf('index in list is %i/%i',ii,n));
disp('Mrr,Mtt,Mpp,Mrt,Mrp,Mtp (N-m):');
for kk=1:6, disp(sprintf('%18.8e',M(kk))); end
disp(sprintf('thetas: %5.1f, %5.1f, %5.1f, %5.1f',thetas));
disp(sprintf('sigmas: %5.1f, %5.1f, %5.1f, %5.1f',sigmas));
disp(sprintf('kappas: %5.1f, %5.1f, %5.1f, %5.1f',kappas));
            
%--------------------------------------------------------------------------

function ipick = pickP1(thetaA,sigmaA,kappaA,thetaB,sigmaB,kappaB)
% choose between two moment tensor orientations based on Figure B1 of TT2012
% NOTE THAT NOT ALL FEATURES OF FIGURE B1 ARE IMPLEMENTED HERE
global EPSVAL

% these choices are based on the strike angle
if abs(thetaA - 90) < EPSVAL
    ipick = find([kappaA kappaB] < 180); return
end
if abs(sigmaA - 90) < EPSVAL
    ipick = find([kappaA kappaB] < 180); return
end
if abs(sigmaA + 90) < EPSVAL
    ipick = find([kappaA kappaB] < 180); return
end 

thetaA,sigmaA,kappaA,thetaB,sigmaB,kappaB
error('no selection criterion was met');

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
