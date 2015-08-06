function [kappaP,thetaP,sigmaP,kappap,thetap,sigmap] = dcfaultpar2Pbox(kappa,theta,sigma)
%DCFAULTPAR2PBOX convert strike-dip-rake angles to P box of Tape and Tape (2012)
%
% The box U contains all possible strike-dip-rake triples;
% it can be thought of as giving all possible orientations.
% The box P is 1/4 the size of U and is the box of all possible
% double couple moment tensors. The four-fold difference arises because you
% can rotate any moment tensor by 180 about three different axes and still
% have the same beachball. See Tape and Tape (2012).
% In other words, there are four possible frames for any moment tensor.
% 
% INPUT:
%   kappa   strike (any angle)
%   theta   dip (any angle)
%   sigma   rake (any angle)
%
% OUTPUT:
%   kappa   strike (0 to 360)
%   theta   dip (0 to 90)
%   sigma   rake (-90 to 90)
%
% See also dcfaultpar2Ubox.m
% Note: The underlying operations are more transparent within CMT2TT.m and
%       TT2CMT.m. For example, one can take any strike-dip-rake triple,
%       feed it into TT2CMT to get M, then CMT2TT to get the strike-dip-rake
%       triple that is in P.
%
% Carl Tape, 10/29/2013
%

n = length(kappa);
deg = 180/pi;

[kappa,theta,sigma] = dcfaultpar2Ubox(kappa,theta,sigma);

cossig = cos(sigma/deg);
sinsig = sin(sigma/deg);
costh  = cos(theta/deg);
sinth  = sin(theta/deg);

% get phi coordinate (for kappaprime)
% note: cart2sph returns [-pi,pi] for azimuth
[az,ele,rho] = cart2sph(-costh.*sinsig, cossig, zeros(size(kappa)));
ph = deg*az;
%th = deg*(pi/2 - ele);

% from extra notes
kappap = wrapTo360(kappa - ph);
thetap = deg * acos(sinsig .* sinth);
sigmap = deg * sign(costh) .* ...
    acos( (-cossig.*sinth) ./ sqrt(cossig.^2 + costh.^2.*sinsig.^2) );

[kappap,thetap,sigmap] = dcfaultpar2Ubox(kappap,thetap,sigmap);

kappaP = NaN(size(kappa));
thetaP = NaN(size(theta));
sigmaP = NaN(size(sigma));
for ii=1:n
   sigx = sigma(ii);
   thtx = theta(ii);
   
   if and(abs(sigx) <= 90,thtx <= 90)
       kappaP(ii) = kappa(ii);
       sigmaP(ii) = sigma(ii);
       thetaP(ii) = theta(ii);
       
   elseif and(abs(sigx) <= 90,thtx >= 90)
       kappaP(ii) = kappa(ii) + 180;
       sigmaP(ii) = -sigma(ii);
       thetaP(ii) = 180 - theta(ii);
       
   elseif sigx >= 90
       kappaP(ii) = kappap(ii);
       sigmaP(ii) = sigmap(ii);
       thetaP(ii) = thetap(ii);
       
   elseif sigx <= -90
       kappaP(ii) = kappap(ii) + 180;
       sigmaP(ii) = -sigmap(ii);
       thetaP(ii) = 180 - thetap(ii);
   else
      ii, sigx, thtx
      error('check conditions'); 
   end
end

[kappaP,thetaP,sigmaP] = dcfaultpar2Ubox(kappaP,thetaP,sigmaP);

%==========================================================================

if 0==1
    % using dcfaultpar2Pbox.m
    %n = 10; kappa = 360*rand(n,1); theta = 180*rand(n,1); sigma = -180 + 360*rand(n,1);
    kappa = 85; theta = 40; sigma = 265;
    %kappa = 285; theta = 161; sigma = 29;
    [kappa,theta,sigma] = dcfaultpar2Pbox(kappa,theta,sigma);
    [kappaP,thetaP,sigmaP,kappap,thetap,sigmap] = dcfaultpar2Pbox(kappa,theta,sigma);
    disp('   kappa     kappap    kappaP     theta    thetap     thetaP    sigma     sigmap   sigmaP');
    disp([kappa kappap kappaP theta thetap thetaP sigma sigmap sigmaP]);
    % check that all three triples give the same moment tensor
    M = TT2CMT(0,0,1,kappa,theta,sigma);
    Mp = TT2CMT(0,0,1,kappap,thetap,sigmap);
    MP = TT2CMT(0,0,1,kappaP,thetaP,sigmaP);
    disp([M Mp MP]);
    
    % alternatively you can use other functions
    M = TT2CMT(0,0,1,kappa,theta,sigma);
    [gamma,delta,M0,kappaP,thetaP,sigmaP] = CMT2TT(M);
    disp([kappa kappaP theta thetaP sigma sigmaP]);
end

%==========================================================================