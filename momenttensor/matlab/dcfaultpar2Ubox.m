function [kappa,theta,sigma] = dcfaultpar2Ubox(kappa,theta,sigma)
%DCFAULTPAR2UBOX convert strike-dip-rake angles to 'normal' ranges
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
%   theta   dip (0 to 180)
%   sigma   rake (-180 to 180)
%
% See also dcfaultpar2Pbox.m
%
% Carl Tape, 10/29/2013
%

% mapping theta (dip) is more complicated than sigma and kappa
% the key identity is: V(kappa,sigma,-theta) = V(kappa+pi,sigma+pi,theta)
theta = wrapTo180(theta);
ineg = find(theta < 0);
if ~isempty(ineg)
    warning('%i/%i theta values are negative',length(ineg),length(theta));
    sigma(ineg) = sigma(ineg) + 180;
    kappa(ineg) = kappa(ineg) + 180;
    theta(ineg) = -theta(ineg);
end
    
% sigma and kappa have simple periodicity of 360
sigma = wrapTo180(sigma);
kappa = wrapTo360(kappa);

%==========================================================================

if 0==1
    % one example
    kappa = 1000; theta = 170; sigma = -1000;
    [kappaU,thetaU,sigmaU] = dcfaultpar2Ubox(kappa,theta,sigma);
    disp([kappa kappaU theta thetaU sigma sigmaU]);
    % check that the moment tensors are the same
    M = TT2CMT(0,0,1,kappa,theta,sigma);
    MU = TT2CMT(0,0,1,kappaU,thetaU,sigmaU);
    disp([M MU]);
    
    % a set or random values
    n = 1000;
    kappa = -500 + 1000*rand(n,1);
    theta = -500 + 1000*rand(n,1);
    sigma = -500 + 1000*rand(n,1);
    [kappaU,thetaU,sigmaU] = dcfaultpar2Ubox(kappa,theta,sigma);
    figure; nr=2; nc=2; 
    subplot(nr,nc,1); plot(kappa,kappaU,'b.'); axis equal; grid on; title('strike');
    subplot(nr,nc,2); plot(theta,thetaU,'b.'); axis equal; grid on; title('dip');
    subplot(nr,nc,3); plot(sigma,sigmaU,'b.'); axis equal; grid on; title('rake');
    % check that the sets of moment tensors are the same
    M = TT2CMT(0,0,1,kappa,theta,sigma);
    MU = TT2CMT(0,0,1,kappaU,thetaU,sigmaU);
    Mdiff = NaN(n,1);
    for ii=1:n, Mdiff(ii) = norm(M(:,ii) - MU(:,ii)); end
    figure; plot(Mdiff,'.');
end

%==========================================================================