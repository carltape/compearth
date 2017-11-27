function [phi,sigma,theta] = U2pst(U)
%U2PST converts basis U to phi-sigma-theta
%
% INPUT
%   U       3 x 3 x n set of orientation matrices
%               
% OUTPUT
%   phi     n x 1 vector of azimuth spherical coordinate angle, degrees
%   sigma   n x 1 vector of longitudinal angle, degrees
%   theta   n x 1 vector of polar spherical angle, degrees
%
% See WTape and CTape (2017) "Volume in moment tensor space in terms of distance"
%
% Carl Tape, 2017/11
%

deg = 180/pi;

% U is assumed to be 3 x 3 x n
[~,~,n] = size(U);

v = NaN(3,n);
for ii=1:n
    V = U(:,:,ii);
    v(:,ii) = V(:,3);   % U*e3
end

% spherical coordinates (in radians)
[az,ele] = cart2sph(v(1,:),v(2,:),v(3,:));
% in degrees
phi = az*deg;
theta = (pi/2 - ele)*deg;

w = NaN(3,n);
for ii=1:n    
    Zp = rotmat(phi(ii),3);
    Yt = rotmat(theta(ii),2);
    A = Zp * Yt;
    W = A' * U(:,:,ii);
    w(:,ii) = W(:,1);   % A'*U*e1
end

% spherical coordinates (in radians)
[sigma,~] = cart2sph(w(1,:),w(2,:),w(3,:));
sigma = sigma*deg;

% column vectors
phi = phi(:);
sigma = sigma(:);
theta = theta(:);

%==========================================================================
% EXAMPLE

if 0==1
    % single set of angles
    deg = 180/pi;
    phi = 0.2*deg;
    sigma = 0.4*deg;
    z = 0.3;
    theta = acos(z)*deg;
    
    U = pst2U(phi,sigma,theta)
    [phicheck,sigmacheck,thetacheck] = U2pst(U);
    phi, phicheck, sigma, sigmacheck, theta, thetacheck
    
    % large set of angles
    deg = 180/pi;
    n = 1e4;
    phi = randomvec(-90,90,n);
    sigma = randomvec(-90,90,n);
    z = randomvec(-1,1,n);
    theta = deg*acos(z);
    
    U = pst2U(phi,sigma,theta);
    [phicheck,sigmacheck,thetacheck] = U2pst(U);
    norm(phi(:) - phicheck(:)) / norm(phi)
    norm(sigma(:) - sigmacheck(:)) / norm(sigma)
    norm(theta(:) - thetacheck(:)) / norm(theta)
    
    % compare with strike/dip/rake
    [kappaX,thetaX,sigmaX] = U2sdr(U);
    % scatterplots
    figure; nr=3; nc=3; kk=0;
    d1 = {kappaX,thetaX,sigmaX}; l1 = {'strike','dip','rake'};
    d2 = {phi,sigma,theta}; l2 = {'\phi','\sigma','\theta'};
    for ii=1:3
        for jj=1:3
            kk = kk+1;
            subplot(nr,nc,kk);
            plot(d1{ii},d2{jj},'.','markersize',2);
            axis equal;
            axis tight;
            xlabel(l1{ii}); ylabel(l2{jj});
        end
    end

end

%==========================================================================
