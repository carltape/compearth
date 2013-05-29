function lam = lune2lam(gamma,delta,M0)
%LUNE2LAM convert lune coordinates (gamma, delta, M0) to eigenvalues
%
% INPUT
%   gamma   n x 1 vector of gamma angles, degrees
%   delta   n x 1 vector of delta angles, degrees
%   M0      1 x n vector of seismic moments
%
% OUTPUT
%   lam     3 x n set of diagonalized moment tensors in GCMT convention
%           note: normalized such that each MT has moment M0
%
% Reverse program for lam2lune.m
%
% See TapeTape2012 "A geometric setting for moment tensors".
%
% Carl Tape, 08-July-2011
%

disp('entering lune2lam.m');

deg = 180/pi;

delta = delta(:);   % (column vector)
gamma = gamma(:);   % (column vector)
beta = 90 - delta;  % colatitude
M0 = M0(:)';        % (row vector)
n = length(gamma);

% magnitude of lambda vectors (TT2012, p. 490 text)
rho = M0*sqrt(2);

% convert to eigenvalues (TT2012, Eq. 20)
% matrix to rotate points such that delta = 90 is (1,1,1) and delta = -90 is (-1,-1,-1)
R = 1/sqrt(6) * [sqrt(3) 0 -sqrt(3) ; -1 2 -1 ; sqrt(2) sqrt(2) sqrt(2)];

% Cartesian points as 3 x n unit vectors (TT2012, Eq. 20)
%Pxyz = latlon2xyz(delta,gamma,ones(n,1));
Pxyz = [cos(gamma/deg).*sin(beta/deg) sin(gamma/deg).*sin(beta/deg) cos(beta/deg)]';

% rotate points and apply magnitudes
lamhat = R' * Pxyz;
lam = lamhat .* repmat(rho,3,1);

idisplay = 0;
icheck = 1;
if idisplay==1
    disp(sprintf('lune2lam.m with %i sets of gamma-delta angles:',n));
    if icheck==1
        % check by performing the reverse operation with lam2lune(M)
        [gammacheck,deltacheck,M0check] = lam2lune(lam);
        for ii=1:n
           disp(sprintf('%3i gamma %6.2f delta %6.2f M0 %7.2e lams: %7.2e %7.2e %7.2e',...
                ii,gamma(ii),delta(ii),M0(ii),lam(:,ii)'));
           disp(sprintf('%3i gamma %6.2f delta %6.2f M0 %7.2e',...
                ii,gammacheck(ii),deltacheck(ii),M0check(ii))); 
        end
    else
        for ii=1:n
           disp(sprintf('%3i gamma %6.2f delta %6.2f M0 %7.2e lams: %7.2e %7.2e %7.2e',...
                ii,gamma(ii),delta(ii),M0(ii),lam(:,ii)'));
        end
    end
end

%==========================================================================
% EXAMPLE

if 0==1
    clear, clc, close all
    
    nx = 100;
    gvec = linspace(-30,30,nx);
    bvec = linspace(-89,89,nx);
    [G,B] = meshgrid(gvec,bvec);
    gamma0 = G(:);
    delta0 = B(:);
    M00 = 1e16*ones(length(gamma0),1);
    D = lune2lam(gamma0,delta0,M00);
    
    [gamma,delta,M0,mu] = lam2lune(D);
    
    figure; nr=2; nc=2;
    subplot(nr,nc,1); plot(gamma-gamma0,'.'); title('gamma residual');
    subplot(nr,nc,2); plot(delta-delta0,'.'); title('delta residual');
    subplot(nr,nc,3); plot(M0-M00,'.'); title('M0 residual');
    
%     % get a set of points evenly distributed on the wedge
%     q = 4;
%     [gamma,delta] = getspheregrid([-30 30 -90 90],q,q,1);
%     figure; plot(gamma,delta,'.'); axis equal;
%     xlabel('gamma angle (CLVD component)');
%     ylabel('delta angle (isotroic combinent)');
%     
%     n = length(gamma);
%     M0 = randi([1 10],1,n) * 1e15;
%     %gamma = [0 0 0]; delta = [-90 0 90]; M0 = [1 1 1];
%     %gamma = [-30:5:30]'; delta = 45*ones(13,1); M0 = ones(13,1);
%     lam = lune2lam(gamma,delta,M0);
%     
%     % check M0
%     M0check = CMT2m0(1,[lam ; zeros(3,n)]);
%     figure; plot(M0-M0check,'.'); title('residuals for M0');
    
    %---------------------
    
    deg = 180/pi;
    n = 100;
    gamma = linspace(-30,30,n);
    delta = zeros(1,n);              % DEVIATORIC
    M0 = 1/sqrt(2) * ones(1,n);     % KEY NORMALIZATION
    lam = lune2lam(gamma,delta,M0);
    lam1 = lam(1,:);
    lam2 = lam(2,:);
    lam3 = lam(3,:);
    
    % NOTE: lambdas must be taken as unit vectors (M0 = 1/sqrt(2))
    epsilonhat = lam2 ./ max(abs(lam([1 3],:)));
    rbailey = sqrt(6)/2 * lam2;
    %rkagan = -3/2*sqrt(6)*lam1.*lam2.*lam3;
    I2 = -1/2*(lam1.^2 + lam2.^2 + lam3.^2);    % = -0.5 here
    rkagan = 1/2*lam1.*lam2.*lam3 .* (-1/3*I2).^(-3/2);
    rkaganmod = (-1/2) * rkagan;  
    %rkaganmod = -3/2*sqrt(6)*lam1.*lam2.*lam3;     % Bailey2009 Eq. (C2)
    %rkagan = 3*sqrt(6)*lam1.*lam2.*lam3;
    
    % check
    %kk = 47; det(diag(lam(:,kk))), lam1(kk)*lam2(kk)*lam3(kk)
    
    figure; nr=3; nc=1; 
    
    subplot(nr,nc,1); hold on; %grid on;
    plot(rbailey,rkaganmod,'r-');
    plot(rbailey,-epsilonhat,'b-');
    legend('Kagan1985mod','epsilon','location','east');
    xlabel('rCLVD of Bailey2009 (M0 = 1/\sqrt{2})');
    ylabel('measure of CLVD component');
    title('Bailey et al. (2009), Figure C1');
    
    subplot(nr,nc,2);
    plot(epsilonhat,rbailey-epsilonhat,'b-'); grid on;
    xlabel('epsilonhat');
    ylabel('Bailey2009 - epsilonhat');
    title('Bailey et al. (2009) is not identical to epsilonhat');
    
    subplot(nr,nc,3); hold on; grid on;
    plot(gamma,epsilonhat,'b-',gamma,rbailey,'r-',gamma,rkaganmod,'k-');
    legend('epsilon-hat','Bailey2009','Kagan1985mod','location','northwest');
    xlabel('gamma of TapeTape2012, deg');
    ylabel('measure of CLVD component');
%     subplot(nr,nc,2); hold on; grid on;
%     plot(gamma,atan(epsilonhat)*deg - gamma,'b-');
%     plot(gamma,atan(rbailey)*deg - gamma,'r-');
%     plot(gamma,atan(rkagan)*deg - gamma,'k-');
%     legend('epsilon-hat','Bailey2009','Kagan1985','location','northwest');
%     xlabel('angular difference');
%     ylabel('measure of CLVD component');
end

%==========================================================================
