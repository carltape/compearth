function beta = u2beta(u)
%U2BETA u(beta) for lune colatitude
%
% INPUT
%   u       n x 1 vector (u = [0, 3*pi/4])
%
% OUTPUT
%   beta    n x 1 vector of lune colatitudes, radians (beta = [0, pi])
%
% From Tape and Tape (2015 GJI) "A uniform parameterization for moment tensors"
%

bdisplayinfo = false;

% u(beta) is not an analytical solution. The accuracy of the numerical
% solution is controlled by two steps below:
%   1. the number of interpolation points
%   2. the parameters set in optimset for fzero
bfzero = false;
if bfzero
    n0 = 10;
    tol = 1e-12;
    emethod = 'linear';
    if bdisplayinfo, disp(sprintf('u2beta.m: using fzero with tol = %.1e',tol)); end
else
    % For some reason, it is NOT possible to simply increase n0 to achieve
    % arbitrarily good accuracy. I run into this error:
    %    Error using griddedInterpolant
    %    The grid vectors are not strictly monotonic increasing.
    n0 = 1000;
    emethod = 'pchip';
    if bdisplayinfo, disp(sprintf('u2beta.m: using %s interpolation with %i points',emethod,n0)); end
end

u = u(:);
%n = length(u);

% we could define an in-line function instead of using beta2u.m
%f = @(x)( (3/4)*x - (1/2)*sin(2*x) + (1/16)*sin(4*x) );

% interpolate
beta0 = linspace(0,pi,n0)';
u0 = beta2u(beta0);
%u0 = f(beta0);
%beta = NaN(n,1);
%for ii=1:n, beta(ii) = interp1(u0,beta0,u(ii)); end
beta = interp1(u0,beta0,u,emethod);

% in some cases there can be NaN appearing
%sum(isnan(beta))

if bfzero
    displaytype = 'off';  % off, iter, notify, final
    options = optimset('TolX',tol,'TolFun',tol,...
        'MaxIter',20,'MaxFunEvals',20,'Display',displaytype);
    for ii=1:n
        x0 = beta(ii);
        beta(ii) = fzero( @(x)(beta2u(x)-u(ii)), x0, options);
        %beta(ii) = fzero( @(x)(f(x)-u(ii)), x0, options);
    end 
end

%==========================================================================

if 0==1
    % how long does it take to interpolate a million points
    tic
    n = 1e6;
    u = 3*pi/4 * rand(n,1);
    beta = u2beta(u);
    disp(sprintf('%.4f s seconds to get %.0e values of beta (lune colatitude) from u',toc,n));
    
    % try a linearly spaced input vector
    n = 191;
    n = 1e6;
    % from beta to u and back
    beta = linspace(0,pi,n)';
    u = beta2u(beta);
    figure; nr=2; nc=1;
    subplot(nr,nc,1); plot(beta,u); xlabel('\beta, radians'); ylabel('u');
    % check
    beta_check = u2beta(u);
    norm(beta - beta_check)
    subplot(nr,nc,2); plot(beta,beta-beta_check,'.-'); grid on;
    xlabel('\beta, radians'); ylabel('beta residual, radians');
    title(sprintf('mean(abs(betadiff)) = %.2e',mean(abs(beta-beta_check))));
    
    % from u to beta and back
    u = linspace(0,3*pi/4,n)';
    tic
    beta = u2beta(u);
    toc
    disp(sprintf('%.3f s to calculate %.0e values of beta from u',toc,n));
    figure; nr=2; nc=1;
    subplot(nr,nc,1); plot(u,beta); ylabel('\beta, radians'); xlabel('u');
    % check
    u_check = beta2u(beta);
    ures = u - u_check;
    subplot(nr,nc,2); plot(u,ures,'.-'); grid on;
    xlabel('u'); ylabel('u residual');
    title(sprintf('mean(abs(ures)) = %.2e',mean(abs(ures)) ));
    
%     % experimenting with some functions that DO have analytical inverses;
%     % these could be useful for picking an initial value for the numerical
%     % solver, but using basic interpolation (as in u2beta.m) is simpler and
%     % even better
%     figure; nr=2; nc=1;
%     subplot(nr,nc,1); plot(beta,u);
%     xlabel('\beta, radians'); ylabel('u');
%     % guess a best-fitting function that has a simple inverse
%     hold on;
%     f1 = pi/8;
%     umax = 3*pi/4;
%     ufit = umax/2 + atan((beta - pi/2)/f1)
%     betafit = pi/2 + f1*tan(u - umax/2);
%     plot(beta,ufit,'r--');
%     axis equal; axis tight;
%     
%     subplot(nr,nc,2); plot(u,beta);
%     xlabel('u'); ylabel('\beta, radians');
%     hold on;
%     plot(u,betafit,'r--');
%     %axis equal; axis tight;
    
end

%==========================================================================
