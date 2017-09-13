function [V,tocv,Vp,Vn,tocvp,tocvn] = Vgammaomega(gamma,omega,atol)
%VGAMMAOMEGA fractional volume as a function of lune longitude and omega
%
% INPUT
%   gamma   in radians
%   omega   in radians
%   atol    absolute tolerance for integral2
%
% OUTPUT
%   V       
%   tocv    
%   Vp      
%   Vn      
%   tocvp   
%   tocvn   
%
% calls zhatpgammaomega.m, zhatngammaomega.m
% called by Vgammaomega.m, run_Vhat.m
%

% LOWER-THAN-DEFAULT TOLERANCES ARE CRITICAL
%rtol = 1e-6; atol = 1e-10; smethod = 'auto';     % default
rtol = 1e-18;   % forces absolute tolerance  (integral2 will use the least restrictive between atol and rtol)
smethod = 'iterated';
if nargin==2, atol = 1e-10; end
if isempty(atol), atol = 1e-10; end

deg = 180/pi;

V=[]; tocv=[]; Vp=[]; Vn=[]; tocvp=[]; tocvn=[];

% modify input gamma and omega, if needed
n = max([numel(gamma) numel(omega)]);
if ~prod(size(gamma) == size(omega))
    if numel(gamma)==n
        if numel(omega)==1
            warning('assigning all %i omega values to be %.1f deg',n,omega*deg);
            omega = omega * ones(size(gamma));
        else
            whos gamma omega
            error('dimension mismatch for gamma and omega');
        end
    else
        if numel(gamma)==1
            warning('assigning all %i gamma values to be %.1f deg',n,gamma*deg);
            gamma = gamma * ones(size(omega));
        else
            whos gamma omega
            error('dimension mismatch for gamma and omega');
        end
    end
end
if ndims(gamma) > 2, error('gamma and omega must have dimension <= 2'); end
[na,nb] = size(gamma);
gamma = gamma(:);
omega = omega(:);

% Tape and Tape, 2016, Eq 46b
% integration is most efficient when using positive gamma
%gamma = abs(gamma);

tocv = NaN(n,1);

% integration limits
phimin = 0;
phimax = pi/2;
sigmamin = 0;
sigmamax = pi/2;

Vp = NaN(n,1);
Vn = NaN(n,1);
tocvp = NaN(n,1);
tocvn = NaN(n,1);
for ii=1:n
    % define functions
    Vposint = @(phi,sigma) ( 1 - zhatpgammaomega(phi,sigma,gamma(ii),omega(ii)) );
    Vnegint = @(phi,sigma) ( 1 + zhatngammaomega(phi,sigma,gamma(ii),omega(ii)) );

    % perform numerical integration (and time it)
    tic
    %Vp(ii) = integral2(Vposint,phimin,phimax,sigmamin,sigmamax);
    Vp(ii) = integral2(Vposint,phimin,phimax,sigmamin,sigmamax,'method',smethod,'RelTol',rtol,'AbsTol',atol);
    tocp = toc;
    tic
    %Vn(ii) = integral2(Vnegint,phimin,phimax,sigmamin,sigmamax);
    Vn(ii) = integral2(Vnegint,phimin,phimax,sigmamin,sigmamax,'method',smethod,'RelTol',rtol,'AbsTol',atol);
    tocn = toc;

    tocvp(ii) = tocp;
    tocvn(ii) = tocn;
    tocv(ii) = tocp + tocn;
    disp(sprintf('Vgammaomega.m [atol %.1e]: gamma = %6.2f omega [%4i/%4i] = %5.2f // %.2f s = %.2f(p) + %.2f(n)',...
        atol,gamma(ii)*deg,ii,n,omega(ii)*deg,tocv(ii),tocvp(ii),tocvn(ii)));
end

Vn = (2/pi^2)*Vn;
Vp = (2/pi^2)*Vp;
V  = Vn + Vp;

%==========================================================================
% EXAMPLE

if 0==1
    % check for a few omega values for V_gamma and V_-gamma
    n = 10;
    deg = 180/pi;
    gamma = 10/deg;
    omega = linspace(0,pi,n);
    [V,tocv,Vp,Vn,tocvp,tocvn] = Vgammaomega(gamma,omega);
    V, sum(tocv)/60  
%  V =
%    0.000000000000000
%    0.005332564539437
%    0.048099363958304
%    0.167252332027672
%    0.351393815328008
%    0.610938166660600
%    0.858630111833154
%    0.972353227459042
%    1.000000000000000
%    1.000000000000000
% ans =
%    0.737186683333333 
%
    [V,tocv,Vp,Vn,tocvp,tocvn] = Vgammaomega(-gamma,omega);
    V, sum(tocv)/60
% V =
%    0.000000000000000
%    0.005332564528060
%    0.048099362516572
%    0.167252355753277
%    0.351393829043454
%    0.610938143938500
%    0.858630111827973
%    0.972353227241735
%    1.000000000000000
%    1.000000000000000
% ans =
%    2.038369466666667

    % plot
    figure; plot(omega*deg,V,'o');
    title(sprintf('\\gamma = %.1f deg',gamma*deg));
    xlabel('\omega, deg'); ylabel('Vhat_\gamma(\omega)');
end

%==========================================================================
