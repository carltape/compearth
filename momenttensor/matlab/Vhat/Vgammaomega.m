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
% calls mugammaomega.m, rgammaomega.m, shatpgammaomega.m, shatngammaomega.m
% called by Vgammaomega.m, run_Vhat.m
%

% LOWER-THAN-DEFAULT TOLERANCES ARE CRITICAL
%rtol = 1e-6; atol = 1e-10; smethod = 'auto';     % default
rtol = 1e-18;   % forces absolute tolerance  (integral2 will use the least restrictive between atol and rtol)
smethod = 'iterated';
if nargin==2, atol = 1e-10; end
if isempty(atol), atol = 1e-10; end

bportion = false;

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
if bportion
    % mu_gamma(omega)
    mugamma = mugammaomega(gamma,omega);
    
    % integration limies
    mumin = 0;
    numin = 0;
    numax = pi;

    V = NaN*ones(n,1);
    for ii=1:n
        tic
        if omega(ii) > pi/3 + 2*gamma(ii)
            continue
        else
            mumax = mugamma(ii);
            Vtemp = @(mu,nu) (1 - rgammaomega(mu,nu,gamma(ii),omega(ii)) );
            V(ii) = integral2(Vtemp,mumin,mumax,numin,numax);
        end
        tocv(ii) = toc;
    end
    V = (2/pi^2)*V;

else
    % integration limits
    mumin = 0;
    mumax = pi/2;
    numin = 0;
    numax = pi;

    Vp = NaN(n,1);
    Vn = NaN(n,1);
    tocvp = NaN(n,1);
    tocvn = NaN(n,1);
    for ii=1:n
        % define functions
        Vposint = @(mu,nu) ( 1 - shatpgammaomega(mu,nu,gamma(ii),omega(ii)) );
        Vnegint = @(mu,nu) shatngammaomega(mu,nu,gamma(ii),omega(ii));
        
        % perform numerical integration (and time it)
        tic
        %Vp(ii) = integral2(Vposint,mumin,mumax,numin,numax);
        Vp(ii) = integral2(Vposint,mumin,mumax,numin,numax,'method',smethod,'RelTol',rtol,'AbsTol',atol);
        tocp = toc;
        tic
        %Vn(ii) = integral2(Vnegint,mumin,mumax,numin,numax);
        Vn(ii) = integral2(Vnegint,mumin,mumax,numin,numax,'method',smethod,'RelTol',rtol,'AbsTol',atol);
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
end

%==========================================================================
