function [M,v,w,kappa,sigma,h] = uniformMT(n,gamma0,delta0)
%UNIFORMMT generate a grid of uniformly spaced (full) moment tensors
%
% INPUT
%   n       length(1): number of moment tensors, randomly generated
%           length(5): number of increments in v,w,kappa,sigma,h
%   gamma0  OPTIONAL: lune longitude, degrees (-30 to 30)
%   delta0  OPTIONAL: lune latitude, degrees (-90 to 90)
%
% See Tape and Tape, 2015, GJI, "A uniform parameterization of moment tensors"
%
% See examples in run_uniformMT.m
%
% Carl Tape, 22-July-2015
%

deg = 180/pi;

FIXEDLAMBDA = false;
if nargin==3
    FIXEDLAMBDA = true;
    disp(sprintf('uniformMT.m: generate uniform moment tensors for lune point (\\gamma = %.2f, \\delta = %.2f)',gamma0,delta0));
end

% calculate the number of points in the set
switch length(n)
    case 1, nq = n;
    case 5
        if FIXEDLAMBDA, nq = prod(n(3:5)); else nq = prod(n); end
    otherwise
        error('n must be length 1 or 5'); 
end
disp(sprintf('%i points in the set',nq));

% min and max limits for each parameter
v1 = -1/3;      v2 = 1/3;       % similar to lune longitude (CLVD)
w1 = -3*pi/8;   w2 = 3*pi/8;    % similar to lune latitude (ISO)
%u1 = 0;         u2 = 3*pi/4;    % similar to lune colatitude (ISO)
k1 = 0;         k2 = 360;       % strike
s1 = -90;       s2 = 90;        % slip (rake)
h1 = 0;         h2 = 1;         % cos(dip)

switch length(n)
    case 1      % random grid, full moment tensors
        disp(sprintf('uniformMT.m: %i random points',nq));
        v       = randomvec(v1,v2,nq);
        w       = randomvec(w1,w2,nq);
        kappa   = randomvec(k1,k2,nq);
        sigma   = randomvec(s1,s2,nq);
        h       = randomvec(h1,h2,nq);

    case 5      % uniform grid, full moment tensors
        nv = n(1);
        nw = n(2);
        nk = n(3);
        ns = n(4);
        nh = n(5);
        disp(sprintf('uniformMT.m: uniform grid: nv = %i, nw = %i, nk = %i, ns = %i, nh = %i',...
            nv,nw,nk,ns,nh));
        dk = (k2-k1)/nk;
        ds = (s2-s1)/ns;
        dh = (h2-h1)/nh;
        % note that with this sampling you NEVER allow for the endpoints of
        % the interval to be used in the grid search
        kvec  = [ (k1+dk/2) : dk : (k2-dk/2) ]';
        svec  = [ (s1+ds/2) : ds : (s2-ds/2) ]';
        hvec  = [ (h1+dh/2) : dh : (h2-dh/2) ]';
        if and(n(1)~=0,n(2)~=0)
            % full moment tensor
            dv = (v2-v1)/nv;
            dw = (w2-w1)/nw;
            vvec  = [ (v1+dv/2) : dv : (v2-dv/2) ]';
            wvec  = [ (w1+dw/2) : dw : (w2-dw/2) ]';
            [v,w,k,s,h] = ndgrid(vvec,wvec,kvec,svec,hvec);
            v = v(:);
            w = w(:);
        else
            % orientation only
            [k,s,h] = ndgrid(kvec,svec,hvec);
        end
        kappa = k(:);
        sigma = s(:);
        h = h(:);
    
    otherwise
        error('n must be length 1 or 5');
end

% calculate moment tensors
if FIXEDLAMBDA
    % fixed lune point
    gamma = gamma0*ones(nq,1);
    delta = delta0*ones(nq,1);
    v0 = gamma2v(gamma0/deg);
    u0 = beta2u((90-delta0)/deg);
    w0 = 3*pi/8 - u0;
    v = v0*ones(nq,1);
    w = w0*ones(nq,1);
else
    [gamma,delta] = rect2lune(v,w);
end
M0 = 1;
theta = acos(h) * deg;
M = TT2CMT(gamma,delta,M0,kappa,theta,sigma);

if length(n)==5
    disp(sprintf('uniformMT.m: %i points in the regular grid',length(kappa)));
    uv = unique(v); uw = unique(w); uk = unique(kappa); us = unique(sigma); uh = unique(h);
    % display lune lon/lat angles in addition to u and h
    if ~FIXEDLAMBDA
        uu = 3*pi/8 - uw;
        ub = u2beta(uu)*deg;
        ud = 90 - ub;
        ug = v2gamma(uv)*deg;
        disp(sprintf('%i unique(v) (dv = %.4f):',length(uv),uv(2)-uv(1))); disp([uv ug]);
        disp(sprintf('%i unique(w) (dw = %.4f):',length(uw),uw(2)-uw(1))); disp([uw ud]);
    end
    disp(sprintf('%i unique(kappa) (dk = %.1f):',length(uk),uk(2)-uk(1))); disp(uk);
    disp(sprintf('%i unique(sigma) (ds = %.1f):',length(us),us(2)-us(1))); disp(us);
    % display dip angle in addition to h
    disp(sprintf('%i unique(h) (dh = %.4f):',length(uh),uh(2)-uh(1))); disp([uh acos(uh)*180/pi]);
else
    disp(sprintf('uniformMT.m: %i points in the random grid',length(kappa)));
end

%==========================================================================

if 0==1
    % an alternative algorithm to generating uniform moment tensors;
    % these are efficient but involve rand (or randn) and a projection
    close all, clc, clear
    if 0==1
        % use rand for the hyperbox, then toss out points outside the hypersphere
        n = 1e6;
        a = -1; b = 1; 
        M = -1 + (b-a)*rand(6,n);

        % throw out the points outside unit norm
        mnorm = sqrt( sum(M.^2) );
        ioutside = find(mnorm > 1);
        M(:,ioutside) = [];
        mnorm(ioutside) = [];
        nq = n-length(ioutside);
        disp(sprintf('%i/%i (%.2f) kept',nq,n,nq/n));

    else
        % use randn [PREFERRED ALGORITHM]
        nq = 1e5;
        sigma = 1;
        M = randn(1,nq*6);
        M = reshape(M,6,nq);
    end

    % project to the hypersphere
    mnorm = sqrt( sum(M.^2) );
    M = M ./ repmat(mnorm,6,1);

    % divide off-diagonal entries by sqrt(2)
    M(4:6,:) = M(4:6,:) / sqrt(2);
    
    % plot the omega distribution
    iref = randi(nq); Mref = M(:,iref);
    omega = CMT2omega(Mref,M);
    plot_omega(omega); 
end

%==========================================================================
