function [gamma,delta] = uv2lune(u,v)
%UV2LUNE convert u-v coordinates to lune coordinates (gamma, delta)
%
% INPUT
%   u       n x 1 vector
%   v       n x 1 vector
%
% OUTPUT
%   gamma   n x 1 vector of gamma angles, degrees
%   delta   n x 1 vector of delta angles, degrees
%
% From Tape and Tape (2015 GJI) "A uniform parameterization for moment tensors"
%
% calls v2gamma.m, u2beta.m
%

disp('entering uv2lune.m');

deg = 180/pi;

u = u(:);
v = v(:);

% output in radians
gamma = v2gamma(v);
beta = u2beta(u);

% convert to degrees
gamma = gamma*deg;
beta = beta*deg;
delta = 90 - beta;

%==========================================================================

if 0==1
    clear, clc, close all
    n = 100;
    u = linspace(0,3*pi/4,n)';
    v = linspace(-1/3,1/3,n)';
    [gamma,delta] = uv2lune(u,v);
    beta = 90 - delta;
    
    figure; nr=2; nc=1;
    subplot(nr,nc,1); hold on; plot(u,beta,'b.-');
    ylabel('\beta, degrees'); xlabel('u'); grid on;
    subplot(nr,nc,2); hold on; plot(v,gamma,'b.-');
    ylabel('\gamma, degrees'); xlabel('v'); grid on;
end

%==========================================================================
