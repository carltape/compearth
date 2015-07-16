function [u,v] = lune2uv(gamma,delta)
%LUNE2UV convert lune coordinates (gamma, delta) to u-v
%
% INPUT
%   gamma   n x 1 vector of gamma angles, degrees
%   delta   n x 1 vector of delta angles, degrees
%
% OUTPUT
%   u       n x 1 vector
%   v       n x 1 vector
%
% From Tape and Tape (2015 GJI) "A uniform parameterization for moment tensors"
%
% calls gamma2u.m, beta2v.m
%

disp('entering lune2uv.m');

deg = 180/pi;

% in radians
delta = delta(:) / deg;   % (column vector)
gamma = gamma(:) / deg;   % (column vector)
beta = pi/2 - delta;  % colatitude

v = gamma2v(gamma);
u = beta2u(beta);

%==========================================================================

if 0==1
    %%
    clear, clc, close all
    n = 100;
    delta = linspace(90,-90,n)';
    gamma = linspace(-30,30,n)';
    beta = 90 - delta;
    [u,v] = lune2uv(gamma,delta);
    
    figure; nr=2; nc=1;
    subplot(nr,nc,1); hold on; plot(beta,v,'b.-');
    xlabel('\beta, degrees'); ylabel('v'); grid on;
    subplot(nr,nc,2); hold on; plot(gamma,u,'b.-');
    xlabel('\gamma, degrees'); ylabel('u'); grid on;
end

%==========================================================================
