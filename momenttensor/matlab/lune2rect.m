function [v,w] = lune2rect(gamma,delta)
%LUNE2RECT convert lune coordinates (gamma, delta) to v-w
%
% INPUT
%   gamma   n x 1 vector of gamma angles, degrees
%   delta   n x 1 vector of delta angles, degrees
%
% OUTPUT
%   v       n x 1 vector (like gamma)
%   w       n x 1 vector (like delta)
%
% From Tape and Tape (2015 GJI) "A uniform parameterization for moment tensors"
%
% calls gamma2v.m, beta2u.m
%

disp('entering lune2rect.m');

deg = 180/pi;

% convert to radians
delta = delta(:) / deg;     % (column vector)
gamma = gamma(:) / deg;     % (column vector)
beta = pi/2 - delta;        % colatitude

v = gamma2v(gamma);
u = beta2u(beta);
w = 3*pi/8 - u;

%==========================================================================

if 0==1
    clear, clc, close all
    n = 100;
    delta = linspace(90,-90,n)';
    gamma = linspace(-30,30,n)';
    [v,w] = lune2rect(gamma,delta);
    
    figure; nr=2; nc=1;
    subplot(nr,nc,1); hold on; plot(gamma,v,'b.-');
    xlabel('\gamma, degrees'); ylabel('v'); grid on; xlim([-30 30]);
    subplot(nr,nc,2); hold on; plot(delta,w,'b.-');
    xlabel('\delta, degrees'); ylabel('w'); grid on; xlim([-90 90]);
end

%==========================================================================
