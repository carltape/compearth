function [gamma,delta] = rect2lune(v,w)
%RECT2LUNE convert v-w coordinates to lune coordinates (gamma, delta)
%
% INPUT
%   v       n x 1 vector (like gamma)
%   w       n x 1 vector (like delta)
%
% OUTPUT
%   gamma   n x 1 vector of gamma angles, degrees
%   delta   n x 1 vector of delta angles, degrees
%
% From Tape and Tape (2015 GJI) "A uniform parameterization for moment tensors"
%
% calls v2gamma.m, u2beta.m
%

disp('entering rect2lune.m');

deg = 180/pi;

v = v(:);
w = w(:);
u = 3*pi/8 - w;

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
    v = linspace(-1/3,1/3,n)';
    w = linspace(-3*pi/8,3*pi/8,n)';
    [gamma,delta] = rect2lune(v,w);
    
    figure; nr=2; nc=1;
    subplot(nr,nc,1); hold on; plot(v,gamma,'b.-');
    ylabel('\gamma, degrees'); xlabel('v'); grid on;
    subplot(nr,nc,2); hold on; plot(w,delta,'b.-');
    ylabel('\delta, degrees'); xlabel('w'); grid on;
end

%==========================================================================
