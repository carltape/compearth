%
% TT2015AppA.m
%
% This script reproduces the calculations in Appendix A of
% Tape and Tape (2015), "A uniform parameterization of moment tensors"
%
% It may be helpful in understanding the equations within the paper.
%
% It uses several functions stored within the github repository 'compearth'.
%
% Carl Tape, 7/14/2015
%

clear, clc, close all

deg = 180/pi;

% example coordinates
% note: the scripts assume that angles are in degrees (not radians)
u     = 3*pi/8;
v     = -1/9;
kappa = 4*pi/5 * deg;
sigma = -pi/2 * deg;
h     = 3/4;
theta = acos(h) * deg;

% w is like lune latitude delta, whereas u is like lune colatitude beta
w = 3*pi/8 - u;

% lune coordinates
[gamma,delta] = rect2lune(v,w);
beta = 90 - delta;

% TWO OPTIONS to compute Lambda, U, M
if 0==1     % OPTION 1: use existing functions
    % note: U is in south-east-up basis, M is in up-south-east basis
    M0 = 1/sqrt(2);     % to give |Lambda| = 1
    [Muse,Lambda,Useu] = TT2CMT(gamma,delta,M0,kappa,theta,sigma);
    % convert U from south-east-up to north-west-up
    U = convertv(5,3,Useu);
    % convert moment tensor from up-south-east to north-west-up
    M = convert_MT(1,3,Muse);
    
    M = Mvec2Mmat(M,1);

else        % OPTION 2: implement equations from TapeTape2015 directly
    % Eq 7
    R = 1/sqrt(6) * [sqrt(3) 0 -sqrt(3) ; -1 2 -1 ; sqrt(2) sqrt(2) sqrt(2)]';
    sb = sin(beta/deg);
    cb = cos(beta/deg);
    Lambda = R * [ sin(beta/deg)*cos(gamma/deg) sin(beta/deg)*sin(gamma/deg) cos(beta/deg) ]';
    
    % Eq 9-10
    Zsigma = rotmat(sigma,3);
    Xtheta = rotmat(theta,1);
    Zkappa = rotmat(-kappa,3);
    V = Zkappa*Xtheta*Zsigma;
    Yrot = rotmat(-45,2);
    U = V*Yrot;
    
    % Eq 4
    M = U * diag(Lambda) * inv(U);
end

disp('  ');
disp('compare output with TapeTape2015, Appendix A:');
disp(sprintf('(beta, gamma) = (%.3f, %.3f) = (%.1f deg, %.1f deg)',beta/deg, gamma/deg, beta, gamma));
disp(sprintf('Lambda = (%.3f, %.3f, %.3f)',Lambda));
disp(sprintf('theta = %.3f = %.3f deg',theta/deg,theta));
U, M

%==========================================================================
