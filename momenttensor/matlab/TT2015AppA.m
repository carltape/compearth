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
% Carl Tape, 2015-07-14
%

clear, close all

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
    % note: M is a matrix of numbers, NOT a moment tensor with an assumed basis.
    %       But if you want to recover the input values of kappa, sigma,
    %       and theta from M, then M needs to be assumed to have a basis of north-west-up.
    M = U * diag(Lambda) * inv(U);
end

disp('  ');
disp('compare output with TapeTape2015, Appendix A:');
disp(sprintf('(beta, gamma) = (%.3f, %.3f) = (%.1f deg, %.1f deg)',beta/deg, gamma/deg, beta, gamma));
disp(sprintf('Lambda = (%.3f, %.3f, %.3f)',Lambda));
disp(sprintf('theta = %.3f = %.3f deg',theta/deg,theta));
disp('note: the basis for U and M is north-west-up');
U, M

%--------------------------------------------------------------------------
% Note about moment tensor basis
% The choice of basis will matter when checking the entries in the moment
% tensor. The example calculations in TT2015 and in TT2012 use north-west-up.
% This is not stated explicitly in TT2015 Appendix A, but here are some passages that do:
%   TT2015, Section S4.2: "(The vector e1 is assumed north and e3 is up.)"
%   TT2015, Figure S3 caption: "If the vector e1 points north and e3 points up,
%           then κ, σ, θ are the strike, slip, and dip angles for
%           one of the two fault planes of [Λ]_VY−π/4 ."
%   TT2012, p. 10: "...where the standard basis is e1 (north), e2 (west), e3 (zenith)"
%   TT2012, Figure 14.
% see convert_MT.m for referenced publications that use various basis conventions

% comment this break out to run the rest
break

disp('moment tensor in south-east-up (GCMT convention):');
Mseu = convert_MT(3,1,Mvec2Mmat(M,0));
Mvec2Mmat(Mseu,1)

%==========================================================================
