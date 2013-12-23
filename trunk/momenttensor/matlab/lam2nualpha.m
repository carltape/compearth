function [nu,alpha] = lam2nualpha(lam)
%LAM2NUALPHA converts eigenvalues to nu and alpha
%
% INPUT
%   lam     3 x n set of eigenvalue triples
%
% OUTPUT
%   nu      n x 1 vector of Poisson parameters (unitless)
%   alpha   n x 1 vector of angles between fault normal and slip vector, degrees
%
% Reverse function is nualpha2lam.m
% See Tape and Tape (2013), "The classical model for moment tensors"
% 
% Carl Tape, 08-Jan-2013
%

lam = lamsort(lam);

lam1 = lam(1,:);
lam2 = lam(2,:);
lam3 = lam(3,:);

% TT2013, Eqs 32ab
alpha = 180/pi* acos( (lam1 - 2*lam2 + lam3) ./ (lam1 - lam3) );
nu = lam2 ./ (lam1 + lam3);

alpha = alpha(:);
nu = nu(:);