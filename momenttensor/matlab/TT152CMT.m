function [M,lam,U] = TT152CMT(rho,u,v,kappa,sigma,h)
%TT152CMT convert parameters of TapeTape2015 into moment tensors
%
% INPUT
%   rho         norm of the moment tensor (note: rho = sqrt(2)*M0)
%   u           coordinate proportional to lune colatitude beta
%   v           coordinate proportional to lune longitude gamma
%   kappa       strike angle, degrees
%   sigma       slip (or rake) angle, degrees
%   h           cos(dip)
%
% OUTPUT
%   M           6 x n set of moment tensors in CMT convention (UP-SOUTH-EAST)
%               M = [Mrr Mtt Mpp Mrt Mrp Mtp]; r=up, theta=south, phi=east
%   lam         3 x n set of eigenvalues
%   U           3 x 3 x n set of bases in SOUTH-EAST-UP convention
%
% Note that the basis for M and U are different.
%
% Reverse program for CMT2TT15.m
% See WTape and CTape (2015) "A uniform parameterization for moment tensors"
%
% Carl Tape, 8/7/2015
%

theta = acos(h)*180/pi;
M0 = rho/sqrt(2);
[gamma,delta] = uv2lune(u,v);
[M,lam,U] = TT2CMT(gamma,delta,M0,kappa,theta,sigma);
    
%==========================================================================
% EXAMPLES

if 0==1
    % using TT152CMT.m
    rho = sqrt(2); u = 3*pi/8; v = 0;
    kappa = 320; sigma = 20; h = cos(10*pi/180);
    M = TT152CMT(rho,u,v,kappa,sigma,h)
    [rhox,ux,vx,kappax,sigmax,hx] = CMT2TT15(M);
    disp([rho rhox u ux v vx]);
    disp([kappa kappax sigma sigmax h hx]);
    
    % using TT2CMT.m
    M0 = 1; gamma = 0; delta = 0;
    kappa = 320; theta = 10; sigma = 20;
    M = TT2CMT(gamma,delta,M0,kappa,theta,sigma)
    [gammax,deltax,M0x,kappax,thetax,sigmax] = CMT2TT(M);
    disp([M0 M0x gamma gammax delta deltax]);
    disp([kappa kappax sigma sigmax theta thetax]);
end

%==========================================================================
