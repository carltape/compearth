function omega_mod = omegacritmod(omega)
%OMEGACRITMOD set of unique critical angles, excluding w0
%
% see lam2omegacrit.m

tol = 1e-6;

% uniquetol appeared in Matlab R2015
omega_mod = uniquetol(omega,tol);

% ensure that omega0 = 0 is included
omega_mod(omega_mod < tol) = [];
omega_mod = [0 omega_mod(:)'];
